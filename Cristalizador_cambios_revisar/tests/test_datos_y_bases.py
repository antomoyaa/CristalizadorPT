"""Pruebas del conversor de datos y de la alineación de bases (perfiles de params)."""
import contextlib
import io
import json
import os
import re

import numpy as np
import pandas as pd
import pytest

from msmpr_core import params as P
from msmpr_core.datos import wide_to_long
from msmpr_core.plant import Plant

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_ANT = os.path.join(_REPO, "Anteriores", "Datos_antiguos")


def _default_nb(env, nombre):
    """Ruta por variable de entorno o, si no esta definida, la copia versionada en Anteriores/Datos_antiguos (si existe)."""
    ruta = os.environ.get(env, "") or os.path.join(_ANT, nombre)
    return ruta if os.path.isfile(ruta) else ""


CAL_NB = _default_nb("CALIBRACION_NB", "calibracion_do_crystal.ipynb")
PARETO_NB = _default_nb("PARETO_NB", "MSMPR_KNO3_DPBE_sensibilidad_Pareto.ipynb")


def _wide(n=6):
    return pd.DataFrame({"run": [f"R{i:02d}" for i in range(n)], "d": np.linspace(0, .1, n), "T_C": 35 + np.arange(n),
                         "tau_h": 2.0 + .1 * np.arange(n), "crystal_flow_kg_h": 1e4 + 100 * np.arange(n), "sd_crystal_flow_kg_h": 150.0,
                         "D50_um": 170.0 + np.arange(n), "sd_D50_um": 5.0, "C_kg_m3": 176.0, "sd_C_kg_m3": 0.5,
                         "impurity_ppm": 0.04, "sd_impurity_ppm": 0.002})


def test_converter_requires_explicit_declarations():
    with pytest.raises(ValueError):
        wide_to_long(_wide(), "A9999", "synthetic", 660.0, 220.0)
    with pytest.raises(ValueError):
        wide_to_long(_wide(), "A3110", "quizas", 660.0, 220.0)


def test_converter_omits_impurity_unless_basis_confirmed():
    with pytest.warns(UserWarning):
        med, ope, _ = wide_to_long(_wide(), "A3110", "synthetic", 660.0, 220.0)
    assert "impurity_ppm" not in set(med.variable)
    med, ope, _ = wide_to_long(_wide(), "A3110", "synthetic", 660.0, 220.0, 3.8e-3, 1040.0, impurity_basis_confirmed=True)
    assert set(med[med.variable == "impurity_ppm"].basis) == {"dry_ClO4"}
    assert {"feed_ClO4_kg_m3", "rho_liquid_kg_m3"} <= set(ope.variable)


def test_converter_operation_starts_at_zero_and_covers_every_observation():
    med, ope, mapa = wide_to_long(_wide(), "A3110", "synthetic", 660.0, 220.0, shuffle_seed=3)
    for var in ("temperature_C", "residence_time", "feed_dilution_fraction"):
        t = ope[ope.variable == var].time.values
        assert t.min() == 0.0 and set(med.time) <= set(t)
    assert sorted(mapa.run) == sorted(_wide().run)


@pytest.mark.skipif(not (CAL_NB and PARETO_NB), reason="falta calibracion_do_crystal.ipynb: definir CALIBRACION_NB (o copiarlo a Anteriores/Datos_antiguos) y PARETO_NB")
def test_converted_file_passes_calibration_notebook_validators():
    import matplotlib
    matplotlib.use("Agg")
    nb = json.load(open(CAL_NB, encoding="utf-8"))
    ns = {"display": lambda *a, **k: None, "__name__": "__main__"}
    for i, c in enumerate(nb["cells"][:7]):
        if c["cell_type"] != "code":
            continue
        s = "".join(c["source"])
        s = re.sub(r"PROJECT\s*=\s*Path\(.*\)", f"PROJECT = Path(r'{os.path.dirname(PARETO_NB)}')", s)
        s = re.sub(r"MODE\s*=\s*'demo'", "MODE = 'experimental'", s)
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(s, f"cell{i}", "exec"), ns)
    med, ope, _ = wide_to_long(_wide(), "A3110", "synthetic", 660.0, 220.0, 3.8e-3, 1040.0, impurity_basis_confirmed=True)
    blocks = [ns["tidy_sheet"](med, "observation"), ns["tidy_sheet"](ope, "input")]
    D = ns["normalize_data"](pd.concat(blocks, ignore_index=True), "demo")
    ns["validate_campaign"](D, "A3110", "estable")
    with pytest.raises(ValueError):                       # el modo experimental debe rechazar sintéticos
        ns["normalize_data"](pd.concat(blocks, ignore_index=True), "experimental")


def test_calibration_profile_matches_pareto_notebook():
    if not PARETO_NB:
        pytest.skip("definir PARETO_NB")
    assert P.load_base_from_notebook(PARETO_NB) == P.BASE_CALIBRACION


def test_diff_bases_flags_what_conditions_calibration():
    diff = P.diff_bases().set_index("parametro")
    assert diff.loc["solubility_slope_K_inv", "condiciona_calibracion"]
    assert diff.loc["rho_liquid_kg_m3", "condiciona_calibracion"]
    assert not diff.loc["Kd", "condiciona_calibracion"]


def test_apply_calibration_refuses_mismatched_basis():
    cal = {"kG_m_s": 2.5e-7, "kB_m3_s_inv": 8e9, "Kd": 0.014}
    with pytest.raises(ValueError):
        P.apply_calibration(P.make_params("literatura"), cal)
    q = P.apply_calibration(P.make_params("calibracion"), cal, "A3110")
    assert q["kG_m_s"] == 2.5e-7 and q["cake_moisture_fraction"] == 0.0 and q["calibrated_for"] == "A3110"
    with pytest.raises(KeyError):
        P.apply_calibration(P.make_params("calibracion"), {"Cstar_ref_kg_m3": 1.0})


def test_calibration_profile_plant_converges_and_closes():
    pl = Plant(P.make_params("calibracion"))
    x, u, d = pl.steady_state()
    assert abs(pl.metrics(x, u, d)["Y_cierre"]) < 1e-8

def test_planta_sqm_profile_has_mg_and_meets_recovery_criterion():
    """Verifica que el perfil 'planta_sqm' incluya Mg=32.1115 g/L y cumpla 0.30 <= Y <= 0.70 (AGENTS.md)."""
    p_sqm = P.make_params("planta_sqm")
    assert p_sqm["feed_stock_kg_m3"]["Mg"] == pytest.approx(32.1115, rel=1e-4)
    assert "Mg" in P.PROVENANCE["feed_stock_kg_m3"]
    pl = Plant(p_sqm)
    x, u, d = pl.steady_state(alpha=0.25)
    m = pl.metrics(x, u, d)
    assert 0.30 <= m["Y_sol"] <= 0.70
    assert abs(m["Y_cierre"]) < 1e-8

def test_convert_file_rechaza_datos_sinteticos_declarados_como_planta(tmp_path):
    """Un Excel con meta.origen = 'sintetico' no puede convertirse con origin distinto de 'synthetic'."""
    import pandas as pd
    from msmpr_core import datos
    wide = pd.DataFrame({"run": ["R01"], "d": [0.0], "T_C": [40.0], "tau_h": [0.6], "crystal_flow_kg_h": [1.0],
                         "D50_um": [200.0], "C_kg_m3": [180.0], "sd_crystal_flow_kg_h": [0.1], "sd_D50_um": [5.0], "sd_C_kg_m3": [0.5]})
    src = tmp_path / "d.xlsx"
    with pd.ExcelWriter(src) as xw:
        wide.to_excel(xw, sheet_name="msmpr", index=False)
        pd.DataFrame({"campo": ["origen"], "valor": ["sintetico"]}).to_excel(xw, sheet_name="meta", index=False)
    with pytest.raises(ValueError, match="sintetico"):
        datos.convert_file(str(src), str(tmp_path / "o.xlsx"), crystallizer="A3110", origin="plant", volume_m3=660.0, feed_kno3_kg_m3=525.0)
