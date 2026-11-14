"""Conversor: Excel ANCHO (datos_cristalizador.xlsx del notebook de estimación) -> formato LARGO de
calibracion_do_crystal.ipynb (hojas `mediciones` y `operacion`).

Qué NO hace (a propósito):
  * No asigna el equipo ni el origen por su cuenta: `crystallizer` ('A3110'/'A3210') y `origin`
    ('measured'/'synthetic') son obligatorios.
  * No inventa la operación que el Excel no trae: volumen, alimentación de KNO3 (y, si hay impureza, ClO4 de
    alimentación y densidad del licor) se piden como argumentos.
  * No incluye `impurity_ppm` salvo que declares la base con `impurity_basis_confirmed=True`
    (el notebook exige basis='dry_ClO4': ppm de masa de producto seco).
  * La hoja de solubilidad no tiene destino en el calibrador (la solubilidad queda fija); se copia como referencia.

Cada corrida del Excel es un estado estacionario independiente. Se pone como UNA campaña `steady` con una
fila de tiempo por corrida (t = 1, 2, ... h, en el orden dado). El calibrador reserva como validación las últimas
fechas: con un diseño factorial ordenado eso equivale a validar en los niveles más altos del último factor.
Usa `shuffle_seed` para mezclar el orden, o pasa los datos en orden cronológico real.

Uso:
    python -m msmpr_core.datos entrada.xlsx salida.xlsx --crystallizer A3110 --origin synthetic \
        --volume-m3 660 --feed-kno3 220 [--feed-clo4 3.799e-3 --rho-liquid 1040 --impurity-basis-confirmed]
"""
from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd

OBS_MAP = {   # columna ancha -> (variable, unidad, basis)
    "crystal_flow_kg_h": ("crystal_flow_kg_h", "kg/h", "model"),
    "D20_um": ("D20_um", "um", "volume"),
    "D50_um": ("D50_um", "um", "volume"),
    "D80_um": ("D80_um", "um", "volume"),
    "C_kg_m3": ("C_kg_m3", "kg/m3", "model"),
    "impurity_ppm": ("impurity_ppm", "ppm", "dry_ClO4"),
}
VARYING_INPUTS = {   # columna ancha -> (variable, unidad)
    "d": ("feed_dilution_fraction", "1"),
    "T_C": ("temperature_C", "degC"),
    "tau_h": ("residence_time", "h"),
}


def wide_to_long(wide: pd.DataFrame, crystallizer: str, origin: str, volume_m3: float, feed_kno3_kg_m3: float,
                 feed_clo4_kg_m3: float | None = None, rho_liquid_kg_m3: float | None = None,
                 impurity_basis_confirmed: bool = False, experiment: str = "estable", shuffle_seed: int | None = None,
                 time_step_h: float = 1.0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Devuelve (mediciones, operacion, mapa_corridas)."""
    if crystallizer not in ("A3110", "A3210"):
        raise ValueError("crystallizer debe ser 'A3110' o 'A3210'")
    if origin not in ("measured", "synthetic"):
        raise ValueError("origin debe ser 'measured' o 'synthetic'")
    need = ["run", "d", "T_C", "tau_h"]
    if [c for c in need if c not in wide.columns]:
        raise ValueError(f"Faltan columnas: {[c for c in need if c not in wide.columns]}")
    w = wide.reset_index(drop=True).copy()
    if shuffle_seed is not None:
        w = w.sample(frac=1.0, random_state=shuffle_seed).reset_index(drop=True)
    w["time"] = (np.arange(len(w)) + 1) * time_step_h
    common = dict(crystallizer=crystallizer, experiment=experiment, regime="steady", time_unit="h", origin=origin)

    # --- observaciones
    obs = []
    include_imp = "impurity_ppm" in w.columns
    if include_imp and not impurity_basis_confirmed:
        warnings.warn("impurity_ppm omitida: declara impurity_basis_confirmed=True si es ppm de masa de producto seco.")
        include_imp = False
    for col, (var, unit, basis) in OBS_MAP.items():
        if col not in w.columns or (col == "impurity_ppm" and not include_imp):
            continue
        sd_col = f"sd_{col}"
        for _, r in w.iterrows():
            if pd.isna(r[col]):
                continue
            obs.append({**common, "time": r["time"], "variable": var, "value": float(r[col]), "unit": unit,
                        "basis": basis, "sd": float(r[sd_col]) if sd_col in w.columns and pd.notna(r[sd_col]) else np.nan})
    mediciones = pd.DataFrame(obs)
    if mediciones.empty:
        raise ValueError("No quedó ninguna observación.")

    # --- operación: constantes en t=0; las que varían por corrida, en t=0 (= corrida 1) y en cada t_k (escalón causal)
    inputs = [{**common, "time": 0.0, "variable": "volume_m3", "value": volume_m3, "unit": "m3", "basis": "model", "sd": np.nan},
              {**common, "time": 0.0, "variable": "feed_KNO3_kg_m3", "value": feed_kno3_kg_m3, "unit": "kg/m3", "basis": "model", "sd": np.nan}]
    if include_imp:
        if feed_clo4_kg_m3 is None or rho_liquid_kg_m3 is None:
            raise ValueError("Con impureza se exigen feed_clo4_kg_m3 y rho_liquid_kg_m3.")
        inputs += [{**common, "time": 0.0, "variable": "feed_ClO4_kg_m3", "value": feed_clo4_kg_m3, "unit": "kg/m3", "basis": "model", "sd": np.nan},
                   {**common, "time": 0.0, "variable": "rho_liquid_kg_m3", "value": rho_liquid_kg_m3, "unit": "kg/m3", "basis": "model", "sd": np.nan}]
    for col, (var, unit) in VARYING_INPUTS.items():
        inputs.append({**common, "time": 0.0, "variable": var, "value": float(w.loc[0, col]), "unit": unit, "basis": "model", "sd": np.nan})
        for _, r in w.iterrows():
            inputs.append({**common, "time": r["time"], "variable": var, "value": float(r[col]), "unit": unit, "basis": "model", "sd": np.nan})
    operacion = pd.DataFrame(inputs)
    mapa = w[["run", "time", "d", "T_C", "tau_h"]].rename(columns={"time": "time_h"})
    return mediciones, operacion, mapa


def convert_file(src: str, dst: str, **kw) -> dict:
    xl = pd.ExcelFile(src)
    sheet = "msmpr" if "msmpr" in xl.sheet_names else xl.sheet_names[0]
    wide = xl.parse(sheet)
    if "meta" in xl.sheet_names:                       # el Excel declara su origen (hoja 'meta': campo/valor)
        meta = xl.parse("meta")
        declarado = str(meta.loc[meta["campo"] == "origen", "valor"].iloc[0]).strip().lower() if (meta["campo"] == "origen").any() else None
        if declarado == "sintetico" and kw.get("origin") != "synthetic":
            raise ValueError(f"{src} declara origen 'sintetico' en su hoja 'meta': usa origin='synthetic' (no se puede presentar como dato de planta).")
    med, ope, mapa = wide_to_long(wide, **kw)
    notes = pd.DataFrame({"nota": [
        f"Convertido de {src} (hoja '{sheet}') con msmpr_core.datos.",
        f"crystallizer={kw['crystallizer']} y origin={kw['origin']} los declaró el usuario; el archivo original no los trae.",
        "volume_m3, feed_KNO3_kg_m3 (y ClO4/densidad si hay impureza) NO vienen del Excel: se pasaron como argumentos.",
        "Cada corrida = un estado estacionario independiente; tiempo = orden de la corrida (h). Ver hoja 'corridas'.",
        "La hoja 'solubilidad' es solo de referencia: el calibrador no la lee (la solubilidad queda fija).",
        "Si origin=synthetic, el calibrador en modo 'experimental' lo rechaza (prohíbe mezclar sintéticos).",
    ]})
    with pd.ExcelWriter(dst) as xw:
        med.to_excel(xw, sheet_name="mediciones", index=False)
        ope.to_excel(xw, sheet_name="operacion", index=False)
        mapa.to_excel(xw, sheet_name="corridas", index=False)
        if "solubilidad" in xl.sheet_names:
            xl.parse("solubilidad").to_excel(xw, sheet_name="solubilidad", index=False)
        notes.to_excel(xw, sheet_name="LEEME", index=False)
    return {"mediciones": len(med), "operacion": len(ope), "corridas": len(mapa)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("entrada"); ap.add_argument("salida")
    ap.add_argument("--crystallizer", required=True, choices=["A3110", "A3210"])
    ap.add_argument("--origin", required=True, choices=["measured", "synthetic"])
    ap.add_argument("--volume-m3", type=float, required=True)
    ap.add_argument("--feed-kno3", type=float, required=True, help="kg/m3")
    ap.add_argument("--feed-clo4", type=float, help="kg/m3")
    ap.add_argument("--rho-liquid", type=float, help="kg/m3")
    ap.add_argument("--impurity-basis-confirmed", action="store_true")
    ap.add_argument("--experiment", default="estable")
    ap.add_argument("--shuffle-seed", type=int)
    a = ap.parse_args(argv)
    info = convert_file(a.entrada, a.salida, crystallizer=a.crystallizer, origin=a.origin, volume_m3=a.volume_m3,
                        feed_kno3_kg_m3=a.feed_kno3, feed_clo4_kg_m3=a.feed_clo4, rho_liquid_kg_m3=a.rho_liquid,
                        impurity_basis_confirmed=a.impurity_basis_confirmed, experiment=a.experiment,
                        shuffle_seed=a.shuffle_seed)
    print(info)


if __name__ == "__main__":
    main()
