"""Pruebas de la revisión: vacío-T, holgura del PID, trim de AIT 3105, balance de masa y procedencia de constantes."""
import ast
import inspect
import re
import warnings

import numpy as np
import pytest

from msmpr_core import params, plant
from msmpr_core.plant import VAC, simulate_plant_npt3_calibrated as sim_npt3

warnings.filterwarnings("ignore")


def psat_kpa(T):  # Antoine agua [mmHg, C]
    return 10 ** (8.07131 - 1730.63 / (233.426 + T)) * 0.133322


def test_vacio_nominal_coherente_con_dwsim():
    assert (VAC["P1_nom"] * 100 + 101.325) == pytest.approx(11.59, abs=0.05)   # A3110, 55 C
    assert (VAC["P2_nom"] * 100 + 101.325) == pytest.approx(5.43, abs=0.05)    # A3210, 40 C


def test_vacio_no_corresponde_a_temperatura_equivocada():
    # P absoluta nominal de agua pura debe quedar por encima de la Psat a (T - 12 C) y por debajo de la Psat a T
    for P, T in ((VAC["P1_nom"], 55.0), (VAC["P2_nom"], 40.0)):
        p_abs = P * 100 + 101.325
        assert psat_kpa(T - 12) < p_abs < psat_kpa(T)


@pytest.mark.parametrize("stage,T", [(1, 55.0), (2, 40.0)])
def test_ganancia_dT_dP_vs_antoine_con_bpe(stage, T):
    # la salmuera hierve BPE grados por encima del agua a la misma P: la pendiente de Antoine se evalua en T - BPE
    K = VAC["K%d" % stage]
    Tw = T - plant._BPE[stage]
    dP = (psat_kpa(Tw + 1) - psat_kpa(Tw - 1)) / 2 / 100.0
    assert K == pytest.approx(1.0 / dP, rel=0.03)
    assert K > {1: 133.0, 2: 255.0}[stage]          # mayor que la estimacion sin BPE


@pytest.mark.parametrize("n", ["P1", "P2"])
def test_pid_vacio_tiene_holgura_en_ambos_sentidos(n):
    lo, hi = VAC[f"{n}_lim"]
    P = VAC[f"{n}_nom"]
    assert hi - lo >= 0.03
    assert P - lo >= 0.01 and hi - P >= 0.01
    # holgura equivalente en temperatura: al menos +-2 C
    K = VAC["K" + n[1]]
    assert K * (hi - P) >= 2.0 and K * (P - lo) >= 2.0


def test_ganancia_de_lazo_conservada():
    assert VAC["K1"] * VAC["Kc_scale1"] == pytest.approx(375.0)
    assert VAC["K2"] * VAC["Kc_scale2"] == pytest.approx(375.0)


def test_sin_constantes_antiguas_en_plant():
    src = inspect.getsource(plant)
    for viejo in ("-0.932", "-0.962", "375.0 * (P_vac", "5.2e-4"):
        assert viejo not in src.replace('"Kc_scale1": 375.0 / 133.0, "Kc_scale2": 375.0 / 255.0', ""), viejo


def test_trim_ait_no_satura_la_valvula_en_regimen():
    # Antes (K=5.2e-4) la consigna quedaba en el tope ~73 % del tiempo. Con el escenario combinado (+8 % de flujo y +650 ppm)
    # un tope transitorio es esperable; el estado final y la fraccion de tiempo saturado deben ser bajos.
    d = sim_npt3(use_AIT3105=True, use_AIC3255=True)
    sat = (d.F_303_sp_t_h >= 42.0 - 0.01).mean()
    assert sat < 0.10
    assert d.F_303_sp_t_h.iloc[-1] < 42.0 - 1.0


def test_trim_ait_dimensionado_para_pulso():
    # +650 ppm en SSK -> exceso en la alimentacion previa; la consigna debe quedar dentro de [12, 42] sin recorte
    dw_pre = 650.0 * 356.85 / plant.F_PRE_NOM
    F_sp = plant.F_PRE_NOM * (plant.R_0_NOM + plant.K_TRIM_AIT * dw_pre)
    assert 12.0 < F_sp < 42.0


def test_razon_dilucion_documentada():
    assert plant.R_0_NOM == pytest.approx(0.07176, abs=1e-4)
    p = params.make_params("literatura")
    d_total = p["feed_dilution_fraction"]
    # razon dilucion/alimentacion equivalente a la fraccion sobre el total: d/(1-d); debe ser distinta de R_0 (otra base)
    assert d_total / (1 - d_total) != pytest.approx(plant.R_0_NOM, rel=0.05)


def test_balance_masa_nominal_cierra():
    d = sim_npt3(use_AIT3105=False, use_AIC3255=False, duration_min=5.0)
    assert abs(d["dMass_A3110_t_h"].iloc[0]) < 0.05           # F1 - vapor - perdidas - F6 = 0 en el nominal


@pytest.fixture(scope="module")
def sims():
    return {m: sim_npt3(tuning_mode=m) for m in ("TL_damped", "ZN_oscillatory", "Hopf_limit_cycle")}


@pytest.fixture(scope="module")
def sim_cascada():
    return sim_npt3(density_cascade=True)


def test_nivel_inicia_en_reposo_y_cierra_balance_nominal():
    d = sim_npt3(use_AIT3105=False, use_AIC3255=False, duration_min=3.0)
    for c in ("dMass_A3110_t_h", "dMass_A3210_t_h"):
        assert abs(d[c].iloc[0]) < 0.05
    assert d.L1.iloc[0] == pytest.approx(1.0) and d.L2.iloc[0] == pytest.approx(1.0)


def test_inventario_coherente_con_los_flujos():
    # M(t) debe ser la integral de dMass: el nivel es un balance de masa real, no un valor impuesto
    d = sim_npt3()
    dt_h = (d.t_min.iloc[1] - d.t_min.iloc[0]) / 60.0
    for M, dm in (("M1_t", "dMass_A3110_t_h"), ("M2_t", "dMass_A3210_t_h")):
        est = d[M].iloc[0] + np.cumsum(d[dm].values) * dt_h
        assert np.max(np.abs(est - d[M].values)) < 0.02 * d[M].iloc[0]


@pytest.mark.parametrize("modo", ["TL_damped", "ZN_oscillatory", "Hopf_limit_cycle"])
def test_nivel_acotado_ante_perturbaciones(sims, modo):
    d = sims[modo]
    assert d.L1.between(0.85, 1.15).all() and d.L2.between(0.85, 1.15).all()
    assert d.F_6_t_h.between(250.0, 520.0).all() and d.F_11_t_h.between(250.0, 520.0).all()
    assert abs(d.dMass_A3110_t_h.iloc[-1]) < 5.0      # sin deriva sostenida al final (t/h)


def test_perturbacion_de_flujo_mueve_la_descarga_y_el_nivel():
    # +8 % de F_208 (t = 25-140 min) debe aumentar F_6 y el nivel; sin lazo de nivel no habia reaccion de F_6 por balance
    d = sim_npt3()
    pico = d[(d.t_min > 25) & (d.t_min < 140)]
    assert pico.F_6_t_h.max() > 390.23 + 5.0
    assert pico.L1.max() > 1.001                 # el nivel se desvia poco (lazo rapido), pero se mueve


def test_nivel_consigna_fija_en_60_por_ciento(sims):
    d = sims["TL_damped"]
    assert (d.L1_sp == 1.0).all() and (d.L2_sp == 1.0).all()
    assert plant.LEVEL["level_sp_pct"] == 60.0
    assert d.level1_pct.iloc[0] == pytest.approx(60.0) and d.level2_pct.iloc[0] == pytest.approx(60.0)
    for c in ("level1_pct", "level2_pct"):
        assert d[c].between(60.0 * 0.95, 60.0 * 1.05).all()      # +-5 % del SP ante las perturbaciones


def test_cascada_de_densidad_opcional_mueve_la_consigna(sim_cascada):
    d = sim_cascada
    assert d.L1_sp.min() >= 0.90 - 1e-9 and d.L1_sp.max() <= 1.10 + 1e-9
    assert d.L1_sp.max() - d.L1_sp.min() > 1e-3


def test_ppm_con_nivel_sigue_cumpliendo_criterios(sims):
    assert sims["TL_damped"].ppm_ClO4_final_6000.max() < 40.0
    sin = sim_npt3(use_AIT3105=False, use_AIC3255=False)
    assert sin.ppm_ClO4_final_6000.max() > 40.0


def test_nivel_placeholders_etiquetados():
    for k in ("tauc_min", "xS_level_gain", "L_sp_lim", "Kc_dens"):
        assert k in plant.LEVEL


def test_coeficientes_empiricos_etiquetados():
    ph = params.placeholders()
    assert any(k.startswith("empirico.") for k in ph)
    for k in params.EMPIRICAL_COEFFS:
        assert f"empirico.{k}" in params.PROVENANCE


def test_referencia_temperatura_arrhenius_usa_parametro():
    from msmpr_core import dpbe
    p = params.make_params("literatura", Eg_J_mol=40000.0)
    p2 = dict(p, reference_temperature_C=50.0, Cstar_ref_kg_m3=p["Cstar_ref_kg_m3"])
    # a T = T_ref el factor de Arrhenius vale 1 -> G = kG sigma^g
    C = dpbe.cstar(50.0, p2) * 1.05
    _, G, _ = dpbe.kinetics(C, 50.0, p2)
    assert G == pytest.approx(p2["kG_m_s"] * 0.05 ** p2["g"], rel=1e-9)


# ---------------------------------------------------------------- vinculo vacio -> T -> granulometria
def test_vac_T_nominal_y_monotonia():
    assert plant.vac_T(VAC["P1_nom"], 1) == pytest.approx(55.0, abs=1e-6)
    assert plant.vac_T(VAC["P2_nom"], 2) == pytest.approx(40.0, abs=1e-6)
    P = np.linspace(-0.97, -0.85, 25)
    for st in (1, 2):
        T = np.array([plant.vac_T(p, st) for p in P])
        assert (np.diff(T) > 0).all()                 # mas vacio (P menor) -> T menor


@pytest.mark.parametrize("stage,T", [(1, 55.0), (2, 40.0), (1, 50.0), (2, 36.0)])
def test_vac_P_es_inversa_de_vac_T(stage, T):
    assert plant.vac_T(plant.vac_P(T, stage), stage) == pytest.approx(T, abs=1e-6)


def test_granulometria_nominal_igual_a_d50_de_planta():
    g = plant.granulometria_estacionaria()
    assert g["d50_1_um"] == pytest.approx(plant.GRANULO["D50_nom1_um"], rel=1e-9)
    assert g["d50_2_um"] == pytest.approx(plant.GRANULO["D50_nom2_um"], rel=1e-9)
    assert g["T1_C"] == pytest.approx(55.0, abs=1e-6) and g["T2_C"] == pytest.approx(40.0, abs=1e-6)


def test_mas_vacio_baja_T_y_cambia_d50_en_el_sentido_del_modelo():
    # Regresion del vinculo (cinetica placeholder): mas vacio -> T menor -> mayor sobresaturacion absoluta de C -> d50 mas grande
    nom = plant.granulometria_estacionaria()
    mas = plant.granulometria_estacionaria(VAC["P1_nom"] - 0.01, VAC["P2_nom"] - 0.01)
    menos = plant.granulometria_estacionaria(VAC["P1_nom"] + 0.01, VAC["P2_nom"] + 0.01)
    assert mas["T1_C"] < nom["T1_C"] < menos["T1_C"] and mas["T2_C"] < nom["T2_C"] < menos["T2_C"]
    assert mas["d50_2_um"] > nom["d50_2_um"] > menos["d50_2_um"]
    assert mas["d50_1_um"] > nom["d50_1_um"] > menos["d50_1_um"]


def test_sensibilidad_granulometrica_acotada():
    # con cinetica placeholder la respuesta debe ser finita y moderada (+-0.02 bar < 10 % de d50)
    for dP in (-0.02, 0.02):
        g = plant.granulometria_estacionaria(VAC["P1_nom"] + dP, VAC["P2_nom"] + dP)
        assert abs(g["d50_1_um"] / plant.GRANULO["D50_nom1_um"] - 1) < 0.10
        assert abs(g["d50_2_um"] / plant.GRANULO["D50_nom2_um"] - 1) < 0.10


def test_mas_agua_de_dilucion_reduce_sobresaturacion():
    base = plant.granulometria_estacionaria()
    dil = plant.granulometria_estacionaria(F_303=42.0)
    assert dil["sigma1"] < base["sigma1"]


def test_dinamico_arranca_en_d50_de_planta_y_columnas_nuevas():
    d = sim_npt3()
    assert d.d50_1_um.iloc[0] == pytest.approx(plant.GRANULO["D50_nom1_um"])
    assert d.d50_2_um.iloc[0] == pytest.approx(plant.GRANULO["D50_nom2_um"])
    for c in ("sigma1", "sigma2", "level1_pct", "level2_pct"):
        assert c in d.columns
    assert d.d50_2_um.between(0.7 * 514, 1.15 * 514).all() and d.d50_1_um.between(0.7 * 541, 1.15 * 541).all()


def test_dinamico_T_sigue_a_la_presion_via_vac_T():
    # en regimen: T de la simulacion ~ vac_T(P) + acoples (0.45*(Tin-70), 0.40*(T1-55)); tolerancia de acoples/dinamica
    d = sim_npt3(duration_min=300.0)
    ult = d.iloc[-1]
    assert abs(ult.T1_C - (plant.vac_T(ult.P_vac1_bar, 1) + 0.45 * 3.5)) < 1.0


# ---------------------------------------------------------------- identificacion de Kcu y modos de sintonia
@pytest.fixture(scope="module")
def modos():
    return {m: sim_npt3(tuning_mode=m) for m in ("ZN_oscillatory", "Hopf_limit_cycle", "TL_damped")}


def test_kcu_identificado_reproduce_las_constantes():
    Kcu, Pu = plant.identificar_Kcu()
    assert Kcu == pytest.approx(plant.KCU_SIM["Kcu_bar_C"], rel=0.02)
    assert Pu == pytest.approx(plant.KCU_SIM["Pu_min"], rel=0.05)


def test_por_debajo_y_por_encima_de_kcu():
    Kcu = plant.KCU_SIM["Kcu_bar_C"]
    est = plant._sim_P(0.9 * Kcu)
    ines = plant._sim_P(1.1 * Kcu)
    assert plant.amplitud_oscilacion(est, 300, 400) < 0.05           # decae
    assert plant.amplitud_oscilacion(ines, 300, 400) > 0.5           # oscilacion sostenida


def test_zn_decae_con_razon_cercana_a_un_cuarto(modos):
    dr = plant.razon_decaimiento(modos["ZN_oscillatory"])
    assert 0.15 <= dr <= 0.5
    assert plant.amplitud_oscilacion(modos["ZN_oscillatory"], 200, 300) < 0.05


def test_hopf_no_decae(modos):
    assert plant.razon_decaimiento(modos["Hopf_limit_cycle"]) > 1.0
    assert plant.amplitud_oscilacion(modos["Hopf_limit_cycle"], 200, 300) > 0.5


def test_tl_amortiguado_sin_oscilacion(modos):
    assert plant.amplitud_oscilacion(modos["TL_damped"], 200, 300) < 0.05


def test_los_tres_modos_cumplen_el_orden_de_agresividad(modos):
    a = {m: plant.amplitud_oscilacion(d, 100, 140) for m, d in modos.items()}
    assert a["Hopf_limit_cycle"] > a["ZN_oscillatory"] and a["Hopf_limit_cycle"] > a["TL_damped"]


def test_override_de_sintonia_cambia_el_resultado():
    base = sim_npt3(duration_min=120.0)
    sub = sim_npt3(duration_min=120.0, Kc_T_override=0.0005)
    assert not np.allclose(base.T2_C.values, sub.T2_C.values)


def test_razon_decaimiento_sin_picos_es_nan(modos):
    assert np.isnan(plant.razon_decaimiento(modos["TL_damped"]))


def test_filtro_pasa_alto_sin_efecto_de_borde():
    x = np.full(6000, 40.0)
    assert np.abs(plant._hp(x, 0.05)).max() < 1e-9


# ---------------------------------------------------------------- procedencia de los datos del perfil planta_sqm
def test_planta_sqm_solubilidad_marcada_como_placeholder():
    for k in ("planta_sqm.Cstar_ref_kg_m3", "planta_sqm.solubility_slope_K_inv"):
        assert params.PROVENANCE[k].startswith("placeholder")
        assert k in params.placeholders()
    src = inspect.getsource(params)
    assert "interno: datos_cristalizador.xlsx" not in src      # ya no se presenta el dato sintetico como interno


def test_planta_sqm_avisa_una_vez_que_la_solubilidad_es_sintetica():
    params._AVISO_PLANTA_SQM["emitido"] = False
    with pytest.warns(UserWarning, match="SINTETICOS"):
        params.make_params("planta_sqm")
    with warnings.catch_warnings(record=True) as w:           # segunda llamada: silencio
        warnings.simplefilter("always")
        params.make_params("planta_sqm")
    assert not [x for x in w if "SINTETICOS" in str(x.message)]


def test_otros_perfiles_no_avisan():
    params._AVISO_PLANTA_SQM["emitido"] = False
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        params.make_params("literatura"); params.make_params("calibracion")
    assert not [x for x in w if "SINTETICOS" in str(x.message)]


# ---------------------------------------------------------------- integridad de plant.py (Mg, recirculacion, sin duplicados)
def test_plant_sin_funciones_duplicadas():
    arbol = ast.parse(inspect.getsource(plant))
    nombres = [n.name for n in arbol.body if isinstance(n, ast.FunctionDef)]
    dup = sorted({n for n in nombres if nombres.count(n) > 1})
    assert not dup, f"definiciones repetidas (la ultima tapa a la anterior): {dup}"


def test_pfd_mg_clo4_tiene_todas_las_claves_usadas():
    usadas = set(re.findall(r'PFD_MG_CLO4\["([A-Za-z0-9_]+)"\]', inspect.getsource(plant)))
    faltan = usadas - set(plant.PFD_MG_CLO4)
    assert not faltan, f"claves usadas y ausentes en PFD_MG_CLO4: {sorted(faltan)}"


@pytest.mark.parametrize("modo", ["sin_reciclo_actual", "reciclo_fijo", "reciclo_inteligente"])
def test_recirculacion_y_magnesio_corre_en_los_tres_modos(modo):
    d = plant.simular_recirculacion_y_magnesio(modo, use_AIT3105=True, duration_min=120.0)
    assert len(d) > 100 and np.isfinite(d.ppm_ClO4_final_6000).all()


def test_tren_integrado_incluye_magnesio_y_sintonias_identificadas():
    d = sim_npt3(duration_min=120.0)
    for c in ("w_Mg_pre_ppm", "w_Mg_1_ppm", "C_Mg_B3520_g_L", "ppm_Mg_crudo_17", "ppm_Mg_final_6000", "d50_1_um", "level1_pct"):
        assert c in d.columns
    assert np.isfinite(d.ppm_Mg_final_6000).all()
    assert "Kc_T_override" in inspect.signature(sim_npt3).parameters


# ---------------------------------------------------------------- agua minima, estrategias y temperatura (B6-B8)
@pytest.fixture(scope="module")
def estrategias():
    return plant.comparar_estrategias_excursion()


def test_saturacion_kclo4_vs_temperatura():
    assert plant.sat_factor_T(40.0) == pytest.approx(1.0)
    assert plant.sat_factor_T(25.0) == pytest.approx(3.48 / 6.54, rel=1e-6)       # dato de SQM bajo 40 C
    Ts = np.linspace(25, 70, 40)
    assert (np.diff([plant.sat_factor_T(t) for t in Ts]) > 0).all()               # mas T -> mas solubilidad
    assert plant.c_sat_clo4_gL(40.0) == pytest.approx(plant.C_SAT_NOM_ClO4)
    assert plant.c_sat_clo4_gL(40.0, d_Mg=0.35) < plant.c_sat_clo4_gL(40.0)       # el Mg2+ saca el KClO4 de solucion


def test_tope_ambiental_de_agua_se_respeta(estrategias):
    sims, tabla = estrategias
    assert plant.F303_MAX_AMBIENTAL <= 36.0
    assert (tabla.agua_max_t_h <= plant.F303_MAX_AMBIENTAL + 1e-6).all()
    d = plant.simular_recirculacion_y_magnesio("reciclo_inteligente", use_AIT3105=True, F303_max=80.0)   # un tope mayor se recorta
    assert d.F_303_t_h.max() <= plant.F303_MAX_AMBIENTAL + 1e-6
    assert sim_npt3(duration_min=120.0).F_303_actual_t_h.max() <= plant.F303_MAX_AMBIENTAL + 1e-6


def test_propuesta_cumple_la_norma_y_las_demas_no_bajo_la_excursion(estrategias):
    _, tabla = estrategias
    assert bool(tabla.loc["4. AIT + AIC + supervisor T2 (propuesta)", "cumple_norma"])
    assert tabla.loc["4. AIT + AIC + supervisor T2 (propuesta)", "ppm_max"] < plant.PPM_NORMA
    for n in ("1. Sin reciclo (actual)", "2. Reciclo fijo", "3. AIT 3105 + AIC 3255"):
        assert not bool(tabla.loc[n, "cumple_norma"]), n        # el veredicto sale de la simulacion, no de un texto fijo


def test_el_supervisor_de_T2_es_la_palanca_no_el_agua():
    sin = plant.simular_recirculacion_y_magnesio("reciclo_inteligente", use_AIT3105=True, F303_max=36.0)
    con = plant.simular_recirculacion_y_magnesio("reciclo_inteligente", use_AIT3105=True, trim_T2=True, F303_max=36.0)
    assert sin.ppm_ClO4_final_6000.max() > plant.PPM_NORMA > con.ppm_ClO4_final_6000.max()
    assert con.T2_C.max() > 40.5 and sin.T2_C.max() < 40.5


def test_optimo_de_agua_bajo_el_tope_y_con_margen():
    df, opt = plant.barrido_agua_optimo()
    assert opt is not None
    assert opt["agua_media_t_h"] <= plant.F303_MAX_AMBIENTAL and opt["agua_media_t_h"] < 36.8
    assert opt["ppm_max"] <= 0.9 * plant.PPM_NORMA + 1e-9
    fact = df[df.factible]
    assert opt["agua_media_t_h"] == pytest.approx(fact.agua_media_t_h.min())      # es el minimo de los factibles


def test_analisis_temperatura_resultados_cualitativos():
    t = plant.analisis_temperatura()
    base, caliente = t.loc["Alimentacion +0 C (sin correccion)"], t.loc["Alimentacion +10 C (sin correccion)"]
    assert caliente.vacio_A3110_saturado_pct > 50.0 > base.vacio_A3110_saturado_pct       # alimentacion caliente satura el vacio de A3110
    assert t.loc["+10 C con agua de dilucion a 35 C", "vacio_A3110_saturado_pct"] < 5.0        # agua fria devuelve el margen
    assert t.loc["+10 C con agua a 35 C + prealimentacion", "T1_final_C"] == pytest.approx(55.0, abs=0.05)
    c40, c42 = t.loc["Cristalizador 2 en 40 C"], t.loc["Cristalizador 2 en 42 C"]
    assert c42.ppm_max < c40.ppm_max and c42.KNO3_medio_t_h < c40.KNO3_medio_t_h          # caliente: menos ClO4-, menos KNO3
    assert c42.perdida_KNO3_por_T_final_t_h / 2.0 == pytest.approx(plant.KNO3_PERDIDA_T_t_h_por_C, rel=0.02)


def test_calor_para_enfriar_el_agua_de_dilucion():
    assert plant.calor_enfriar_agua_dilucion_kW() == pytest.approx(27.14 * 1000 / 3600 * 4.172 * 25.0, rel=1e-9)


def test_notebook_06_sin_optimo_ni_veredictos_escritos_a_mano():
    import json, os
    ruta = os.path.join(os.path.dirname(__file__), "..", "notebooks", "06_Control_por_etapas_dilucion.ipynb")
    src = "\n".join("".join(c["source"]) for c in json.load(open(ruta, encoding="utf-8"))["cells"] if c["cell_type"] == "code")
    for prohibido in ("R* = 0.0973", "25.3 ppm", "np.log10(9.66)", "F_303* = 36.8", "(< 40 ppm)", "100 * np.exp", "ÓPTIMO GLOBAL (PARETO)", "+ 5.8", "Óptimo Económico"):
        assert prohibido not in src, prohibido


def test_supervisor_de_T2_en_el_modelo_de_agua_de_dilucion():
    sin = plant.simular_planta_agua_dilucion_magnesio("AIT3105_and_Mg_trim", "smart_AIC3255", True)
    con = plant.simular_planta_agua_dilucion_magnesio("AIT3105_and_Mg_trim", "smart_AIC3255", True, trim_T2=True)
    assert sin.F_303.max() <= plant.F303_MAX_AMBIENTAL + 1e-6 and con.F_303.max() <= plant.F303_MAX_AMBIENTAL + 1e-6
    assert sin.ppm_ClO4_final_6000.max() > plant.PPM_NORMA > con.ppm_ClO4_final_6000.max()
    assert con.dT2_trim_C.max() > 0.5 and sin.dT2_trim_C.max() == 0.0
