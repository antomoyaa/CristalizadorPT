"""Parámetros base del MSMPR de KNO3 con su procedencia.

Regla del proyecto (ver AGENTS.md): ningún valor se inventa en silencio. Cada clave de BASE
tiene una etiqueta en PROVENANCE:

  'respaldado'  : valor con fuente (Perry, IAPWS, Rolfs, Selzer...). Ver parametros/parametros_DPBE.csv
  'interno'     : dato de SQM/GEA o de tus balances, sin confirmar con la planta
  'placeholder' : valor ilustrativo; hay que calibrarlo o ajustarlo con datos

Unidades: SI. Concentraciones en kg/m3 de pulpa (~licor), temperaturas en degC en la interfaz.
"""
from __future__ import annotations

import math
import warnings
from copy import deepcopy

M_KNO3 = 0.1011032          # kg/mol
M_H2O = 0.0180153           # kg/mol
KCAL = 4184.0               # J/kcal
R_GAS = 8.314462618         # J/(mol K)

BASE = {
    # --- Geometría y flujos (interno: hojas GEA / balance de masa; actualizar con datos de planta)
    "volume_m3": 660.0,
    "feed_flow_m3_s": 280.83 / 3600.0,       # caudal total de alimentación (stock + agua de dilución)
    "feed_dilution_fraction": 0.05,          # d = Q_dil / Q_total
    "crystallizer_temperature_C": 40.0,
    "feed_temperature_C": 55.0,
    "dilution_water_temperature_C": 20.0,    # agua de proceso (P&ID: process water 20 degC)
    "tank_volume_m3": 150.0,                 # tanque de licor madre (B3520); placeholder
    # --- Alimentación
    "feed_stock_kg_m3": {
        "KNO3": 220.0,       # placeholder (ver nota: SSK-d da ~525 kg/m3 a 70 degC)
        "ClO4": 5.56379,     # interno: informe §2, base iónica SSK-d (antes 3.8e-3: base distinta)
        "Mg": 32.1115,       # interno: Evaporative_MSMPR_DPBE-2.ipynb feed_composition['Mg_g_L'] = 32.1114745 g/L
        "OTHER": 386.5,      # interno: demás sales disueltas (1-0.3994-0.3460)*1518 kg/m3, SSK-d; trazador inerte
    },
    # --- Termodinámica
    "reference_temperature_C": 40.0,
    "Cstar_ref_kg_m3": 180.0,                # placeholder: solubilidad en salmuera NO medida
    "solubility_slope_K_inv": 0.0312,        # Perry Tabla 2-122 (pendiente de q, agua pura, 30-50 degC)
    "deltaH_crystallization_J_kg": 7.73 * KCAL / M_KNO3,   # interno: -7.73 kcal/mol (319.9 kJ/kg)
    "latent_heat_J_kg": 2.406081e6,          # respaldado: IAPWS SR1, agua pura a 40 degC
    "solution_cp_J_kg_K": 2178.0,            # interno/estimado (alternativa: 2248, corriente 334)
    "rho_liquid_kg_m3": 1470.0,              # interno: PFD GEA (densidad nominal del licor)
    "density_k_KNO3": 0.57,                  # respaldado: Perry Tabla 2-79, d(rho)/dC a 40 C entre 16 y 24 % en masa
    "density_k_other": 0.57,                 # supuesto: igual al del KNO3 para las demás sales disueltas
    "rho_water_kg_m3": 998.2,                # respaldado: IAPWS SR6 a 20 degC
    "cp_water_J_kg_K": 4184.19,              # respaldado: IAPWS SR6 a 20 degC
    "water_mass_fraction_liquor": 0.3994,    # interno: corriente 334 SSK-d
    # --- Cristal
    "rho_crystal_kg_m3": 2110.0,             # respaldado: Perry Tabla 2-1
    "kv": math.pi / 6.0,                     # exacto por definición (diámetro de esfera equivalente)
    # --- Cinética (placeholders: calibrar con datos; Rolfs/Selzer no son directamente transferibles)
    "kG_m_s": 2.0e-7,
    "g": 1.0,
    "Eg_J_mol": 0.0,                         # energía de activación del crecimiento (sin dato)
    "kB_m3_s_inv": 1.0e10,
    "b": 2.0,
    # --- Impureza ClO4 (placeholders; el informe dice: no hay Kd validado)
    "Kd": 0.002,                             # w_cristal / w_licor
    "cake_moisture_fraction": 0.03,
    "centrifuge_wash_eff": 0.75,
    # --- Recirculación
    "purge_fraction_alpha": 0.25,            # fracción del licor madre que se purga
    # --- Calidad (interno, provisional)
    "spec_ClO4_ppm_max": 40.0,
}

PROVENANCE = {k: "placeholder" for k in BASE}
PROVENANCE.update({
    "volume_m3": "interno", "feed_flow_m3_s": "interno", "feed_temperature_C": "interno",
    "crystallizer_temperature_C": "interno", "tank_volume_m3": "placeholder",
    "dilution_water_temperature_C": "interno",
    "feed_stock_kg_m3": "placeholder (KNO3) / interno (ClO4, Mg, OTHER)",
    "solubility_slope_K_inv": "respaldado", "deltaH_crystallization_J_kg": "interno",
    "latent_heat_J_kg": "respaldado", "solution_cp_J_kg_K": "interno",
    "rho_liquid_kg_m3": "interno", "rho_water_kg_m3": "respaldado", "cp_water_J_kg_K": "respaldado",
    "water_mass_fraction_liquor": "interno", "rho_crystal_kg_m3": "respaldado", "kv": "respaldado",
    "density_k_KNO3": "respaldado", "density_k_other": "placeholder",
    "spec_ClO4_ppm_max": "interno",
})


# Coeficientes EMPÍRICOS del tren integrado (plant.simulate_plant_npt3_calibrated y afines). No son mecanísticos ni están
# calibrados con datos de planta: etiqueta 'placeholder' (AGENTS.md regla 1). Entran en placeholders().
EMPIRICAL_COEFFS = {
    "vac_P1_P2_nom": "interno: Entrega 2, DWSIM (11.59 / 5.43 kPa abs)",
    "vac_K1_K2": "respaldado: Antoine (agua) 133 / 255 K/bar; no incluye elevación de punto de ebullición de la salmuera",
    "vac_limits": "placeholder: rango de la válvula de vacío",
    "K_TRIM_AIT": "placeholder: 5.2e-5 1/ppm (compensa ~40 % del exceso; decisión económica)",
    "coprecip_threshold_ClO4_g_L": "placeholder: 4.0-4.5 g/L ClO4 (= 5.6-6.3 g/L KClO4; ensayo 6.54 g/L a 40 C, x1.393)",
    "coprecip_slope_ppm_per_g_L": "placeholder: 82-125 ppm por g/L",
    "moisture_d50_exponent": "placeholder: 0.55-0.65", "rem_eff_d50_exponent": "placeholder: 0.65-0.85",
    "rem_eff_6000_base": "interno: 0.899 (Compilado pruebas A6000)",
    "xS_d50_gains": "placeholder: ganancias de xS, d50, T",
    "level_loop": "placeholder: tauc 6 min, xS_level_gain 1.0, L_sp_lim 0.9-1.1, Kc_dens 3.0; tau1/tau2 0.7/1.4 h son del protocolo (Entrega 2)",
}
# Claves propias del perfil 'planta_sqm' cuya fuente NO es planta (ver BASE_PLANTA_SQM): entran en placeholders().
PROVENANCE.update({
    "planta_sqm.Cstar_ref_kg_m3": "placeholder: dato sintetico (datos_cristalizador.xlsx, hoja solubilidad); sin medicion de SQM",
    "planta_sqm.solubility_slope_K_inv": "placeholder: pendiente de la base sintetica del notebook de sensibilidad",
})
PROVENANCE.update({f"empirico.{k}": ("placeholder" if v.startswith("placeholder") else v.split(":")[0]) for k, v in EMPIRICAL_COEFFS.items()})

# ---------------------------------------------------------------------------------------------
# Perfil 'calibracion': la base EXACTA del notebook MSMPR_KNO3_DPBE_sensibilidad_Pareto.ipynb, que es el modelo
# que usa calibracion_do_crystal.ipynb. Los kG, kB y Kd que se calibren allí solo valen sobre ESTA base
# (solubilidad, densidades, factor de forma): al cambiarla hay que recalibrar.
# ---------------------------------------------------------------------------------------------
BASE_CALIBRACION = {
    "volume_m3": 660.0,
    "feed_dilution_fraction": 0.05,
    "crystallizer_temperature_C": 40.0,
    "residence_time": 660.0 / (280.83 / 3600.0),          # s; en este notebook Q se deriva de V/tau
    "feed_stock_kg_m3": {"KNO3": 220.0, "ClO4": 3.799262673e-3},
    "reference_temperature_C": 40.0,
    "Cstar_ref_kg_m3": 180.0,
    "solubility_slope_K_inv": 0.010,
    "kG_m_s": 2.0e-7, "g": 1.0, "kB_m3_s_inv": 1.0e10, "b": 2.0,
    "rho_crystal_kg_m3": 2100.0, "kv": math.pi / 6.0,
    "Kd": 0.01, "rho_liquid_kg_m3": 1040.0,
}

# ---------------------------------------------------------------------------------------------
# Perfil 'planta_sqm': datos reconciliados de la corriente industrial 334 SSK-d (informe_parametros.md §2)
# y, PROVISORIAMENTE, solubilidad de datos_cristalizador.xlsx (hoja 'solubilidad'), que es un conjunto SINTETICO generado con la base del
# notebook de sensibilidad (C*_ref = 180, pendiente 0.010 + ruido), NO una medicion de SQM. Permite simular Plant con
# recuperación industrial realista (Y ~ 50-65 %, cumpliendo AGENTS.md: Y entre 0.3 y 0.7).
# ---------------------------------------------------------------------------------------------
BASE_PLANTA_SQM = {
    "feed_stock_kg_m3": {
        "KNO3": 525.228,      # interno: corriente 334 SSK-d (0.3460 * 1518 kg/m3 a 70 degC, informe_parametros.md §2)
        "ClO4": 5.56379,      # interno: base iónica SSK-d (g/L = kg/m3)
        "Mg": 32.1115,        # interno: Evaporative_MSMPR_DPBE-2.ipynb feed_composition['Mg_g_L'] = 32.1114745 g/L
        "OTHER": 386.5,       # interno: demás sales disueltas (1-0.3994-0.3460)*1518 kg/m3, SSK-d
    },
    "feed_temperature_C": 70.0,      # interno: corriente 334 SSK-d a 70 degC (informe_parametros.md §2)
    "Cstar_ref_kg_m3": 177.940756,   # PLACEHOLDER: dato SINTETICO (datos_cristalizador.xlsx, hoja 'solubilidad' a 40 degC); reemplazar por solubilidad real
    "solubility_slope_K_inv": 0.010, # PLACEHOLDER: pendiente de la base sintetica del notebook de sensibilidad; reemplazar por curva real
    "rho_liquid_kg_m3": 1518.0,      # interno: densidad líquida corriente 334 SSK-d (informe_parametros.md §2)
    "solution_cp_J_kg_K": 2248.0,    # interno: Cp líquido corriente 334 SSK-d (informe_parametros.md §2)
}

# claves que condicionan los parámetros calibrados (si difieren, kG/kB/Kd no son transferibles)
BASIS_KEYS = ["reference_temperature_C", "Cstar_ref_kg_m3", "solubility_slope_K_inv",
              "rho_crystal_kg_m3", "kv", "rho_liquid_kg_m3"]
CALIBRATED_KEYS = ["kG_m_s", "kB_m3_s_inv", "Kd", "g", "b"]


def _deep_update(p: dict, src: dict) -> None:
    for k, v in src.items():
        if isinstance(v, dict) and isinstance(p.get(k), dict):
            p[k].update(v)
        else:
            p[k] = deepcopy(v)


_AVISO_PLANTA_SQM = {"emitido": False}


def _avisar_planta_sqm():
    """Avisa una vez por sesion que C*_ref y la pendiente de solubilidad del perfil 'planta_sqm' son placeholders sinteticos."""
    if not _AVISO_PLANTA_SQM["emitido"]:
        _AVISO_PLANTA_SQM["emitido"] = True
        warnings.warn("perfil 'planta_sqm': Cstar_ref_kg_m3 y solubility_slope_K_inv vienen de datos SINTETICOS (placeholder); "
                      "los resultados de sobresaturacion/D50 son indicativos hasta tener la solubilidad real.", UserWarning, stacklevel=3)


def make_params(profile: str = "literatura", **overrides) -> dict:
    """Copia de BASE con sobrescrituras (admite 'feed_stock_kg_m3' parcial: {'KNO3': 525.0}).

    profile='literatura'  : valores respaldados con fuente donde existen (Perry, IAPWS...). Por defecto.
    profile='calibracion' : base del notebook de Pareto/calibración (ver BASE_CALIBRACION). Con ella, el
                            término de humedad de torta se anula porque su Kd ya absorbe todo el ClO4 del producto seco.
    profile='planta_sqm'  : alimentación industrial corriente 334 SSK-d (informe_parametros.md §2, KNO3 = 525.23 kg/m3,
                            Mg = 32.1115 kg/m3) con recuperación estacionaria Y ~ 64.3 % (cumple 0.3 <= Y <= 0.7).
    """
    if profile not in ("literatura", "calibracion", "planta_sqm"):
        raise ValueError("profile debe ser 'literatura', 'calibracion' o 'planta_sqm'")
    p = deepcopy(BASE)
    p["profile"] = profile
    if profile == "calibracion":
        src = {k: v for k, v in BASE_CALIBRACION.items() if k != "residence_time"}
        _deep_update(p, src)
        p["feed_flow_m3_s"] = BASE_CALIBRACION["volume_m3"] / BASE_CALIBRACION["residence_time"]
        p["cake_moisture_fraction"] = 0.0
    elif profile == "planta_sqm":
        _deep_update(p, BASE_PLANTA_SQM)
        _avisar_planta_sqm()
    for k, v in overrides.items():
        if k not in p:
            raise KeyError(f"Parámetro desconocido: {k}")
        if isinstance(p[k], dict) and isinstance(v, dict):
            p[k].update(v)
        else:
            p[k] = v
    return p


def diff_bases():
    """Tabla de diferencias entre el perfil 'literatura' y el de 'calibracion' (lo que hay que alinear)."""
    import pandas as pd
    a, b = make_params("literatura"), make_params("calibracion")
    rows = []
    for k in list(BASE_CALIBRACION) + ["cake_moisture_fraction"]:
        if k == "residence_time":
            continue
        va, vb = a[k], b[k]
        if isinstance(va, dict):
            for kk in BASE_CALIBRACION[k]:
                if va[kk] != vb[kk]:
                    rows.append({"parametro": f"{k}.{kk}", "literatura": va[kk], "calibracion": vb[kk],
                                 "condiciona_calibracion": False})
        elif va != vb:
            rows.append({"parametro": k, "literatura": va, "calibracion": vb,
                         "condiciona_calibracion": k in BASIS_KEYS})
    return pd.DataFrame(rows)


def load_base_from_notebook(path: str) -> dict:
    """Extrae BASE del notebook de Pareto (sin ejecutarlo) para comprobar que BASE_CALIBRACION no se desfasó."""
    import ast
    import json
    nb = json.load(open(path, encoding="utf-8"))
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        src = "".join(c["source"])
        for node in ast.parse(src).body:
            if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "BASE" for t in node.targets):
                return eval(compile(ast.Expression(node.value), "BASE", "eval"), {"math": math})
    raise ValueError("No se encontró BASE en el notebook")


def apply_calibration(p: dict, calibration, equipment: str = "", force: bool = False) -> dict:
    """Inserta kG, kB, Kd (y g, b si vienen) calibrados en `p`, comprobando que la base coincide."""
    if hasattr(calibration, "columns"):
        calibration = dict(zip(calibration["parameter"], calibration["optimized"].astype(float)))
    bad = [k for k in BASIS_KEYS if p[k] != BASE_CALIBRACION[k]]
    if bad and not force:
        raise ValueError(f"La base difiere de la de calibración en {bad}: kG/kB/Kd no son transferibles. "
                         "Usa make_params('calibracion') o recalibra.")
    q = deepcopy(p)
    for k, v in calibration.items():
        if k not in CALIBRATED_KEYS:
            raise KeyError(f"Solo se aceptan {CALIBRATED_KEYS}; llegó {k}")
        q[k] = float(v)
    q["cake_moisture_fraction"] = 0.0          # el Kd calibrado ya incluye el arrastre de licor
    q["calibrated_for"] = equipment or "sin especificar"
    return q


def placeholders(p: dict | None = None) -> list[str]:
    """Claves cuyo valor sigue siendo un supuesto: lo que hay que calibrar antes de sacar conclusiones."""
    return [k for k, v in PROVENANCE.items() if v.startswith("placeholder")]


def literature_table(csv_path: str):
    """Lee parametros/parametros_DPBE.csv (informe de parámetros) como DataFrame para consulta."""
    import pandas as pd
    return pd.read_csv(csv_path, encoding="utf-8-sig")

# ---------------------------------------------------------------------------------------------
# Datos conciliados de Magnesio (Mg2+) y Perclorato (ClO4-) desde las 5 planillas de planta:
# 1) Balance_de_masa_segun_PFD_completo.xlsx (hojas 'Propiedades' y 'Corrientes')
# 2) Balances de Masa detallado v6.xlsx (hojas 'Balance Químico' y 'Balance NPT3')
# 3) Base de datos.xlsx (histórico 2025-2026: SSK, MLK, ML Refi, Cristal Crudo y Refinado)
# 4) AQ SSKd.xlsx (campaña oct-nov 2026: SSK, SSKd y Flujos)
# 5) AG NPT_3 (3000 y 6000)_.xlsx (granulometría histórica D50 y D80 en Áreas 3000 y 6000)
# ---------------------------------------------------------------------------------------------
MG_CLO4_PLANT_DATA = {
    # --- Flujos másicos nominales PFD / v6 [t/h] ---
    "F_208_SSK_t_h": 356.85,
    "F_234_ML_ref_t_h": 21.33,
    "F_303_agua_dil_t_h": 27.14,
    "F_313_agua_rep_t_h": 11.97,
    "F_327_agua_lav_t_h": 4.61,
    "F_1_SSKd_t_h": 405.32,
    "F_6_A3110_out_t_h": 390.23,
    "F_11_A3210_out_t_h": 390.23,
    "F_17_queque_crudo_t_h": 76.48508,
    "F_prod_refinado_v6_t_h": 65.30817,
    "F_MLK_recirc_v6_t_h": 223.670822,
    "F_MLK_purga_v6_t_h": 109.196982,
    # --- Fracción másica total por corriente [ppm = mg/kg] (PFD / v6) ---
    "w_Mg_208_ppm": 10206.0,       # 0.010206 en SSK total (10858 ppm en fase líquida)
    "w_ClO4_208_ppm": 3892.0,      # 0.003892 en SSK total (4140 ppm en fase líquida)
    "w_Mg_234_ppm": 2397.0,        # 0.002397 en ML refinación (2445 ppm en fase líquida)
    "w_ClO4_234_ppm": 738.0,       # 0.000738 en ML refinación (753 ppm en fase líquida)
    "w_Mg_303_ppm": 329.0,         # 0.000329 en agua industrial (303, 313, 327)
    "w_ClO4_303_ppm": 36.0,        # 0.000036 en agua industrial (303, 313, 327)
    "w_Mg_1_SSKd_ppm": 9134.0,     # 0.009134 en corriente 1/334 SSK-d (9654 ppm en líquido)
    "w_ClO4_1_SSKd_ppm": 3468.0,   # 0.003468 en corriente 1/334 SSK-d (3665 ppm en líquido)
    # --- Fracción másica en el líquido por etapa [ppm] (PFD Corrientes) ---
    "wL_Mg_A3110_ppm": 10767.0,    # 0.010767 en licor salida A3110 (corriente 6)
    "wL_ClO4_A3110_ppm": 4088.0,   # 0.004088 en licor salida A3110 (corriente 6)
    "wL_Mg_A3210_ppm": 11720.0,    # 0.011720 en licor salida A3210 / S3540 (corriente 11/29/12)
    "wL_ClO4_A3210_ppm": 4446.0,   # 0.004446 en licor salida A3210 / S3540 (corriente 11/29/12)
    "wL_Mg_lav_S3550_ppm": 10748.0,# 0.010748 en licor retenido y filtrado S3550 (corrientes 15-24)
    "wL_ClO4_lav_S3550_ppm": 4070.0,# 0.004070 en licor retenido y filtrado S3550 (corrientes 15-24)
    "w_Mg_MLK_v6_ppm": 10578.0,    # 0.010578 en MLK purga/recirculación (v6)
    "w_ClO4_MLK_v6_ppm": 4012.0,   # 0.004012 en MLK purga/recirculación (v6)
    # --- Sales reconstituidas en fase líquida A3210 [fracción másica] (PFD) ---
    "x_KClO4_A3210": 0.006194,
    "x_MgSO4_A3210": 0.022750,
    "x_MgCl2_A3210": 0.027914,
    "h_licor_por_crudo_v6": 0.070012,  # w_H2O_crudo (0.0339) / w_H2O_MLK (0.4842)
    "moisture_post_wash_PFD": 0.04997, # 1.912127 / 38.24254 (5.0 % licor retenido de lavado en queque 17)
    # --- Calidad en sólido crudo (corriente 17) y refinado (Área 6000) [ppm] ---
    "ppm_Mg_crudo_PFD": 537.0,
    "ppm_Mg_crudo_v6": 740.0,
    "ppm_ClO4_crudo_PFD": 203.0,
    "ppm_ClO4_crudo_v6": 281.0,
    "ppm_Mg_refinado_v6": 144.0,
    "ppm_Mg_refinado_espec_v6": 72.0,
    "ppm_ClO4_refinado_v6": 44.0,
    "ppm_ClO4_refinado_espec_v6": 22.1486,
    # --- Estadísticas históricas de laboratorio IC/AA (Base de datos.xlsx y AQ SSKd.xlsx) ---
    "lab_SSK_Mg_g_L": {"mean": 32.1115, "median": 32.0, "std": 4.3406, "p25": 29.3, "p75": 35.0, "max": 50.0, "aq_2026_mean": 28.6820},
    "lab_SSK_ClO4_g_L": {"mean": 3.7993, "median": 3.90, "std": 0.4890, "p25": 3.50, "p75": 4.20, "max": 6.00, "aq_2026_mean": 4.4735},
    "lab_SSK_KClO4_g_L": {"mean": 5.2909, "median": 5.40, "std": 0.6805, "p25": 4.90, "p75": 5.80, "max": 8.40, "aq_2026_mean": 6.2304},
    "lab_MLK_Mg_g_L": {"mean": 31.0317, "median": 31.0, "std": 4.2458, "p25": 28.6, "p75": 34.0, "max": 47.0},
    "lab_MLK_ClO4_g_L": {"mean": 3.5862, "median": 3.60, "std": 0.4970, "p25": 3.30, "p75": 4.00, "max": 6.30},
    "lab_ML_refi_Mg_g_L": {"mean": 5.7698, "median": 5.60, "std": 2.4145, "p25": 4.10, "p75": 7.30, "max": 55.0},
    "lab_ML_refi_ClO4_g_L": {"mean": 1.5876, "median": 1.30, "std": 1.0116, "p25": 0.88, "p75": 2.10, "max": 6.60},
    "lab_crudo_Mg_ppm": {"mean": 991.9, "median": 1000.0, "std": 425.3, "p25": 690.0, "p75": 1300.0, "max": 2900.0},
    "lab_crudo_ClO4_ppm": {"mean": 314.7, "median": 202.0, "std": 291.8, "p25": 136.0, "p75": 432.75, "max": 3545.0},
    "lab_refinado_Mg_ppm": {"mean": 61.1, "median": 60.0, "std": 28.6, "p25": 40.0, "p75": 70.0, "max": 350.0},
    "lab_refinado_ClO4_ppm": {"mean": 30.86, "median": 21.0, "std": 19.19, "p25": 20.0, "p75": 36.0, "max": 263.0},
    # --- Granulometría histórica AG NPT_3 (3000 y 6000)_.xlsx [um] ---
    "granulo_AG_A3000_D50_um": {"mean": 380.4, "min": 295.7, "max": 470.0},
    "granulo_AG_A6000_D50_um": {"mean": 330.1, "min": 272.0, "max": 422.2},
}

PROVENANCE["MG_CLO4_PLANT_DATA"] = (
    "interno: conciliado desde Balance_de_masa_segun_PFD_completo.xlsx, "
    "Balances de Masa detallado v6.xlsx, Base de datos.xlsx, AQ SSKd.xlsx y AG NPT_3 (3000 y 6000)_.xlsx"
)
