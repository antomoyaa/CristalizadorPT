"""Cristalizador al vacío (A3110) con reciclo de licor madre y purga: modelo dinámico y estacionario.

UNA sola fuente de ecuaciones: `rhs` es el modelo dinámico y el estado estacionario se obtiene resolviendo rhs = 0
(así, en lazo abierto y en el estacionario, el modelo no se mueve; ver tests/test_core.py).

Simplificaciones (declaradas):
  * C y CP son concentraciones por m3 de pulpa; el volumen de cristales se desprecia en los caudales. OJO: en A3110 es ~13 % másico (~9 % en volumen) y en A3210 ~20 %
    (Entrega 2), no <3 %; supuesto optimista a revisar (afecta Q_out, tau y la recuperación).
  * Todo el caudal de descarga Q_out pasa por la separación: licor madre Q_ML = Q_out; purga alpha*Q_ML, reciclo (1-alpha)*Q_ML.
  * Tanque de licor madre de volumen constante, mezcla perfecta.
  * La humedad de la torta solo entra en la calidad del producto (ppm), no en los balances.
  * El calor sensible usa rho_liquid y cp de solución para todas las corrientes, salvo el agua de dilución.
  * d50 dinámico estimado como 1.224 d32 (MSMPR ideal); solo vale cerca del estacionario.

Estados: V, T, C, CP, mu0..mu3, Cm, CPm, Qd, S, Sm  (Qd = caudal real de agua de dilución, con dinámica de válvula;
S = demás sales disueltas, trazador inerte que solo se diluye y recircula; Sm = su valor en el tanque de licor madre).
Entradas manipulables: Q_out [m3/s], m_evap [kg/s], Q_dil_cmd [m3/s], alpha [-].
Perturbaciones: Q_stock [m3/s], T_feed [C], C_stock, CP_stock, S_stock [kg/m3], T_dil [C],
valve_gain [-] (opcional, 1 = nominal; <1 modela caída de presión en la línea de agua de dilución).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import math

import numpy as np
from scipy.optimize import root

from . import dpbe

NAMES = ["V", "T", "C", "CP", "mu0", "mu1", "mu2", "mu3", "Cm", "CPm", "Qd", "S", "Sm"]
IDX = {n: i for i, n in enumerate(NAMES)}

# --- Constantes del tren integrado NPT-3 (todas PLACEHOLDER salvo que se indique; ver params.EMPIRICAL_COEFFS) ---
# Vacío: presión manométrica [bar] que da la T de cada cristalizador. Fuente: Entrega 2, DWSIM (11.59 kPa a 55 C; 5.43 kPa a 40 C,
# salmuera con elevación del punto de ebullición). K = dT/dP [K/bar] de Antoine (agua) a 55 y 40 C: 133 y 255.
VAC = {
    "P1_nom": (11.59 - 101.325) / 100.0, "P2_nom": (5.43 - 101.325) / 100.0,   # -0.897 / -0.959 bar man
    "K1": 133.0, "K2": 255.0,
    "P1_lim": (-0.920, -0.870), "P2_lim": (-0.975, -0.940),                      # (min, max) supuestos de diseño
    "Kc_scale1": 375.0 / 133.0, "Kc_scale2": 375.0 / 255.0,                      # conserva la ganancia de lazo con que se sintonizó
}
# Ganancia ultima y periodo de los lazos de temperatura (TIR -> PIRCAL), identificados SOBRE ESTA SIMULACION con accion P pura y valvula
# sin saturar (identificar_Kcu). Kc en [bar/C] antes del reescalado por etapa. Reemplazan a los de la tabla FOPDT (Kp=375, tau_p=5.5, theta=1).
KCU_SIM = {"Kcu_bar_C": 0.012407, "Pu_min": 9.70}
K_TRIM_AIT = 5.2e-5          # [1/ppm] razón extra de dilución por ppm de ClO4 sobre 3714.1 (placeholder; decisión económica)
F303_NOM, F_PRE_NOM = 27.14, 378.18
# Tope AMBIENTAL del agua de dilucion F_303 [t/h]: ninguna estrategia puede pedir mas (valvula FV 3103 y consigna FIRC 3103).
# Decision de proyecto: el agua es un recurso escaso para la comunidad; se busca el MINIMO que cumple la norma, nunca por sobre este tope.
F303_MAX_AMBIENTAL = 36.0
# KNO3 cristalizado por grado de enfriamiento en A3210 [t/h/C] = 27.26 t/h / 15 C (Entrega 2, Tabla 8): cada C de T2 por sobre 40 C pierde esa produccion.
KNO3_PERDIDA_T_t_h_por_C = 27.26 / 15.0
R_0_NOM = F303_NOM / F_PRE_NOM          # 0.07176: razón dilución / alimentación previa (no es feed_dilution_fraction)
F_VAPOR_A3110, F_PERDIDAS_A3110 = 10.03, 5.06   # t/h, balance de masa detallado
# Lazo de nivel (LIRCSAHL 3101 / 3201): inventario M [t] = tau * F nominal; nivel normalizado L = M / M_nom (SP = 1).
# tau de 0.7 h (A3110) y 1.4 h (A3210): protocolo experimental, Entrega 2. F_25 = F_21 (Supuesto 10) -> A3210 no acumula por reposición.
LEVEL = {"tau1_h": 0.7, "tau2_h": 1.4, "F_nom": 390.23, "tauc_min": 6.0, "tauI_factor": 4.0,   # tauc, tauI: placeholder
         "xS_level_gain": 1.0,       # placeholder: d(xS)/xS = gain * d(L)/L (mas retencion -> mas cristal); sin dato de planta
         "L_sp_lim": (0.90, 1.10),   # placeholder: rango de la consigna de nivel que puede dar el lazo de densidad
         "Kc_dens": 3.0,             # placeholder: [L / fraccion de solidos] del maestro de densidad (cascada DIC -> LIC)
         "level_sp_pct": 60.0}       # SP de nivel de A3110 y A3210 [% de la escala LIRCSAHL] (dato del usuario); L = 1 <-> 60 %

# Granulometria (vinculo vacio -> T -> sobresaturacion -> d50). Anclas = D50 medidos en planta (Reporte seguimiento cristalizadores NPT3,
# agosto 2026: promedios del under de A3110 y A3210). La SENSIBILIDAD viene del modelo MSMPR (dpbe) con cinetica placeholder.
GRANULO = {"D50_nom1_um": 541.0, "D50_nom2_um": 514.0, "perfil": "planta_sqm",
           "C_in1_nom_kg_m3": 525.228 * 405.32 / (405.32 - 10.03),   # KNO3 en la alimentacion a A3110 ya concentrada por el flash (placeholder)
           "update_every_min": 1.0, "tau_d50_min": 18.0, "d50_floor_frac": 0.49}


def tsat_water_C(P_abs_kPa):
    """Temperatura de saturacion del agua [C] a P absoluta [kPa] (Antoine, 1-100 C)."""
    return 1730.63 / (8.07131 - math.log10(P_abs_kPa / 0.133322)) - 233.426


# Elevacion del punto de ebullicion de la salmuera [C]: se ajusta para que las presiones de DWSIM den las T de operacion.
_BPE = {1: 55.0 - tsat_water_C(11.59), 2: 40.0 - tsat_water_C(5.43)}


def vac_T(P_bar_man, stage):
    """Temperatura de equilibrio del cristalizador `stage` (1: A3110, 2: A3210) para una presion de vacio [bar man]: Tsat(P) + BPE."""
    return tsat_water_C(max(101.325 + 100.0 * P_bar_man, 0.5)) + _BPE[stage]


def vac_P(T_C, stage):
    """Inversa de vac_T: presion de vacio [bar man] que da la temperatura T_C."""
    Tw = T_C - _BPE[stage]
    return (10 ** (8.07131 - 1730.63 / (Tw + 233.426)) * 0.133322 - 101.325) / 100.0


# K = dT/dP [K/bar] en el nominal, derivada numerica de vac_T (incluye la elevacion del punto de ebullicion: el agua equivalente
# esta ~6 C mas fria que la salmuera, donde la curva P-T es mas empinada). Reemplaza los valores 133/255 calculados sin BPE.
_h = 1e-4
for _n, _k in ((1, "K1"), (2, "K2")):
    _P = VAC["P%d_nom" % _n]
    VAC[_k] = (vac_T(_P + _h, _n) - vac_T(_P - _h, _n)) / (2 * _h)
    VAC["Kc_scale%d" % _n] = 375.0 / VAC[_k]

_PERFILES = {}


def _d50_mech(C_in, tau_h, T_C):
    """(d50 [um], sigma, C [kg/m3]) del MSMPR estacionario a T_C. Devuelve d50 = 0 si no hay sobresaturacion."""
    from . import params as _params
    p = _PERFILES.setdefault(GRANULO["perfil"], _params.make_params(GRANULO["perfil"]))
    r = dpbe.stationary_msmpr(C_in, tau_h * 3600.0, T_C, p)
    return float(r["d50"]) * 1e6, float(r["sigma"]), float(r["C"])


def _C_in1(F_pre, F_303):
    return GRANULO["C_in1_nom_kg_m3"] * (F_pre / (F_pre + F_303)) / (F_PRE_NOM / (F_PRE_NOM + F303_NOM))


def _cadena_d50(T1, T2, tau1_h, tau2_h, F_pre, F_303):
    d1, s1, C1 = _d50_mech(_C_in1(F_pre, F_303), tau1_h, T1)
    d2, s2, _ = _d50_mech(C1, tau2_h, T2)           # la etapa 2 recibe el licor de la etapa 1 (sin semillas: simplificacion)
    return d1, d2, s1, s2


_REF_CACHE = {}


def _ref_nominal():
    if "v" not in _REF_CACHE:
        _REF_CACHE["v"] = _cadena_d50(55.0, 40.0, LEVEL["tau1_h"], LEVEL["tau2_h"], F_PRE_NOM, F303_NOM)
    return _REF_CACHE["v"]


def granulometria_estacionaria(P_vac1_bar=None, P_vac2_bar=None, F_pre=F_PRE_NOM, F_303=F303_NOM,
                               tau1_h=None, tau2_h=None):
    """Vinculo estacionario P de vacio -> T -> sobresaturacion -> D50 de A3110 y A3210.

    T2 incluye el acople con la etapa 1 (0.40 K/K) usado en el modelo dinamico. D50 = D50_nominal * (d50_mec / d50_mec_nominal).
    La magnitud de la sensibilidad depende de la cinetica placeholder (kG, kB, g, b); no es un resultado de planta
    hasta calibrar con los protocolos de nucleacion y solubilidad.
    """
    P1 = VAC["P1_nom"] if P_vac1_bar is None else P_vac1_bar
    P2 = VAC["P2_nom"] if P_vac2_bar is None else P_vac2_bar
    tau1_h = LEVEL["tau1_h"] if tau1_h is None else tau1_h
    tau2_h = LEVEL["tau2_h"] if tau2_h is None else tau2_h
    T1 = vac_T(P1, 1)
    T2 = vac_T(P2, 2) + 0.40 * (T1 - 55.0)
    d1, d2, s1, s2 = _cadena_d50(T1, T2, tau1_h, tau2_h, F_pre, F_303)
    r1, r2, _, _ = _ref_nominal()
    return {"T1_C": T1, "T2_C": T2, "sigma1": s1, "sigma2": s2,
            "d50_1_um": GRANULO["D50_nom1_um"] * d1 / r1, "d50_2_um": GRANULO["D50_nom2_um"] * d2 / r2}

TAU_VALVE_S = 8.0          # constante de tiempo de la válvula/medidor del lazo de dilución (placeholder)


@dataclass
class Plant:
    p: dict
    tau_valve: float = TAU_VALVE_S
    d_nom: dict = field(default_factory=dict)
    anchor: tuple | None = None        # (C, CP) a los que rho_liquid_kg_m3 es la densidad del licor; se fija en el primer estacionario

    def __post_init__(self):
        p = self.p
        Q = p["feed_flow_m3_s"]
        d = p["feed_dilution_fraction"]
        self.d_nom = {
            "Q_stock": (1 - d) * Q,
            "T_feed": p["feed_temperature_C"],
            "C_stock": p["feed_stock_kg_m3"]["KNO3"],
            "CP_stock": p["feed_stock_kg_m3"]["ClO4"],
            "Mg_stock": p["feed_stock_kg_m3"].get("Mg", 32.1115),
            "S_stock": p["feed_stock_kg_m3"]["OTHER"],
            "T_dil": p["dilution_water_temperature_C"],
        }
    # ------------------------------------------------------------------ ecuaciones
    def rhs(self, x, u, d):
        p = self.p
        V, T, C, CP, m0, m1, m2, m3, Cm, CPm, Qd, S, Sm = x
        sigma, G, B = dpbe.kinetics(C, T, p)
        rho_c, kv = p["rho_crystal_kg_m3"], p["kv"]
        Rc = 3.0 * rho_c * kv * G * m2                      # consumo de soluto por crecimiento [kg/m3/s]
        Qout, mev, a = u["Q_out"], u["m_evap"], u["alpha"]
        Qrec = (1.0 - a) * Qout
        dV = d["Q_stock"] + Qd + Qrec - Qout - mev / p["rho_water_kg_m3"]
        dC = (d["Q_stock"] * d["C_stock"] + Qrec * Cm - Qout * C - Rc * V - C * dV) / V
        w_s = p["Kd"] * CP / p["rho_liquid_kg_m3"]          # kg ClO4 / kg cristal
        dCP = (d["Q_stock"] * d["CP_stock"] + Qrec * CPm - Qout * CP - Rc * w_s * V - CP * dV) / V
        f = Qout / V + dV / V
        dm0 = B - m0 * f
        dm1 = G * m0 - m1 * f
        dm2 = 2.0 * G * m1 - m2 * f
        dm3 = 3.0 * G * m2 - m3 * f
        dCm = Qout * (C - Cm) / p["tank_volume_m3"]
        dCPm = Qout * (CP - CPm) / p["tank_volume_m3"]
        rcp = p["rho_liquid_kg_m3"] * p["solution_cp_J_kg_K"]
        dT = (d["Q_stock"] * rcp * (d["T_feed"] - T)
              + Qd * p["rho_water_kg_m3"] * p["cp_water_J_kg_K"] * (d["T_dil"] - T)
              + p["deltaH_crystallization_J_kg"] * Rc * V
              - mev * p["latent_heat_J_kg"]) / (V * rcp)
        dQd = (u["Q_dil_cmd"] * d.get("valve_gain", 1.0) - Qd) / self.tau_valve
        dS = (d["Q_stock"] * d["S_stock"] + Qrec * Sm - Qout * S - S * dV) / V
        dSm = Qout * (S - Sm) / p["tank_volume_m3"]
        return np.array([dV, dT, dC, dCP, dm0, dm1, dm2, dm3, dCm, dCPm, dQd, dS, dSm])

    def step(self, x, u, d, dt):
        """RK4 con entradas constantes en el paso; moments con piso positivo."""
        k1 = self.rhs(x, u, d)
        k2 = self.rhs(x + 0.5 * dt * k1, u, d)
        k3 = self.rhs(x + 0.5 * dt * k2, u, d)
        k4 = self.rhs(x + dt * k3, u, d)
        xn = x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        xn[4:8] = np.maximum(xn[4:8], 0.0)
        return xn

    # ------------------------------------------------------------------ estacionario
    def steady_state(self, alpha=None, T_set=None, V_set=None, d=None, guess=None, Q_dil=None):
        """Resuelve rhs=0 para [C, CP, m_evap, Q_out] con T y V en consigna. Devuelve (x, u, d).
        Q_dil [m3/s] fija el agua de dilución (por defecto d_fraction * feed_flow; con Q_stock fijo, más agua = más caudal total)."""
        p = self.p
        alpha = p["purge_fraction_alpha"] if alpha is None else alpha
        T = p["crystallizer_temperature_C"] if T_set is None else T_set
        V = p["volume_m3"] if V_set is None else V_set
        d = dict(self.d_nom if d is None else d)
        Qd = p["feed_dilution_fraction"] * p["feed_flow_m3_s"] if Q_dil is None else Q_dil
        Qtot = d["Q_stock"] + Qd

        def build(z):
            C, CP, mev, Qout = z
            _, G, B = dpbe.kinetics(C, T, p)
            mu = dpbe.moments_steady(B, G, V / Qout)
            S_ss = d["Q_stock"] * d["S_stock"] / (alpha * Qout)      # trazador inerte: se va solo por la purga
            x = np.array([V, T, C, CP, *mu, C, CP, Qd, S_ss, S_ss])
            u = {"Q_out": Qout, "m_evap": mev, "Q_dil_cmd": Qd, "alpha": alpha}
            return x, u

        scale = np.array([Qtot, Qtot * d["C_stock"] / V, Qtot * d["CP_stock"] / V, 1.0 / 3600.0])

        def resid(z):
            x, u = build(z)
            r = self.rhs(x, u, d)
            return np.array([r[0], r[2], r[3], r[1]]) / np.array([scale[0] / V, scale[1], scale[2], scale[3]])

        Cs = dpbe.cstar(T, p)
        z0 = guess if guess is not None else np.array(
            [Cs * 1.03, d["CP_stock"] * Qtot / max(alpha, 0.05) / Qtot, 2.0, Qtot * 0.98])
        sol = root(resid, z0, method="hybr", tol=1e-13)
        if not sol.success or np.max(np.abs(resid(sol.x))) > 1e-7:
            raise RuntimeError(f"Estacionario sin convergencia (alpha={alpha}): {sol.message}")
        x, u = build(sol.x)
        if self.anchor is None:
            self.anchor = (float(x[2]), float(x[11]))
        return x, u, d

    # ------------------------------------------------------------------ densidad (sensor DIRCAHL 3101)
    def liquor_density(self, C, S):
        """Densidad del licor [kg/m3]: nominal + efecto de la concentración de KNO3 y de las demás sales.

        S = demás sales disueltas (trazador inerte con su propio balance). El agua de dilución baja la densidad del licor
        directamente; una perturbación de ClO4 no la mueve (el ClO4 es traza). Solo valen las diferencias respecto del
        ancla (C0, S0); la densidad absoluta depende de supuestos."""
        p = self.p
        C0, S0 = self.anchor if self.anchor is not None else (C, S)
        return (p["rho_liquid_kg_m3"] + p["density_k_KNO3"] * (C - C0) + p["density_k_other"] * (S - S0))

    # ------------------------------------------------------------------ métricas
    def metrics(self, x, u, d):
        p = self.p
        V, T, C, CP = x[0], x[1], x[2], x[3]
        sigma, G, B = dpbe.kinetics(C, T, p)
        MT = dpbe.crystal_mass(x[7], p)
        Qout, a = u["Q_out"], u["alpha"]
        fresh = d["Q_stock"] * d["C_stock"]
        crystal_kg_s = Qout * MT
        Y_sol = crystal_kg_s / fresh
        Y_liq = 1.0 - a * Qout * C / fresh                  # licor que sale del sistema (purga)
        rho_L = self.liquor_density(C, x[11])
        w_cryst = p["Kd"] * CP / p["rho_liquid_kg_m3"]
        wm, eta = p["cake_moisture_fraction"], p["centrifuge_wash_eff"]
        w_cake = wm / (1.0 - wm) * (1.0 - eta) * CP / p["rho_liquid_kg_m3"]
        mg_ratio = d.get("Mg_stock", 0.0) / max(d["S_stock"], 1e-12)
        C_Mg_liquor = x[11] * mg_ratio
        d32 = x[7] / x[6] if x[6] > 0 else 0.0
        return {
            "sigma": sigma, "G_m_s": G, "B_m3_s": B, "M_T_kg_m3": MT,
            "crystal_t_h": 3.6 * crystal_kg_s, "Y_sol": Y_sol, "Y_liq": Y_liq,
            "Y_cierre": Y_sol - Y_liq, "recovery_pct": 100.0 * Y_sol,
            "tau_h": V / Qout / 3600.0, "d32_um": 1e6 * d32, "d50_um": 1e6 * dpbe.Q50 / 3.0 * d32,
            "N_per_kg_water": x[4] / (p["rho_liquid_kg_m3"] * p["water_mass_fraction_liquor"]),
            "liquor_density_kg_m3": rho_L,
            "slurry_density_kg_m3": rho_L + MT * (1.0 - rho_L / p["rho_crystal_kg_m3"]),
            "ppm_ClO4_crystal": 1e6 * w_cryst, "ppm_ClO4_cake": 1e6 * w_cake,
            "ppm_ClO4": 1e6 * (w_cryst + w_cake),
            "ClO4_liquor_g_L": CP, "evap_t_h": 3.6 * u["m_evap"], "purge_m3_h": 3600.0 * a * Qout,
            "recycle_m3_h": 3600.0 * (1 - a) * Qout, "Mg_liquor_g_L": C_Mg_liquor
        }


def simulate(plant: Plant, x0, t_end_s, dt_s, controller, disturbance=None, log_every=6):
    """Simula con `controller(t, x, plant) -> dict u` y `disturbance(t) -> dict d` (None = nominal)."""
    import pandas as pd
    x = np.array(x0, dtype=float)
    rows = []
    n = int(round(t_end_s / dt_s))
    for i in range(n + 1):
        t = i * dt_s
        d = plant.d_nom if disturbance is None else disturbance(t)
        u = controller(t, x, plant)
        if i % log_every == 0:
            m = plant.metrics(x, u, d)
            rows.append({"t_h": t / 3600.0, **{k: x[j] for k, j in IDX.items()},
                         **{f"u_{k}": v for k, v in u.items()}, **m})
        if i < n:
            x = plant.step(x, u, d, dt_s)
    return pd.DataFrame(rows)


#Se incorpora al final de **msmpr_core/plant.py** la función `simulate_plant_npt3_calibrated`, que extiende el modelo de **A3110** hacia las etapas acopladas de **B2510** (Área 2000), **A3210**, **B3520** (con el lazo de purga **AIC 3255**) y el Área 6000 de Refinación ($w_{\text{KNO}_3,\text{in}} = 0.95678$, especificación $< 40\text{ ppm}$):
# === ACOPLAMIENTO INTEGRADO PLANTA NPT-3 (ÁREAS 2000 - 3000 - 6000) ===
from msmpr_core.control import DelayedSensor, CommercialPID
import pandas as pd

# Constantes conciliadas desde Balance_de_masa_segun_PFD_completo.xlsx, v6 y Base de datos.xlsx
PFD_MG_CLO4 = {
    "w_Mg_208_ppm": 10206.0,      # Corriente 208 SSK (0.010206 total; 10858 ppm líquido)
    "w_ClO4_208_ppm": 3892.0,     # Corriente 208 SSK (0.003892 total; 4140 ppm líquido)
    "w_Mg_234_ppm": 2397.0,       # Corriente 234 ML Refinación (0.002397 total; 5.77 g/L BD)
    "w_ClO4_234_ppm": 738.0,      # Corriente 234 ML Refinación (0.000738 total; 1.59 g/L BD)
    "w_Mg_303_ppm": 329.0,        # Corriente 303 Agua dilución (0.000329)
    "w_ClO4_303_ppm": 36.0,       # Corriente 303 Agua dilución (0.000036)
    "w_Mg_pre_nom_ppm": (356.85 * 10206.0 + 21.33 * 2397.0) / 378.18,   # 9765.5 ppm antes de agua 303
    "w_ClO4_pre_nom_ppm": (356.85 * 3892.0 + 21.33 * 738.0) / 378.18,   # 3714.1 ppm antes de agua 303
    "w_Mg_1_nom_ppm": 9134.0,     # Corriente 1/334 SSK-d entrada A3110
    "w_ClO4_1_nom_ppm": 3468.0,   # Corriente 1/334 SSK-d entrada A3110
    "wL_Mg_A3110_ppm": 10767.0,   # Licor salida A3110 (corriente 6)
    "wL_Mg_A3210_ppm": 11720.0,   # Licor salida A3210 / S3540 (corriente 11)
    "wL_Mg_lav_3000_ppm": 10748.0,# Licor lavado S3550 (corrientes 15-24)
    "wL_ClO4_lav_3000_ppm": 4070.0,# Licor lavado S3550 (corrientes 15-24)
    "C_Mg_SSK_lab_g_L": 32.1115,  # Base de datos.xlsx hoja SSK (AQ SSKd 2026: 28.68 g/L)
    "C_Mg_MLK_lab_g_L": 31.0317,  # Base de datos.xlsx hoja Purga-Recirculación (MLK)
    "C_ClO4_MLK_lab_g_L": 3.5862, # Base de datos.xlsx hoja Purga-Recirculación (KClO4: 5.00 g/L)
    "C_Mg_MLref_lab_g_L": 5.7698, # Base de datos.xlsx hoja ML Refi
    "w_Mg_MLK_v6_ppm": 10578.0,   # MLK purga/recirculacion v6 (0.010578)  <- faltaba: KeyError en simular_recirculacion_y_magnesio
    "w_ClO4_MLK_v6_ppm": 4012.0,  # MLK purga/recirculacion v6 (0.004012)
}



# Solubilidad de KClO4 en funcion de la temperatura, expresada como ClO4- disuelto [g/L] (KClO4 [g/L] / 1.393).
# T <= 40 C: ensayo de SQM con SSKd (Mg 27.1 g/L, KClO4 5.7 g/L), Informaciones_Varias.docx: 25 C 3.48, 30 C 3.90, 35 C 4.59, 40 C 6.54 g/L KClO4
#   (el punto de inflexion experimental esta en ~40 C; bajo esa T la concentracion disuelta cae por precipitacion de KClO4).
# T > 40 C: NO hay dato de SQM. Se usa la razon S(T)/S(40) de la solubilidad de KClO4 en AGUA PURA (valores de referencia tipo CRC Handbook,
#   g/100 g agua: 40 C 3.73; 50 C 5.44; 60 C 7.3; 70 C 10.3; 80 C 13.4) -> PLACEHOLDER hasta medir la curva en salmuera real.
KCLO4_SOL_EXP = {25.0: 3.48 / 1.393, 30.0: 3.90 / 1.393, 35.0: 4.59 / 1.393, 40.0: 6.54 / 1.393}
KCLO4_SOL_AGUA_REF = {40.0: 3.73, 50.0: 5.44, 60.0: 7.3, 70.0: 10.3, 80.0: 13.4}
C_SAT_NOM_ClO4 = 4.50       # g/L de ClO4- (equivale a ~6.3 g/L de KClO4) a 40 C sin efecto de Mg: valor historico del modelo


def sat_factor_T(T_C):
    """Razon S(T)/S(40 C) de la solubilidad de KClO4: tabla experimental bajo 40 C; razon de agua pura sobre 40 C (placeholder)."""
    if T_C <= 40.0:
        xs = sorted(KCLO4_SOL_EXP)
        return float(np.interp(T_C, xs, [KCLO4_SOL_EXP[x] for x in xs])) / KCLO4_SOL_EXP[40.0]
    xs = sorted(KCLO4_SOL_AGUA_REF)
    return float(np.interp(T_C, xs, [KCLO4_SOL_AGUA_REF[x] for x in xs])) / KCLO4_SOL_AGUA_REF[40.0]


def c_sat_clo4_gL(T_C, d_Mg=0.0, salting_coef=0.85):
    """ClO4- disuelto de saturacion de KClO4 [g/L] a T_C, con salting-out por Mg (d_Mg = fraccion sobre el nominal)."""
    return C_SAT_NOM_ClO4 * sat_factor_T(T_C) - salting_coef * d_Mg



def simulate_plant_npt3_calibrated(
    use_AIT3105=True,
    use_AIC3255=True,
    tuning_mode="TL_damped",
    duration_min=300.0,
    dt_min=0.05,
    theta_AIT_min=10.0,
    density_cascade=False,
    Kc_T_override=None,
    tauI_T_override=None,
    tauD_T_override=None,
    vac_limits_override=None,
):
    """Simula el tren integrado Área 2000 (B2510) -> Área 3000 (A3110/A3210/B3520) -> Área 6000.

    Integra las composiciones de Mg2+ y ClO4- de Balance_de_masa_segun_PFD_completo.xlsx, v6 y Base de datos.xlsx.
    """
    n_steps = int(duration_min / dt_min)
    w_ClO4_pre = PFD_MG_CLO4["w_ClO4_pre_nom_ppm"]
    w_Mg_pre = PFD_MG_CLO4["w_Mg_pre_nom_ppm"]
    T1, T2 = 55.0, 40.0
    v_T1, v_T2 = 0.0, 0.0
    xS1, xS2 = 0.130, 0.200
    D50_1N, D50_2N = GRANULO["D50_nom1_um"], GRANULO["D50_nom2_um"]
    d50_1, d50_2 = D50_1N, D50_2N
    d50_k1 = d50_k2 = 1.0
    sigma1 = sigma2 = float("nan")
    n_upd = max(int(round(GRANULO["update_every_min"] / dt_min)), 1)
    C_ClO4_B3520 = 4.00
    C_Mg_B3520 = PFD_MG_CLO4["C_Mg_MLK_lab_g_L"]
    w_ClO4_ML_ref = PFD_MG_CLO4["w_ClO4_234_ppm"]
    w_Mg_ML_ref = PFD_MG_CLO4["w_Mg_234_ppm"]
    F_303_actual = F303_NOM

    is_osc = tuning_mode == "ZN_oscillatory"
    is_hopf = tuning_mode == "Hopf_limit_cycle"

    sens_TIR3101 = DelayedSensor(1.20, 0.5, dt_min, T1)
    sens_TIR3201 = DelayedSensor(1.20, 0.5, dt_min, T2)
    sens_DIR3101 = DelayedSensor(0.80, 0.4, dt_min, xS1, resolution=0.0002)
    sens_DIR3201 = DelayedSensor(0.80, 0.4, dt_min, xS2, resolution=0.0002)
    sens_AIT3105 = DelayedSensor(theta_AIT_min, 2.0, dt_min, w_ClO4_pre)
    sens_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520)
    sens_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual)

    # Sintonias de T: Hopf = Kcu (limite de estabilidad); Z-N PI = 0.45 Kcu, tauI = Pu/1.2 (decaimiento ~1/4 en la simulacion);
    # T-L amortiguada = valores historicos (Kc = 0.42 Kcu, con accion derivativa).
    if is_hopf:
        Kc_T, tauI_T, tauD_T = KCU_SIM["Kcu_bar_C"], KCU_SIM["Pu_min"] / 1.2, 0.0
        Kc_D, tauI_D = -2100.0, 8.5
    elif is_osc:
        Kc_T, tauI_T, tauD_T = 0.45 * KCU_SIM["Kcu_bar_C"], KCU_SIM["Pu_min"] / 1.2, 0.0
        Kc_D, tauI_D = -1350.0, 11.0
    else:
        Kc_T, tauI_T, tauD_T = 0.0052, 9.5, 0.60
        Kc_D, tauI_D = -650.0, 28.0

    # Sobrescritura de la sintonia de T (identificacion de Kcu, analisis de sensibilidad). Kc_T en [bar/C] antes del reescalado por etapa.
    Kc_T = Kc_T if Kc_T_override is None else Kc_T_override
    tauI_T = tauI_T if tauI_T_override is None else tauI_T_override
    tauD_T = tauD_T if tauD_T_override is None else tauD_T_override
    lim1, lim2 = (VAC["P1_lim"], VAC["P2_lim"]) if vac_limits_override is None else (vac_limits_override, vac_limits_override)
    pid_TP1 = CommercialPID(Kc_T * VAC["Kc_scale1"], tauI_T, tauD_T, 55.0, VAC["P1_nom"], *lim1)
    pid_TP2 = CommercialPID(Kc_T * VAC["Kc_scale2"], tauI_T, tauD_T, 40.0, VAC["P2_nom"], *lim2)
    # Cascada DIC (densidad, maestro lento) -> LIC (nivel, esclavo rapido) -> caudal de descarga (SIC 3101/3201).
    # El maestro de densidad ya no manda el caudal: ajusta la consigna de nivel (retencion) dentro de LEVEL["L_sp_lim"].
    # La severidad del modo (ZN / Hopf) se conserva escalando su ganancia respecto de la sintonia amortiguada (-650).
    k_dens = LEVEL["Kc_dens"] * (Kc_D / -650.0)
    M1_nom, M2_nom = LEVEL["tau1_h"] * LEVEL["F_nom"], LEVEL["tau2_h"] * LEVEL["F_nom"]
    M1, M2 = M1_nom, M2_nom
    L1, L2 = 1.0, 1.0
    pid_DIR1 = CommercialPID(k_dens, tauI_D, 0.0, 0.130, 1.0, *LEVEL["L_sp_lim"])
    pid_DIR2 = CommercialPID(k_dens, tauI_D, 0.0, 0.200, 1.0, *LEVEL["L_sp_lim"])
    kc_L1 = -60.0 * M1_nom / LEVEL["tauc_min"]
    kc_L2 = -60.0 * M2_nom / LEVEL["tauc_min"]
    ti_L = LEVEL["tauI_factor"] * LEVEL["tauc_min"]
    pid_LIC1 = CommercialPID(kc_L1, ti_L, 0.0, 1.0, LEVEL["F_nom"], 250.0, 520.0)
    pid_LIC2 = CommercialPID(kc_L2, ti_L, 0.0, 1.0, LEVEL["F_nom"], 250.0, 520.0)
    sens_LIC1 = DelayedSensor(0.20, 0.20, dt_min, 1.0)
    sens_LIC2 = DelayedSensor(0.20, 0.20, dt_min, 1.0)
    pid_AIC3255 = CommercialPID(-0.48 if is_osc else -0.28, 18.0 if is_osc else 45.0, 0.0, 4.00, 0.22, 0.08, 0.65)
    pid_FIRC3103 = CommercialPID(Kc=1.25, tau_I=0.45, tau_D=0.0, sp=27.14, u_bias=27.14, u_min=12.0, u_max=F303_MAX_AMBIENTAL)

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208 = 356.85 * (1.08 if 25.0 <= t < 140.0 else 1.0)
        T_in_1 = 70.0 + (3.5 if t >= 25.0 else 0.0)
        w_ClO4_208 = PFD_MG_CLO4["w_ClO4_208_ppm"] + (650.0 if t >= 80.0 else 0.0)
        w_Mg_208 = PFD_MG_CLO4["w_Mg_208_ppm"]
        F_234 = 21.33
        F_pre = F_208 + F_234

        w_ClO4_pre_target = (F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref) / F_pre
        w_Mg_pre_target = (F_208 * w_Mg_208 + F_234 * w_Mg_ML_ref) / F_pre
        w_ClO4_pre += (dt_min / 8.0) * (w_ClO4_pre_target - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * (w_Mg_pre_target - w_Mg_pre)

        R_0 = R_0_NOM
        if use_AIT3105:
            ait3105_pv = sens_AIT3105.read(w_ClO4_pre)
            R_eff = R_0 + K_TRIM_AIT * (ait3105_pv - PFD_MG_CLO4["w_ClO4_pre_nom_ppm"])
        else:
            R_eff = R_0

        F_303_sp = np.clip(F_pre * R_eff, 12.0, F303_MAX_AMBIENTAL)
        pid_FIRC3103.sp = F_303_sp
        valve_cmd = pid_FIRC3103.step(sens_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (valve_cmd - F_303_actual)
        F_303_actual = min(F_303_actual, F303_MAX_AMBIENTAL)   # la valvula FV 3103 no puede entregar mas que el tope

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * PFD_MG_CLO4["w_ClO4_303_ppm"]) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * PFD_MG_CLO4["w_Mg_303_ppm"]) / F_1

        P_vac1 = pid_TP1.step(sens_TIR3101.read(T1), dt_min)
        P_vac2 = pid_TP2.step(sens_TIR3201.read(T2), dt_min)
        if density_cascade:
            pid_LIC1.sp = pid_DIR1.step(sens_DIR3101.read(xS1), dt_min)
            pid_LIC2.sp = pid_DIR2.step(sens_DIR3201.read(xS2), dt_min)
        F_6 = pid_LIC1.step(sens_LIC1.read(L1), dt_min)
        F_11 = pid_LIC2.step(sens_LIC2.read(L2), dt_min)
        M1 += (F_1 - F_VAPOR_A3110 - F_PERDIDAS_A3110 - F_6) * dt_min / 60.0
        M2 += (F_6 - F_11) * dt_min / 60.0
        L1, L2 = M1 / M1_nom, M2 / M2_nom
        alpha_purge = pid_AIC3255.step(sens_AIC3255.read(C_ClO4_B3520), dt_min) if use_AIC3255 else 0.22

        T1_eq = vac_T(P_vac1, 1) + 0.45 * (T_in_1 - 70.0)
        T2_eq = vac_T(P_vac2, 2) + 0.40 * (T1 - 55.0) + (1.1 if t >= 25.0 else 0.0)
        tau1_T, tau2_T = 4.5, 1.8
        v_T1 += ((T1_eq - T1 - (tau1_T + tau2_T) * v_T1) / (tau1_T * tau2_T)) * dt_min
        v_T2 += ((T2_eq - T2 - (tau1_T + tau2_T) * v_T2) / (tau1_T * tau2_T)) * dt_min
        T1 += v_T1 * dt_min
        T2 += v_T2 * dt_min

        dil_effect = -0.0004 * (F_303_actual - 27.14)
        xS1_eq = 0.130 - 4.5e-4 * (F_6 - 390.23) - 0.003 * (T1 - 55.0) + dil_effect + (0.009 if t >= 25.0 else 0.0)
        xS1_eq += LEVEL["xS_level_gain"] * 0.130 * (L1 - 1.0)
        xS2_eq = 0.200 + 0.75 * (xS1 - 0.130) - 4.5e-4 * (F_11 - 390.23) - 0.0035 * (T2 - 40.0) + (0.008 if t >= 25.0 else 0.0)
        xS2_eq += LEVEL["xS_level_gain"] * 0.200 * (L2 - 1.0)
        xS1 += (dt_min / 9.0) * (xS1_eq - xS1)
        xS2 += (dt_min / 9.0) * (xS2_eq - xS2)

        if k % n_upd == 0:
            try:
                r1n, r2n, _, _ = _ref_nominal()
                m1, m2, sigma1, sigma2 = _cadena_d50(T1, T2, M1 / max(F_6, 1.0), M2 / max(F_11, 1.0), F_pre, F_303_actual)
                if m1 > 0.0 and m2 > 0.0:
                    d50_k1, d50_k2 = m1 / r1n, m2 / r2n
            except (ValueError, RuntimeError):
                pass
        s2 = D50_2N / 185.0
        d50_1_eq = D50_1N * d50_k1
        d50_2_eq = D50_2N * d50_k2 + s2 * (180.0 * (xS2 - 0.200) + 0.35 * (390.23 - F_11))
        d50_1 += (dt_min / GRANULO["tau_d50_min"]) * (d50_1_eq - d50_1)
        d50_2 += (dt_min / GRANULO["tau_d50_min"]) * (d50_2_eq - d50_2)
        d50_floor = GRANULO["d50_floor_frac"] * D50_2N

        C_ClO4_in_eq = 4.00 * (w_ClO4_1 / 3468.0) ** 1.45 * (0.22 / max(alpha_purge, 0.05)) ** 0.50
        C_Mg_in_eq = PFD_MG_CLO4["C_Mg_MLK_lab_g_L"] * (w_Mg_1 / 9134.0) * (0.22 / max(alpha_purge, 0.05)) ** 0.35
        C_ClO4_B3520 += (dt_min / 24.0) * (C_ClO4_in_eq - C_ClO4_B3520)
        C_Mg_B3520 += (dt_min / 24.0) * (C_Mg_in_eq - C_Mg_B3520)

        coprecip_ppm = max(C_ClO4_B3520 - 4.50, 0.0) * 82.0
        moisture_17 = 0.050 * (D50_2N / max(d50_2, d50_floor)) ** 0.65
        w_ClO4_lav_3000 = PFD_MG_CLO4["wL_ClO4_lav_3000_ppm"] * (C_ClO4_B3520 / 4.00)
        w_Mg_lav_3000 = PFD_MG_CLO4["wL_Mg_lav_3000_ppm"] * (C_Mg_B3520 / PFD_MG_CLO4["C_Mg_MLK_lab_g_L"])
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm
        ppm_Mg_17 = moisture_17 * w_Mg_lav_3000

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (D50_2N / max(d50_2, d50_floor)) ** 0.85
        w_ClO4_ML_ref_eq = PFD_MG_CLO4["w_ClO4_234_ppm"] * (ppm_ClO4_17 / 203.5)
        w_Mg_ML_ref_eq = PFD_MG_CLO4["w_Mg_234_ppm"] * (ppm_Mg_17 / 537.4)
        w_ClO4_ML_ref += (dt_min / 30.0) * (w_ClO4_ML_ref_eq - w_ClO4_ML_ref)
        w_Mg_ML_ref += (dt_min / 30.0) * (w_Mg_ML_ref_eq - w_Mg_ML_ref)
        ppm_superficial_ref = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.4
        ppm_ClO4_final_6000 = ppm_superficial_ref + 0.80 * coprecip_ppm
        ppm_Mg_final_6000 = ppm_Mg_17 * (1.0 - rem_eff_6000) * (w_Mg_ML_ref / 2397.0) ** 0.4

        hist.append({
            "t_min": t, "T1_C": T1, "T2_C": T2, "P_vac1_bar": P_vac1, "P_vac2_bar": P_vac2,
            "xS1": xS1, "xS2": xS2, "F_6_t_h": F_6, "F_11_t_h": F_11,
            "F_pre_t_h": F_pre, "w_ClO4_pre_ppm": w_ClO4_pre, "w_ClO4_1_ppm": w_ClO4_1,
            "w_Mg_pre_ppm": w_Mg_pre, "w_Mg_1_ppm": w_Mg_1, "C_Mg_B3520_g_L": C_Mg_B3520,
            "dMass_A3110_t_h": F_1 - F_VAPOR_A3110 - F_PERDIDAS_A3110 - F_6,
            "dMass_A3210_t_h": F_6 - F_11, "L1": L1, "L2": L2, "M1_t": M1, "M2_t": M2,
            "L1_sp": pid_LIC1.sp, "L2_sp": pid_LIC2.sp,
            "R_eff": R_eff, "F_303_sp_t_h": F_303_sp, "F_303_actual_t_h": F_303_actual,
            "C_ClO4_B3520_g_L": C_ClO4_B3520, "alpha_purge_pct": 100.0 * alpha_purge,
            "d50_1_um": d50_1, "d50_2_um": d50_2, "sigma1": sigma1, "sigma2": sigma2,
            "level1_pct": LEVEL["level_sp_pct"] * L1, "level2_pct": LEVEL["level_sp_pct"] * L2,
            "ppm_ClO4_crudo_17": ppm_ClO4_17, "ppm_ClO4_final_6000": ppm_ClO4_final_6000,
            "ppm_Mg_crudo_17": ppm_Mg_17, "ppm_Mg_final_6000": ppm_Mg_final_6000,
        })
    return pd.DataFrame(hist)

# === FUNCIONES CENTRALIZADAS: MAGNESIO (Mg2+), RECIRCULACIÓN B3520, FLUJO ÓPTIMO Y ROBUSTEZ ===


def simular_planta_agua_dilucion_magnesio(
    strategy="AIT3105_and_Mg_trim",   # "ratio_FFIC3102", "high_water_manual", "AIT3105_only", "AIT3105_and_Mg_trim"
    recycle_mode="no_recycle_actual", # "no_recycle_actual", "fixed_recycle", "smart_AIC3255"
    include_noise=True,
    theta_AIT_min=10.0,
    duration_min=360.0,
    dt_min=0.05,
    trim_T2=False,       # supervisor: sube T2_sp (minimo necesario) si ClO4/C_sat supera margen_sat (+solubilidad de KClO4, -KNO3)
    margen_sat=0.93,
    dT2_max=4.0,
):
    """Simula el aporte del agua de dilución (F_303) y trim de Mg2+ bajo ruido coloreado de Ornstein-Uhlenbeck."""
    n_steps = int(duration_min / dt_min)
    rng = np.random.default_rng(101)
    phi = np.exp(-dt_min / 8.0)

    def gen_ou(sigma):
        if not include_noise:
            return np.zeros(n_steps)
        e = rng.standard_normal(n_steps)
        x = np.zeros(n_steps)
        x[0] = e[0]
        for i in range(1, n_steps):
            x[i] = phi * x[i - 1] + np.sqrt(1.0 - phi**2) * e[i]
        return sigma * x

    n_F, n_T, n_Cl, n_Mg, n_V = gen_ou(0.015), gen_ou(0.45), gen_ou(0.020), gen_ou(0.025), gen_ou(0.008)

    w_ClO4_pre, w_Mg_pre = 3714.1, 8200.0
    w_ClO4_1, w_Mg_1 = 3468.0, 7650.0
    C_ClO4_A3110, C_Mg_A3110 = 3.65, 8.10
    C_ClO4_A3210, C_Mg_A3210 = 3.95, 8.95
    C_ClO4_B3520, C_Mg_B3520 = 4.00, 9.10
    T1, T2 = 55.0, 40.0
    dT2_trim = 0.0
    P_vac1, P_vac2 = VAC["P1_nom"], VAC["P2_nom"]
    xS1, xS2 = 0.130, 0.200
    d50_1, d50_2 = 145.0, 185.0
    w_ClO4_ML_ref = 738.0
    F_303_actual = 27.14
    M_KNO3_prod = 67.5 if recycle_mode == "no_recycle_actual" else 74.5

    ns = 1.0 if include_noise else 0.0
    s_TIR1 = DelayedSensor(1.0, 0.4, dt_min, T1, noise_std=0.05 * ns, seed=11)
    s_TIR2 = DelayedSensor(1.0, 0.4, dt_min, T2, noise_std=0.05 * ns, seed=12)
    s_DIR1 = DelayedSensor(0.6, 0.3, dt_min, xS1, resolution=0.0002, noise_std=0.0003 * ns, seed=13)
    s_DIR2 = DelayedSensor(0.6, 0.3, dt_min, xS2, resolution=0.0002, noise_std=0.0003 * ns, seed=14)
    s_AIT3105 = DelayedSensor(theta_AIT_min, 2.0, dt_min, w_ClO4_pre, noise_std=18.0 * ns, seed=15)
    s_AIT_Mg = DelayedSensor(theta_AIT_min, 2.0, dt_min, w_Mg_pre, noise_std=35.0 * ns, seed=16)
    s_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520, noise_std=0.025 * ns, seed=17)
    s_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual, noise_std=0.15 * ns, seed=18)

    pid_TIR1 = CommercialPID(0.0055 * VAC["Kc_scale1"], 9.0, 0.4, 55.0, VAC["P1_nom"], *VAC["P1_lim"], beta=0.5)
    pid_TIR2 = CommercialPID(0.0055 * VAC["Kc_scale2"], 9.0, 0.4, 40.0, VAC["P2_nom"], *VAC["P2_lim"], beta=0.5)
    pid_PIR1 = CommercialPID(1.5, 0.3, 0.0, VAC["P1_nom"], VAC["P1_nom"], *VAC["P1_lim"])
    pid_PIR2 = CommercialPID(1.5, 0.3, 0.0, VAC["P2_nom"], VAC["P2_nom"], *VAC["P2_lim"])
    pid_DIR1 = CommercialPID(-650.0, 28.0, 0.0, 0.130, 390.23, 250.0, 520.0)
    pid_DIR2 = CommercialPID(-650.0, 28.0, 0.0, 0.200, 390.23, 250.0, 520.0)
    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, F303_MAX_AMBIENTAL)
    pid_AIC3255 = CommercialPID(-0.42, 20.0, 0.0, 3.92, 0.25, 0.12, 1.00)

    R_0 = R_0_NOM
    hist = []

    for k in range(n_steps):
        t = k * dt_min
        F_208_base = 356.85 * (1.10 if 25.0 <= t < 105.0 else (0.93 if 250.0 <= t < 295.0 else 1.0))
        F_208 = F_208_base * (1.0 + n_F[k])
        T_in_1 = 70.0 + (3.5 if 25.0 <= t < 125.0 else 0.0) + n_T[k]
        w_ClO4_208 = (3892.0 + (520.0 if t >= 75.0 else 0.0)) * (1.0 + n_Cl[k])
        w_Mg_208 = (8200.0 * (1.42 if t >= 155.0 else 1.0)) * (1.0 + n_Mg[k])
        valve_gain = (0.82 if 230.0 <= t < 285.0 else 1.0) * (1.0 + n_V[k])
        F_234 = 21.33

        if recycle_mode == "no_recycle_actual":
            alpha_purge, F_rec = 1.00, 0.0
        elif recycle_mode == "fixed_recycle":
            alpha_purge, F_rec = 0.25, 80.0 * 0.75
        else:
            d_Mg_b3520 = max((C_Mg_B3520 - 9.10) / 9.10, 0.0)
            pid_AIC3255.sp = 3.92 - 0.45 * d_Mg_b3520
            alpha_purge = pid_AIC3255.step(s_AIC3255.read(C_ClO4_B3520), dt_min)
            F_rec = 80.0 * (1.0 - alpha_purge)

        F_pre = F_208 + F_234 + F_rec
        w_ClO4_pre += (dt_min / 8.0) * ((F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref + F_rec * (C_ClO4_B3520 * 815.0)) / F_pre - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * ((F_208 * w_Mg_208 + F_234 * 1200.0 + F_rec * (C_Mg_B3520 * 815.0)) / F_pre - w_Mg_pre)

        ait_cl_pv = s_AIT3105.read(w_ClO4_pre)
        ait_mg_pv = s_AIT_Mg.read(w_Mg_pre)

        if strategy == "ratio_FFIC3102":
            R_eff = R_0
        elif strategy == "high_water_manual":
            R_eff = 0.115
        elif strategy == "AIT3105_only":
            R_eff = R_0 + K_TRIM_AIT * max(ait_cl_pv - 3714.1, 0.0)
        else:
            d_mg_in = max((ait_mg_pv - 8200.0) / 8200.0, 0.0)
            R_eff = R_0 + K_TRIM_AIT * max(ait_cl_pv - 3714.1, 0.0) + 0.065 * d_mg_in

        F_303_sp = np.clip(F_pre * R_eff, 12.0, F303_MAX_AMBIENTAL)
        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (valve_gain * v_cmd - F_303_actual)
        F_303_actual = min(F_303_actual, F303_MAX_AMBIENTAL)   # la valvula FV 3103 no puede entregar mas que el tope

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * 10.0) / F_1

        pid_PIR1.sp = pid_TIR1.step(s_TIR1.read(T1), dt_min)
        pid_TIR2.sp = 40.0 + dT2_trim
        pid_PIR2.sp = pid_TIR2.step(s_TIR2.read(T2), dt_min)
        P_vac1 += (dt_min / 0.25) * (pid_PIR1.step(P_vac1, dt_min) - P_vac1)
        P_vac2 += (dt_min / 0.25) * (pid_PIR2.step(P_vac2, dt_min) - P_vac2)

        T1 += (dt_min / 5.5) * (vac_T(P_vac1, 1) + 0.40 * (T_in_1 - 70.0) - 0.025 * (F_303_actual - 27.14) - T1)
        T2 += (dt_min / 5.5) * (vac_T(P_vac2, 2) + 0.35 * (T1 - 55.0) - T2)

        F_6 = pid_DIR1.step(s_DIR1.read(xS1), dt_min)
        F_11 = pid_DIR2.step(s_DIR2.read(xS2), dt_min)

        xS1 += (dt_min / 9.0) * (0.130 - 4.5e-4 * (F_6 - 390.23) - 0.003 * (T1 - 55.0) - 0.00040 * (F_303_actual - 27.14) - xS1)
        xS2 += (dt_min / 9.0) * (0.200 + 0.70 * (xS1 - 0.130) - 4.5e-4 * (F_11 - 390.23) - 0.0035 * (T2 - 40.0) - xS2)

        rec_factor = (1.0 / max(alpha_purge, 0.12)) ** 0.18
        C_Mg_A3110 += (dt_min / 11.0) * (8.10 * (w_Mg_1 / 7650.0) - C_Mg_A3110)
        C_ClO4_A3110 += (dt_min / 11.0) * (3.65 * (w_ClO4_1 / 3468.0) - C_ClO4_A3110)
        C_Mg_A3210 += (dt_min / 13.0) * (8.95 * (C_Mg_A3110 / 8.10) - C_Mg_A3210)
        C_ClO4_A3210 += (dt_min / 13.0) * (3.95 * (C_ClO4_A3110 / 3.65) ** 1.18 - C_ClO4_A3210)
        C_Mg_B3520 += (dt_min / 16.0) * (9.10 * (C_Mg_A3210 / 8.95) * rec_factor - C_Mg_B3520)
        C_ClO4_B3520 += (dt_min / 16.0) * (4.00 * (C_ClO4_A3210 / 3.95) * rec_factor - C_ClO4_B3520)

        d_Mg_proc = max((C_Mg_A3210 - 8.95) / 8.95, 0.0)
        water_added = max(F_303_actual - 27.14, 0.0)
        C_sat_KClO4_40C = (4.38 - 1.15 * d_Mg_proc + 0.042 * water_added) * sat_factor_T(T2)   # satura a la T real de A3210
        if trim_T2:
            exceso = C_ClO4_B3520 / max(C_sat_KClO4_40C, 1e-6) - margen_sat
            dT2_trim += (dt_min / 10.0) * (float(np.clip(exceso * 60.0, 0.0, dT2_max)) - dT2_trim)

        d50_1 += (dt_min / 14.0) * (145.0 - 25.0 * max((C_Mg_A3110 - 8.10) / 8.10, 0.0) + 0.55 * water_added - d50_1)
        d50_2 += (dt_min / 16.0) * (185.0 * (d50_1 / 145.0) - 28.0 * d_Mg_proc + 0.65 * water_added - d50_2)

        keff_mult = 1.0 + 0.38 * d_Mg_proc
        coprecip_ppm = max(C_ClO4_B3520 - C_sat_KClO4_40C, 0.0) * 85.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 100.0)) ** 0.65 * (1.0 + 0.12 * d_Mg_proc)
        w_ClO4_lav_3000 = 4070.0 * (C_ClO4_B3520 / 4.00) * keff_mult
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm

        excess_water = max(F_303_actual - (F_pre * R_0), 0.0)
        M_KNO3_eq = (67.5 + 8.8 * (1.0 - alpha_purge)) * (F_208 / 356.85) - 0.36 * excess_water - KNO3_PERDIDA_T_t_h_por_C * max(T2 - 40.0, 0.0)
        M_KNO3_prod += (dt_min / 15.0) * (M_KNO3_eq - M_KNO3_prod)

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 100.0)) ** 0.85
        w_ClO4_ML_ref += (dt_min / 24.0) * (738.0 * (ppm_ClO4_17 / 203.5) - w_ClO4_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.38 + 0.82 * coprecip_ppm

        hist.append({
            "t_min": t, "F_303": F_303_actual, "w_Mg_1": w_Mg_1,
            "C_Mg_A3110": C_Mg_A3110, "C_Mg_A3210": C_Mg_A3210, "C_Mg_B3520": C_Mg_B3520,
            "C_ClO4_B3520": C_ClO4_B3520, "C_sat_KClO4_40C": C_sat_KClO4_40C, "d50_2": d50_2,
            "M_KNO3_t_h": M_KNO3_prod, "dT2_trim_C": dT2_trim, "T2_C": T2, "ppm_ClO4_final_6000": ppm_final_6000,
        })
    return pd.DataFrame(hist)


def simular_AIT_AIC_por_flujo_agua(
    F_303_target=F303_MAX_AMBIENTAL,
    modo_control="AIT_and_AIC",
    delta_ClO4_ssk=650.0,
    include_noise=True,
    duration_min=320.0,
    dt_min=0.05,
):
    """Evalúa los lazos AIT 3105 y AIC 3255 en función del flujo objetivo de agua F_303 (el flujo se limita a F303_MAX_AMBIENTAL)."""
    n_steps = int(duration_min / dt_min)
    rng = np.random.default_rng(105)
    phi = np.exp(-dt_min / 8.0)

    def gen_ou(sigma):
        if not include_noise:
            return np.zeros(n_steps)
        e = rng.standard_normal(n_steps)
        x = np.zeros(n_steps)
        x[0] = e[0]
        for i in range(1, n_steps):
            x[i] = phi * x[i - 1] + np.sqrt(1.0 - phi**2) * e[i]
        return sigma * x

    n_F, n_Cl, n_Mg = gen_ou(0.012), gen_ou(0.018), gen_ou(0.020)

    w_ClO4_pre, w_Mg_pre = 3714.1, 8400.0
    C_ClO4_A3110, C_Mg_A3110 = 3.65, 8.15
    C_ClO4_A3210, C_Mg_A3210 = 3.95, 8.95
    C_ClO4_B3520, C_Mg_B3520 = 4.00, 9.10
    v_B3520 = 0.0
    d50_2 = 185.0
    w_ClO4_ML_ref = 738.0
    F_303_actual = 27.14
    alpha_purge = 0.25
    M_KNO3_prod = 72.8

    ns = 1.0 if include_noise else 0.0
    s_AIT3105 = DelayedSensor(10.0, 2.0, dt_min, w_ClO4_pre, noise_std=15.0 * ns, seed=21)
    s_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520, noise_std=0.02 * ns, seed=22)
    s_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual, noise_std=0.12 * ns, seed=23)

    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, F303_MAX_AMBIENTAL)
    pid_AIC3255 = CommercialPID(-0.35, 22.0, 0.0, 4.02, 0.25, 0.12, 0.80)

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208 = 356.85 * (1.0 + n_F[k])
        F_234 = 21.33
        w_ClO4_208 = (3892.0 + (delta_ClO4_ssk if t >= 45.0 else 0.0)) * (1.0 + n_Cl[k])
        w_Mg_208 = (8400.0 * (1.32 if t >= 115.0 else 1.0)) * (1.0 + n_Mg[k])

        aic_pv = s_AIC3255.read(C_ClO4_B3520)
        alpha_purge = pid_AIC3255.step(aic_pv, dt_min) if modo_control == "AIT_and_AIC" else 0.25

        F_rec = 80.0 * (1.0 - alpha_purge)
        F_pre = F_208 + F_234 + F_rec

        w_ClO4_pre += (dt_min / 8.0) * ((F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref + F_rec * (C_ClO4_B3520 * 815.0)) / F_pre - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * ((F_208 * w_Mg_208 + F_234 * 1200.0 + F_rec * (C_Mg_B3520 * 815.0)) / F_pre - w_Mg_pre)

        ait_pv = s_AIT3105.read(w_ClO4_pre)

        if t < 45.0:
            F_303_sp = 27.14
        else:
            factor_ait = np.clip((ait_pv - 3714.1) / max(delta_ClO4_ssk * 0.72, 1.0), 0.0, 1.15)
            F_303_sp = np.clip(27.14 + (F_303_target - 27.14) * factor_ait, 12.0, F303_MAX_AMBIENTAL)

        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)
        F_303_actual = min(F_303_actual, F303_MAX_AMBIENTAL)   # la valvula FV 3103 no puede entregar mas que el tope

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * 10.0) / F_1

        tau_scale = 405.3 / max(F_1, 300.0)
        C_Mg_A3110 += (dt_min / (10.0 * tau_scale)) * (8.15 * (w_Mg_1 / 7650.0) - C_Mg_A3110)
        C_ClO4_A3110 += (dt_min / (10.0 * tau_scale)) * (3.65 * (w_ClO4_1 / 3468.0) - C_ClO4_A3110)

        C_Mg_A3210 += (dt_min / (12.0 * tau_scale)) * (8.95 * (C_Mg_A3110 / 8.15) - C_Mg_A3210)
        C_ClO4_A3210 += (dt_min / (12.0 * tau_scale)) * (3.95 * (C_ClO4_A3110 / 3.65) ** 1.25 - C_ClO4_A3210)

        rec_factor = (0.25 / max(alpha_purge, 0.10)) ** 0.24
        C_ClO4_B3520_eq = 4.00 * (C_ClO4_A3210 / 3.95) * rec_factor
        tau1_b, tau2_b = 12.0 * tau_scale, 4.0 * tau_scale
        accel_b = (C_ClO4_B3520_eq - C_ClO4_B3520 - (tau1_b + tau2_b) * v_B3520) / (tau1_b * tau2_b)
        v_B3520 += accel_b * dt_min
        C_ClO4_B3520 += v_B3520 * dt_min
        C_Mg_B3520 += (dt_min / tau1_b) * (9.10 * (C_Mg_A3210 / 8.95) * (rec_factor ** 0.6) - C_Mg_B3520)

        d_Mg_proc = max((C_Mg_A3210 - 8.95) / 8.95, 0.0)
        water_delta = F_303_actual - 27.14
        C_sat_KClO4_40C = 4.02 - 1.12 * d_Mg_proc + 0.036 * water_delta

        d50_2_eq = 185.0 - 26.0 * d_Mg_proc + 0.60 * water_delta
        d50_2 += (dt_min / 15.0) * (d50_2_eq - d50_2)

        keff_mult = 1.0 + 0.38 * d_Mg_proc
        coprecip_ppm = max(C_ClO4_B3520 - C_sat_KClO4_40C, 0.0) * 125.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 105.0)) ** 0.65 * (1.0 + 0.12 * d_Mg_proc)
        w_ClO4_lav_3000 = 4070.0 * (C_ClO4_B3520 / 4.00) * keff_mult
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm

        excess_water = max(F_303_actual - 27.14, 0.0)
        purge_loss = 10.5 * max(alpha_purge - 0.25, 0.0)
        M_KNO3_eq = 73.5 - 0.36 * excess_water - purge_loss
        M_KNO3_prod += (dt_min / 14.0) * (M_KNO3_eq - M_KNO3_prod)

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 105.0)) ** 0.85
        w_ClO4_ML_ref += (dt_min / 22.0) * (738.0 * (ppm_ClO4_17 / 203.5) - w_ClO4_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.38 + 0.85 * coprecip_ppm

        hist.append({
            "t_min": t, "F_303": F_303_actual, "w_ClO4_1": w_ClO4_1, "ait_pv": ait_pv,
            "C_ClO4_B3520": C_ClO4_B3520, "C_sat_KClO4": C_sat_KClO4_40C,
            "alpha_purge_pct": 100.0 * alpha_purge, "M_KNO3_t_h": M_KNO3_prod,
            "ppm_ClO4_final_6000": ppm_final_6000, "d50_2": d50_2,
        })
    return pd.DataFrame(hist)


def simular_robustez_parametros(
    k_g_mult=1.0,
    k_b_mult=1.0,
    Keff_mult=1.0,
    tau_3000_mult=1.0,
    control_active=True,
    duration_min=320.0,
    dt_min=0.05,
    seed=101,
):
    """Evalúa la robustez del control PID frente a incertidumbre paramétrica (L-Alanina vs KNO3 real)."""
    n_steps = int(duration_min / dt_min)
    rng = np.random.default_rng(seed)
    phi = np.exp(-dt_min / 8.0)

    def gen_ou(sigma):
        e = rng.standard_normal(n_steps)
        x = np.zeros(n_steps)
        x[0] = e[0]
        for i in range(1, n_steps):
            x[i] = phi * x[i - 1] + np.sqrt(1.0 - phi**2) * e[i]
        return sigma * x

    n_F, n_Cl, n_Mg = gen_ou(0.012), gen_ou(0.018), gen_ou(0.020)
    w_ClO4_pre, w_Mg_pre = 3714.1, 8400.0
    C_ClO4_A3110, C_Mg_A3110 = 3.65, 8.15
    C_ClO4_A3210, C_Mg_A3210 = 3.95, 8.95
    C_ClO4_B3520, C_Mg_B3520 = 4.00, 9.10
    d50_base = 185.0 * (k_g_mult ** 0.28) / (k_b_mult ** 0.14)
    d50_2 = d50_base
    w_ClO4_ML_ref, F_303_actual, alpha_purge = 738.0, 27.14, 0.25

    s_AIT3105 = DelayedSensor(10.0, 2.0, dt_min, w_ClO4_pre, noise_std=14.0, seed=seed + 1)
    s_AIT_Mg = DelayedSensor(10.0, 2.0, dt_min, w_Mg_pre, noise_std=30.0, seed=seed + 2)
    s_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520, noise_std=0.018, seed=seed + 3)
    s_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual, noise_std=0.10, seed=seed + 4)

    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, F303_MAX_AMBIENTAL)
    pid_AIC3255 = CommercialPID(-0.48, 16.0, 0.0, 4.00, 0.25, 0.12, 0.85)
    R_0 = R_0_NOM

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208, F_234 = 356.85 * (1.0 + n_F[k]), 21.33
        w_ClO4_208 = (3892.0 + (620.0 if t >= 45.0 else 0.0)) * (1.0 + n_Cl[k])
        w_Mg_208 = (8400.0 * (1.32 if t >= 115.0 else 1.0)) * (1.0 + n_Mg[k])

        aic_pv = s_AIC3255.read(C_ClO4_B3520)
        if control_active:
            d_Mg_b3520 = max((C_Mg_B3520 - 9.10) / 9.10, 0.0)
            pid_AIC3255.sp = 4.00 - 0.32 * d_Mg_b3520
            alpha_purge = pid_AIC3255.step(aic_pv, dt_min)
        else:
            alpha_purge = 0.25

        F_rec = 80.0 * (1.0 - alpha_purge)
        F_pre = F_208 + F_234 + F_rec
        w_ClO4_pre += (dt_min / 8.0) * ((F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref + F_rec * (C_ClO4_B3520 * 815.0)) / F_pre - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * ((F_208 * w_Mg_208 + F_234 * 1200.0 + F_rec * (C_Mg_B3520 * 815.0)) / F_pre - w_Mg_pre)

        ait_pv, ait_mg_pv = s_AIT3105.read(w_ClO4_pre), s_AIT_Mg.read(w_Mg_pre)
        if control_active:
            d_mg_in = max((ait_mg_pv - 8400.0) / 8400.0, 0.0)
            aic_trim = 0.015 * max(aic_pv - 4.00, 0.0)
            R_eff = R_0 + K_TRIM_AIT * max(ait_pv - 3714.1, 0.0) + 0.055 * d_mg_in + aic_trim
            F_303_sp = np.clip(F_pre * R_eff, 12.0, F303_MAX_AMBIENTAL)
        else:
            F_303_sp = 27.14

        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)
        F_303_actual = min(F_303_actual, F303_MAX_AMBIENTAL)   # la valvula FV 3103 no puede entregar mas que el tope

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * 10.0) / F_1

        tau_scale = (405.3 / max(F_1, 300.0)) * tau_3000_mult
        C_Mg_A3110 += (dt_min / (10.0 * tau_scale)) * (8.15 * (w_Mg_1 / 7650.0) - C_Mg_A3110)
        C_ClO4_A3110 += (dt_min / (10.0 * tau_scale)) * (3.65 * (w_ClO4_1 / 3468.0) - C_ClO4_A3110)
        C_Mg_A3210 += (dt_min / (12.0 * tau_scale)) * (8.95 * (C_Mg_A3110 / 8.15) - C_Mg_A3210)
        C_ClO4_A3210 += (dt_min / (12.0 * tau_scale)) * (3.95 * (C_ClO4_A3110 / 3.65) ** 1.22 - C_ClO4_A3210)

        rec_factor = (0.25 / max(alpha_purge, 0.10)) ** 0.28
        C_ClO4_B3520_eq = 4.00 * ((C_ClO4_A3210 / 3.95) ** (1.05 * (Keff_mult ** 0.3))) * rec_factor
        C_ClO4_B3520 += (dt_min / (14.0 * tau_scale)) * (C_ClO4_B3520_eq - C_ClO4_B3520)
        C_Mg_B3520 += (dt_min / (14.0 * tau_scale)) * (9.10 * (C_Mg_A3210 / 8.95) * (rec_factor ** 0.6) - C_Mg_B3520)

        d_Mg_proc = max((C_Mg_A3210 - 8.95) / 8.95, 0.0)
        water_delta = F_303_actual - 27.14
        C_sat_KClO4_40C = 4.18 - 1.05 * d_Mg_proc + 0.044 * water_delta
        d50_2 += (dt_min / (15.0 * tau_3000_mult)) * (d50_base - 22.0 * d_Mg_proc + 0.55 * water_delta - d50_2)

        keff_total = (1.0 + 0.32 * d_Mg_proc) * (Keff_mult ** 0.45)
        coprecip_ppm = max(C_ClO4_B3520 - C_sat_KClO4_40C, 0.0) * 115.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 110.0)) ** 0.55 * (1.0 + 0.08 * d_Mg_proc)
        w_ClO4_lav_3000 = 4070.0 * (C_ClO4_B3520 / 4.00) * keff_total
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 110.0)) ** 0.65
        w_ClO4_ML_ref += (dt_min / 22.0) * (738.0 * (ppm_ClO4_17 / 203.5) - w_ClO4_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.35 + 0.85 * coprecip_ppm

        hist.append({
            "t_min": t, "F_303": F_303_actual, "error_AIC": pid_AIC3255.last_error,
            "I_term_AIC": pid_AIC3255.I, "P_term_AIC": pid_AIC3255.last_P,
            "alpha_purge_pct": 100.0 * alpha_purge, "ppm_ClO4_final_6000": ppm_final_6000,
        })
    return pd.DataFrame(hist)

def simular_recirculacion_y_magnesio(
    modo_reciclo="sin_reciclo_actual",  # "sin_reciclo_actual", "reciclo_fijo", "reciclo_inteligente"
    use_AIT3105=False,
    hipotesis_mg="combinada",           # "salting_out", "cinetica_particion", "combinada"
    duration_min=360.0,
    dt_min=0.05,
    aic_sp0=3.95,            # consigna base de ClO4 en B3520 para AIC 3255 [g/L]
    aic_gain_Mg=0.55,        # descenso de la consigna por fraccion de Mg sobre el nominal
    k_trim=K_TRIM_AIT,       # trim de AIT 3105 sobre la razon de dilucion [1/ppm]
    trim_Mg=0.018,           # trim por Mg sobre la razon de dilucion
    F303_max=None,           # tope de agua de dilucion [t/h]; None = F303_MAX_AMBIENTAL
    T_in_extra=0.0,          # sobretemperatura de la alimentacion previa a la dilucion [C] (escalon en t = 40 min)
    T2_sp=40.0,              # consigna de temperatura del cristalizador 2 [C]
    T1_sp=55.0,              # consigna de temperatura del cristalizador 1 [C]
    T_agua_dilucion=60.0,    # temperatura del agua de dilucion F_303 [C] (60 C en el balance de masa)
    ff_temperatura=False,    # prealimentacion: la consigna de vacio descuenta la sobretemperatura medida de la alimentacion
    trim_T2=False,           # supervisor: sube T2_sp en la cantidad minima que mantiene el margen de saturacion de KClO4
    margen_sat=0.93,         # relacion ClO4/C_sat a partir de la cual actua el supervisor de T2
    dT2_max=4.0,             # maxima elevacion de T2_sp que puede pedir el supervisor [C]
):
    """Evalúa la recirculación B3520 -> A3110 y las 3 hipótesis de interacción Mg2+ / ClO4- con datos PFD y Base de datos.xlsx.

    Nuevos argumentos de analisis: sintonia de AIC/AIT, tope de agua y temperaturas (alimentacion y consignas de cristalizadores).
    """
    F303_max = F303_MAX_AMBIENTAL if F303_max is None else min(F303_max, F303_MAX_AMBIENTAL)
    n_steps = int(duration_min / dt_min)
    w_ClO4_pre = PFD_MG_CLO4["w_ClO4_pre_nom_ppm"]
    w_Mg_pre = PFD_MG_CLO4["w_Mg_pre_nom_ppm"]
    T1, T2 = T1_sp, T2_sp
    T2_sp_nom = T2_sp
    dT2_trim = 0.0
    T2_perdida_KNO3 = 0.0
    P_vac1, P_vac2 = VAC["P1_nom"], VAC["P2_nom"]
    d50_2 = 185.0
    C_ClO4_B3520 = 3.85
    C_Mg_B3520 = 9.20                  # Fracción porcentual en peso líquido x10 (equiv. 31.03 g/L en Base de datos.xlsx)
    k_gL_Mg = PFD_MG_CLO4["C_Mg_MLK_lab_g_L"] / 9.20
    w_ClO4_ML_ref = PFD_MG_CLO4["w_ClO4_234_ppm"]
    w_Mg_ML_ref = PFD_MG_CLO4["w_Mg_234_ppm"]
    F_303_actual = F303_NOM
    M_KNO3 = 66.5 if modo_reciclo in ("sin_reciclo_actual", "no_recycle_actual") else 74.0

    s_TIR1 = DelayedSensor(1.0, 0.5, dt_min, T1)
    s_TIR2 = DelayedSensor(1.0, 0.5, dt_min, T2)
    s_AIT3105 = DelayedSensor(10.0, 2.0, dt_min, w_ClO4_pre)
    s_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520)
    s_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual)

    pid_TIR1 = CommercialPID(0.0055 * VAC["Kc_scale1"], 9.0, 0.4, T1_sp, VAC["P1_nom"], *VAC["P1_lim"])
    pid_TIR2 = CommercialPID(0.0055 * VAC["Kc_scale2"], 9.0, 0.4, T2_sp, VAC["P2_nom"], *VAC["P2_lim"])
    pid_PIR1 = CommercialPID(1.5, 0.3, 0.0, VAC["P1_nom"], VAC["P1_nom"], *VAC["P1_lim"])
    pid_PIR2 = CommercialPID(1.5, 0.3, 0.0, VAC["P2_nom"], VAC["P2_nom"], *VAC["P2_lim"])
    pid_FIRC3103 = CommercialPID(1.25, 0.45, 0.0, 27.14, 27.14, 12.0, F303_max)
    # arranque SIN recirculacion (alpha = 1): la planta parte en regimen sin reciclo y AIC solo la abre cuando hay margen
    pid_AIC3255 = CommercialPID(-0.38, 22.0, 0.0, aic_sp0, 1.00 if modo_reciclo in ("reciclo_inteligente", "smart_recycle") else 0.25, 0.12, 1.00)

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208, F_234 = 356.85, 21.33
        w_ClO4_208 = PFD_MG_CLO4["w_ClO4_208_ppm"] + (550.0 if t >= 40.0 else 0.0)
        # Escalon de Mg2+ en t=140 min (+35%, pasa de 32.1 g/L a 43.3 g/L, dentro del rango max=50.0 g/L de Base de datos.xlsx)
        w_Mg_208 = PFD_MG_CLO4["w_Mg_208_ppm"] * (1.35 if t >= 140.0 else 1.0)

        if modo_reciclo in ("sin_reciclo_actual", "no_recycle_actual"):
            alpha_purge, F_rec_B3520 = 1.00, 0.0
        elif modo_reciclo in ("reciclo_fijo", "fixed_recycle"):
            alpha_purge = 0.25
            F_rec_B3520 = 85.0 * (1.0 - alpha_purge)
        else:
            delta_Mg_sp = max((C_Mg_B3520 - 9.20) / 9.20, 0.0)
            pid_AIC3255.sp = aic_sp0 - aic_gain_Mg * delta_Mg_sp
            alpha_purge = pid_AIC3255.step(s_AIC3255.read(C_ClO4_B3520), dt_min)
            F_rec_B3520 = 85.0 * (1.0 - alpha_purge)

        F_pre = F_208 + F_234 + F_rec_B3520
        w_ClO4_rec = C_ClO4_B3520 * (PFD_MG_CLO4["wL_ClO4_lav_3000_ppm"] / 4.00)
        w_Mg_rec = C_Mg_B3520 * (PFD_MG_CLO4["w_Mg_MLK_v6_ppm"] / 9.20)
        w_ClO4_pre_eq = (F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref + F_rec_B3520 * w_ClO4_rec) / F_pre
        w_Mg_pre_eq = (F_208 * w_Mg_208 + F_234 * w_Mg_ML_ref + F_rec_B3520 * w_Mg_rec) / F_pre
        w_ClO4_pre += (dt_min / 8.0) * (w_ClO4_pre_eq - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * (w_Mg_pre_eq - w_Mg_pre)

        R_0 = R_0_NOM
        if use_AIT3105:
            ait_pv = s_AIT3105.read(w_ClO4_pre)
            delta_Mg_trim = max((w_Mg_pre - PFD_MG_CLO4["w_Mg_pre_nom_ppm"]) / PFD_MG_CLO4["w_Mg_pre_nom_ppm"], 0.0)
            R_eff = R_0 + k_trim * (ait_pv - PFD_MG_CLO4["w_ClO4_pre_nom_ppm"]) + trim_Mg * delta_Mg_trim
        else:
            R_eff = R_0

        F_303_sp = np.clip(F_pre * R_eff, 12.0, F303_max)
        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)
        F_303_actual = min(F_303_actual, F303_MAX_AMBIENTAL)   # la valvula FV 3103 no puede entregar mas que el tope

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * PFD_MG_CLO4["w_ClO4_303_ppm"]) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * PFD_MG_CLO4["w_Mg_303_ppm"]) / F_1

        dT_in = T_in_extra if t >= 40.0 else 0.0
        dT_dil = -(F_303_actual / max(F_1, 1.0)) * (60.0 - T_agua_dilucion) * (4.172 / 2.180)
        if ff_temperatura:   # prealimentacion: se corrige la consigna de T1 con la sobretemperatura medida (TIR 2501) y el agua fria
            pid_TIR1.sp = T1_sp
        pid_TIR2.sp = T2_sp_nom + dT2_trim
        pid_PIR1.sp = pid_TIR1.step(s_TIR1.read(T1), dt_min)
        if ff_temperatura:
            pid_PIR1.sp = float(np.clip(pid_PIR1.sp - (vac_P(T1_sp + 0.45 * (dT_in + dT_dil), 1) - vac_P(T1_sp, 1)), *VAC["P1_lim"]))
        pid_PIR2.sp = pid_TIR2.step(s_TIR2.read(T2), dt_min)
        P_vac1 += (dt_min / 0.25) * (pid_PIR1.step(P_vac1, dt_min) - P_vac1)
        P_vac2 += (dt_min / 0.25) * (pid_PIR2.step(P_vac2, dt_min) - P_vac2)
        # Perturbaciones termicas sobre T1 (open loop): sobretemperatura de la alimentacion previa (+0.45 K/K, como en el tren integrado)
        # y agua de dilucion mas fria que 60 C (cada K menos baja la mezcla ~ F_303/F_1 * Cp_agua/Cp_mezcla).
        dT_in = T_in_extra if t >= 40.0 else 0.0
        dT_dil = -(F_303_actual / max(F_1, 1.0)) * (60.0 - T_agua_dilucion) * (4.172 / 2.180)
        T1 += (dt_min / 5.5) * (vac_T(P_vac1, 1) + 0.45 * (dT_in + dT_dil) - T1)
        T2 += (dt_min / 5.5) * (vac_T(P_vac2, 2) + 0.35 * (T1 - T1_sp) - T2)

        rec_factor = (1.0 / max(alpha_purge, 0.12)) ** 0.28
        C_ClO4_eq = 3.85 * (w_ClO4_1 / PFD_MG_CLO4["w_ClO4_1_nom_ppm"]) ** 1.35 * rec_factor
        C_Mg_eq = 9.20 * (w_Mg_1 / PFD_MG_CLO4["w_Mg_1_nom_ppm"]) * rec_factor
        C_ClO4_B3520 += (dt_min / 22.0) * (C_ClO4_eq - C_ClO4_B3520)
        C_Mg_B3520 += (dt_min / 22.0) * (C_Mg_eq - C_Mg_B3520)

        d_Mg = max((C_Mg_B3520 - 9.20) / 9.20, 0.0)
        salting = 0.85 * d_Mg if hipotesis_mg in ("salting_out", "combinada", "combined") else 0.0
        C_sat_eff = c_sat_clo4_gL(T2, 0.0) - salting      # satura a la T real de A3210 (y cae con el Mg2+)
        if trim_T2:       # supervisor de margen: sube T2_sp (si la razon ClO4/Csat supera margen_sat) -> +solubilidad, -KNO3
            exceso = C_ClO4_B3520 / max(C_sat_eff, 1e-6) - margen_sat
            dT2_trim += (dt_min / 10.0) * (float(np.clip(exceso * 60.0, 0.0, dT2_max)) - dT2_trim)

        keff_mult = (1.0 + 0.45 * d_Mg) if hipotesis_mg in ("cinetica_particion", "kinetic_partition", "combinada", "combined") else 1.0
        d50_eq = 185.0 - (28.0 * d_Mg if hipotesis_mg in ("cinetica_particion", "kinetic_partition", "combinada", "combined") else 0.0)
        d50_2 += (dt_min / 18.0) * (d50_eq - d50_2)

        coprecip_ppm = max(C_ClO4_B3520 - C_sat_eff, 0.0) * 85.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 90.0)) ** 0.65
        w_ClO4_lav = PFD_MG_CLO4["wL_ClO4_lav_3000_ppm"] * (C_ClO4_B3520 / 4.00) * keff_mult
        w_Mg_lav = PFD_MG_CLO4["wL_Mg_lav_3000_ppm"] * (C_Mg_B3520 / 9.20)
        ppm_17 = moisture_17 * w_ClO4_lav + coprecip_ppm
        ppm_Mg_17 = moisture_17 * w_Mg_lav

        water_loss = 0.42 * max(F_303_actual - 27.14, 0.0)
        T2_perdida_KNO3 = KNO3_PERDIDA_T_t_h_por_C * max(T2 - 40.0, 0.0) + KNO3_PERDIDA_T_t_h_por_C * max(T1 - 55.0, 0.0) * (2.0 / 1.82)
        M_KNO3 += (dt_min / 15.0) * (66.5 + 10.0 * (1.0 - alpha_purge) - water_loss - T2_perdida_KNO3 - M_KNO3)

        rem_eff = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 90.0)) ** 0.85
        w_ClO4_ML_ref += (dt_min / 30.0) * (738.0 * (ppm_17 / 203.5) - w_ClO4_ML_ref)
        w_Mg_ML_ref += (dt_min / 30.0) * (2397.0 * (ppm_Mg_17 / 537.4) - w_Mg_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav) * (1.0 - rem_eff) * (w_ClO4_ML_ref / 738.0) ** 0.4 + 0.80 * coprecip_ppm
        ppm_Mg_final_6000 = ppm_Mg_17 * (1.0 - rem_eff) * (w_Mg_ML_ref / 2397.0) ** 0.4

        hist.append({
            "t_min": t,
            "C_Mg_B3520": C_Mg_B3520, "C_Mg_B3520_g_L": C_Mg_B3520 * k_gL_Mg,
            "C_ClO4_B3520": C_ClO4_B3520, "C_ClO4_B3520_g_L": C_ClO4_B3520,
            "C_sat_eff": C_sat_eff, "C_sat_KClO4_eff": C_sat_eff,
            "F_rec_t_h": F_rec_B3520, "alpha_purge_pct": 100.0 * alpha_purge,
            "M_KNO3_t_h": M_KNO3, "ppm_ClO4_final_6000": ppm_final_6000,
            "F_303_t_h": F_303_actual, "dT2_trim_C": dT2_trim, "T2_sp_C": T2_sp_nom + dT2_trim, "P_vac1_bar": P_vac1, "P_vac2_bar": P_vac2, "sat_ratio": C_ClO4_B3520 / max(C_sat_eff, 1e-6), "perdida_KNO3_T_t_h": T2_perdida_KNO3, "R_eff": R_eff, "T1_C": T1, "T2_C": T2, "coprecip_ppm": coprecip_ppm,
            "ppm_Mg_crudo_17": ppm_Mg_17, "ppm_Mg_final_6000": ppm_Mg_final_6000,
        })
    return pd.DataFrame(hist)


simulate_plant_npt3_mg_recycle = simular_recirculacion_y_magnesio

# ---------------------------------------------------------------------------------------------
# Identificacion de la ganancia ultima sobre la simulacion
def _hp(x, dt_min, win_min=15.0):
    k = max(int(win_min / dt_min), 1) | 1               # ventana impar
    xp = np.pad(np.asarray(x, float), k // 2, mode="edge")   # sin ceros en los bordes (evita falsos picos al inicio/final)
    return x - np.convolve(xp, np.ones(k) / k, mode="valid")


def amplitud_oscilacion(df, t0, t1, col="T2_C", win_min=15.0):
    """Desviacion estandar de la componente oscilatoria (sin tendencia lenta) de `col` entre t0 y t1 [min]."""
    dt = float(df.t_min.iloc[1] - df.t_min.iloc[0])
    m = ((df.t_min > t0) & (df.t_min < t1)).values
    return float(np.std(_hp(df[col].values, dt, win_min)[m]))


def razon_decaimiento(df, col="T2_C", t0=25.0, t1=260.0, umbral=0.02):
    """Razon de decaimiento = (2do pico)/(1er pico) de la oscilacion de `col` tras la perturbacion. nan si hay < 2 picos."""
    dt = float(df.t_min.iloc[1] - df.t_min.iloc[0])
    y = _hp(df[col].values, dt)
    t = df.t_min.values
    m = (t > t0) & (t < t1)
    yy = y[m]
    pk = [i for i in range(1, len(yy) - 1) if yy[i] > yy[i - 1] and yy[i] >= yy[i + 1] and yy[i] > umbral]
    return float(yy[pk[1]] / yy[pk[0]]) if len(pk) >= 2 else float("nan")


def _sim_P(Kc, **kw):
    return simulate_plant_npt3_calibrated(Kc_T_override=Kc, tauI_T_override=1e6, tauD_T_override=0.0,
                                          vac_limits_override=(-2.0, 0.0), duration_min=420.0, **kw)


def identificar_Kcu(lo=0.012, hi=0.013, n_iter=9):
    """Ganancia ultima Kcu [bar/C] y periodo Pu [min] del lazo de T sobre simulate_plant_npt3_calibrated.

    Accion P pura (tauI muy grande), valvula sin saturar. Kcu = menor Kc con que la oscilacion de T2 no decae
    (amplitud tardia / amplitud intermedia >= 1), por biseccion en [lo, hi] (el intervalo debe encerrar la transicion). Pu: promedio
    del periodo entre cruces por cero a Kc = 1.03 Kcu.
    """
    def crece(Kc):
        d = _sim_P(Kc)
        return amplitud_oscilacion(d, 300, 400) / max(amplitud_oscilacion(d, 110, 170), 1e-9) >= 1.0
    if crece(lo) or not crece(hi):
        raise ValueError("el intervalo [lo, hi] no encierra la transicion estable -> inestable")
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if crece(mid) else (mid, hi)
    Kcu = 0.5 * (lo + hi)
    d = _sim_P(1.03 * Kcu)
    x = _hp(d.T2_C.values, float(d.t_min.iloc[1] - d.t_min.iloc[0]))
    m = ((d.t_min > 250) & (d.t_min < 420)).values
    xs, tt = x[m], d.t_min.values[m]
    z = np.where((xs[:-1] < 0) & (xs[1:] >= 0))[0]
    return Kcu, float(np.mean(np.diff(tt[z])))

# ---------------------------------------------------------------------------------------------
# Analisis de estrategias, agua minima y temperatura sobre simular_recirculacion_y_magnesio
PPM_NORMA = 40.0            # ppm de ClO4- en producto refinado (KPI de SQM)

ESTRATEGIAS_EXCURSION = {
    "1. Sin reciclo (actual)": dict(modo_reciclo="sin_reciclo_actual", use_AIT3105=False),
    "2. Reciclo fijo": dict(modo_reciclo="reciclo_fijo", use_AIT3105=False),
    "3. AIT 3105 + AIC 3255": dict(modo_reciclo="reciclo_inteligente", use_AIT3105=True),
    "4. AIT + AIC + supervisor T2 (propuesta)": dict(modo_reciclo="reciclo_inteligente", use_AIT3105=True, trim_T2=True),
}


def resumen_simulacion(d):
    """Indicadores de una corrida de simular_recirculacion_y_magnesio."""
    sat_P1 = float((d.P_vac1_bar <= VAC["P1_lim"][0] + 1e-4).mean() * 100.0)
    return {"ppm_max": float(d.ppm_ClO4_final_6000.max()), "ppm_final": float(d.ppm_ClO4_final_6000.iloc[-1]),
            "agua_media_t_h": float(d.F_303_t_h.mean()), "agua_max_t_h": float(d.F_303_t_h.max()),
            "KNO3_medio_t_h": float(d.M_KNO3_t_h.mean()), "T1_max_C": float(d.T1_C.max()), "T2_max_C": float(d.T2_C.max()),
            "T1_final_C": float(d.T1_C.iloc[-1]), "T2_final_C": float(d.T2_C.iloc[-1]), "perdida_KNO3_por_T_final_t_h": float(d.perdida_KNO3_T_t_h.iloc[-1]),
            "dT2_max_C": float(d.dT2_trim_C.max()), "vacio_A3110_saturado_pct": sat_P1,
            "perdida_KNO3_por_T_max_t_h": float(d.perdida_KNO3_T_t_h.max()), "razon_ClO4_sat_max": float(d.sat_ratio.max()),
            "cumple_norma": bool(d.ppm_ClO4_final_6000.max() <= PPM_NORMA)}


def comparar_estrategias_excursion(F303_max=None, **kw):
    """Corre las estrategias de ESTRATEGIAS_EXCURSION bajo la excursion de ClO4- (+550 ppm, t = 40 min) y de Mg2+ (+35 %, t = 140 min)."""
    sims = {n: simular_recirculacion_y_magnesio(F303_max=F303_max, **{**a, **kw}) for n, a in ESTRATEGIAS_EXCURSION.items()}
    tabla = pd.DataFrame({n: resumen_simulacion(d) for n, d in sims.items()}).T
    return sims, tabla


def barrido_agua_optimo(caps=(12.0, 16.0, 20.0, 22.0, 24.0, 27.14, 30.0, 33.0, 36.0), margen=0.90, tol_KNO3=1.0, **kw):
    """Barre el tope de agua F_303 con la estrategia propuesta (AIT + AIC + supervisor de T2) y elige el MINIMO que cumple.

    Criterio (parametrizable): ppm_max <= margen * PPM_NORMA y KNO3 medio >= (mejor KNO3 del barrido - tol_KNO3). Entre los factibles se
    toma el de menor agua media. Todos los topes se limitan a F303_MAX_AMBIENTAL.
    """
    filas = []
    for cap in caps:
        d = simular_recirculacion_y_magnesio("reciclo_inteligente", use_AIT3105=True, trim_T2=True, F303_max=cap, **kw)
        filas.append({"tope_agua_t_h": min(cap, F303_MAX_AMBIENTAL), **resumen_simulacion(d)})
    df = pd.DataFrame(filas)
    mejor = df.KNO3_medio_t_h.max()
    df["factible"] = (df.ppm_max <= margen * PPM_NORMA) & (df.KNO3_medio_t_h >= mejor - tol_KNO3)
    fact = df[df.factible]
    opt = None if fact.empty else fact.sort_values("agua_media_t_h").iloc[0].to_dict()
    return df, opt


def analisis_temperatura(**kw):
    """Sobretemperatura de la alimentacion previa a la dilucion y cristalizador operando caliente, con sus correcciones."""
    def run(**a): return resumen_simulacion(simular_recirculacion_y_magnesio("reciclo_inteligente", use_AIT3105=True, **{**kw, **a}))
    filas = {}
    for dT in (0.0, 5.0, 10.0, 15.0):
        filas[f"Alimentacion +{dT:.0f} C (sin correccion)"] = run(T_in_extra=dT)
    for nombre, a in (("+10 C con prealimentacion de vacio", dict(ff_temperatura=True)),
                      ("+10 C con agua de dilucion a 35 C", dict(T_agua_dilucion=35.0)),
                      ("+10 C con agua a 35 C + prealimentacion", dict(T_agua_dilucion=35.0, ff_temperatura=True)),
                      ("+10 C con prealimentacion + supervisor T2", dict(ff_temperatura=True, trim_T2=True))):
        filas[nombre] = run(T_in_extra=10.0, **a)
    for T2 in (40.0, 42.0, 44.0, 46.0, 48.0):
        filas[f"Cristalizador 2 en {T2:.0f} C"] = run(T2_sp=T2)
    filas["Cristalizador 1 en 58 C"] = run(T1_sp=58.0)
    return pd.DataFrame(filas).T


def calor_enfriar_agua_dilucion_kW(F303_t_h=F303_NOM, T_inicial=60.0, T_final=35.0, cp_kJ_kgK=4.172):
    """Carga termica del intercambiador que enfria el agua de dilucion [kW]."""
    return F303_t_h * 1000.0 / 3600.0 * cp_kJ_kgK * (T_inicial - T_final)
