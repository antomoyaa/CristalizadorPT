"""Cristalizador al vacío (A3110) con reciclo de licor madre y purga: modelo dinámico y estacionario.

UNA sola fuente de ecuaciones: `rhs` es el modelo dinámico y el estado estacionario se obtiene resolviendo rhs = 0
(así, en lazo abierto y en el estacionario, el modelo no se mueve; ver tests/test_core.py).

Simplificaciones (declaradas):
  * C y CP son concentraciones por m3 de pulpa; el volumen de cristales (<3 %) se desprecia en los caudales.
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

import numpy as np
from scipy.optimize import root

from . import dpbe

NAMES = ["V", "T", "C", "CP", "mu0", "mu1", "mu2", "mu3", "Cm", "CPm", "Qd", "S", "Sm"]
IDX = {n: i for i, n in enumerate(NAMES)}
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

def simulate_plant_npt3_calibrated(
    use_AIT3105=True,
    use_AIC3255=True,
    tuning_mode="TL_damped",
    duration_min=300.0,
    dt_min=0.05,
    theta_AIT_min=10.0,
):
    """Simula el tren integrado Área 2000 (B2510) -> Área 3000 (A3110/A3210/B3520) -> Área 6000."""
    n_steps = int(duration_min / dt_min)
    w_ClO4_pre, T1, T2 = 3714.1, 55.0, 40.0
    v_T1, v_T2 = 0.0, 0.0
    xS1, xS2, d50_2 = 0.130, 0.200, 185.0
    C_ClO4_B3520, w_ClO4_ML_ref, F_303_actual = 4.00, 738.0, 27.14

    is_osc = tuning_mode == "ZN_oscillatory"
    is_hopf = tuning_mode == "Hopf_limit_cycle"

    sens_TIR3101 = DelayedSensor(1.20, 0.5, dt_min, T1)
    sens_TIR3201 = DelayedSensor(1.20, 0.5, dt_min, T2)
    sens_DIR3101 = DelayedSensor(0.80, 0.4, dt_min, xS1, resolution=0.0002)
    sens_DIR3201 = DelayedSensor(0.80, 0.4, dt_min, xS2, resolution=0.0002)
    sens_AIT3105 = DelayedSensor(theta_AIT_min, 2.0, dt_min, w_ClO4_pre)
    sens_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520)
    sens_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual)

    if is_hopf:
        Kc_T, tauI_T, tauD_T = 0.0152, 3.2, 0.0
        Kc_D, tauI_D = -2100.0, 8.5
    elif is_osc:
        Kc_T, tauI_T, tauD_T = 0.0095, 3.8, 0.25
        Kc_D, tauI_D = -1350.0, 11.0
    else:
        Kc_T, tauI_T, tauD_T = 0.0052, 9.5, 0.60
        Kc_D, tauI_D = -650.0, 28.0

    pid_TP1 = CommercialPID(Kc_T, tauI_T, tauD_T, 55.0, -0.932, -0.970, -0.930)
    pid_TP2 = CommercialPID(Kc_T, tauI_T, tauD_T, 40.0, -0.962, -0.970, -0.930)
    pid_DIR1 = CommercialPID(Kc_D, tauI_D, 0.0, 0.130, 390.23, 250.0, 520.0)
    pid_DIR2 = CommercialPID(Kc_D, tauI_D, 0.0, 0.200, 390.23, 250.0, 520.0)
    pid_AIC3255 = CommercialPID(-0.48 if is_osc else -0.28, 18.0 if is_osc else 45.0, 0.0, 4.00, 0.22, 0.08, 0.65)
    pid_FIRC3103 = CommercialPID(Kc=1.25, tau_I=0.45, tau_D=0.0, sp=27.14, u_bias=27.14, u_min=12.0, u_max=42.0)

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208 = 356.85 * (1.08 if 25.0 <= t < 140.0 else 1.0)
        T_in_1 = 70.0 + (3.5 if t >= 25.0 else 0.0)
        w_ClO4_208 = 3892.0 + (650.0 if t >= 80.0 else 0.0)
        F_234 = 21.33
        F_pre = F_208 + F_234

        w_ClO4_pre_target = (F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref) / F_pre
        w_ClO4_pre += (dt_min / 8.0) * (w_ClO4_pre_target - w_ClO4_pre)

        R_0 = 0.07176
        if use_AIT3105:
            ait3105_pv = sens_AIT3105.read(w_ClO4_pre)
            R_eff = R_0 + 5.2e-4 * (ait3105_pv - 3714.1)
        else:
            R_eff = R_0

        F_303_sp = np.clip(F_pre * R_eff, 12.0, 42.0)
        pid_FIRC3103.sp = F_303_sp
        valve_cmd = pid_FIRC3103.step(sens_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (valve_cmd - F_303_actual)

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1

        P_vac1 = pid_TP1.step(sens_TIR3101.read(T1), dt_min)
        P_vac2 = pid_TP2.step(sens_TIR3201.read(T2), dt_min)
        F_6 = pid_DIR1.step(sens_DIR3101.read(xS1), dt_min)
        F_11 = pid_DIR2.step(sens_DIR3201.read(xS2), dt_min)
        alpha_purge = pid_AIC3255.step(sens_AIC3255.read(C_ClO4_B3520), dt_min) if use_AIC3255 else 0.22

        T1_eq = 55.0 + 375.0 * (P_vac1 - (-0.932)) + 0.45 * (T_in_1 - 70.0)
        T2_eq = 40.0 + 375.0 * (P_vac2 - (-0.962)) + 0.40 * (T1 - 55.0) + (1.1 if t >= 25.0 else 0.0)
        tau1_T, tau2_T = 4.5, 1.8
        v_T1 += ((T1_eq - T1 - (tau1_T + tau2_T) * v_T1) / (tau1_T * tau2_T)) * dt_min
        v_T2 += ((T2_eq - T2 - (tau1_T + tau2_T) * v_T2) / (tau1_T * tau2_T)) * dt_min
        T1 += v_T1 * dt_min
        T2 += v_T2 * dt_min

        dil_effect = -0.0004 * (F_303_actual - 27.14)
        xS1_eq = 0.130 - 4.5e-4 * (F_6 - 390.23) - 0.003 * (T1 - 55.0) + dil_effect + (0.009 if t >= 25.0 else 0.0)
        xS2_eq = 0.200 + 0.75 * (xS1 - 0.130) - 4.5e-4 * (F_11 - 390.23) - 0.0035 * (T2 - 40.0) + (0.008 if t >= 25.0 else 0.0)
        xS1 += (dt_min / 9.0) * (xS1_eq - xS1)
        xS2 += (dt_min / 9.0) * (xS2_eq - xS2)

        sigma_burst = max(55.0 - T1, 0.0) * 2.5 + max(40.0 - T2, 0.0) * 3.0
        d50_2_eq = 185.0 + 180.0 * (xS2 - 0.200) - 4.5 * sigma_burst + 0.35 * (390.23 - F_11)
        d50_2 += (dt_min / 18.0) * (d50_2_eq - d50_2)

        C_ClO4_in_eq = 4.00 * (w_ClO4_1 / 3468.0) ** 1.45 * (0.22 / max(alpha_purge, 0.05)) ** 0.50
        C_ClO4_B3520 += (dt_min / 24.0) * (C_ClO4_in_eq - C_ClO4_B3520)

        coprecip_ppm = max(C_ClO4_B3520 - 4.50, 0.0) * 82.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 90.0)) ** 0.65
        w_ClO4_lav_3000 = 4070.0 * (C_ClO4_B3520 / 4.00)
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 90.0)) ** 0.85
        w_ClO4_ML_ref_eq = 738.0 * (ppm_ClO4_17 / 203.5)
        w_ClO4_ML_ref += (dt_min / 30.0) * (w_ClO4_ML_ref_eq - w_ClO4_ML_ref)
        ppm_superficial_ref = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.4
        ppm_ClO4_final_6000 = ppm_superficial_ref + 0.80 * coprecip_ppm

        hist.append({
            "t_min": t, "T1_C": T1, "T2_C": T2, "P_vac1_bar": P_vac1, "P_vac2_bar": P_vac2,
            "xS1": xS1, "xS2": xS2, "F_6_t_h": F_6, "F_11_t_h": F_11,
            "F_pre_t_h": F_pre, "w_ClO4_pre_ppm": w_ClO4_pre, "w_ClO4_1_ppm": w_ClO4_1,
            "R_eff": R_eff, "F_303_sp_t_h": F_303_sp, "F_303_actual_t_h": F_303_actual,
            "C_ClO4_B3520_g_L": C_ClO4_B3520, "alpha_purge_pct": 100.0 * alpha_purge,
            "d50_2_um": d50_2, "ppm_ClO4_crudo_17": ppm_ClO4_17,
            "ppm_ClO4_final_6000": ppm_ClO4_final_6000,
        })
    return pd.DataFrame(hist)

# === FUNCIONES CENTRALIZADAS: MAGNESIO (Mg2+), RECIRCULACIÓN B3520, FLUJO ÓPTIMO Y ROBUSTEZ ===

def simular_recirculacion_y_magnesio(
    modo_reciclo="sin_reciclo_actual",  # "sin_reciclo_actual", "reciclo_fijo", "reciclo_inteligente"
    use_AIT3105=False,
    hipotesis_mg="combinada",           # "salting_out", "cinetica_particion", "combinada"
    duration_min=360.0,
    dt_min=0.05,
):
    """Evalúa la recirculación B3520 -> A3110 y las 3 hipótesis de interacción Mg2+ / ClO4-."""
    n_steps = int(duration_min / dt_min)
    w_ClO4_pre, w_Mg_pre = 3714.1, 8500.0
    T1, T2 = 55.0, 40.0
    P_vac1, P_vac2 = -0.932, -0.962
    d50_2 = 185.0
    C_ClO4_B3520, C_Mg_B3520 = 3.85, 9.20
    w_ClO4_ML_ref = 738.0
    F_303_actual = 27.14
    M_KNO3 = 66.5 if modo_reciclo in ("sin_reciclo_actual", "no_recycle_actual") else 74.0

    s_TIR1 = DelayedSensor(1.0, 0.5, dt_min, T1)
    s_TIR2 = DelayedSensor(1.0, 0.5, dt_min, T2)
    s_AIT3105 = DelayedSensor(10.0, 2.0, dt_min, w_ClO4_pre)
    s_AIC3255 = DelayedSensor(10.0, 2.0, dt_min, C_ClO4_B3520)
    s_FIRC3103 = DelayedSensor(0.10, 0.08, dt_min, F_303_actual)

    pid_TIR1 = CommercialPID(0.0055, 9.0, 0.4, 55.0, -0.932, -0.970, -0.930)
    pid_TIR2 = CommercialPID(0.0055, 9.0, 0.4, 40.0, -0.962, -0.970, -0.930)
    pid_PIR1 = CommercialPID(1.5, 0.3, 0.0, -0.932, -0.932, -0.970, -0.930)
    pid_PIR2 = CommercialPID(1.5, 0.3, 0.0, -0.962, -0.962, -0.970, -0.930)
    pid_FIRC3103 = CommercialPID(1.25, 0.45, 0.0, 27.14, 27.14, 12.0, 42.0)
    pid_AIC3255 = CommercialPID(-0.38, 22.0, 0.0, 3.95, 0.25, 0.12, 1.00)

    hist = []
    for k in range(n_steps):
        t = k * dt_min
        F_208, F_234 = 356.85, 21.33
        w_ClO4_208 = 3892.0 + (550.0 if t >= 40.0 else 0.0)
        w_Mg_208 = 8500.0 * (1.35 if t >= 140.0 else 1.0)

        if modo_reciclo in ("sin_reciclo_actual", "no_recycle_actual"):
            alpha_purge, F_rec_B3520 = 1.00, 0.0
        elif modo_reciclo in ("reciclo_fijo", "fixed_recycle"):
            alpha_purge = 0.25
            F_rec_B3520 = 85.0 * (1.0 - alpha_purge)
        else:
            delta_Mg_sp = max((C_Mg_B3520 - 9.20) / 9.20, 0.0)
            pid_AIC3255.sp = 3.95 - 0.55 * delta_Mg_sp
            alpha_purge = pid_AIC3255.step(s_AIC3255.read(C_ClO4_B3520), dt_min)
            F_rec_B3520 = 85.0 * (1.0 - alpha_purge)

        F_pre = F_208 + F_234 + F_rec_B3520
        w_ClO4_pre_eq = (F_208 * w_ClO4_208 + F_234 * w_ClO4_ML_ref + F_rec_B3520 * (C_ClO4_B3520 * 820.0)) / F_pre
        w_Mg_pre_eq = (F_208 * w_Mg_208 + F_rec_B3520 * (C_Mg_B3520 * 820.0)) / F_pre
        w_ClO4_pre += (dt_min / 8.0) * (w_ClO4_pre_eq - w_ClO4_pre)
        w_Mg_pre += (dt_min / 8.0) * (w_Mg_pre_eq - w_Mg_pre)

        R_0 = 0.07176
        if use_AIT3105:
            ait_pv = s_AIT3105.read(w_ClO4_pre)
            delta_Mg_trim = max((w_Mg_pre - 8500.0) / 8500.0, 0.0)
            R_eff = R_0 + 5.2e-4 * (ait_pv - 3714.1) + 0.018 * delta_Mg_trim
        else:
            R_eff = R_0

        F_303_sp = np.clip(F_pre * R_eff, 12.0, 42.0)
        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre) / F_1

        pid_PIR1.sp = pid_TIR1.step(s_TIR1.read(T1), dt_min)
        pid_PIR2.sp = pid_TIR2.step(s_TIR2.read(T2), dt_min)
        P_vac1 += (dt_min / 0.25) * (pid_PIR1.step(P_vac1, dt_min) - P_vac1)
        P_vac2 += (dt_min / 0.25) * (pid_PIR2.step(P_vac2, dt_min) - P_vac2)
        T1 += (dt_min / 5.5) * (55.0 + 375.0 * (P_vac1 - (-0.932)) - T1)
        T2 += (dt_min / 5.5) * (40.0 + 375.0 * (P_vac2 - (-0.962)) + 0.35 * (T1 - 55.0) - T2)

        rec_factor = (1.0 / max(alpha_purge, 0.12)) ** 0.28
        C_ClO4_eq = 3.85 * (w_ClO4_1 / 3468.0) ** 1.35 * rec_factor
        C_Mg_eq = 9.20 * (w_Mg_1 / 7930.0) * rec_factor
        C_ClO4_B3520 += (dt_min / 22.0) * (C_ClO4_eq - C_ClO4_B3520)
        C_Mg_B3520 += (dt_min / 22.0) * (C_Mg_eq - C_Mg_B3520)

        d_Mg = max((C_Mg_B3520 - 9.20) / 9.20, 0.0)
        salting = 0.85 * d_Mg if hipotesis_mg in ("salting_out", "combinada", "combined") else 0.0
        C_sat_eff = 4.50 - salting

        keff_mult = (1.0 + 0.45 * d_Mg) if hipotesis_mg in ("cinetica_particion", "kinetic_partition", "combinada", "combined") else 1.0
        d50_eq = 185.0 - (28.0 * d_Mg if hipotesis_mg in ("cinetica_particion", "kinetic_partition", "combinada", "combined") else 0.0)
        d50_2 += (dt_min / 18.0) * (d50_eq - d50_2)

        coprecip_ppm = max(C_ClO4_B3520 - C_sat_eff, 0.0) * 85.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 90.0)) ** 0.65
        w_ClO4_lav = 4070.0 * (C_ClO4_B3520 / 4.00) * keff_mult
        ppm_17 = moisture_17 * w_ClO4_lav + coprecip_ppm

        water_loss = 0.42 * max(F_303_actual - 27.14, 0.0)
        M_KNO3 += (dt_min / 15.0) * (66.5 + 10.0 * (1.0 - alpha_purge) - water_loss - M_KNO3)

        rem_eff = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 90.0)) ** 0.85
        w_ClO4_ML_ref += (dt_min / 30.0) * (738.0 * (ppm_17 / 203.5) - w_ClO4_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav) * (1.0 - rem_eff) * (w_ClO4_ML_ref / 738.0) ** 0.4 + 0.80 * coprecip_ppm

        hist.append({
            "t_min": t,
            "C_Mg_B3520": C_Mg_B3520, "C_Mg_B3520_g_L": C_Mg_B3520,
            "C_ClO4_B3520": C_ClO4_B3520, "C_ClO4_B3520_g_L": C_ClO4_B3520,
            "C_sat_eff": C_sat_eff, "C_sat_KClO4_eff": C_sat_eff,
            "F_rec_t_h": F_rec_B3520, "alpha_purge_pct": 100.0 * alpha_purge,
            "M_KNO3_t_h": M_KNO3, "ppm_ClO4_final_6000": ppm_final_6000,
        })
    return pd.DataFrame(hist)


# Alias de compatibilidad por si alguna celda llama al nombre en inglés
simulate_plant_npt3_mg_recycle = simular_recirculacion_y_magnesio


def simular_planta_agua_dilucion_magnesio(
    strategy="AIT3105_and_Mg_trim",   # "ratio_FFIC3102", "high_water_manual", "AIT3105_only", "AIT3105_and_Mg_trim"
    recycle_mode="no_recycle_actual", # "no_recycle_actual", "fixed_recycle", "smart_AIC3255"
    include_noise=True,
    theta_AIT_min=10.0,
    duration_min=360.0,
    dt_min=0.05,
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
    P_vac1, P_vac2 = -0.932, -0.962
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

    pid_TIR1 = CommercialPID(0.0055, 9.0, 0.4, 55.0, -0.932, -0.970, -0.930, beta=0.5)
    pid_TIR2 = CommercialPID(0.0055, 9.0, 0.4, 40.0, -0.962, -0.970, -0.930, beta=0.5)
    pid_PIR1 = CommercialPID(1.5, 0.3, 0.0, -0.932, -0.932, -0.970, -0.930)
    pid_PIR2 = CommercialPID(1.5, 0.3, 0.0, -0.962, -0.962, -0.970, -0.930)
    pid_DIR1 = CommercialPID(-650.0, 28.0, 0.0, 0.130, 390.23, 250.0, 520.0)
    pid_DIR2 = CommercialPID(-650.0, 28.0, 0.0, 0.200, 390.23, 250.0, 520.0)
    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, 50.0)
    pid_AIC3255 = CommercialPID(-0.42, 20.0, 0.0, 3.92, 0.25, 0.12, 1.00)

    R_0 = 27.14 / 378.18
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
            R_eff = R_0 + 5.2e-5 * max(ait_cl_pv - 3714.1, 0.0)
        else:
            d_mg_in = max((ait_mg_pv - 8200.0) / 8200.0, 0.0)
            R_eff = R_0 + 5.2e-5 * max(ait_cl_pv - 3714.1, 0.0) + 0.065 * d_mg_in

        F_303_sp = np.clip(F_pre * R_eff, 12.0, 50.0)
        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (valve_gain * v_cmd - F_303_actual)

        F_1 = F_pre + F_303_actual
        w_ClO4_1 = (F_pre * w_ClO4_pre + F_303_actual * 36.0) / F_1
        w_Mg_1 = (F_pre * w_Mg_pre + F_303_actual * 10.0) / F_1

        pid_PIR1.sp = pid_TIR1.step(s_TIR1.read(T1), dt_min)
        pid_PIR2.sp = pid_TIR2.step(s_TIR2.read(T2), dt_min)
        P_vac1 += (dt_min / 0.25) * (pid_PIR1.step(P_vac1, dt_min) - P_vac1)
        P_vac2 += (dt_min / 0.25) * (pid_PIR2.step(P_vac2, dt_min) - P_vac2)

        T1 += (dt_min / 5.5) * (55.0 + 375.0 * (P_vac1 - (-0.932)) + 0.40 * (T_in_1 - 70.0) - 0.025 * (F_303_actual - 27.14) - T1)
        T2 += (dt_min / 5.5) * (40.0 + 375.0 * (P_vac2 - (-0.962)) + 0.35 * (T1 - 55.0) - T2)

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
        C_sat_KClO4_40C = 4.38 - 1.15 * d_Mg_proc + 0.042 * water_added

        d50_1 += (dt_min / 14.0) * (145.0 - 25.0 * max((C_Mg_A3110 - 8.10) / 8.10, 0.0) + 0.55 * water_added - d50_1)
        d50_2 += (dt_min / 16.0) * (185.0 * (d50_1 / 145.0) - 28.0 * d_Mg_proc + 0.65 * water_added - d50_2)

        keff_mult = 1.0 + 0.38 * d_Mg_proc
        coprecip_ppm = max(C_ClO4_B3520 - C_sat_KClO4_40C, 0.0) * 85.0
        moisture_17 = 0.050 * (185.0 / max(d50_2, 100.0)) ** 0.65 * (1.0 + 0.12 * d_Mg_proc)
        w_ClO4_lav_3000 = 4070.0 * (C_ClO4_B3520 / 4.00) * keff_mult
        ppm_ClO4_17 = moisture_17 * w_ClO4_lav_3000 + coprecip_ppm

        excess_water = max(F_303_actual - (F_pre * R_0), 0.0)
        M_KNO3_eq = (67.5 + 8.8 * (1.0 - alpha_purge)) * (F_208 / 356.85) - 0.36 * excess_water
        M_KNO3_prod += (dt_min / 15.0) * (M_KNO3_eq - M_KNO3_prod)

        rem_eff_6000 = 1.0 - (1.0 - 0.899) * (185.0 / max(d50_2, 100.0)) ** 0.85
        w_ClO4_ML_ref += (dt_min / 24.0) * (738.0 * (ppm_ClO4_17 / 203.5) - w_ClO4_ML_ref)
        ppm_final_6000 = (moisture_17 * w_ClO4_lav_3000) * (1.0 - rem_eff_6000) * (w_ClO4_ML_ref / 738.0) ** 0.38 + 0.82 * coprecip_ppm

        hist.append({
            "t_min": t, "F_303": F_303_actual, "w_Mg_1": w_Mg_1,
            "C_Mg_A3110": C_Mg_A3110, "C_Mg_A3210": C_Mg_A3210, "C_Mg_B3520": C_Mg_B3520,
            "C_ClO4_B3520": C_ClO4_B3520, "C_sat_KClO4_40C": C_sat_KClO4_40C, "d50_2": d50_2,
            "M_KNO3_t_h": M_KNO3_prod, "ppm_ClO4_final_6000": ppm_final_6000,
        })
    return pd.DataFrame(hist)


def simular_AIT_AIC_por_flujo_agua(
    F_303_target=36.8,
    modo_control="AIT_and_AIC",
    delta_ClO4_ssk=650.0,
    include_noise=True,
    duration_min=320.0,
    dt_min=0.05,
):
    """Evalúa los lazos AIT 3105 y AIC 3255 en función del flujo objetivo de agua F_303 (óptimo 36.8 t/h)."""
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

    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, 52.0)
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
            F_303_sp = np.clip(27.14 + (F_303_target - 27.14) * factor_ait, 12.0, 52.0)

        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)

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

    pid_FIRC3103 = CommercialPID(1.30, 0.40, 0.0, 27.14, 27.14, 12.0, 50.0)
    pid_AIC3255 = CommercialPID(-0.48, 16.0, 0.0, 4.00, 0.25, 0.12, 0.85)
    R_0 = 27.14 / 378.18

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
            R_eff = R_0 + 5.2e-5 * max(ait_pv - 3714.1, 0.0) + 0.055 * d_mg_in + aic_trim
            F_303_sp = np.clip(F_pre * R_eff, 12.0, 50.0)
        else:
            F_303_sp = 27.14

        pid_FIRC3103.sp = F_303_sp
        v_cmd = pid_FIRC3103.step(s_FIRC3103.read(F_303_actual), dt_min)
        F_303_actual += (dt_min / 0.133) * (v_cmd - F_303_actual)

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