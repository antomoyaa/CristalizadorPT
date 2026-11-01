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
        self.d_nom = {"Q_stock": (1 - d) * Q, "T_feed": p["feed_temperature_C"],
                      "C_stock": p["feed_stock_kg_m3"]["KNO3"], "CP_stock": p["feed_stock_kg_m3"]["ClO4"],
                      "S_stock": p["feed_stock_kg_m3"]["OTHER"],
                      "T_dil": p["dilution_water_temperature_C"]}

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
            "recycle_m3_h": 3600.0 * (1 - a) * Qout,
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
