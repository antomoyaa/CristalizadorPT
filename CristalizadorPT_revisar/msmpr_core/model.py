"""API de la guía de ajuste (sección 2): Y = f(P).

    msmpr(tau_h, f_dil, T, Pn) -> dict(d50 [um], span, N_kg [#/kg agua], Y, sigma, prod_tph, ...)

Entradas P: tau_h [h], f_dil [kg agua añadida / kg alimentación] (0-0.2), T [degC].
Parámetros cinéticos Pn: {'B_ref','b','G_ref','g'} (+ 'Eg', 'sigma_ref'); B_ref en #/(m3 s), G_ref en m/s.
Sin reciclo ni evaporación (modelo de la guía). Con reciclo y vacío use msmpr_core.plant.

Base de concentración: kg KNO3 / m3 de pulpa. La recuperación se calcula de DOS formas que deben coincidir
(pizarra del 15-10):  Y_liq = 1 - C/C_in   e   Y_sol = M_T/C_in   -> 'Y_cierre' = Y_sol - Y_liq ~ 1e-12.
"""
from __future__ import annotations

from copy import deepcopy

import numpy as np

from . import dpbe
from .params import make_params

DEFAULT_SIGMA_REF = 0.04


def _kin_params(p: dict, Pn: dict | None) -> dict:
    q = deepcopy(p)
    if Pn:
        q.update(dpbe.from_ref(Pn["B_ref"], Pn["b"], Pn["G_ref"], Pn["g"],
                               Pn.get("sigma_ref", DEFAULT_SIGMA_REF), Pn.get("Eg", 0.0)))
    return q


def msmpr(tau_h: float, f_dil: float, T: float, Pn: dict | None = None, p: dict | None = None) -> dict:
    p = make_params() if p is None else p
    q = _kin_params(p, Pn)
    C_in = q["feed_stock_kg_m3"]["KNO3"] / (1.0 + f_dil)        # dilución con agua: baja C_in
    tau = tau_h * 3600.0
    s = dpbe.stationary_msmpr(C_in, tau, T, q)
    Y_sol = s["M_T"] / C_in
    Y_liq = 1.0 - s["C"] / C_in
    water_per_m3 = q["rho_liquid_kg_m3"] * q["water_mass_fraction_liquor"]
    d50 = 1e6 * s["d50"]
    span = (s["d90"] - s["d10"]) / s["d50"] if s["d50"] > 0 else np.nan
    return {
        "d50_um": d50, "span": span, "N_kg": s["mu"][0] / water_per_m3, "Y": Y_sol,
        "Y_liq": Y_liq, "Y_sol": Y_sol, "Y_cierre": Y_sol - Y_liq, "sigma": s["sigma"],
        "prod_tph": q.get("F_kno3_tph", 150.0) * Y_sol,
        "G_m_s": s["G"], "B_m3_s": s["B"], "M_T_kg_m3": s["M_T"],
        "d10_um": 1e6 * s["d10"], "d90_um": 1e6 * s["d90"],
    }


def calibrate_Fin(p: dict, target_prod_tph: float = 70.0, base=(1.5, 0.10, 40.0), Pn: dict | None = None) -> float:
    """Ajusta F_kno3_tph para que el caso base produzca `target_prod_tph` (guía §7.1, punto 4)."""
    Y = msmpr(base[0], base[1], base[2], Pn, p)["Y"]
    return target_prod_tph / Y
