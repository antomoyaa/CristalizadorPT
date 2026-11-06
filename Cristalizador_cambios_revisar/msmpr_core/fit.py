"""Ajuste ESTACIONARIO rápido de la cinética por mínimos cuadrados ponderados (guía §3.1).

NO es el calibrador de referencia: la calibración de A3110/A3210 (dinámica, validación temporal, un ajuste por
equipo) vive en calibracion_do_crystal.ipynb. Este módulo sirve para chequeos rápidos con datos estacionarios y
sigma medido, sobre el modelo analítico de msmpr_core (no sobre el DPBE del notebook de Pareto).

ajustar_cineticas(df) estima theta = [ln B_ref, b, ln G_ref, g] con la parametrización centrada en sigma_ref
(menos correlacionada) y exige medir sigma en el licor madre; sin sigma los parámetros no son identificables.

Columnas de df: tau_h, f_dil, T, d50_um, N_kg, Y, sigma  y  s_d50_um, s_N_kg, s_Y, s_sigma (error de cada medición;
s_d50_um y s_N_kg relativos si < 1.0, absolutos si no; s_Y y s_sigma absolutos).
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.optimize import least_squares

from .model import msmpr, DEFAULT_SIGMA_REF

LABELS = ["ln B_ref", "b", "ln G_ref", "g"]


def _pn(theta, sigma_ref):
    return {"B_ref": float(np.exp(theta[0])), "b": theta[1], "G_ref": float(np.exp(theta[2])), "g": theta[3],
            "sigma_ref": sigma_ref}


def _rel(s, y):
    return s if s < 1.0 else s / y


def ajustar_cineticas(df, p=None, theta0=None, sigma_ref=DEFAULT_SIGMA_REF, bounds=None) -> dict:
    theta0 = np.array([np.log(1e10 * sigma_ref ** 2.0), 2.0, np.log(2e-7 * sigma_ref), 1.0]) if theta0 is None else np.array(theta0, float)
    lb = np.array([-60.0, 0.5, -60.0, 0.5]) if bounds is None else np.array(bounds[0], float)
    ub = np.array([60.0, 5.0, 0.0, 5.0]) if bounds is None else np.array(bounds[1], float)

    def resid(theta):
        Pn = _pn(theta, sigma_ref)
        r = []
        for row in df.itertuples(index=False):
            m = msmpr(row.tau_h, row.f_dil, row.T, Pn, p)
            r += [np.log(m["d50_um"] / row.d50_um) / _rel(row.s_d50_um, row.d50_um),
                  np.log(m["N_kg"] / row.N_kg) / _rel(row.s_N_kg, row.N_kg),
                  (m["Y"] - row.Y) / row.s_Y,
                  (m["sigma"] - row.sigma) / row.s_sigma]
        return np.array(r)

    sol = least_squares(resid, theta0, bounds=(lb, ub), x_scale=[5.0, 1.0, 5.0, 1.0])
    n, k = len(sol.fun), len(sol.x)
    dof = max(n - k, 1)
    chi2_red = float(2.0 * sol.cost / dof)
    J = sol.jac
    cov = np.linalg.pinv(J.T @ J) * max(chi2_red, 1.0)
    se = np.sqrt(np.diag(cov))
    tcrit = stats.t.ppf(0.975, dof)
    corr = cov / np.outer(se, se)
    return {
        "theta": dict(zip(LABELS, sol.x)), "theta_vec": sol.x,
        "Pn": _pn(sol.x, sigma_ref), "se": dict(zip(LABELS, se)),
        "IC95": {l: (float(v - tcrit * s), float(v + tcrit * s)) for l, v, s in zip(LABELS, sol.x, se)},
        "corr": corr, "chi2_red": chi2_red, "exito": bool(sol.success), "n_resid": n,
    }
