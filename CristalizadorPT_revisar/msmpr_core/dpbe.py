"""Cinética y CSD estacionaria de un MSMPR con G y B independientes del tamaño.

Cierre del balance de soluto (guía, sección 6.4):  C_in - C = M_T = 6 rho_c kv B G^3 tau^4
con n(L) = (B/G) exp(-L/(G tau)); la CSD en masa es Gamma(4) en L/(G tau):
  d10 : d50 : d90 = 1.745 : 3.672 : 6.681  (x G tau),  (d90-d10)/d50 = 1.345 siempre.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.stats import gamma as _gamma

from .params import R_GAS

T_REF_K = 313.15


# ------------------------------------------------------------------ termodinámica y cinética
def cstar(T_C: float, p: dict) -> float:
    """Solubilidad (kg/m3). Exponencial alrededor de T_ref; reemplazar por la curva de la salmuera real."""
    return p["Cstar_ref_kg_m3"] * math.exp(p["solubility_slope_K_inv"] * (T_C - p["reference_temperature_C"]))


def kinetics(C: float, T_C: float, p: dict):
    """Devuelve (sigma, G [m/s], B [#/(m3 s)]). sigma = C/C* - 1 (0 si subsaturado)."""
    sigma = max(C / cstar(T_C, p) - 1.0, 0.0)
    arr = math.exp(-p.get("Eg_J_mol", 0.0) / R_GAS * (1.0 / (T_C + 273.15) - 1.0 / T_REF_K))
    G = p["kG_m_s"] * arr * sigma ** p["g"]
    B = p["kB_m3_s_inv"] * sigma ** p["b"]
    return sigma, G, B


def from_ref(B_ref: float, b: float, G_ref: float, g: float, sigma_ref: float = 0.04, Eg: float = 0.0) -> dict:
    """Parametrización centrada (menos correlacionada): B = B_ref (s/s_ref)^b, G = G_ref (s/s_ref)^g."""
    return {"kB_m3_s_inv": B_ref / sigma_ref ** b, "b": b, "kG_m_s": G_ref / sigma_ref ** g, "g": g, "Eg_J_mol": Eg}


def moments_steady(B: float, G: float, tau: float):
    """mu_k = k! B G^(k-1) tau^(k+1)*G  ->  mu0=B tau, mu1=B G tau^2, mu2=2 B G^2 tau^3, mu3=6 B G^3 tau^4."""
    return np.array([B * tau, B * G * tau ** 2, 2.0 * B * G ** 2 * tau ** 3, 6.0 * B * G ** 3 * tau ** 4])


def crystal_mass(mu3: float, p: dict) -> float:
    """M_T [kg cristal / m3 pulpa]."""
    return p["rho_crystal_kg_m3"] * p["kv"] * mu3


# ------------------------------------------------------------------ percentiles de la CSD en masa
_G4 = _gamma(4.0)
Q10, Q50, Q90 = (float(_G4.ppf(q)) for q in (0.1, 0.5, 0.9))   # 1.745, 3.672, 6.681
SPAN_IDEAL = (Q90 - Q10) / Q50                                  # 1.345


def sizes_from_Gtau(G: float, tau: float) -> dict:
    """d10, d50, d90, d32 [m] de un MSMPR ideal."""
    x = G * tau
    return {"d10": Q10 * x, "d50": Q50 * x, "d90": Q90 * x, "d32": 3.0 * x}


# ------------------------------------------------------------------ estacionario analítico
def stationary_msmpr(C_in: float, tau: float, T_C: float, p: dict) -> dict:
    """Resuelve C_in - C = M_T(sigma) con brentq (monótona en sigma). C_in ya incluye la concentración por evaporación."""
    Cs = cstar(T_C, p)
    if C_in <= Cs:
        return {"sigma": 0.0, "C": C_in, "G": 0.0, "B": 0.0, "M_T": 0.0, "mu": np.zeros(4),
                "d10": 0.0, "d50": 0.0, "d90": 0.0, "d32": 0.0, "tau": tau}

    def f(sigma):
        C = Cs * (1.0 + sigma)
        _, G, B = kinetics(C, T_C, p)
        return (C_in - C) - crystal_mass(moments_steady(B, G, tau)[3], p)

    s_hi = C_in / Cs - 1.0
    sigma = brentq(f, 1e-12, s_hi, xtol=1e-14, rtol=1e-13)
    C = Cs * (1.0 + sigma)
    _, G, B = kinetics(C, T_C, p)
    mu = moments_steady(B, G, tau)
    out = {"sigma": sigma, "C": C, "G": G, "B": B, "mu": mu, "M_T": crystal_mass(mu[3], p), "tau": tau}
    out.update(sizes_from_Gtau(G, tau))
    return out


# ------------------------------------------------------------------ versión en malla (para extender a G(L), aglomeración...)
def make_grid(Lmax: float, bins: int) -> dict:
    edges = np.linspace(0.0, Lmax, bins + 1)
    dl = edges[1] - edges[0]
    q3 = (edges[1:] ** 4 - edges[:-1] ** 4) / (4.0 * dl)      # integral de L^3 sobre cada celda / dl
    return {"edges": edges, "dl": dl, "centers": 0.5 * (edges[1:] + edges[:-1]), "q3": q3}


def stationary_grid(C_in: float, tau: float, T_C: float, p: dict, Lmax: float, bins: int = 1000) -> dict:
    """Upwind de 1er orden: N_i = B/(a+1/tau) (a/(a+1/tau))^i con a = G/dl. N_i en #/m3 por celda."""
    grid = make_grid(Lmax, bins)
    Cs = cstar(T_C, p)
    idx = np.arange(bins)

    def dist(C):
        _, G, B = kinetics(C, T_C, p)
        a = G / grid["dl"]
        return (B / (a + 1.0 / tau)) * (a / (a + 1.0 / tau)) ** idx

    def f(C):
        mass = p["rho_crystal_kg_m3"] * p["kv"] * float(grid["q3"] @ dist(C))
        return (C_in - C) - mass

    C = brentq(f, Cs * (1 + 1e-12), C_in, xtol=1e-12, rtol=1e-13)
    N = dist(C)
    w = grid["q3"] * N
    cum = np.cumsum(w) / w.sum()
    pct = {name: float(np.interp(q, cum, grid["edges"][1:])) for name, q in (("d10", .1), ("d50", .5), ("d90", .9))}
    return {"C": C, "N": N, "grid": grid, "M_T": float(p["rho_crystal_kg_m3"] * p["kv"] * w.sum()), **pct}
