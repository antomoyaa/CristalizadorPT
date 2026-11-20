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



# ------------------------------------------------------------------ termodinámica y cinética
def cstar(T_C: float, p: dict) -> float:
    """Solubilidad (kg/m3). Exponencial alrededor de T_ref; reemplazar por la curva de la salmuera real."""
    return p["Cstar_ref_kg_m3"] * math.exp(p["solubility_slope_K_inv"] * (T_C - p["reference_temperature_C"]))


def kinetics(C: float, T_C: float, p: dict):
    """Devuelve (sigma, G [m/s], B [#/(m3 s)]). sigma = C/C* - 1 (0 si subsaturado)."""
    sigma = max(C / cstar(T_C, p) - 1.0, 0.0)
    arr = math.exp(-p.get("Eg_J_mol", 0.0) / R_GAS * (1.0 / (T_C + 273.15) - 1.0 / (p["reference_temperature_C"] + 273.15)))
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


# ------------------------------------------------------------------------------
# Cinética de Disolución (Sub-saturación: C < C*, Nagy, Fujiwara & Braatz Ec. 11.31)
# ------------------------------------------------------------------------------
def dissolution_rate_m_s(C_kg_m3, C_star_kg_m3, p=None):
    """
    Tasa lineal de disolución D(-S) <= 0 [m/s] cuando C < C*:
        D(-S) = -k_d * ((C* - C) / C*)^d_exp
    """
    C = float(C_kg_m3)
    C_star = max(float(C_star_kg_m3), 1e-9)
    if C >= C_star:
        return 0.0
    p = p or {}
    k_d = float(p.get("kd_dissolution_m_s", 1.2e-5))
    d_exp = float(p.get("d_dissolution_order", 1.0))
    undersat = (C_star - C) / C_star
    return -k_d * (undersat ** d_exp)


def net_growth_or_dissolution_rate(C_kg_m3, C_star_kg_m3, T_K, p, allow_dissolution=True):
    """
    Velocidad lineal neta [m/s]: crecimiento G > 0 si C > C*, 0 en saturación,
    y disolución D < 0 si C < C* (cuando allow_dissolution=True).
    """
    C = float(C_kg_m3)
    C_star = max(float(C_star_kg_m3), 1e-9)
    if C < C_star:
        return dissolution_rate_m_s(C, C_star, p) if allow_dissolution else 0.0
    if C == C_star:
        return 0.0
    sigma = (C - C_star) / C_star
    kg = float(p.get("kg_m_s", p.get("k_g", 1.0e-6)))
    g_ord = float(p.get("g_order", p.get("g", 1.31)))
    Eg = float(p.get("Eg_J_mol", 0.0))
    T_ref = float(p.get("T_ref_K", 313.15))
    R_gas = 8.314462618
    arrh = float(np.exp(-Eg / R_gas * (1.0 / max(float(T_K), 200.0) - 1.0 / T_ref))) if Eg != 0.0 else 1.0
    return kg * arrh * (sigma ** g_ord)


# ------------------------------------------------------------------------------
# Termodinámica Multicomponente KNO3 / KClO4 Calibrada (Modelo_Cristalizador_Solubilidad)
# Ecuaciones de Apelblat (CRC, 0-100 °C) + ion común K+ + ensayo SQM (SSKd, 25-35 °C)
# ------------------------------------------------------------------------------
MOLAR_MASS_G_MOL = {
    "KNO3": 101.10, "KClO4": 138.55, "K": 39.10, "NO3": 62.00,
    "ClO4": 99.45, "Mg": 24.305, "Na": 22.99, "NaNO3": 84.99,
}

# Parámetros Apelblat: m_sat [mol/kg H2O] = exp(A + B/T_K + C*ln(T_K))
# Error máx. vs tabla CRC (0-100 °C): KNO3 = 1.4 %, KClO4 = 2.4 %
APELBLAT_PARAMS = {
    "KNO3":  (110.02522951788562, -7605.497487167708, -14.600400228580432),
    "KClO4": (-11.130388918105526, -2481.900639355448, 3.084914560684994),
}

# Agua por litro de SSKd (x_agua=0.3994, rho=1518 kg/m3) y composición molal base [mol/kg H2O]
W_AGUA_SSKD_KG_L = 0.3994 * 1.518
COMP_SSKD_DEFAULT = {
    "K":    223.5 / MOLAR_MASS_G_MOL["K"] / W_AGUA_SSKD_KG_L,
    "NO3":  (607.5 * MOLAR_MASS_G_MOL["NO3"] / MOLAR_MASS_G_MOL["NaNO3"]) / MOLAR_MASS_G_MOL["NO3"] / W_AGUA_SSKD_KG_L,
    "Mg":   27.1 / MOLAR_MASS_G_MOL["Mg"] / W_AGUA_SSKD_KG_L,
    "ClO4": 4.1 / MOLAR_MASS_G_MOL["ClO4"] / W_AGUA_SSKD_KG_L,
}

# Factores de actividad aparentes en salmuera (lnGamma0 KClO4 calibrado con ensayo SQM 25-35 °C, error 2.3 %)
# k_Mg parametrizable (en 0.0 por defecto hasta contar con segundo ensayo experimental a distinto Mg2+)
TERMO_SQM_DEFAULT = {
    "lnGamma0": {"KNO3": 0.0, "KClO4": 1.57086574},
    "k_Mg":     {"KNO3": 0.0, "KClO4": 0.0},
    "m_Mg_ref": COMP_SSKD_DEFAULT["Mg"],
}


def make_termo_sqm(k_Mg_KClO4=0.0, k_Mg_KNO3=0.0, lnGamma0_KClO4=None, lnGamma0_KNO3=0.0, m_Mg_ref=None):
    """Construye un diccionario termodinámico compatible con Cristalizador_SQM.ipynb."""
    return {
        "lnGamma0": {
            "KNO3": float(lnGamma0_KNO3),
            "KClO4": float(TERMO_SQM_DEFAULT["lnGamma0"]["KClO4"] if lnGamma0_KClO4 is None else lnGamma0_KClO4),
        },
        "k_Mg": {"KNO3": float(k_Mg_KNO3), "KClO4": float(k_Mg_KClO4)},
        "m_Mg_ref": float(TERMO_SQM_DEFAULT["m_Mg_ref"] if m_Mg_ref is None else m_Mg_ref),
    }


def apelblat_molal(sal: str, T_C):
    """Solubilidad en agua pura [mol/kg H2O] por ecuación de Apelblat."""
    A, B, C = APELBLAT_PARAMS[sal]
    T_K = np.asarray(T_C, dtype=float) + 273.15
    val = np.exp(A + B / T_K + C * np.log(T_K))
    return float(val) if np.ndim(T_C) == 0 else val


def ksp_aparente(sal: str, T_C):
    """Producto de solubilidad aparente en agua pura: K'_sp = m_sat^2 [(mol/kg H2O)^2]."""
    m = apelblat_molal(sal, T_C)
    return m * m


def gamma_salmuera(sal: str, m_Mg: float, termo: dict | None = None) -> float:
    """Factor de actividad aparente Gamma(sal, m_Mg) = exp(lnGamma0 + k_Mg * (m_Mg - m_Mg_ref))."""
    t = TERMO_SQM_DEFAULT if termo is None else termo
    return float(np.exp(t["lnGamma0"][sal] + t["k_Mg"][sal] * (float(m_Mg) - float(t["m_Mg_ref"]))))


def sobresaturacion_salmuera(sal: str, T_C: float, comp: dict, termo: dict | None = None) -> float:
    """
    Sobresaturación unificada (Cristalizador_SQM.ipynb §2.3):
        sigma = sqrt(IAP / (Gamma * K'_sp(T))) - 1
    donde sigma < 0 indica subsaturación y sigma > 0 sobresaturación.
    """
    anion = {"KNO3": "NO3", "KClO4": "ClO4"}[sal]
    IAP = max(float(comp["K"]) * float(comp[anion]), 0.0)
    denom = max(gamma_salmuera(sal, float(comp["Mg"]), termo) * float(ksp_aparente(sal, T_C)), 1e-18)
    return float(math.sqrt(IAP / denom) - 1.0)


def m_K_eq_KNO3(T_C: float, comp: dict, termo: dict | None = None) -> float:
    """Molalidad de K+ [mol/kg H2O] en equilibrio con KNO3 sólido conservando exceso d = NO3 - K."""
    d = float(comp["NO3"]) - float(comp["K"])
    Kp = gamma_salmuera("KNO3", float(comp["Mg"]), termo) * float(ksp_aparente("KNO3", T_C))
    return float(0.5 * (-d + math.sqrt(d * d + 4.0 * Kp)))


def m_ClO4_sat(T_C: float, m_K: float, m_Mg: float, termo: dict | None = None) -> float:
    """Solubilidad molal de ClO4- [mol/kg H2O] a temperatura T_C para K+ y Mg2+ dados (ion común)."""
    return float(gamma_salmuera("KClO4", m_Mg, termo) * float(ksp_aparente("KClO4", T_C)) / max(float(m_K), 1e-12))


def c_KClO4_sat_g_L(T_C: float, comp: dict | None = None, termo: dict | None = None,
                    w_agua: float = W_AGUA_SSKD_KG_L, as_ion: bool = False) -> float:
    """Solubilidad de KClO4 (o ion ClO4- si as_ion=True) en [g/L] para licor saturado en KNO3."""
    c = COMP_SSKD_DEFAULT if comp is None else comp
    mK = m_K_eq_KNO3(T_C, c, termo)
    m_sat = m_ClO4_sat(T_C, mK, float(c["Mg"]), termo)
    mm = MOLAR_MASS_G_MOL["ClO4"] if as_ion else MOLAR_MASS_G_MOL["KClO4"]
    return float(m_sat * mm * float(w_agua))


def T_saturacion_KClO4(comp: dict, termo: dict | None = None, T_lo: float = 5.0, T_hi: float = 95.0) -> float:
    """Temperatura [°C] a la cual sigma_KClO4 = 0 para licor en equilibrio con KNO3."""
    def _f(T):
        c = dict(comp)
        c["K"] = min(float(comp["K"]), m_K_eq_KNO3(T, comp, termo))
        c["NO3"] = float(comp["NO3"]) - (float(comp["K"]) - c["K"])
        return sobresaturacion_salmuera("KClO4", T, c, termo)

    if _f(T_lo) < 0.0:
        return float("nan")
    if _f(T_hi) > 0.0:
        return float("inf")
    return float(brentq(_f, T_lo, T_hi))


def holgura_termica_KClO4(T_op_C: float, comp: dict, termo: dict | None = None) -> float:
    """Holgura térmica respecto a saturación de KClO4 [°C]: > 0 seguro, < 0 sobresaturado."""
    return float(T_op_C) - T_saturacion_KClO4(comp, termo)
