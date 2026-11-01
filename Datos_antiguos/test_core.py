"""Pruebas de aceptación (guía §8, criterios): balance, estacionario, identificabilidad y control."""
import itertools

import numpy as np
import pandas as pd
import pytest

from msmpr_core import dpbe, control as C
from msmpr_core.fit import ajustar_cineticas
from msmpr_core.model import msmpr
from msmpr_core.params import make_params, placeholders
from msmpr_core.plant import Plant, simulate

P = make_params()


def test_gamma4_ratios_and_span():
    assert dpbe.Q10 == pytest.approx(1.745, abs=2e-3)
    assert dpbe.Q50 == pytest.approx(3.672, abs=2e-3)
    assert dpbe.Q90 == pytest.approx(6.681, abs=2e-3)
    assert dpbe.SPAN_IDEAL == pytest.approx(1.344, abs=2e-3)


def test_grid_matches_analytic():
    tau = P["volume_m3"] / P["feed_flow_m3_s"]
    Cin = 0.95 * P["feed_stock_kg_m3"]["KNO3"]
    a = dpbe.stationary_msmpr(Cin, tau, 40.0, P)
    g = dpbe.stationary_grid(Cin, tau, 40.0, P, Lmax=1.2e-3, bins=3000)
    assert g["d50"] == pytest.approx(a["d50"], rel=0.02)
    assert g["M_T"] == pytest.approx(a["M_T"], rel=0.02)


def test_d50_exponent_vs_tau():
    """guía §2(b): d50 ~ tau^((b-g)/(b+3g))."""
    q = make_params(Cstar_ref_kg_m3=120.0)
    b, g = q["b"], q["g"]
    d1, d2 = (msmpr(t, 0.1, 40.0, None, q)["d50_um"] for t in (1.0, 2.0))
    # el exponente teórico vale para sigma<<1 (C* ~ constante); tolerancia amplia
    assert np.log(d2 / d1) / np.log(2.0) == pytest.approx((b - g) / (b + 3 * g), abs=0.35)


def test_mass_balance_closes_two_ways():
    for tau, f, T in itertools.product([0.75, 1.5, 2.25], [0.0, 0.1, 0.2], [30, 35, 40]):
        m = msmpr(tau, f, T, None, make_params(Cstar_ref_kg_m3=120.0))
        assert abs(m["Y_cierre"]) < 1e-6
        assert 0.0 <= m["Y"] <= 1.0


def test_steady_state_is_stationary_and_closes():
    pl = Plant(P)
    for alpha in (1.0, 0.5, 0.25):
        x, u, d = pl.steady_state(alpha=alpha)
        assert np.max(np.abs(pl.rhs(x, u, d)[:4])) < 1e-9
        assert abs(pl.metrics(x, u, d)["Y_cierre"]) < 1e-8


def test_open_loop_has_no_drift():
    pl = Plant(P)
    x, u, d = pl.steady_state()
    df = simulate(pl, x, 4 * 3600, 10.0, lambda t, xx, pp: u)
    for c in ("V", "T", "C", "CP", "mu3"):
        assert abs(df[c].iloc[-1] - df[c].iloc[0]) / abs(df[c].iloc[0]) < 1e-9


def test_purge_is_set_by_water_balance():
    """Con volumen fijo la purga NO depende de alpha: la fija (alimentación - evaporación)."""
    pl = Plant(P)
    purge = []
    for a in (0.5, 0.25, 0.1):
        x, u, d = pl.steady_state(alpha=a)
        purge.append(pl.metrics(x, u, d)["purge_m3_h"])
    assert max(purge) / min(purge) < 1.02


def test_pid_antiwindup_and_bumpless():
    pid = C.PID(Kp=1.0, Ki=0.1, Kd=0.0, u_bar=0.5, u_min=0.0, u_max=1.0)
    for k in range(200):
        pid.update(k, pv=0.0, sp=10.0)           # error enorme: salida saturada
    assert pid.u == 1.0 and abs(pid.I) < 1.5       # el integral no se dispara
    pid.update(201, pv=10.0, sp=10.0, track=0.3)   # seguimiento
    assert pid.u == pytest.approx(0.3)


def test_fit_recovers_synthetic_parameters():
    p = make_params(Cstar_ref_kg_m3=120.0)
    true = {"B_ref": 3e6, "b": 2.5, "G_ref": 2e-8, "g": 1.5, "sigma_ref": 0.04}
    rng = np.random.default_rng(1)
    rows = []
    for tau, f, T in itertools.product([0.75, 1.5, 2.25], [0.0, 0.1, 0.2], [30, 40, 50]):
        r = msmpr(tau, f, T, true, p)
        nz = lambda v, s: v * np.exp(rng.normal(0, s))
        rows.append(dict(tau_h=tau, f_dil=f, T=T, d50_um=nz(r["d50_um"], .05), N_kg=nz(r["N_kg"], .05),
                         Y=r["Y"] + rng.normal(0, .01), sigma=r["sigma"] * np.exp(rng.normal(0, .08)),
                         s_d50_um=.05, s_N_kg=.05, s_Y=.01, s_sigma=.08 * r["sigma"]))
    out = ajustar_cineticas(pd.DataFrame(rows), p)
    truth = [np.log(true["B_ref"]), true["b"], np.log(true["G_ref"]), true["g"]]
    assert out["exito"]
    for (lab, (lo, hi)), t in zip(out["IC95"].items(), truth):
        assert lo <= t <= hi, lab


def test_closed_loop_rejects_disturbances():
    pl = Plant(P)
    x, u, d0 = pl.steady_state(alpha=0.25)
    tun, _ = C.autotune(pl, x, u, d0)
    m0 = pl.metrics(x, u, d0)
    sp = {"T": x[1], "V": x[0], "rho": m0["slurry_density_kg_m3"]}

    def dist(t):
        d = dict(d0)
        if t >= 2 * 3600:
            d["T_feed"] += 6.0
            d["Q_stock"] *= 1.12
        return d

    ol = simulate(pl, x, 10 * 3600, 10.0, lambda t, xx, pp: u, disturbance=dist)
    cl = simulate(pl, x, 10 * 3600, 10.0, C.PlantController(pl, u, sp, tun, enabled=("T", "V", "rho")), disturbance=dist)
    assert cl["T"].max() < 41.0 and ol["T"].max() > 43.0
    assert cl["V"].max() < 680.0 and ol["V"].max() > 800.0


def test_placeholders_are_declared():
    assert "Cstar_ref_kg_m3" in placeholders() and "kG_m_s" in placeholders()


def test_dilution_loop_does_not_chatter():
    """Regresión: el esclavo FIRCAL con tau_c < dt oscilaba entre sus límites."""
    pl = Plant(P)
    x, u, d0 = pl.steady_state(alpha=0.25)
    tun, _ = C.autotune(pl, x, u, d0)
    m0 = pl.metrics(x, u, d0)
    sp = {"T": x[1], "V": x[0], "rho": m0["slurry_density_kg_m3"]}

    def dist(t):
        d = dict(d0)
        if t >= 2 * 3600:
            d["T_feed"] += 6.0
            d["Q_stock"] *= 1.12
        return d

    cl = simulate(pl, x, 8 * 3600, 10.0, C.PlantController(pl, u, sp, tun, enabled=("T", "V", "rho")), disturbance=dist)
    q = cl["u_Q_dil_cmd"].values * 3600.0
    assert np.abs(np.diff(q)).max() < 1.0
