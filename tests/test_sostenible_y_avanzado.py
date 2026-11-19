import numpy as np
import pandas as pd

from msmpr_core import control as C
from msmpr_core import dpbe, params, plant, replay


def test_compute_rga_propiedades_bristol():
    """Verifica que filas y columnas de la RGA suman 1.0 y detecta emparejamiento diagonal."""
    G = np.array([
        [ 2.0, -0.15,  0.05,  0.00],
        [ 0.1,  3.50, -0.10,  0.02],
        [ 0.0,  0.10, -1.80,  0.05],
        [ 0.0,  0.02,  0.12, -2.40],
    ])
    lam = C.compute_rga(G)
    assert lam.shape == (4, 4)
    assert np.allclose(lam.sum(axis=1), np.ones(4), atol=1e-10)
    assert np.allclose(lam.sum(axis=0), np.ones(4), atol=1e-10)
    assert np.all(np.diag(lam) > 0.90)


def test_lead_lag_compensator_ganancia_y_dinamica():
    """Verifica ganancia estacionaria exacta, impulso inicial en Lead neto y suavizado en Lag neto."""
    ll_lead = C.LeadLagCompensator(tau_lead=60.0, tau_lag=20.0, gain=2.0, u0=0.0)
    y_first = ll_lead.step(1.0, dt=1.0)
    assert y_first > 2.0
    for _ in range(500):
        y_ss = ll_lead.step(1.0, dt=1.0)
    assert np.isclose(y_ss, 2.0, atol=1e-3)

    ll_lag = C.LeadLagCompensator(tau_lead=5.0, tau_lag=40.0, gain=2.0, u0=0.0)
    y1 = ll_lag.step(1.0, dt=1.0)
    assert 0.0 < y1 < 2.0
    for _ in range(800):
        y_ss2 = ll_lag.step(1.0, dt=1.0)
    assert np.isclose(y_ss2, 2.0, atol=1e-3)


def test_supersaturation_override_protege_cristales():
    """El Override limita R para que C_in / (1 + R) >= (1 + sigma_min) * C_star."""
    C_in, C_star, sigma_min = 660.0, 600.0, 0.02
    out_hi = C.supersaturation_override_limit(C_in, C_star, R_cmd=0.12, sigma_min=sigma_min, R_min=0.001, R_max=0.15)
    assert out_hi["override_active"] is True
    assert out_hi["excess_ratio"] > 0.03
    C_mix = C_in / (1.0 + out_hi["R_applied"])
    assert C_mix >= (1.0 + sigma_min) * C_star - 1e-6

    out_ok = C.supersaturation_override_limit(C_in, C_star, R_cmd=0.005, sigma_min=sigma_min, R_min=0.001, R_max=0.15)
    assert out_ok["override_active"] is False
    assert np.isclose(out_ok["R_applied"], 0.005)


def test_dpbe_cinetica_disolucion_y_crecimiento():
    """Verifica rama negativa de disolución (C < C*), cero en saturación y crecimiento (C > C*)."""
    p = params.make_params()
    assert dpbe.dissolution_rate_m_s(550.0, 600.0, p) < 0.0
    assert dpbe.dissolution_rate_m_s(600.0, 600.0, p) == 0.0
    assert dpbe.dissolution_rate_m_s(650.0, 600.0, p) == 0.0

    r_sub = dpbe.net_growth_or_dissolution_rate(550.0, 600.0, 330.0, p, allow_dissolution=True)
    r_sat = dpbe.net_growth_or_dissolution_rate(600.0, 600.0, 330.0, p, allow_dissolution=True)
    r_sup = dpbe.net_growth_or_dissolution_rate(650.0, 600.0, 330.0, p, allow_dissolution=True)
    r_off = dpbe.net_growth_or_dissolution_rate(550.0, 600.0, 330.0, p, allow_dissolution=False)
    assert r_sub < 0.0
    assert r_sat == 0.0
    assert r_sup > 0.0
    assert r_off == 0.0


def test_magnesium_perchlorate_coupling_monotonia():
    """Al subir Mg2+ debe caer salting_factor y subir viscosidad, Kd_eff y humedad de torta."""
    p = params.make_params(**params.SOSTENIBLE_DILUCION_OVERRIDES)
    low = plant.magnesium_perchlorate_coupling(28.0, 1.5e-7, 550.0, p)
    high = plant.magnesium_perchlorate_coupling(42.0, 1.5e-7, 550.0, p)
    assert high["salting_factor"] < low["salting_factor"]
    assert high["eta_ratio"] > low["eta_ratio"]
    assert high["Kd_eff"] > low["Kd_eff"]
    assert high["w_moist_eff"] > low["w_moist_eff"]


def test_modo_sostenible_05tph_cumple_kpi_y_ahorra_agua():
    """Con F_303_nom=0.50 t/h y K_trim escalado, el caudal se mantiene acotado y Área 6000 < 40 ppm."""
    df = plant.simulate_plant_npt3_calibrated(
        use_AIT3105=True, use_AIC3255=True, tuning_mode="TL_damped",
        F_303_nom=0.50, F_303_min=0.20, F_303_max=2.50,
    )
    assert df["F_303_actual_t_h"].mean() < 2.0
    assert df["F_303_actual_t_h"].max() <= 2.50 + 1e-6
    assert df["ppm_ClO4_final_6000"].max() < 40.0


def test_replay_modulo_contrafactual_sintetico():
    """Verifica el módulo msmpr_core.replay sin depender de archivos Excel externos."""
    idx = pd.date_range("2026-11-11 00:00:00", periods=180, freq="10min")
    t_h = np.arange(len(idx)) / 6.0
    pulse = (t_h >= 6.0) & (t_h <= 18.0)
    df_syn = pd.DataFrame({
        "t_h": t_h,
        "Q_SSKd_m3_h": np.full(len(idx), 243.2),
        "Q_dil_real": np.full(len(idx), 20.0),
        "Q_purga_real": np.where(np.arange(len(idx)) % 2 == 0, 36.5, 4.0),
        "T_A3210_real": np.full(len(idx), 43.0),
        "ClO4_SSK": np.where(pulse, 5.4, 4.7),
        "Mg_SSK": np.where(pulse, 33.5, 30.0),
        "ClO4_MLK_real": np.where(pulse, 5.15, 4.50),
        "Mg_MLK_real": np.where(pulse, 33.5, 30.0),
        "ClO4_MLRefi_real": np.full(len(idx), 1.60),
    }, index=idx)

    r_base = replay.simular_replay_npt3(df_syn, mode="baseline_sqm")
    r_prop = replay.simular_replay_npt3(df_syn, mode="proposed_integrated", lab_sampling_h=2.0)

    assert r_prop["Q_dil_t_h"].mean() < 1.5
    assert r_prop["Q_dil_t_h"].mean() < 0.10 * r_base["Q_dil_t_h"].mean()
    assert r_prop["ClO4_refi_ppm"].max() < 40.0
    assert r_prop["ClO4_refi_ppm"].max() < r_base["ClO4_refi_ppm"].max()
