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


def test_dpbe_solubilidad_kclo4_calibrada_sqm():
    """Valida Apelblat KNO3/KClO4 y ajuste SQM (lnGamma0=1.571, T_sat=41.3 °C, ensayo 25-35 °C)."""
    g_kno3_40 = dpbe.apelblat_molal("KNO3", 40.0) * dpbe.MOLAR_MASS_G_MOL["KNO3"] / 10.0
    g_kclo4_40 = dpbe.apelblat_molal("KClO4", 40.0) * dpbe.MOLAR_MASS_G_MOL["KClO4"] / 10.0
    assert abs(g_kno3_40 - 63.9) / 63.9 < 0.02
    assert abs(g_kclo4_40 - 3.63) / 3.63 < 0.03

    for T_exp, c_exp in [(35.0, 4.59), (30.0, 3.90), (25.0, 3.48)]:
        c_pred = dpbe.c_KClO4_sat_g_L(T_exp, dpbe.COMP_SSKD_DEFAULT)
        assert abs(c_pred - c_exp) / c_exp < 0.03
    T_sat = dpbe.T_saturacion_KClO4(dpbe.COMP_SSKD_DEFAULT)
    assert 40.5 <= T_sat <= 42.0


def test_dpbe_salting_out_k_Mg_parametrizable():
    """Un k_Mg < 0 debe reducir m_ClO4_sat y elevar T_saturacion_KClO4 al subir Mg2+."""
    comp_hi_mg = dict(dpbe.COMP_SSKD_DEFAULT)
    comp_hi_mg["Mg"] = 1.30 * dpbe.COMP_SSKD_DEFAULT["Mg"]

    t_cero = dpbe.make_termo_sqm(k_Mg_KClO4=0.0)
    t_neg = dpbe.make_termo_sqm(k_Mg_KClO4=-0.20)

    mK = dpbe.m_K_eq_KNO3(40.0, comp_hi_mg, t_neg)
    assert dpbe.m_ClO4_sat(40.0, mK, comp_hi_mg["Mg"], t_neg) < dpbe.m_ClO4_sat(40.0, mK, comp_hi_mg["Mg"], t_cero)
    assert dpbe.T_saturacion_KClO4(comp_hi_mg, t_neg) > dpbe.T_saturacion_KClO4(comp_hi_mg, t_cero)


def test_control_restriccion_dura_kclo4_y_derivacion_purga():
    """
    Verifica EDT 4.3:
      - Caso A: en modo sostenible (F_303 <= 2.5 t/h) con SSKd sin purga a 43 °C, la restricción
        dura de KClO4 supera el tope de agua y activa purge_required=True (AIC 3522 -> FV 3522).
      - Caso B: con purga AIC 3522 activa (ClO4 reducido 25 %), cumple holgura >= 0 a 0.5 t/h.
      - Caso C: con techo amplio de agua (40 t/h), el controlador eleva F_303 justo al mínimo
        termodinámico de KClO4 (~16 t/h) sin violar la sobresaturación de KNO3.
    """
    out_sust = C.dilution_kclo4_hard_constraint(
        C_in_kg_m3=660.0, C_star_kg_m3=600.0, R_cmd=0.50 / 378.18,
        T_stage2_C=43.0, R_min=0.20 / 378.18, R_max=2.50 / 378.18,
    )
    assert out_sust["purge_required"] is True
    assert out_sust["F_303_applied_t_h"] <= 2.50 + 1e-6
    assert out_sust["F_303_kclo4_min_t_h"] > 10.0
    assert 0.05 < out_sust["clo4_reduction_needed_frac"] < 0.25

    comp_purged = dict(dpbe.COMP_SSKD_DEFAULT)
    comp_purged["ClO4"] *= 0.75
    out_ok = C.dilution_kclo4_hard_constraint(
        C_in_kg_m3=660.0, C_star_kg_m3=600.0, R_cmd=0.50 / 378.18,
        T_stage2_C=43.0, comp_ref=comp_purged,
        R_min=0.20 / 378.18, R_max=2.50 / 378.18,
    )
    assert out_ok["purge_required"] is False
    assert out_ok["holgura_termica_C"] >= 0.0
    assert np.isclose(out_ok["F_303_applied_t_h"], 0.50, atol=1e-6)

    out_wide = C.dilution_kclo4_hard_constraint(
        C_in_kg_m3=660.0, C_star_kg_m3=600.0, R_cmd=0.50 / 378.18,
        T_stage2_C=43.0, R_min=0.20 / 378.18, R_max=40.0 / 378.18,
    )
    assert out_wide["purge_required"] is False
    assert out_wide["kclo4_constraint_active"] is True
    assert abs(out_wide["holgura_termica_C"]) < 0.05
    assert 14.0 < out_wide["F_303_applied_t_h"] < 18.0
