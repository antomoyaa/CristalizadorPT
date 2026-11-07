"""Pruebas del control integrado (lazo abierto/cerrado, ruido OU, IMC, Routh-Hopf y tren NPT-3).

Ruta: tests/test_control_integrado.py. Criterios en docs/CAMBIOS.md.
"""
import numpy as np
import pytest
from msmpr_core.params import make_params
from msmpr_core.plant import Plant, simulate, simulate_plant_npt3_calibrated
from msmpr_core import control as C
from msmpr_core.control import (
    ConsignaVariable,
    Ruido,
    imc_pi,
    iae,
    metricas_consigna,
    sigma_ventana,
    routh_hopf_analysis,
    DelayedSensor,
    CommercialPID,
)


@pytest.fixture(scope="module")
def setup_a3110():
    """Inicializa el modelo estacionario nominal del cristalizador A3110 y su sintonía base."""
    p = make_params()
    pl = Plant(p)
    x0, u0, d0 = pl.steady_state(alpha=0.25)
    m0 = pl.metrics(x0, u0, d0)
    dt = 5.0
    tun, info = C.autotune(pl, x0, u0, d0, dt_s=dt)
    sp = {"T": x0[1], "V": x0[0], "rho": 0.0}
    R0 = u0["Q_dil_cmd"] / d0["Q_stock"]
    return p, pl, x0, u0, d0, m0, dt, tun, info, sp, R0


def test_open_loop_steady_state_and_balances(setup_a3110):
    """Verifica reposo estacionario, identificación FOPDT (+/-10%) y balances de volumen/válvula."""
    p, pl, x0, u0, d0, m0, dt, tun, info, sp, R0 = setup_a3110
    tp = 3600.0

    # 1. Reposo estacionario nulo en lazo abierto durante 6 h
    df0 = simulate(pl, x0, 6 * 3600, dt, lambda t, x, pl_: u0, disturbance=None, log_every=60)
    for col in ("V", "T", "C", "CP", "ppm_ClO4", "recovery_pct"):
        assert float(np.max(np.abs(df0[col] - df0[col].iloc[0]))) < 1e-6

    # 2. Validación FOPDT en escalón positivo (+10%) y negativo (-10%) de m_evap
    K_, tau_, th_ = info["T"]["K"], info["T"]["tau_s"], info["T"]["theta_s"]
    for sg in (+1, -1):
        dmev = sg * 0.10 * u0["m_evap"]
        t_s, T_sim = C.step_test(pl, x0, u0, d0, "m_evap", dmev, lambda x, u, dd: x[1], 12 * 3600, dt_s=20.0)
        T_fopdt = x0[1] + K_ * dmev * (1.0 - np.exp(-np.clip(t_s - th_, 0.0, None) / tau_))
        rmse = float(np.sqrt(np.mean((T_sim - T_fopdt) ** 2)))
        assert rmse < 0.35
        assert np.sign(T_sim[-1] - x0[1]) == np.sign(K_ * dmev)

    # 3. Balance de volumen integrador (dV/dt inicial = delta Q_stock, sube y baja)
    qs0 = d0["Q_stock"]
    for sg in (+1, -1):
        dist_q = lambda t, s=sg: {**d0, "Q_stock": qs0 * (1.0 + s * 0.10 if t >= tp else 1.0)}
        df_q = simulate(pl, x0, 3 * 3600, dt, lambda t, x, pl_: u0, disturbance=dist_q, log_every=12)
        win = df_q[(df_q.t_h * 3600 >= tp) & (df_q.t_h * 3600 <= tp + 600)]
        pend = (win.V.iloc[-1] - win.V.iloc[0]) / (3600.0 * (win.t_h.iloc[-1] - win.t_h.iloc[0]))
        assert abs(pend / (sg * 0.10 * qs0) - 1.0) < 0.03

    # 4. Ganancia efectiva de válvula de dilución (0.80 y 1.20)
    q0 = u0["Q_dil_cmd"]
    for gain in (0.8, 1.2):
        dist_v = lambda t, g=gain: {**d0, "valve_gain": g if t >= tp else 1.0}
        df_v = simulate(pl, x0, 3 * 3600, dt, lambda t, x, pl_: u0, disturbance=dist_v, log_every=12)
        assert abs(df_v.Qd.iloc[-1] / q0 - gain) < 1e-3


def test_closed_loop_zero_offset_consistency(setup_a3110):
    """Comprueba que el estado final en lazo cerrado converge al estacionario teórico (offset nulo)."""
    p, pl, x0, u0, d0, m0, dt, tun, info, sp, R0 = setup_a3110
    tp = 3600.0
    escenarios = [
        lambda t: {**d0, "T_feed": d0["T_feed"] + (5.0 if t >= tp else 0.0)},
        lambda t: {**d0, "T_feed": d0["T_feed"] - (5.0 if t >= tp else 0.0)},
        lambda t: {**d0, "Q_stock": d0["Q_stock"] * (1.10 if t >= tp else 1.0)},
    ]
    for dist_fn in escenarios:
        qsf = lambda t, f=dist_fn: f(t)["Q_stock"]
        ctrl = C.PlantController(pl, u0, sp, tun, enabled=("T", "V", "ratio"), q_stock_meas=qsf, ratio_nom=R0)
        dfc = simulate(pl, x0, 24 * 3600, dt, ctrl, disturbance=dist_fn, log_every=60)
        d_end = dist_fn(24 * 3600)
        xs, us, ds = pl.steady_state(alpha=0.25, d=d_end, Q_dil=R0 * d_end["Q_stock"])
        ms = pl.metrics(xs, us, ds)
        assert abs(dfc["T"].iloc[-1] - sp["T"]) < 0.05
        assert abs(dfc["V"].iloc[-1] - sp["V"]) < 0.5
        assert 100.0 * abs(dfc.recovery_pct.iloc[-1] / ms["recovery_pct"] - 1.0) < 0.1


def test_colored_noise_attenuation_and_invariance(setup_a3110):
    """Verifica que el lazo cerrado atenúa sigma_T y sigma_V bajo ruido, pero no sigma_ppm sin AIT."""
    p, pl, x0, u0, d0, m0, dt, tun, info, sp, R0 = setup_a3110
    ruido = Ruido(seed=7, dt=dt, tc=600.0)
    dist_ruido = lambda t: ruido(dict(d0), t)
    H = 18 * 3600

    ol_r = simulate(pl, x0, H, dt, lambda t, x, pl_: u0, disturbance=dist_ruido, log_every=12)
    ctrl = C.PlantController(
        pl, u0, sp, tun, enabled=("T", "V", "ratio"), q_stock_meas=lambda t: dist_ruido(t)["Q_stock"], ratio_nom=R0
    )
    cl_r = simulate(pl, x0, H, dt, ctrl, disturbance=dist_ruido, log_every=12)

    assert sigma_ventana(cl_r, "T") < 0.65 * sigma_ventana(ol_r, "T")
    assert sigma_ventana(cl_r, "V") < 0.55 * sigma_ventana(ol_r, "V")
    # Sin AIT 3105, la variabilidad de ppm de ClO4- permanece prácticamente igual
    assert abs(sigma_ventana(cl_r, "ppm_ClO4") - sigma_ventana(ol_r, "ppm_ClO4")) < 0.01


def test_setpoint_weighting_beta_and_imc_pareto(setup_a3110):
    """Verifica arranque sin salto con precarga beta != 1, reducción de sobrepaso y monotonía IMC."""
    p, pl, x0, u0, d0, m0, dt, tun, info, sp, R0 = setup_a3110
    tp, H_sp = 3600.0, 8 * 3600
    fT_down = lambda t: x0[1] - (2.0 if t >= tp else 0.0)

    resultados_beta = {}
    for beta in (1.0, 0.5):
        ctrl = C.PlantController(
            pl, u0, dict(sp), tun, enabled=("T", "V", "ratio"), q_stock_meas=lambda t: d0["Q_stock"], ratio_nom=R0
        )
        ctrl.pid_T.beta = beta
        ctrl.pid_T.I = -ctrl.pid_T.Kp * ctrl.pid_T.s * (beta * x0[1] - x0[1])
        df_b = simulate(pl, x0, H_sp, dt, ConsignaVariable(ctrl, {"T": fT_down}), log_every=12)
        # Antes del escalón (t < 1 h), el estacionario no debe perturbarse (arranque sin salto)
        pre_step = df_b[df_b.t_h < 0.9]
        assert float(np.max(np.abs(pre_step["T"] - x0[1]))) < 1e-5
        over, tset, iae_v = metricas_consigna(df_b.t_h.values * 3600, df_b["T"].values, x0[1], x0[1] - 2.0, tp)
        resultados_beta[beta] = (over, df_b.u_m_evap.max())

    # Con beta = 0.5 se reduce el golpe proporcional (menor pico de m_evap y menor sobrepaso)
    assert resultados_beta[0.5][0] < resultados_beta[1.0][0]
    assert resultados_beta[0.5][1] < resultados_beta[1.0][1]

    # Monotonía de la ganancia IMC respecto a epsilon
    K_, tau_, th_ = info["T"]["K"], info["T"]["tau_s"], max(info["T"]["theta_s"], 60.0)
    t_fast = imc_pi(K_, tau_, th_, 0.1 * tau_)
    t_slow = imc_pi(K_, tau_, th_, 1.0 * tau_)
    assert t_fast.Kp > t_slow.Kp


def test_routh_hopf_and_npt3_plantwide_kpi():
    """Valida el criterio de Routh-Hurwitz/Hopf y el cumplimiento de < 40 ppm en Área 6000."""
    # 1. Criterio de estabilidad Routh-Hurwitz con aproximación de Padé
    rh = routh_hopf_analysis(Kp=375.0, tau_p=6.3, theta=1.2, tau_I=3.8)
    assert 0.0 < abs(rh["Kc_TL"]) < abs(rh["Kc_ZN"]) < abs(rh["Kc_hopf"])
    assert rh["omega_u"] > 0.0
    assert np.isclose(rh["tau_u"], 2.0 * np.pi / rh["omega_u"])
    assert rh["tau_I_ZN"] < rh["tau_I_TL"]

    # 2. Comparación de sintonía Z-N vs Tyreus-Luyben en el tren integrado NPT-3
    df_zn = simulate_plant_npt3_calibrated(use_AIT3105=True, use_AIC3255=True, tuning_mode="ZN_oscillatory")
    df_tl = simulate_plant_npt3_calibrated(use_AIT3105=True, use_AIC3255=True, tuning_mode="TL_damped")
    tv_zn = float(np.sum(np.abs(np.diff(df_zn["P_vac2_bar"].values))))
    tv_tl = float(np.sum(np.abs(np.diff(df_tl["P_vac2_bar"].values))))
    assert tv_tl < tv_zn  # Tyreus-Luyben exige menor variación total (TV) en la válvula de vacío

    # 3. Verificación del lazo FFIC 3102 -> FIRC 3103 SIN vs CON AIT 3105 + AIC 3255
    df_sin_ait = simulate_plant_npt3_calibrated(use_AIT3105=False, use_AIC3255=False, tuning_mode="TL_damped")
    df_con_ambos = simulate_plant_npt3_calibrated(use_AIT3105=True, use_AIC3255=True, tuning_mode="TL_damped")

    max_ppm_sin = float(df_sin_ait["ppm_ClO4_final_6000"].max())
    max_ppm_con = float(df_con_ambos["ppm_ClO4_final_6000"].max())

    assert max_ppm_sin > 40.0  # Sin AIT 3105 viola la especificación de 40 ppm en Área 6000
    assert max_ppm_con < 40.0  # Con AIT 3105 + AIC 3255 cumple holgadamente (< 40 ppm)