"""Control PID para el modelo de planta (A3110), alineado con el P&ID 172164 hoja 10 (GEA/ABB).

Inspirado en los apuntes CBE30338 de J. Kantor (MIT, código): ponderación de consigna (beta, gamma), derivada
sobre la medición, anti-windup por retrocálculo, seguimiento (bumpless transfer) y sintonía a partir de FOPDT.
Implementación propia; no copia código de esos repositorios.

Estructura confirmada por el equipo para el agua de dilución (A3110):
  AIT 3105   analizador de ClO4 en la corriente de entrada; si viene alto, manda su señal a ...
  FFIC 3102  controlador de razón (dilución/alimentación): fija la consigna de ...
  FIRC 3103  caudal de agua de dilución -> válvula.  (el retardo del analizador aún NO está definido: parámetro libre)
Otros lazos del P&ID:
  PIRCSAHL 3101 vacío -> T (m_evap);  LIRCSAHL 3101 nivel -> descarga Q_out;
  DIRCAHL 3101 (densidad) -> FIRCAL 3101: actuador por identificar; NO es el agua de dilución. El modo 'rho' (maestro de densidad
  sobre la dilución) queda como alternativa de estudio, no como P&ID.
  Lazo D propuesto (ppm de laboratorio recortando consigna de densidad): alternativo, no está en el plano.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


# ------------------------------------------------------------------ PID
class PID:
    """PID paralelo, con ponderación de consigna, derivada filtrada sobre la medición y anti-windup.

    direct=False (acción inversa): si PV > SP, la salida baja.  direct=True: si PV > SP, la salida sube.
    Ki en 1/s (= Kp/Ti).  Todas las magnitudes de salida en las unidades de la variable manipulada.
    """

    def __init__(self, Kp, Ki, Kd, u_bar, u_min, u_max, direct=False, beta=1.0, gamma=0.0, N=10.0, Tt=None):
        self.Kp, self.Ki, self.Kd = float(Kp), float(Ki), float(Kd)
        self.u_bar, self.u_min, self.u_max = float(u_bar), float(u_min), float(u_max)
        self.s = -1.0 if direct else 1.0
        self.beta, self.gamma, self.N = beta, gamma, N
        Ti = self.Kp / self.Ki if self.Ki > 0 else np.inf
        Td = self.Kd / self.Kp if self.Kp > 0 else 0.0
        self.Tt = Tt if Tt is not None else (np.sqrt(Ti * Td) if np.isfinite(Ti) and Td > 0 else (Ti if np.isfinite(Ti) else 1.0))
        self.I = 0.0
        self.D = 0.0
        self.u = self.u_bar
        self._t = None
        self._ed = None

    def reset(self, u=None):
        self.I, self.D, self._t, self._ed = 0.0, 0.0, None, None
        self.u = self.u_bar if u is None else u

    def update(self, t, pv, sp, track=None):
        """Calcula la salida. `track` (opcional) fuerza seguimiento de una salida externa (modo manual/cascada abierta)."""
        dt = 0.0 if self._t is None else max(t - self._t, 0.0)
        s = self.s
        P = self.Kp * s * (self.beta * sp - pv)
        e_i = s * (sp - pv)
        e_d = s * (self.gamma * sp - pv)
        if dt > 0 and self.Kd > 0:
            Tf = self.Kd / (self.Kp * self.N) if self.Kp > 0 else dt
            self.D = Tf / (Tf + dt) * self.D + self.Kd / (Tf + dt) * (e_d - self._ed)
        if track is not None:
            self.I = track - self.u_bar - P - self.D                       # bumpless
            self.u = float(np.clip(track, self.u_min, self.u_max))
        else:
            u_raw = self.u_bar + P + self.I + self.D
            u_sat = float(np.clip(u_raw, self.u_min, self.u_max))
            self.I += self.Ki * e_i * dt + dt / self.Tt * (u_sat - u_raw)      # retrocálculo anti-windup
            self.u = float(np.clip(self.u_bar + P + self.I + self.D, self.u_min, self.u_max))
        self._t, self._ed = t, e_d
        return self.u


class SampledSensor:
    """Analizador de laboratorio: toma una muestra cada `period_s` y la entrega `delay_s` después (retención de orden cero)."""

    def __init__(self, period_s: float, delay_s: float, initial: float):
        self.period, self.delay = period_s, delay_s
        self.value = initial
        self._queue = deque()
        self._next = 0.0

    def read(self, t, true_value):
        if t >= self._next:
            self._queue.append((t + self.delay, true_value))
            self._next += self.period
        while self._queue and self._queue[0][0] <= t:
            self.value = self._queue.popleft()[1]
        return self.value


class DelayedSensor:
    """Medición continua con tiempo muerto compatible con PlantController (2 args) y modelo NPT-3 (4 a 7 args).

    ATENCIÓN UNIDADES: en modo PlantController (2 args) `delay_s` está en SEGUNDOS y el 2º arg es el valor inicial.
    En modo NPT-3 (si se pasa dt_min o initial_val) el 1er arg es el retardo en MINUTOS y el 2º es la constante de
    tiempo del filtro en minutos (tau_m). Preferir argumentos con nombre. TODO: separar en dos clases.
    """

    def __init__(
        self,
        delay_s: float,
        initial: float = 0.0,
        dt_min: float | None = None,
        initial_val: float | None = None,
        resolution: float = 0.0,
        noise_std: float = 0.0,
        seed: int = 42,
    ):
        if dt_min is None and initial_val is None:
            # Modo original de PlantController: DelayedSensor(delay_s, initial)
            self.delay = float(delay_s)
            self.value = float(initial)
            self._queue = deque()
            self._npt3_mode = False
        else:
            # Modo discreto para modelo NPT-3: DelayedSensor(theta_min, tau_m_min, dt_min, initial_val, resolution, noise_std, seed)
            theta_min = float(delay_s)
            tau_m_min = float(initial)
            self.delay_steps = max(int(round(theta_min / float(dt_min))), 1)
            self.buffer = deque([initial_val] * self.delay_steps, maxlen=self.delay_steps)
            self.tau_m = tau_m_min
            self.dt = float(dt_min)
            self.filtered_val = float(initial_val)
            self.res = float(resolution)
            self.noise_std = float(noise_std)
            self.rng = np.random.default_rng(seed)
            self._npt3_mode = True

    def read(self, t_or_val, true_value=None):
        if not self._npt3_mode and true_value is not None:
            # Llamado por PlantController: sensor.read(t, true_value)
            t = float(t_or_val)
            self._queue.append((t + self.delay, true_value))
            while self._queue and self._queue[0][0] <= t + 1e-9:
                self.value = self._queue.popleft()[1]
            return self.value
        else:
            # Llamado por modelo NPT-3: sensor.read(true_val)
            val = float(t_or_val)
            if self.noise_std > 0.0:
                val += float(self.rng.normal(0.0, self.noise_std))
            self.buffer.append(val)
            delayed = self.buffer[0]
            if self.tau_m > 0.0:
                self.filtered_val += (self.dt / self.tau_m) * (delayed - self.filtered_val)
            else:
                self.filtered_val = delayed
            if self.res > 0.0:
                return round(self.filtered_val / self.res) * self.res
            return self.filtered_val


# ------------------------------------------------------------------ identificación y sintonía
def fopdt_fit(t, y, du):
    """Ajusta y(t) = y0 + K du (1 - exp(-(t-theta)/tau)) para t > theta. Devuelve (K, tau, theta)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    y0 = y[0]
    span = max(t[-1] - t[0], 1e-9)

    def model(th):
        K, tau, theta = th
        z = np.clip(t - theta, 0.0, None)
        return y0 + K * du * (1.0 - np.exp(-z / max(tau, 1e-9)))

    K0 = (y[-1] - y0) / du
    best = None
    for tau0 in (0.05 * span, 0.2 * span, 0.5 * span):
        sol = least_squares(lambda th: model(th) - y, [K0, tau0, 0.01 * span],
                            bounds=([-np.inf, 1e-6, 0.0], [np.inf, 10 * span, span / 2]))
        if best is None or sol.cost < best.cost:
            best = sol
    return tuple(float(v) for v in best.x)


def simc_pi(K, tau, theta, tau_c=None):
    """Sintonía SIMC (Skogestad) para PI. Devuelve (Kp>0, Ti, direct) con direct=True si K<0."""
    tau_c = max(theta, 0.1 * tau) if tau_c is None else tau_c
    Kp = tau / (abs(K) * (tau_c + theta))
    Ti = min(tau, 4.0 * (tau_c + theta))
    return Kp, Ti, bool(K < 0)


def ziegler_nichols_pi(K, tau, theta):
    """Regla de reacción de curva de Ziegler-Nichols (PI); más agresiva que SIMC."""
    theta = max(theta, 1e-6)
    return 0.9 * tau / (abs(K) * theta), 3.33 * theta, bool(K < 0)


def step_test(plant, x_ss, u_ss, d, mv, delta, pv_fn, t_end_s, dt_s=10.0):
    """Escalón en lazo abierto sobre la variable manipulada `mv`; devuelve (t, pv)."""
    x = np.array(x_ss, float)
    u = dict(u_ss)
    ts, pvs = [], []
    for i in range(int(t_end_s / dt_s) + 1):
        t = i * dt_s
        uu = dict(u)
        if t >= dt_s:
            uu[mv] = u[mv] + delta
        ts.append(t)
        pvs.append(pv_fn(x, uu, d))
        x = plant.step(x, uu, d, dt_s)
    return np.array(ts), np.array(pvs)


# ------------------------------------------------------------------ controlador de planta (lazos A, B, C, D)
@dataclass
class LoopTuning:
    Kp: float
    Ki: float
    Kd: float = 0.0
    direct: bool = False


class PlantController:
    """Implementa los lazos A (T), B (nivel), C (cascada densidad -> flujo de dilución) y D (ClO4, opcional).

    sp: dict con consignas 'T', 'V', 'rho' (kg/m3). tun: dict {'T','V','rho','flow','ppm'} de LoopTuning.
    enabled: lazos activos; los inactivos mantienen la salida nominal (lazo abierto).
    """

    def __init__(self, plant, u_nom, sp, tun, enabled=("T", "V", "rho"), ppm_sp=None, lab=None, mv_limits=None,
                 flow_sp=None, noise_rel=None, seed=0, q_stock_meas=None, ratio_nom=None, ait=None):
        """enabled puede incluir 'flow': lazo FIRC solo (consigna `flow_sp`, número o función de t), sin maestro de densidad.
        noise_rel: desviación estándar relativa de la medición, p. ej. {'flow': 0.01}."""
        self.plant, self.u_nom, self.sp = plant, dict(u_nom), dict(sp)
        self.flow_sp = flow_sp
        self.q_stock_meas = q_stock_meas
        self.ratio_mode = "ratio" in set(enabled)
        if self.ratio_mode and (q_stock_meas is None or ratio_nom is None):
            raise ValueError("El modo 'ratio' exige q_stock_meas (función de t) y ratio_nom")
        self.ratio_nom = ratio_nom
        self.ait = None
        if "ait" in set(enabled):
            if not self.ratio_mode or ait is None:
                raise ValueError("'ait' exige modo 'ratio' y el dict ait={'cp_true': f(t), 'delay_s': ..., 'gain': 0-1}")
            cp_nom = plant.d_nom["CP_stock"]
            self.ait = {"cp_true": ait["cp_true"], "R_min": ait.get("R_min", 0.0), "R_max": ait.get("R_max", 3.0 * ratio_nom),
                        "CP_target": ait.get("CP_target", cp_nom / (1.0 + ratio_nom)),
                        "noise": ait.get("noise_rel", 0.0), "gain": ait.get("gain", 1.0), "R_nom": ratio_nom,
                        "sensor": DelayedSensor(ait.get("delay_s", 0.0), cp_nom)}
            self.last_R = ratio_nom
        self.noise_rel = noise_rel or {}
        self.rng = np.random.default_rng(seed)
        self.enabled = set(enabled)
        lim = {"m_evap": (0.2 * u_nom["m_evap"], 3.0 * u_nom["m_evap"]),
               "Q_out": (0.5 * u_nom["Q_out"], 1.6 * u_nom["Q_out"]),
               "Q_dil": (0.0, max(3.0 * u_nom["Q_dil_cmd"], 1e-4))}
        if mv_limits:
            lim.update(mv_limits)
        t = tun
        self.pid_T = PID(t["T"].Kp, t["T"].Ki, t["T"].Kd, u_nom["m_evap"], *lim["m_evap"], direct=t["T"].direct)
        self.pid_V = PID(t["V"].Kp, t["V"].Ki, t["V"].Kd, u_nom["Q_out"], *lim["Q_out"], direct=t["V"].direct)
        if self.ratio_mode:
            qs0 = plant.d_nom["Q_stock"]
            self.pid_rho = PID(t["rho"].Kp / qs0, t["rho"].Ki / qs0, t["rho"].Kd / qs0, ratio_nom, 0.0, 3.0 * ratio_nom,
                               direct=t["rho"].direct)
        else:
            self.pid_rho = PID(t["rho"].Kp, t["rho"].Ki, t["rho"].Kd, u_nom["Q_dil_cmd"], *lim["Q_dil"], direct=t["rho"].direct)
        self.pid_flow = PID(t["flow"].Kp, t["flow"].Ki, 0.0, u_nom["Q_dil_cmd"], *lim["Q_dil"], direct=False)
        self.pid_ppm = None
        if "ppm" in self.enabled:
            tp = t["ppm"]
            self.pid_ppm = PID(tp.Kp, tp.Ki, 0.0, sp["rho"], 0.97 * sp["rho"], 1.03 * sp["rho"], direct=tp.direct)
            self.ppm_sp = ppm_sp
            self.lab = lab
        self.last = dict(u_nom)

    def __call__(self, t, x, plant):
        u = dict(self.u_nom)
        d = plant.d_nom
        m = plant.metrics(x, u, d)
        V, T, Qd = x[0], x[1], x[10]
        Qd_meas = Qd * (1.0 + self.noise_rel.get("flow", 0.0) * self.rng.standard_normal())
        if "T" in self.enabled:
            u["m_evap"] = self.pid_T.update(t, T, self.sp["T"])
        if "V" in self.enabled:
            u["Q_out"] = self.pid_V.update(t, V, self.sp["V"])
        rho_sp = self.sp["rho"]
        if self.pid_ppm is not None:
            ppm_meas = self.lab.read(t, m["ppm_ClO4"])
            rho_sp = self.pid_ppm.update(t, ppm_meas, self.ppm_sp)
        if "rho" in self.enabled:
            out = self.pid_rho.update(t, m["slurry_density_kg_m3"], rho_sp)        # maestro DIRCAHL 3101
            q_sp = out * self.q_stock_meas(t) if self.ratio_mode else out            # razón x alimentación, o caudal directo
            u["Q_dil_cmd"] = self.pid_flow.update(t, Qd_meas, q_sp)                   # esclavo FIRCAL 3101
        elif self.ait is not None:                                                   # AIT 3105 -> FFIC 3102 -> FIRC 3103
            a = self.ait
            cp = a["sensor"].read(t, a["cp_true"](t) * (1.0 + a["noise"] * self.rng.standard_normal()))
            full = cp / a["CP_target"] - 1.0                                                # razón que deja la mezcla en CP_target
            self.last_R = float(np.clip(a["R_nom"] + a["gain"] * (full - a["R_nom"]), a["R_min"], a["R_max"]))
            u["Q_dil_cmd"] = self.pid_flow.update(t, Qd_meas, self.last_R * self.q_stock_meas(t))
        elif self.ratio_mode:                                                        # solo prealimentación (razón fija)
            u["Q_dil_cmd"] = self.pid_flow.update(t, Qd_meas, self.ratio_nom * self.q_stock_meas(t))
        elif "flow" in self.enabled:                                                 # Etapa 1: solo FIRCAL 3101
            q_sp = self.flow_sp(t) if callable(self.flow_sp) else (self.u_nom["Q_dil_cmd"] if self.flow_sp is None else self.flow_sp)
            u["Q_dil_cmd"] = self.pid_flow.update(t, Qd_meas, q_sp)
        self.last = u
        return u


# ------------------------------------------------------------------ sintonía automática por ensayos de escalón
def simc_pi_integrating(k_prime, theta, tau_c):
    """SIMC para proceso integrador: Kp = 1/(|k'| (tau_c+theta)), Ti = 4 (tau_c+theta). k' = d(dPV/dt)/dMV."""
    Kp = 1.0 / (abs(k_prime) * (tau_c + theta))
    return Kp, 4.0 * (tau_c + theta), bool(k_prime < 0)


def _to_tuning(Kp, Ti, direct):
    return LoopTuning(Kp=Kp, Ki=Kp / Ti, Kd=0.0, direct=direct)


def autotune(plant, x_ss, u_ss, d, aggressiveness=1.0, t_test_h=6.0, dt_s=10.0, verbose=False):
    """Sintoniza los lazos A, B y C (maestro y esclavo) con ensayos de escalón y reglas SIMC.

    Orden: A y B primero (lazo abierto); C se identifica con A y B ya cerrados. `aggressiveness` >1 acorta tau_c.
    Devuelve (tun: dict de LoopTuning, info: dict con los modelos identificados).
    """
    info, tun = {}, {}
    t_end = t_test_h * 3600.0

    # A: T <- m_evap (de primer orden)
    dmev = 0.10 * u_ss["m_evap"]
    t, y = step_test(plant, x_ss, u_ss, d, "m_evap", dmev, lambda x, u, dd: x[1], t_end, dt_s)
    K, tau, th = fopdt_fit(t, y, dmev)
    tun["T"] = _to_tuning(*simc_pi(K, tau, th, tau_c=max(th, 0.1 * tau) / aggressiveness))
    info["T"] = {"K": K, "tau_s": tau, "theta_s": th, "unid": "K [C/(kg/s)]"}

    # B: V <- Q_out (integrador): pendiente de V tras el escalón
    dq = 0.02 * u_ss["Q_out"]
    t, y = step_test(plant, x_ss, u_ss, d, "Q_out", dq, lambda x, u, dd: x[0], 1800.0, dt_s)
    kprime = np.polyfit(t[len(t) // 2:], y[len(t) // 2:], 1)[0] / dq
    tau_c = 600.0 / aggressiveness
    tun["V"] = _to_tuning(*simc_pi_integrating(kprime, 0.0, tau_c))
    info["V"] = {"k_prime": kprime, "tau_c_s": tau_c, "unid": "dV/dt por Q_out [m3/s por m3/s]"}

    # esclavo FIRCAL 3101: primer orden de constante plant.tau_valve (K=1)
    # tau_c >= 2 dt: un lazo muestreado cada dt con tau_c < dt es inestable en tiempo discreto
    tun["flow"] = _to_tuning(*simc_pi(1.0, plant.tau_valve, 0.0, tau_c=max(plant.tau_valve / 3.0, 2.0 * dt_s)))
    info["flow"] = {"K": 1.0, "tau_s": plant.tau_valve}

    # C maestro: densidad <- Q_dil_cmd con A y B cerrados
    sp = {"T": x_ss[1], "V": x_ss[0], "rho": plant.metrics(x_ss, u_ss, d)["slurry_density_kg_m3"]}
    base = PlantController(plant, u_ss, sp, {**tun, "rho": LoopTuning(0, 0), "flow": tun["flow"]}, enabled=("T", "V"))
    dqd = max(0.15 * u_ss["Q_dil_cmd"], 1e-4)
    x = np.array(x_ss, float)
    ts, ys = [], []
    for i in range(int(t_end / dt_s) + 1):
        tt = i * dt_s
        u = base(tt, x, plant)
        u["Q_dil_cmd"] = u_ss["Q_dil_cmd"] + (dqd if tt >= dt_s else 0.0)
        ts.append(tt)
        ys.append(plant.metrics(x, u, d)["slurry_density_kg_m3"])
        x = plant.step(x, u, d, dt_s)
    K, tau, th = fopdt_fit(np.array(ts), np.array(ys), dqd)
    tun["rho"] = _to_tuning(*simc_pi(K, tau, th, tau_c=max(th, 0.25 * tau) / aggressiveness))
    info["rho"] = {"K": K, "tau_s": tau, "theta_s": th, "unid": "K [kg/m3 por m3/s]"}
    if verbose:
        for k, v in info.items():
            print(k, v)
    return tun, info


def tune_ppm_loop(plant, u_ss_alpha=None, delta_d=0.02, lab_delay_s=3600.0, lab_period_s=7200.0, tau_s=None):
    """Lazo D (propuesto): ppm ClO4 -> consigna de densidad del lazo C."""
    from .plant import Plant
    p = plant.p
    x0, u0, d0 = plant.steady_state(alpha=u_ss_alpha)
    m0 = plant.metrics(x0, u0, d0)
    p2 = dict(p)
    p2["feed_dilution_fraction"] = p["feed_dilution_fraction"] + delta_d
    pl2 = Plant(p2, anchor=plant.anchor)
    x1, u1, d1 = pl2.steady_state(alpha=u_ss_alpha)
    m1 = pl2.metrics(x1, u1, d1)
    drho = m1["slurry_density_kg_m3"] - m0["slurry_density_kg_m3"]
    K = (m1["ppm_ClO4"] - m0["ppm_ClO4"]) / drho if abs(drho) > 1e-9 else np.nan
    tau = tau_s if tau_s is not None else plant.p["volume_m3"] / (d0["Q_stock"] + u0["Q_dil_cmd"] - u0["m_evap"] / p["rho_water_kg_m3"])
    theta = lab_delay_s + 0.5 * lab_period_s
    Kp, Ti, direct = simc_pi(K, tau, theta, tau_c=theta)
    lim = max(abs(drho) * 5.0, 0.5)
    return (LoopTuning(Kp=Kp, Ki=Kp / Ti, direct=direct),
            {"K_ppm_per_kg_m3": K, "tau_s": tau, "theta_s": theta, "drho_kg_m3": drho}, lim)


# === EXTENSIONES INTEGRADAS: LAZO ABIERTO/CERRADO, RUIDO, IMC Y ROUTH-HOPF ===
class ConsignaVariable:
    """Envoltorio de PlantController para cambios de consigna dinámicos (escalón/rampa)."""
    def __init__(self, ctrl, funcs):
        self.ctrl = ctrl
        self.funcs = funcs

    def __call__(self, t, x, plant):
        for k, f in self.funcs.items():
            self.ctrl.sp[k] = f(t)
        return self.ctrl(t, x, plant)


class Ruido:
    """Ruido coloreado de Ornstein-Uhlenbeck (AR(1)) para perturbaciones de planta."""
    SIG = {
        "Q_stock": ("rel", 0.015),
        "T_feed": ("add", 0.4),
        "C_stock": ("rel", 0.005),
        "CP_stock": ("rel", 0.02),
        "T_dil": ("add", 0.5),
        "valve_gain": ("rel", 0.005),
    }

    def __init__(self, seed=7, escala=1.0, H=30 * 3600, dt=5.0, tc=600.0):
        rng = np.random.default_rng(seed)
        n = int(H / dt) + 2
        phi = np.exp(-dt / tc)
        self.dt, self.n, self.s = dt, n, {}
        for k, (kind, sig) in self.SIG.items():
            e = rng.standard_normal(n)
            x = np.empty(n)
            x[0] = e[0]
            for i in range(1, n):
                x[i] = phi * x[i - 1] + np.sqrt(1.0 - phi**2) * e[i]
            self.s[k] = escala * sig * x

    def __call__(self, d, t):
        i = min(int(t / self.dt), self.n - 1)
        out = dict(d)
        for k, (kind, _) in self.SIG.items():
            base = out.get(k, 1.0) if k == "valve_gain" else out[k]
            out[k] = base * (1.0 + self.s[k][i]) if kind == "rel" else base + self.s[k][i]
        return out


def imc_pi(K, tau, theta, eps):
    """Sintonía PI por método IMC (Rivera et al., 1986; Seborg Cap. 12)."""
    Kp = (tau + 0.5 * theta) / (abs(K) * (eps + 0.5 * theta))
    Ti = tau + 0.5 * theta
    return LoopTuning(Kp=Kp, Ki=Kp / Ti, Kd=0.0, direct=bool(K < 0))


_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def iae(t_h, e):
    """Integral del error absoluto (IAE) en [unid * h]."""
    return float(_trapz(np.abs(e), t_h))


def metricas_consigna(t_s, y, y0, y1, t_step, band=0.05):
    """Calcula sobrepaso [%], tiempo de asentamiento [h] e IAE [unid*h]."""
    s = (y - y0) / (y1 - y0)
    m = t_s >= t_step
    over = 100.0 * max(s[m].max() - 1.0, 0.0)
    fuera = np.abs(s[m] - 1.0) > band
    tset = (t_s[m][fuera][-1] - t_step) / 3600.0 if fuera.any() else 0.0
    iae_val = float(_trapz(np.abs(y[m] - y1), t_s[m]) / 3600.0)
    return over, tset, iae_val


def sigma_ventana(df, col, t0_h=12.0):
    """Desviación estándar descontando deriva lineal en la ventana estacionaria final."""
    w = df[df.t_h >= t0_h]
    y = w[col].values
    t = w.t_h.values
    return float(np.std(y - np.polyval(np.polyfit(t, y, 1), t)))


def routh_hopf_analysis(Kp, tau_p, theta, tau_I):
    """Criterio de Routh-Hurwitz con Padé 1er orden y detección de Bifurcación de Hopf."""
    theta_eff = max(float(theta), 1e-3)
    a0 = 0.5 * tau_p * tau_I * theta_eff
    K_abs = abs(Kp)
    Kc_max_a1 = (1.0 / K_abs) * (2.0 * tau_p / theta_eff + 1.0)
    kc_grid = np.linspace(1e-4, Kc_max_a1 * 0.999, 10000)
    a1 = tau_I * (tau_p + 0.5 * theta_eff - 0.5 * kc_grid * K_abs * theta_eff)
    a2 = tau_I + kc_grid * K_abs * tau_I - 0.5 * kc_grid * K_abs * theta_eff
    a3 = kc_grid * K_abs
    A1 = (a1 * a2 - a0 * a3) / a1
    idx_hopf = np.where(A1 <= 0)[0]
    Kc_hopf = kc_grid[idx_hopf[0]] if len(idx_hopf) > 0 else Kc_max_a1
    a1_h = tau_I * (tau_p + 0.5 * theta_eff - 0.5 * Kc_hopf * K_abs * theta_eff)
    a3_h = Kc_hopf * K_abs
    omega_u = np.sqrt(max(a3_h / max(a1_h, 1e-9), 1e-9))
    tau_u = 2.0 * np.pi / omega_u
    Kc_overshoot = (1.0 / K_abs) * max((2.0 * tau_p / tau_I) - 1.0, 0.05)
    sgn = -1.0 if Kp < 0 else 1.0
    return {
        "Kc_overshoot": sgn * Kc_overshoot,
        "Kc_hopf": sgn * Kc_hopf,
        "omega_u": float(omega_u),
        "tau_u": float(tau_u),
        "Kc_ZN": sgn * (Kc_hopf / 2.2),
        "tau_I_ZN": float(tau_u / 1.2),
        "Kc_TL": sgn * (Kc_hopf / 3.2),
        "tau_I_TL": float(2.2 * tau_u),
    }


class CommercialPID:
    """PID industrial 2-DOF con derivada sobre PV, filtro derivativo, anti-windup y trazabilidad de P e I."""

    def __init__(self, Kc, tau_I, tau_D, sp, u_bias, u_min, u_max, alpha=0.1, beta=1.0):
        self.Kc, self.tau_I, self.tau_D = float(Kc), float(tau_I), float(tau_D)
        self.sp, self.u_bias = float(sp), float(u_bias)
        self.u_min, self.u_max = float(u_min), float(u_max)
        self.alpha, self.beta = float(alpha), float(beta)
        self.integral = -self.Kc * (self.beta * self.sp - self.sp)
        self.I = self.integral
        self.prev_pv = None
        self.d_filt = 0.0
        self.last_error = 0.0
        self.last_P = 0.0

    def step(self, pv, dt):
        error = self.sp - pv
        self.last_error = error
        P = self.Kc * (self.beta * self.sp - pv)
        self.last_P = P
        if self.prev_pv is None or self.tau_D <= 0 or dt <= 0:
            D = 0.0
        else:
            raw_d = -self.Kc * self.tau_D * (pv - self.prev_pv) / dt
            tau_f = max(self.alpha * self.tau_D, 1e-6)
            self.d_filt += (dt / (tau_f + dt)) * (raw_d - self.d_filt)
            D = self.d_filt
        tentative_I = self.integral + (self.Kc / max(self.tau_I, 1e-6)) * error * dt if self.tau_I > 0 else 0.0
        u_unsat = self.u_bias + P + tentative_I + D
        if self.u_min <= u_unsat <= self.u_max:
            self.integral = tentative_I
            self.I = tentative_I
        self.prev_pv = pv
        return float(np.clip(u_unsat, self.u_min, self.u_max))