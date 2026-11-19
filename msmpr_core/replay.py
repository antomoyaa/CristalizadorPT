"""
Módulo de reproducción contrafactual (disturbance replay) sobre campañas reales SQM (NPT-3).
Sincroniza sensores físicos (10 min) y química (2 h) y evalúa la estrategia sostenible
(F_303 = 0.50 t/h nominal + AIT 3105 + purga continua AIC 3255 + trim térmico T2).
"""
from pathlib import Path
import unicodedata
import numpy as np
import pandas as pd


def find_sqm_excel(keyword, search_bases=(Path("."), Path(".."))):
    """Busca archivos Excel en data/ tolerando normalización NFC/NFD de acentos en macOS."""
    def _norm(s):
        return unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode("ascii").lower()
    kw = _norm(keyword)
    for base in search_bases:
        b = Path(base)
        for folder in (b / "data", b):
            if folder.exists():
                for p in sorted(folder.glob("*.xlsx")):
                    if kw in _norm(p.name) and not p.name.startswith("~$"):
                        return p
    return None


def cargar_ventana_real_sqm(
    t_inicio="2026-11-10 22:00:00",
    t_fin="2026-11-15 05:00:00",
    usar_sensores_4d=True,
    ruta_fis=None,
    ruta_v2=None,
):
    """Sincroniza sensores de planta (10 min) y química (2 h) desde los Excel de SQM."""
    ruta_fis = Path(ruta_fis) if ruta_fis is not None else find_sqm_excel("analisis fisicos")
    ruta_v2 = Path(ruta_v2) if ruta_v2 is not None else find_sqm_excel("base de datos_v2")
    if ruta_v2 is None or not ruta_v2.exists():
        raise FileNotFoundError("Copia 'Base de datos_v2.xlsx' en la carpeta data/ del repositorio.")

    xl = pd.ExcelFile(ruta_v2)

    def leer_quim(hoja, cols_map):
        d = xl.parse(hoja)
        ts = pd.to_datetime(
            pd.to_datetime(d["Fecha"]).dt.strftime("%Y-%m-%d") + " " + d["Hora"].astype(str),
            errors="coerce",
        )
        d = d.assign(TS=ts).dropna(subset=["TS"]).set_index("TS").sort_index()
        d = d[~d.index.duplicated(keep="last")].loc[t_inicio:t_fin]
        out = pd.DataFrame(index=d.index)
        for orig, nuevo in cols_map.items():
            out[nuevo] = pd.to_numeric(d[orig], errors="coerce")
        return out

    if usar_sensores_4d and ruta_fis is not None and ruta_fis.exists():
        sens = pd.read_excel(ruta_fis, sheet_name="Data sensores", header=1)
        sens.columns = [str(c).strip() for c in sens.columns]
        sens["TS"] = pd.to_datetime(sens["Unnamed: 1"], errors="coerce")
        sens = sens.dropna(subset=["TS"]).set_index("TS").sort_index().loc[t_inicio:t_fin]

        df = pd.DataFrame(index=sens.index)
        df["Q_SSKd_m3_h"]  = pd.to_numeric(sens["SSKd"], errors="coerce").interpolate().bfill().ffill()
        df["Q_dil_real"]   = (pd.to_numeric(sens["Agua de dilución"], errors="coerce") / 1000.0).interpolate().bfill().ffill()
        df["Q_purga_real"] = pd.to_numeric(sens["Flujo Ml purga"], errors="coerce").clip(lower=0.0).interpolate().bfill().ffill()
        df["T_A3110_real"] = pd.to_numeric(sens["Temperatura CRY A3110"], errors="coerce").interpolate().bfill().ffill()
        df["T_A3210_real"] = pd.to_numeric(sens["Temperatura CRY A3210"], errors="coerce").interpolate().bfill().ffill()
        df["Nivel_A3110"]  = pd.to_numeric(sens["Nivel A3110"], errors="coerce").interpolate().bfill().ffill()
        df["LV3101_real"]  = pd.to_numeric(sens["Valvula LV 3101"], errors="coerce").interpolate().bfill().ffill()
    else:
        d_ssk_raw = xl.parse("SSK")
        df = pd.DataFrame({
            "TS": pd.to_datetime(d_ssk_raw["Unnamed: 13"], errors="coerce").dt.round("10min"),
            "Q_SSKd_m3_h": pd.to_numeric(d_ssk_raw["Unnamed: 15"], errors="coerce"),
            "Q_dil_real": pd.to_numeric(d_ssk_raw["Unnamed: 16"], errors="coerce") / 1000.0,
        }).dropna(subset=["TS"]).set_index("TS").sort_index().loc[t_inicio:t_fin]
        df = df[~df.index.duplicated(keep="last")]
        df["Q_SSKd_m3_h"]  = df["Q_SSKd_m3_h"].clip(lower=120.0).interpolate().bfill().ffill()
        df["Q_dil_real"]   = df["Q_dil_real"].clip(lower=0.0, upper=35.0).interpolate().bfill().ffill()
        df["Q_purga_real"] = 20.43
        df["T_A3110_real"] = 57.48
        df["T_A3210_real"] = 43.04

    ssk_q   = leer_quim("SSK", {"ClO4 IC": "ClO4_SSK", "Mg": "Mg_SSK", "K": "K_SSK"})
    mlk_q   = leer_quim("Purga-Recirculación", {"ClO4 IC": "ClO4_MLK_real", "Mg": "Mg_MLK_real"})
    ref_q   = leer_quim("Cristal Refinado", {"ClO4": "ClO4_Refi_real"})
    mlref_q = leer_quim("ML Refi", {"ClO4 IC": "ClO4_MLRefi_real"})

    df = df.join([ssk_q, mlk_q, ref_q, mlref_q], how="left")
    for col in ["ClO4_SSK", "Mg_SSK", "K_SSK", "ClO4_MLK_real", "Mg_MLK_real", "ClO4_MLRefi_real"]:
        df[col] = df[col].interpolate("time").ffill().bfill()

    df["t_h"] = (df.index - df.index[0]).total_seconds() / 3600.0
    return df


class PIControllerDiscrete:
    """Controlador PI discreto con anti-windup por retrocálculo."""
    def __init__(self, Kc, tau_I_h, sp, u_bias, u_min, u_max):
        self.Kc = float(Kc)
        self.tau_I = float(tau_I_h)
        self.sp = float(sp)
        self.u_bias = float(u_bias)
        self.u_min = float(u_min)
        self.u_max = float(u_max)
        self.I = 0.0

    def step(self, pv, dt_h):
        e = self.sp - float(pv)
        u_unsat = self.u_bias + self.Kc * e + self.I
        u_sat = float(np.clip(u_unsat, self.u_min, self.u_max))
        aw = (u_sat - u_unsat) / max(self.tau_I, 0.1)
        self.I += (self.Kc / self.tau_I * e + aw) * float(dt_h)
        return u_sat


def simular_replay_npt3(df_in, mode="baseline_sqm", theta_ait_h=0.25, lab_sampling_h=0.0, mg_offset=0.0):
    """
    Simula el tren acoplado Área 2000 -> Área 3000 -> Área 6000 sobre una ventana de planta.
    Modos: 'baseline_sqm', 'dil_05_only', 'proposed_integrated'.
    """
    dt_h = 10.0 / 60.0
    n = len(df_in)

    C_MLK    = float(df_in["ClO4_MLK_real"].iloc[0])
    Mg_MLK   = float(df_in["Mg_MLK_real"].iloc[0]) + float(mg_offset)
    C_MLRefi = float(df_in["ClO4_MLRefi_real"].iloc[0])
    T2       = float(df_in["T_A3210_real"].iloc[0])
    Q_dil    = float(df_in["Q_dil_real"].iloc[0]) if mode == "baseline_sqm" else 0.50

    a_valve = 1.0 - np.exp(-dt_h / 0.05)
    a_T2    = 1.0 - np.exp(-dt_h / 0.50)
    a_B3520 = 1.0 - np.exp(-dt_h / 2.20)
    a_refi  = 1.0 - np.exp(-dt_h / 3.50)

    d_steps = max(int(round(theta_ait_h / dt_h)), 1)
    ait_buf = [float(df_in["ClO4_SSK"].iloc[0])] * d_steps

    pid_aic = PIControllerDiscrete(Kc=-35.0, tau_I_h=1.0, sp=4.25, u_bias=28.0, u_min=18.0, u_max=65.0)
    pid_t2  = PIControllerDiscrete(Kc=-2.2,  tau_I_h=1.2, sp=4.30, u_bias=43.0, u_min=43.0, u_max=44.3)

    C_MLK_meas, Mg_MLK_meas = C_MLK, Mg_MLK
    rows = []

    for k in range(n):
        t_h      = float(df_in["t_h"].iloc[k])
        Q_sskd   = max(float(df_in["Q_SSKd_m3_h"].iloc[k]), 50.0)
        C_ssk    = float(df_in["ClO4_SSK"].iloc[k])
        T2_real  = float(df_in["T_A3210_real"].iloc[k])
        Q_dil_r  = float(df_in["Q_dil_real"].iloc[k])
        Q_purg_r = float(df_in["Q_purga_real"].iloc[k])

        ait_buf.append(C_ssk)
        C_ait = ait_buf.pop(0)

        if lab_sampling_h > 0:
            if k % max(int(round(lab_sampling_h / dt_h)), 1) == 0:
                C_MLK_meas, Mg_MLK_meas = C_MLK, Mg_MLK
        else:
            C_MLK_meas, Mg_MLK_meas = C_MLK, Mg_MLK

        C_sat = (4.92 + 0.14 * (T2 - 43.0)) * np.exp(-0.022 * (Mg_MLK - 30.5))

        if mode == "baseline_sqm":
            Q_dil, Q_purga, T2 = Q_dil_r, Q_purg_r, T2_real
        elif mode == "dil_05_only":
            R0 = 0.50 / 243.2
            R_cmd = R0 + 0.006 * max(C_ait - 4.40, 0.0)
            Q_dil += a_valve * (float(np.clip(Q_sskd * R_cmd, 0.20, 2.50)) - Q_dil)
            Q_purga, T2 = Q_purg_r, T2_real
        else:
            R0 = 0.50 / 243.2
            R_cmd = R0 + 0.005 * max(C_ait - 4.60, 0.0) + 0.0005 * max(Mg_MLK_meas - 31.0, 0.0)
            Q_dil += a_valve * (float(np.clip(Q_sskd * R_cmd, 0.20, 2.50)) - Q_dil)

            C_sat_est = (4.92 + 0.14 * (T2 - 43.0)) * np.exp(-0.022 * (Mg_MLK_meas - 30.5))
            pid_aic.sp = min(4.45, C_sat_est - 0.35)
            ff_purge = 12.0 * max(C_ait - 4.80, 0.0) + 2.5 * max(Mg_MLK_meas - 31.2, 0.0)
            Q_purga = float(np.clip(pid_aic.step(C_MLK_meas, dt_h) + ff_purge, 18.0, 65.0))

            pid_t2.sp = min(4.50, C_sat_est - 0.28)
            T2 += a_T2 * (pid_t2.step(C_MLK_meas, dt_h) - T2)

        dil_ratio = (Q_sskd + 20.0) / (Q_sskd + Q_dil)
        purge_factor = (20.4 + 95.0) / (Q_purga + 95.0)

        if mode == "baseline_sqm":
            C_MLK_eq  = float(df_in["ClO4_MLK_real"].iloc[k])
            Mg_MLK_eq = float(df_in["Mg_MLK_real"].iloc[k]) + float(mg_offset)
        else:
            C_MLK_eq  = float(df_in["ClO4_MLK_real"].iloc[k]) * dil_ratio * purge_factor
            Mg_MLK_eq = (float(df_in["Mg_MLK_real"].iloc[k]) + float(mg_offset)) * dil_ratio * (0.35 + 0.65 * purge_factor)

        C_MLK  += a_B3520 * (C_MLK_eq - C_MLK)
        Mg_MLK += a_B3520 * (Mg_MLK_eq - Mg_MLK)

        eta_ratio = max(1.0 + 0.025 * (Mg_MLK - 30.5), 0.78)
        w_moist   = 0.0286 * (eta_ratio ** 0.30)

        excess = max(C_MLK - C_sat, 0.0)
        coprecip_crudo_ppm = 2600.0 * (excess ** 1.15)
        brine_crudo_ppm    = (w_moist * C_MLK * 1000.0 / 1.40) * 1.45
        ClO4_crudo_ppm     = brine_crudo_ppm + coprecip_crudo_ppm

        C_MLRefi += a_refi * ((0.90 + 0.0028 * ClO4_crudo_ppm) - C_MLRefi)
        ppm_refi_raw = (
            11.0
            + 0.095 * brine_crudo_ppm * (C_MLRefi / 1.50) ** 0.35
            + 0.058 * coprecip_crudo_ppm * (C_MLRefi / 1.50) ** 0.45
        )
        ppm_refi = max(20.0, ppm_refi_raw)

        dKNO3_net = 0.38 * (Q_dil_r - Q_dil) - 0.055 * (Q_purga - Q_purg_r) - 0.85 * (T2 - T2_real)

        rows.append({
            "t_h": t_h,
            "Q_dil_t_h": Q_dil,
            "Q_purga_m3_h": Q_purga,
            "T2_C": T2,
            "C_MLK_g_L": C_MLK,
            "C_sat_g_L": C_sat,
            "Mg_MLK_g_L": Mg_MLK,
            "ClO4_refi_ppm": ppm_refi,
            "dKNO3_net_t_h": dKNO3_net,
        })

    return pd.DataFrame(rows, index=df_in.index)
