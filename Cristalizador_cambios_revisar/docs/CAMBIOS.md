# Cambios tras incorporar el trabajo de calibración del equipo

## 1. Bases alineadas (`msmpr_core/params.py`)
* `make_params("calibracion")` reproduce EXACTAMENTE la `BASE` del notebook de Pareto (que usa `calibracion_do_crystal.ipynb`).
  Una prueba la compara con el notebook (variable `PARETO_NB`) para detectar desfases.
* `make_params("literatura")` (por defecto) conserva los valores con fuente (Perry, IAPWS).
* `diff_bases()` lista qué difiere y cuáles diferencias condicionan los parámetros calibrados:

| Parámetro | literatura | calibracion | condiciona kG/kB/Kd |
|---|---|---|---|
| `solubility_slope_K_inv` | 0.0312 | 0.010 | sí |
| `rho_liquid_kg_m3` | 1470 | 1040 | sí |
| `rho_crystal_kg_m3` | 2110 | 2100 | sí |
| `Kd` | 0.002 | 0.01 | no (es lo que se calibra) |
| `feed_stock_kg_m3.ClO4` | 5.56 | 3.8e-3 | no (base de ppm sin resolver) |
| `cake_moisture_fraction` | 0.03 | 0 | no |

* `apply_calibration(p, calibracion, equipo)` inserta kG, kB, Kd (y g, b) calibrados; **rechaza** si la base de `p` no es la de
  calibración (salvo `force=True`) y anula la humedad de torta, porque el Kd calibrado ya absorbe ese ClO4 (si no, se cuenta dos veces).

## 2. Conversor de datos (`msmpr_core/datos.py`)
Excel ancho (`datos_cristalizador.xlsx`) → formato largo del calibrador. Validado con las propias funciones del notebook
(`tidy_sheet`, `normalize_data`, `validate_campaign`, `simulate_campaign`). Ejemplo: `data/EJEMPLO_datos_largo_SINTETICO_A3110.xlsx`.
Obliga a declarar equipo y origen, y a pasar volumen y alimentación (no vienen en el Excel).

## 3. Qué quedó como estaba
`plant.py`, `control.py` y los notebooks 04–05 no cambian: no tienen calibración que los reemplace. Siguen sin calibrar.
`fit.py`/`model.py` se mantienen solo como chequeo estacionario rápido; se eliminó `data/experimentos_plantilla.csv` (formato redundante).

## 4. Control por etapas (nuevo)
* `PlantController(..., enabled=("flow",), flow_sp=..., noise_rel={"flow":0.01})`: etapa 1 (solo FIRCAL 3101).
* `plant.rhs` acepta la perturbación opcional `valve_gain` (caída de presión del agua de dilución).
* Notebook `06_Control_por_etapas_dilucion` y `docs/ROADMAP_control_por_etapas.md`.

## 5. Modelado del P&ID en torno al agua de dilución
* `DIRCAHL 3101` (no OIRCAHL) mide densidad. Se renombró en código y documentos.
* `plant.py`: nuevo trazador inerte `S` (demás sales disueltas) y `liquor_density`: el agua de dilución baja la densidad del licor
  (antes constante). La densidad usa un ancla (primer estacionario) y valen sus diferencias, no su valor absoluto.
* `PlantController`: modo `ratio` (FFIC 3102: consigna = R × caudal de alimentación; el maestro puede ajustar R).
* Notebook 05 reescrito: seis estrategias de dilución (S0–S5) frente a las mismas perturbaciones.

## 6. Estructura confirmada del agua de dilución (AIT 3105 → FFIC 3102 → FIRC 3103)
* El agua de dilución la controla `FIRC 3103`, con consigna de `FFIC 3102` (razón) y señal de `AIT 3105` (ClO4 de entrada).
  La cascada `DIRCAHL 3101 → FIRCAL 3101` queda sin modelar (actuador por identificar).
* `control.DelayedSensor` y `PlantController(..., enabled=(..., "ratio", "ait"), ait={"cp_true":..., "delay_s":..., "gain":...})`.
* `plant.steady_state(Q_dil=...)` para evaluar el estacionario con agua adicional.
* Notebook 05 reescrito: tabla estática de dilución necesaria, FIRC 3103 solo, FFIC 3102 y barrido del retardo del analizador.



**NUEVO**
# Registro de Cambios — Integración de Control Dinámico y Planta NPT-3

## [v2.0.0] - Integración Lazo Abierto/Cerrado, Ruido, Routh-Hopf y Tren NPT-3

### Módulos Core (msmpr_core)
- **msmpr_core/control.py**: Se incorporan `ConsignaVariable`, `Ruido` (proceso coloreado de Ornstein-Uhlenbeck), `imc_pi`, `metricas_consigna`, `sigma_ventana`, `routh_hopf_analysis`, `DelayedSensor` y `CommercialPID`.
- **msmpr_core/plant.py**: Se añade `simulate_plant_npt3_calibrated` acoplando el Área 2000 (B2510), Área 3000 (A3110, A3210 y B3520 con purga AIC 3255) y Área 6000 de Refinación ($w_{\text{KNO}_3,\text{in}} = 0.95678$, KPI $< 40\text{ ppm}$ de $\text{ClO}_4^-$).

### Notebooks Integrados (notebooks/)
- **notebooks/05_Control_PID_dinamico.ipynb**: Unifica `05_Control_PID_Lazo_AbiertoCerrado.ipynb` y `05_Control_PID_Perturbaciones.ipynb` con validación FOPDT ($\pm 10\%$), 12 escenarios de perturbación ($\pm$), ruido coloreado, seguimiento 2-DOF ($\beta$, $\varepsilon$), Frente de Pareto y estabilidad Routh-Hurwitz/Hopf.
- **notebooks/06_Control_por_etapas_dilucion.ipynb**: Evalúa la cascada `AIT 3105 -> FFIC 3102 -> FIRC 3103` frente a razón fija sin analizador y el impacto del retardo $\theta_{\text{AIT}}$ junto con la purga `AIC 3255 -> FV 3522`.

### Pruebas Automatizadas (tests/)
- **tests/test_control_integrado.py**: Suite en `pytest` que verifica reposo estacionario, balances integradores, offset nulo a 24 h, precarga sin salto para $\beta \neq 1$, límites de Hopf y cumplimiento de $< 40\text{ ppm}$ en el producto final.

## [v2.1.0] - Revisión: vacío, trim de AIT 3105, lazo de nivel
- `plant.VAC`: consignas y ganancias de vacío coherentes con DWSIM/Antoine (antes -0.932 bar para 55 C y K = 375 K/bar).
- `plant.K_TRIM_AIT = 5.2e-5` unificado (antes 5.2e-4 en dos funciones: saturaba FIRC 3103).
- **Lazo de nivel** en `simulate_plant_npt3_calibrated`: balance de inventario de A3110 y A3210 (`LEVEL`), LIC manda el caudal de descarga (F_6, F_11) y el maestro de densidad (DIC) ajusta la consigna de nivel (cascada DIC -> LIC). Nuevas columnas: `L1, L2, M1_t, M2_t, L1_sp, L2_sp, dMass_A3110_t_h, dMass_A3210_t_h`.
- `params.EMPIRICAL_COEFFS`: coeficientes empíricos etiquetados como placeholder. `.gitignore` agregado.
- Pruebas nuevas: `tests/test_revision_vacio_trim_balance.py`.

## [v2.2.0] - Nivel 60 %, vínculo vacío -> T -> granulometría
- Nivel de A3110/A3210 con consigna fija de 60 % (`LEVEL["level_sp_pct"]`); la cascada densidad -> nivel queda opcional (`density_cascade=True`). Columnas `level1_pct`, `level2_pct`.
- `plant.vac_T(P, etapa)` / `vac_P(T, etapa)`: T = Tsat_agua(P) + elevación del punto de ebullición (ajustada a DWSIM: 11,59 kPa -> 55 C; 5,43 kPa -> 40 C). K1, K2 se calculan de esta función (≈161 y ≈320 K/bar; los 133/255 de v2.1.0 ignoraban la elevación del punto de ebullición).
- `plant.granulometria_estacionaria(P1, P2, ...)` y D50 dinámico en `simulate_plant_npt3_calibrated`: P -> T -> sobresaturación (MSMPR `dpbe`) -> D50, anclado a los D50 de planta (541 / 514 um, reporte de seguimiento, agosto 2026). Se reemplaza el término empírico `sigma_burst`. Columnas `d50_1_um`, `d50_2_um`, `sigma1`, `sigma2`.
- OJO: la sensibilidad usa la cinética placeholder del perfil `planta_sqm`; con ella +-0,02 bar mueve el D50 ~3-5 %.

## [v2.3.0] - Sintonías de T identificadas sobre la simulación
- `plant.identificar_Kcu()`: Kcu y Pu del lazo de T (acción P pura, válvula sin saturar) por bisección sobre `simulate_plant_npt3_calibrated`. Resultado: Kcu = 0,01241 bar/°C, Pu = 9,70 min (`plant.KCU_SIM`). La tabla FOPDT del notebook 05 (Kp = 375, tau_p = 5,5, theta = 1) daba Kcu = 0,028 y Pu = 3,3 min: no representa la planta simulada.
- Modos: `Hopf_limit_cycle` = Kcu (PI, tauI = Pu/1,2; no decae, razón de decaimiento 1,74); `ZN_oscillatory` = 0,45 Kcu, tauI = Pu/1,2 (razón de decaimiento 0,31 ≈ 1/4; antes oscilaba sin decaer); `TL_damped` sin cambios (Kc = 0,42 Kcu).
- Nuevos kwargs de `simulate_plant_npt3_calibrated`: `Kc_T_override`, `tauI_T_override`, `tauD_T_override`, `vac_limits_override`. Nuevas utilidades: `razon_decaimiento`, `amplitud_oscilacion`.
- Notebook 05: tabla 1b (sintonías identificadas); la figura 3 usa simulaciones.

## [v2.4.0] - Procedencia de datos y documentación de bases
- `params`: `Cstar_ref_kg_m3 = 177,94` y `solubility_slope_K_inv = 0,010` del perfil `planta_sqm` pasan de "interno" a **placeholder (dato sintético)**; `make_params("planta_sqm")` avisa una vez por sesión. Entradas `planta_sqm.*` en `PROVENANCE`.
- `docs/BASES_DE_PARAMETROS.md`: tabla de las tres bases y qué no se puede mezclar.

## [v2.5.0] - Fusión con los commits "revisandou" de PID
- `plant.py`: el bloque agregado en "revisandou" duplicaba `simulate_plant_npt3_calibrated` y `simular_recirculacion_y_magnesio` (la segunda definición tapaba a la primera, incluidas las mejoras de Kcu/nivel/granulometría). Queda **una sola** definición de cada una: el tren integrado con Mg2+ del equipo + sintonías identificadas.
- Corrige `KeyError: 'w_Mg_MLK_v6_ppm'` en `simular_recirculacion_y_magnesio` (faltaban `w_Mg_MLK_v6_ppm` y `w_ClO4_MLK_v6_ppm` en `PFD_MG_CLO4`).
- Notebooks 04-06: la primera celda que importa `msmpr_core` fija la ruta y el directorio de trabajo (VS Code puede abrir el kernel en la raíz del repo); se agrega `display` de IPython; la celda de conciliación Mg/ClO4 del 06 queda ejecutable.
- Pruebas: sin funciones duplicadas, claves de `PFD_MG_CLO4`, los tres modos de recirculación, columnas de Mg.

## [v2.6.0] - Tope de agua, supervisor de T2 y análisis térmico
- `plant.F303_MAX_AMBIENTAL = 36 t/h`: ninguna simulación pide más agua de dilución (antes topes de 42, 50 y 52 t/h). `F_303_target` por defecto ya no es 36,8.
- Saturación de KClO4 dependiente de T (`sat_factor_T`, `c_sat_clo4_gL`): dato de SQM bajo 40 °C; sobre 40 °C razón del agua pura (**placeholder**).
- Supervisor de T2 (`trim_T2`) en `simular_recirculacion_y_magnesio` y `simular_planta_agua_dilucion_magnesio`: sube la consigna de T2 lo mínimo para mantener el margen de saturación. Pérdida de KNO3: 1,82 t/h por °C (Entrega 2, Tabla 8).
- Nuevas funciones: `comparar_estrategias_excursion`, `barrido_agua_optimo`, `analisis_temperatura`, `calor_enfriar_agua_dilucion_kW`; nuevos argumentos de análisis (`T_in_extra`, `T_agua_dilucion`, `ff_temperatura`, `T1_sp`, `T2_sp`, ...).
- Arranque del reciclo inteligente sin recirculación (α = 1): elimina el pico espurio de ~67 ppm de los primeros 20 min.
- Notebook 06: B6 recalculada con 4 estrategias (la propuesta cumple; las demás no, con veredicto calculado); B7 óptimo de agua; B8 análisis de temperatura; tablas sin veredictos escritos a mano ni el ajuste `+ 5.8` de KNO3; Pareto sintético reemplazado por simulaciones. Se eliminan las figuras PPT_Fig5/6 antiguas (curvas escritas a mano).
