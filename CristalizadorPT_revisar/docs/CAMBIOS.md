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
