# Control por etapas: de la corriente de dilución al sistema completo

Regla: **un lazo nuevo por etapa; antes de cerrarlo, reidentificar el proceso con los lazos anteriores ya cerrados**
(la ganancia depende de qué otros lazos están cerrados: ver notebook 06, etapa 0). Cada etapa termina con una prueba en `tests/`.

| Etapa | Se agrega | Qué se hace | Criterio de salida (se vuelve prueba) | Estado |
|---|---|---|---|---|
| 0 | — | Estacionario `rhs=0`, sin deriva en lazo abierto; escalones de identificación | cierre de balance < 1e-8; deriva < 1e-9 | hecho |
| 1 | FIRC 3103 (caudal de dilución) | PI por SIMC con tau_c >= 2·dt; consigna externa; ruido de medición; perturbación de ganancia de la válvula; límites | sigue la consigna (<2 %), rechaza la perturbación (<3 %), no satura, no oscila | hecho (`test_stage1_*`, notebook 06) |
| 2 | FFIC 3102 (razón dilución/alimentación) | La consigna de FIRC 3103 pasa a ser R × caudal de alimentación; la dilución sigue las variaciones de alimentación | la razón dilución/alimentación se mantiene dentro de 2 % ante +12 % de alimentación | hecho (`test_ratio_mode_*`, notebook 05) |
| 2b | AIT 3105 (ClO4 de entrada) | Ley de razón sobre el ClO4 medido, con ganancia y límites; barrer el retardo del analizador (aún sin definir) | ppm del pulso y recuperación mínima aceptables para la ganancia y el retardo definidos | hecho en simulación (`test_ait_*`); falta definir retardo, ganancia y capacidad de la válvula |
| 3 | nivel (LIRCSAHL 3101) | La dilución agrega volumen: cerrar V y reidentificar el maestro | V dentro de ±1 % ante el escalón de dilución | implementado en el modelo; falta ensayo por etapas |
| 4 | T/vacío (PIRCSAHL 3101) | Calcular la matriz de ganancias y la RGA de los pares candidatos antes de cerrar | pareamiento justificado; sin interacción destructiva | pendiente |
| 5 | ppm de ClO4 (propuesto) y prealimentación | Retardo de laboratorio (muestreo 2 h, entrega 1 h); prealimentación con Q_stock, T_feed y ClO4 de alimentación | pico de ppm menor sin saturar la dilución | modelo del lazo hecho; prealimentación pendiente |
| 6 | opcional | MPC/NMPC (do-mpc en `examples/two_stage_msmpr`) con el modelo calibrado | mejora frente a PID con las mismas perturbaciones | pendiente |

## Al activar un lazo en la planta real o en simulación
* Cierre en secuencia con seguimiento (bumpless): `PID.update(..., track=u_manual)`.
* Periodo de muestreo: un lazo muestreado cada dt con tau_c < dt es inestable. El caudal necesita dt de pocos segundos; el nivel y la temperatura toleran decenas.
* Anti-windup por retrocálculo ya incluido; probar siempre con la válvula saturada.
* Los valores del modelo siguen siendo placeholders (solubilidad, cinética, Kd): las conclusiones son de estructura, no de planta.

## Nota
La cascada `DIRCAHL 3101 -> FIRCAL 3101` del plano **no** se modela como agua de dilución (el equipo confirmó que la dilución la controla
`FIRC 3103`). El modo `rho` de `PlantController` (maestro de densidad sobre la dilución) queda como alternativa de estudio.
