# Control PID alineado con el P&ID del cristalizador A3110

Fuente: plano GEA/ABB 172164 (cajetín: hoja 10), "Crystallization", leído a mano desde `A3110_PID.pdf` (imagen sin texto).
**No se revisó `A3210_PID.pdf`**: los tags 32xx del bloque de código anterior (p. ej. `TIRCSAHL 3202`) no están verificados.

## Instrumentos legibles en A3110
| Tag | Lectura probable | Uso en el modelo |
|---|---|---|
| `PIRCSAHL 3101` | presión/vacío de la cabeza del cristalizador (con interruptores de alarma) | Lazo A: en un cristalizador al vacío la T queda fijada por la presión; el modelo lo representa como lazo de T con vapor extraído `m_evap` |
| `LIRCSAHL 3101` | nivel del cristalizador | Lazo B: volumen V -> caudal de descarga `Q_out` |
| `DIRCAHL 3101 -> FIRCAL 3101` | `DIRCAHL` mide la **densidad** (confirmado) y fija la consigna de `FIRCAL 3101`; ese actuador **no** es el agua de dilución (por identificar) | solo se mide la densidad; el lazo no está modelado |
| `AIT 3105 -> FFIC 3102 -> FIRC 3103` | analizador de ClO4 en la entrada -> razón dilución/alimentación -> caudal de agua de dilución (confirmado por el equipo) | modelado: anticipación sobre el ClO4 de alimentación; **retardo del analizador sin definir** (parámetro libre) |
| `TIR 3101`, `TIRAL 3102` | indicadores de temperatura (descarga hacia B3113) | solo medición |
| `SIC 3151`, `SIC 3102`, `EIRAH 3101` | velocidad del agitador P3111 y de la bomba de circulación; corriente del agitador | no modelados |
| `HS 3101` | válvula de tres vías en la línea de circulación | no modelada |

## Supuestos que hay que confirmar con la leyenda y el proceso
1. Dónde mide `AIT 3105`: el modelo supone la alimentación *antes* de mezclar el agua (anticipación). Si mide después de la dilución, sería realimentación.
2. Que la T se regule vía presión (`PIRCSAHL`); no hay `TIRC` en la hoja leída.
3. El lazo D (ppm de ClO4 -> consigna de densidad) **no existe en el P&ID**: es una propuesta. Necesita un análisis
   de laboratorio (CI) con retardo; el código lo modela con muestreo cada 2 h y entrega 1 h después.

## Qué enseñó el modelo (con los valores actuales, que son placeholders)
- Con volumen fijo, la **purga queda fijada por el balance de agua** (alimentación − evaporación). La fracción `alpha`
  solo reparte `Q_out` entre purga y reciclo: cambia el tiempo de residencia efectivo y el tamaño, no la recuperación ni el ClO4 de equilibrio.
- Las palancas del ClO4 son la dilución, la evaporación (T) y el ClO4 de alimentación.
- La densidad de pulpa sí tiene señal frente a la dilución porque el agua baja la densidad del *licor* (el cambio por cristales es mucho menor). Esto requiere modelar las demás sales disueltas como trazador inerte (`S` en `plant.py`); con densidad de licor constante la señal parecía despreciable, lo que era un artefacto del modelo.
- Compensar el ClO4 de alimentación con agua es caro: mantener la mezcla en su concentración nominal exige varias veces la dilución nominal (6x para +25 % de ClO4 con estos placeholders) y la recuperación de KNO3 cae hacia cero. La capacidad de la válvula de FIRC 3103 y la ganancia de FFIC 3102 son decisiones de diseño.

## Sintonía
`control.autotune` identifica por escalón (FOPDT; integrador para el nivel) y aplica SIMC. Es un punto de partida:
la planta real tiene tiempos muertos de transporte, ruido y no linealidades que este modelo no incluye.
Referencias de método: apuntes CBE30338 (J. Kantor; código MIT, texto CC BY-NC-ND) y Lima et al. (NMPC experimental
de cristalización por lotes, MIT) para el uso de sensores blandos (IR/FBRM) como PV de lazos de calidad.

---

## 3. Transición Ambiental a Dilución Sostenible ($0{,}5\text{ t/h}$) y Validación sobre Campaña Real SQM

El análisis de los registros operacionales reales de SQM (`10–15 nov 2026`, sensores cada $10\text{ min}$ y química cada $2\text{ h}$) demostró que la planta opera actualmente con un promedio de **$19{,}96\text{ t/h}$** de agua de dilución y una purga **on/off** ($0\text{–}36{,}5\text{ m}^3/\text{h}$, media $20{,}4\text{ m}^3/\text{h}$), presentando aun así **$26{,}8\%$ de muestras fuera de especificación ($> 40\text{ ppm}$)** en Área 6000 (pico de $99\text{ ppm}$), gatilladas principalmente por alzas de **$\text{Mg}^{2+}$ en licor madre** ($r = +0{,}523$, efecto *salting-out* y mayor viscosidad/retención en centrífugas `CF 3550`)[cite: 8].

### 3.1 Estrategia Integrada Propuesta (Modificaciones P&ID en Rojo)
* **Prealimentación sostenible (`AIT 3105 → FFIC 3102 → FIRC 3103`):** Opera con consigna base de **$0{,}50\text{ t/h}$** (rango acotado $[0{,}20,\; 2{,}50]\text{ t/h}$), compensación *Lead-Lag* y selector *Override* de sobresaturación ($\sigma_{\min} = 0{,}02$).
* **Purga continua compensada por $\text{Mg}^{2+}$ (`AIC 3255 → FV 3522`):** Reemplaza la conmutación *on/off* actual por modulación continua ($18\text{–}65\text{ m}^3/\text{h}$) con *setpoint* dinámico bajo el límite de saturación $C_{\text{sat}}(T_2, \text{Mg}^{2+})$, absorbiendo los picos de perclorato sin malgastar agua fresca.
* **Trim térmico suave en `A3210` (`TIR 3201 → PIRCAL 3201`):** Ajuste acotado en $[43{,}0,\; 44{,}3]\text{ }^\circ\text{C}$ para ampliar el margen a co-precipitación de $\text{KClO}_4$ ante eventos extremos.

### 3.2 Resultados del *Replay* Contrafactual (`msmpr_core/replay.py` y Notebook `06`)

| Escenario Evaluado sobre Datos Reales SQM | Agua Media [t/h] | Purga Media [$\text{m}^3/\text{h}$] | $\text{ClO}_4^-$ Refinado Medio [ppm] | $\text{ClO}_4^-$ Refinado Máximo [ppm] | Muestras $> 40\text{ ppm}$ [%] | $\Delta\text{KNO}_3$ Neto [t/h] |
| :--- | :-: | :-: | :-: | :-: | :-: | :-: |
| **Campaña 4 días — Planta Real SQM ($n=56$)** | $19{,}96$ | $20{,}4$ *(on/off)* | $34{,}4$ | $99{,}0$ | $26{,}8\%$ | Base ($0{,}00$) |
| **Campaña 4 días — Solo bajar agua a $0{,}5\text{ t/h}$ SIN `AIC 3255`** | $1{,}28$ | $20{,}4$ *(on/off)* | $226{,}1$ | $536{,}5$ | $89{,}8\%$ | $+7{,}10$ *(fuera de norma)* |
| **Campaña 4 días — Propuesta Integrada (En línea)** | **$0{,}98$** | **$38{,}2$** *(continua)* | **$22{,}8$** | **$23{,}3$** | **$0{,}0\%$** | **$+6{,}23$** |
| **Campaña 4 días — Propuesta Integrada (Lab c/2 h)** | **$0{,}98$** | **$38{,}4$** *(continua)* | **$22{,}8$** | **$23{,}7$** | **$0{,}0\%$** | **$+6{,}13$** |
| **Robustez 30 días — Planta Real SQM ($n=393$)** | $20{,}44$ | $20{,}4$ | $31{,}3$ | $170{,}0$ | $16{,}0\%$ | Base ($0{,}00$) |
| **Robustez 30 días — Propuesta Integrada (Lab c/2 h)** | **$0{,}83$** | **$33{,}4$** | **$22{,}6$** | **$27{,}4$** | **$0{,}0\%$** | **$+6{,}64$** |
| **Estrés Extremo $\text{Mg}^{2+}$ ($43\text{–}46\text{ g/L}$, Lab c/2 h)** | $1{,}85$ | $61{,}3$ | $21{,}6$ | $59{,}5$ | $0{,}2\%$ | $+4{,}53$ |

> **Validación de robustez y economía:** Debido a la inercia del inventario `A3110 + A3210 + B3520` ($\tau \approx 2{,}2\text{ h}$), realimentar `AIC 3255` con el **muestreo de cromatografía iónica de laboratorio cada 2 horas** mantiene $0{,}0\%$ de violaciones tanto en la campaña de 4 días como en el mes completo (ahorrando **$1\,958\text{ m}^3$ de agua en 4 días**, $-95{,}1\%$, y evitando la inversión de un cromatógrafo industrial en línea). El escenario de $\text{Mg}^{2+} = 43\text{–}46\text{ g/L}$ establece el límite operativo superior de la válvula de purga ($61{,}3\text{ m}^3/\text{h}$).

### 3.3 Validación Automatizada del Repositorio
* **Suite `pytest`:** `93 passed, 2 skipped` (`tests/test_core.py`, `tests/test_datos_y_bases.py`, `tests/test_control_integrado.py`, `tests/test_revision_vacio_trim_balance.py` y `tests/test_sostenible_y_avanzado.py`).
* **Figura de síntesis:** `docs/fig/10_replay_campana_real_SQM.png`.