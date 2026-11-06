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
