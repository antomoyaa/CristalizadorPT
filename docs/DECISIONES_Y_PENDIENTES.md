# Decisiones confirmadas y pendientes del modelo de control (A3110)

Última actualización: 2026-11-01. Fuente de lo "confirmado": indicaciones del equipo durante el desarrollo del modelo. Lo no confirmado queda como supuesto.

## Confirmado por el equipo
- `DIRCAHL 3101` es un sensor de **densidad** (se había leído como `OIRCAHL`).
- El agua de dilución la controla **`FIRC 3103`**; su consigna la entrega **`FFIC 3102`** (razón dilución/alimentación).
- **`AIT 3105`** es un analizador de **perclorato en la corriente de entrada al cristalizador**; si la concentración viene alta, envía la señal a `FFIC 3102`, que la entrega a `FIRC 3103` para abrir o cerrar la válvula.
- El **retardo del analizador aún no está definido**.
- Se modela el P&ID por etapas, partiendo del agua de dilución y sumando los lazos propuestos.

## Supuestos del modelo (por confirmar)
- `AIT 3105` mide la alimentación **antes** de mezclar el agua (acción anticipada). Si mide después de la dilución, es realimentación y la ley cambia.
- La cascada `DIRCAHL 3101 -> FIRCAL 3101` **no** es el agua de dilución; su actuador está por identificar y no se modela.
- Constante de tiempo de la válvula de dilución: 8 s (placeholder).
- Resolución de `DIRCAHL`: 0.5 kg/m3 (supuesto, solo en el notebook 06).
- Todos los parámetros de planta (solubilidad, cinética, Kd, alimentación de KNO3) siguen siendo placeholders: ver `msmpr_core.params.placeholders()`.

## Pendiente de definir
| Qué | Quién/cómo | Dónde se usa |
|---|---|---|
| Retardo de `AIT 3105` | equipo de proceso / hoja del analizador | `ait={"delay_s": ...}` en `PlantController`; barrido en el notebook 05 |
| Ganancia y límites de `FFIC 3102` | evaluación económica (ppm fuera de norma frente a producción perdida) | `ait={"gain":..., "R_max":...}` |
| Capacidad real de la válvula de `FIRC 3103` | hoja de datos de la válvula | `mv_limits={"Q_dil": (0, ...)}` |
| Punto de muestreo exacto de `AIT 3105` | P&ID / planta | decide anticipación o realimentación |
| Actuador de `DIRCAHL 3101 -> FIRCAL 3101` | leyenda del plano | hoy no modelado |
| Solubilidad del licor real y calibración de kG, kB, Kd | calibrador del equipo (`calibracion_do_crystal.ipynb`) con datos medidos | `params.apply_calibration` |
| Especificación de producto (wKNO3 mín., ClO4 máx.) | SQM, por escrito | `config/specs.yaml` |

## Resultados de referencia (parámetros placeholder; solo tendencias)
- Compensar +25 % de ClO4 de entrada con agua exige ~6 veces la dilución nominal y lleva la recuperación de KNO3 a ~0 %.
- Con un retardo de `AIT` comparable a la duración del pulso, el agua llega tarde: se gasta igual y el ppm pico no mejora.



**NUEVO**

# Decisiones de Diseño de Control y Pendientes de Planta (SQM)

## 1. Decisiones Técnicas Adoptadas
- **Control de Temperatura vía Vacío (`TIR 3101/3201 -> PIRCAL 3101/3201`):** Se manipula la presión de vacío entre $-0.93\text{ bar}$ y $-0.97\text{ bar man.}$ (o $m_{\text{evap}} \in [0.2\times, 3.0\times]$ nominal) con sintonía amortiguada Tyreus-Luyben e IMC para evitar ciclos límite de Hopf.
- **Control de Densidad en Cascada (`DIRCAHL 3101/3201 -> FIRCAL -> SIC 3101/3201`):** Regula la fracción de sólidos hacia A3210 ($x_{S,1} = 0.13$) y hacia el hidrociclón S3540 ($x_{S,2} = 0.20$) variando el caudal de las bombas P3112 y P3212 ($[250, 520]\text{ t/h}$).
- **Ponderación de Consigna ($\beta = 0.5$) y Anti-Windup:** Se adopta $\beta = 0.5$ con precarga integral $I_0 = -K_p s (\beta y_{sp,0} - y_0)$ para suprimir la saturación del actuador ante escalones de consigna.
- **Mitigación Dual de Perclorato (`AIT 3105` + `AIC 3255/3522`):** La razón variable en `FFIC 3102` ($R_0 = 0.07176$, $K_{\text{trim}} = 2.914\times 10^{-4}\text{–}5.2\times 10^{-4}\text{ ppm}^{-1}$) combinada con la purga activa post-B3520 mantiene el producto refinado del Área 6000 en $23.0\text{–}24.5\text{ ppm}$ ($< 40\text{ ppm}$).

## 2. Parámetros Pendientes de Confirmación en Planta
- **Tiempo muerto de analizadores en línea (`AIT 3105` y `AIC 3255/3522`):** Modelado como supuesto de trabajo en $10\text{–}15\text{ min}$ (con barrido de sensibilidad de $0$ a $120\text{ min}$); requiere validación frente a la duración real de las excursiones de $\text{ClO}_4^-$.
- **Capacidad máxima y constante de tiempo de la válvula `FIRC 3103`:** Se asume $\tau_v = 8\text{ s}$ ($0.133\text{ min}$) y rango operativo de $[12, 42]\text{ t/h}$ (comparado frente a topes de $3\times$ y $8\times$ el caudal nominal).
- **Resolución de densímetros `DIRCAHL 3101/3201`:** Evaluada bajo supuestos de cuantización de $0.0002$ a $0.0015$ en fracción másica de sólidos.
- **Caracterización del ruido industrial (Historiador PI):** Las desviaciones estándar del proceso de Ornstein-Uhlenbeck ($\tau_c = 10\text{ min}$) deben calibrarse con series temporales de planta.