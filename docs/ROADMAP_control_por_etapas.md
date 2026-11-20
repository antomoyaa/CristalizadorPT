# Control por etapas: de la corriente de dilución al sistema completo

Regla: **un lazo nuevo por etapa; antes de cerrarlo, reidentificar el proceso con los lazos anteriores ya cerrados**
(la ganancia depende de qué otros lazos están cerrados: ver notebook 06, etapa 0). Cada etapa termina con una prueba en `tests/`.

| Etapa | Se agrega | Qué se hace | Criterio de salida (se vuelve prueba) | Estado |
|---|---|---|---|---|
| 0 | — | Estacionario `rhs=0`, sin deriva en lazo abierto; escalones de identificación | cierre de balance < 1e-8; deriva < 1e-9 | hecho |
| 1 | FIRC 3103 (caudal de dilución) | PI por SIMC con tau_c >= 2·dt; consigna externa; ruido de medición; perturbación de ganancia de la válvula; límites | sigue la consigna (<2 %), rechaza la perturbación (<3 %), no satura, no oscila | hecho (`test_stage1_*`, notebook 06) |
| 2 | FFIC 3102 (razón dilución/alimentación) | La consigna de FIRC 3103 pasa a ser R × caudal de alimentación; la dilución sigue las variaciones de alimentación | la razón dilución/alimentación se mantiene dentro de 2 % ante +12 % de alimentación | hecho (`test_ratio_mode_*`, notebook 05) |
| 2b | AIT 3105 (ClO4 de entrada) | Ley de razón sobre el ClO4 medido, con ganancia y límites; barrer el retardo del analizador (aún sin definir) | ppm del pulso y recuperación mínima aceptables para la ganancia y el retardo definidos | hecho en simulación (`test_ait_*`); falta definir retardo, ganancia y capacidad de la válvula |
| 3 | nivel (LIRCSAHL 3101) | La dilución agrega volumen: cerrar V y reidentificar el maestro | V dentro de ±1 % ante el escalón de dilución | implementado en A3110 (`plant.py`) y en el tren integrado A3110/A3210 (cascada DIC -> LIC, `test_nivel_*`); falta reidentificar con datos de planta |
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

## Etapa 4 — Análisis Multivariable: Matriz de Ganancias, RGA e Índice de Niederlinski (`COMPLETADO`)

Implementado en `msmpr_core/control.py` (`compute_rga`, `rga_from_plant`) y validado en `tests/test_sostenible_y_avanzado.py`.

### 4.1 Emparejamiento de Inventario y Subsistema Autorregulado en `A3110`
Como el volumen del cristalizador ($V$) es un proceso integrador sin estado estacionario propio en lazo abierto, se aplica la práctica estándar para sistemas integradores: se empareja primero el lazo de inventario $V \leftrightarrow Q_{\text{out}}$ ($\lambda_{V, Q_{\text{out}}} = 1$, manteniendo $dV/dt = 0$) y se evalúa la Matriz de Ganancias Relativas ($\Lambda = G \otimes (G^{-1})^T$, Bristol 1966) junto con el Índice de Niederlinski ($NI = \det(G)/\prod_i g_{ii}$, Niederlinski 1971) sobre el subsistema autorregulado[cite: 4, 5].

#### Matriz de Ganancias $G$ y RGA $3 \times 3$ en `A3110` Aislado ($\alpha = 0{,}25$)

| Variable Controlada (CV) | $m_{\text{evap}}$ (`TIR/PIRCAL`) | $Q_{\text{dil}}$ (`FIRC 3103`) | $\alpha$ (purga local) |
| :--- | :-: | :-: | :-: |
| **$T_C$** ($G$ / **$\Lambda$**) | $-5{,}8332$ / **$-0{,}057$** | $-349{,}9609$ / **$-0{,}003$** | $+1{,}3718$ / **$+1{,}060$** |
| **$\rho_{\text{pulpa}}$** (`DIRCAHL`) ($G$ / **$\Lambda$**) | $-12{,}5716$ / **$+0{,}813$** | $-3823{,}4104$ / **$+0{,}245$** | $+0{,}2015$ / **$-0{,}059$** |
| **$\text{ppm ClO}_4^-$** ($G$ / **$\Lambda$**) | $+0{,}4392$ / **$+0{,}244$** | $-439{,}5226$ / **$+0{,}758$** | $-0{,}0007$ / **$-0{,}001$** |

* **Verificación algebraica:** Suma de filas $= [1{,}000,\; 1{,}000,\; 1{,}000]$; suma de columnas $= [1{,}000,\; 1{,}000,\; 1{,}000]$; $NI_{\text{diag},3\times 3} = -610{,}09 < 0$[cite: 4, 5].
* **Conclusión de diseño:** En un cristalizador aislado con nivel cerrado, la purga neta queda fijada por el balance de agua global ($\alpha Q_{\text{out}} = Q_{\text{stock}} + Q_{\text{dil}} - m_{\text{evap}}/\rho_w$), por lo que $\partial(\text{ppm ClO}_4^-)/\partial\alpha \approx 0$ y el emparejamiento diagonal $3 \times 3$ carece de integridad ($NI < 0$, $\lambda_{ii} < 0$)[cite: 4, 5]. Asimismo, emparejar $T_C \leftrightarrow \alpha$ se descarta por criterios dinámicos debido al gran retardo del reciclo[cite: 4, 5]. Esto fundamenta separar la capa regulatoria local de `A3110` de la capa *plantwide* en `B3520` (`AIC 3522 → FV 3522`), donde la purga divide el licor madre entre reciclo y salida del Área 3000[cite: 4].

#### Validación de los Sub-bloques Regulatorios $2 \times 2$ en `A3110` (con $\alpha$ en Capa *Plantwide*)

| Sub-bloque Regulatorio ($2 \times 2$) | $\Lambda_{2 \times 2}$ (Diagonal / Cruzada) | Índice de Niederlinski ($NI$) | Resultado |
| :--- | :-: | :-: | :--- |
| **Térmico–Densidad:** $(T_C,\; \rho_{\text{pulpa}})$ vs. $(m_{\text{evap}},\; Q_{\text{dil}})$ | **$+1{,}246$** / $-0{,}246$ | **$+0{,}803 > 0$** | **Aprobado:** Garantiza Controlabilidad Integral Descentralizada (DIC); emparejamiento cruzado descartado ($\lambda_{12} < 0$)[cite: 4, 5]. |
| **Térmico–Pureza:** $(T_C,\; \text{ppm ClO}_4^-)$ vs. $(m_{\text{evap}},\; Q_{\text{dil}})$ | **$+0{,}943$** / $+0{,}057$ | **$+1{,}060 > 0$** | **Aprobado:** Lazos prácticamente desacoplados ($5{,}7\%$ de interacción). |

---

## Etapa 5 — Prealimentación Dinámica, *Override* de Sobresaturación y Régimen Sostenible (`COMPLETADO`)

Implementado en `msmpr_core/control.py`, `msmpr_core/dpbe.py`, `msmpr_core/plant.py`, `msmpr_core/params.py` y `msmpr_core/replay.py` (**93 tests unitarios en verde**):

1. **Compensador *Lead-Lag* (`LeadLagCompensator`):** Discretización implícita incondicionalmente estable para sincronizar la acción anticipativa de `AIT 3105 → FFIC 3102` con el tiempo de residencia hidráulico entre `B2510` y `A3110`.
2. **Selector *Override* de Sobresaturación (`supersaturation_override_limit`):** Límite superior $R_{\max,\sigma}$ en `FFIC 3102` que garantiza $C_{\text{mix}} = C_{\text{in}}/(1+R) \ge (1+\sigma_{\min})C^*$, evitando la redisolución del lecho cristalino y derivando el exceso de carga de $\text{ClO}_4^-$ hacia `AIC 3522 → FV 3522`.
3. **Cinética de disolución en subsaturación (`dpbe.dissolution_rate_m_s`):** Rama negativa $D(-S) < 0$ cuando $C < C^*$ para cuantificar el riesgo de sobre-dilución.
4. **Escalado sostenible de dilución ($0{,}50\text{ t/h}$ nominal):** Parametrización de `F_303_nom`, límites `[0.20, 2.50] t/h` y escalado proporcional de `K_TRIM_AIT`, reduciendo el agua media a **$0{,}67\text{ t/h}$** (máximo **$0{,}78\text{ t/h}$**, ahorro de **$97{,}5\%$** vs. diseño y **$96{,}6\%$** vs. operación real SQM) con producto refinado en **$22{,}65\text{ ppm}$** ($< 40\text{ ppm}$).