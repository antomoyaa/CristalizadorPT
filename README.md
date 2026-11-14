# Diseño y Simulacion de Lazos de Control PID — Cristalizadores de KNO₃ (Planta NPT-3, SQM)

Este repositorio contiene el modelo fenomenologico (`msmpr_core`) y los cuadernos interactivos de diseño, sintonia y validacion dinamica de los **lazos de control regulatorio y supervisor** para los cristalizadores evaporativos al vacio **A3110** y **A3210** (Area 3000), integrados con el **Area 2000** (alimentacion SSK) y el **Area 6000** (refinacion final, meta de calidad **< 40 ppm de ClO₄⁻**).

---

## 1. Arquitectura de Lazos de Control e Instrumentacion (P&ID GEA/SQM 172164)

| Tag P&ID | Clasificacion | Variable Controlada (CV) | Variable Manipulada (MV) | Funcion Operacional |
|---|---|---|---|---|
| **TIR 3101/3201 -> PIRCAL** | Existente (Base) | Temperatura (55 °C / 40 °C) | Presion de vacio (-0.93 a -0.97 bar) / Evaporacion | Sostiene la sobresaturacion nominal en cada etapa |
| **LIRCSAHL 3101** | Existente (Base) | Nivel / Volumen (V = 660 m³) | Caudal de descarga (Q_out) | Estabiliza el inventario liquido del cristalizador |
| **DIRCAHL 3101/3201 -> FIRCAL -> SIC** | Existente (Base) | Densidad de pulpa / Solidos (13% y 20%) | Velocidad bombas P3112 y P3212 | Regula la consistencia de pulpa enviada a separacion |
| **FIRC 3103** | Existente (Base) | Caudal de agua de dilucion (corriente 303) | Apertura de valvula de agua (tau_v = 8 s) | Lazo esclavo rapido que sigue la consigna de FFIC 3102 |
| **FFIC 3102** | Existente / Mejorado | Razon agua / alimentacion (R = F_303 / F_alim) | Consigna de caudal para FIRC 3103 | Mantiene la proporcion de dilucion ante cambios de carga (FIRCAHL 3101) |
| **AIT 3105** | **Propuesto (Nuevo)** | ClO₄⁻ en alimentacion pre-A3110 | Ajuste (trim) de la razon R(t) en FFIC 3102 | Anticipa alzas de perclorato en SSK antes de entrar a A3110 |
| **AIC 3255 -> FV 3522** | **Propuesto (Nuevo)** | ClO₄⁻ en licor madre (B3520) | Fraccion de purga alpha(t) post-B3520 | Evita co-precipitacion de KClO₄ y asegura < 40 ppm en Area 6000 |

---

## 2. Guia de Lectura de los Notebooks (`notebooks/`)

1. **`notebooks/04_Recirculacion_purga.ipynb` — Analisis Estacionario y Frente de Compromiso:**
   * Verifica el cierre exacto de balances (Y_liq = Y_sol) y cuantifica el conflicto entre recuperacion de KNO₃ y pureza (ppm de ClO₄⁻) frente a dilucion, temperatura y purga alpha.
2. **`notebooks/05_Control_PID_dinamico.ipynb` — Lazo Abierto, Lazo Cerrado, Ruido y Estabilidad (Routh-Hopf):**
   * Valida la identificacion FOPDT (±10%), evalua 12 escenarios de perturbacion (±) con ruido coloreado (Ornstein-Uhlenbeck), comprueba offset nulo a 24 h y determina la ganancia ultima (K_cu) y **Bifurcacion de Hopf** por Routh-Hurwitz (Ziegler-Nichols vs. Tyreus-Luyben).
3. **`notebooks/06_Control_por_etapas_dilucion.ipynb` — Control por Etapas (`AIT 3105 -> FFIC 3102 -> FIRC 3103` + `AIC 3255`):**
   * Despliega en secuencia el cierre de FIRC 3103 (Etapa 1), FFIC 3102 (Etapa 2), barrido de retardo theta_AIT en [0, 2 h] (Etapa 3) y la simulacion acoplada del tren **NPT-3 (Areas 2000-3000-6000)** incluyendo el efecto del magnesio (Mg²⁺) y la comparacion **SIN AIT 3105** vs. **CON AIT 3105 + AIC 3255**.

---

## 3. Ejecucion y Pruebas Automatizadas

```bash
pip install -r requirements.txt
python3 -m pytest tests/ -v
```
