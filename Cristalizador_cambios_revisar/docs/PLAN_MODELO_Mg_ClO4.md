# Plan: balances de Mg2+ y su interacción con ClO4- en el modelo del cristalizador
Ramas: `estimacion-parametros` (ajuste de parámetros) y `modelo-cristalizador2-kno3` (MSMPR evaporativo, notebook `Evaporative_MSMPR_DPBE-2`).

## 1. Estados y balances (ambas ramas)
Líquido con volumen V constante, evaporación `m_ev`, caudales de alimentación y producto Q_f, Q_p:
- Mg2+ (no cristaliza ni entra a la red de KNO3): `V dC_Mg/dt = Q_f C_Mg,f - Q_p C_Mg  (+ efecto de concentración por m_ev/rho)`. Sale solo con el licor adherido al cristal (h ~ 7 % del producto crudo) y con el licor madre.
- ClO4-: igual, más un sumidero por precipitación de KClO4: `- R_KClO4`, con `R_KClO4 = k_p * max(C_ClO4 * C_K / Ksp_eff - 1, 0)` (primer orden en la sobresaturación; k_p a estimar).
- Producto: `ppm ClO4 = h * C_ClO4,adherido + fraccion_coprecipitada` (la parte ocluida/coprecipitada no sale con el lavado).

## 2. Interacción Mg2+ - ClO4-
La saturación de KClO4 depende de T, Mg2+ y NaNO3: `ln(C*_KClO4 / C*_0(T)) = - beta_Mg * C_Mg - beta_Na * C_NaNO3` (tipo Setschenow) o, mejor, con coeficientes de Pitzer. Datos SQM disponibles: con Mg 27,1 g/L, saturación de 3,5 a 6,5 g/L de KClO4 entre 25 y 40 °C; cualitativos: más Mg o menos NaNO3 => menos solubilidad.
- En `estimacion-parametros`: agregar `beta_Mg`, `beta_Na`, `k_p` a `BASIS_KEYS`/`fit.py` y ajustar con `Base de datos.xlsx`.
- En `modelo-cristalizador2-kno3`: definir `steady_crystal_perchlorate_ppm(run)` (hoy no existe y el heatmap queda en NaN) con el balance de la sección 1.

## 2.1 Datos de literatura / SQM que necesito
1. Solubilidad de KClO4 vs T (25-60 °C) en salmuera saturada de KNO3 con al menos 3 niveles de Mg (p. ej. 10, 30, 50 g/L) y 2 de NaNO3. Si no hay ensayos: parámetros de Pitzer K+/ClO4-, Mg2+/ClO4-, Na+/ClO4-, Mg2+/NO3- (Pitzer; Kim & Frederick; Pitzer-Mayorga) y la solubilidad de KClO4 en agua pura 0-80 °C.
2. Solubilidad de KNO3 en la salmuera (T y Mg) para el rendimiento.
3. Distribución de Mg y ClO4- entre cristal y licor adherido (composición del producto crudo y del licor madre por lote) para fijar `h` y la fracción ocluida.
4. Cinética de nucleación de KClO4 (tiempos de inducción, protocolo 12.10.2) para `k_p`.
5. Serie de tiempo de Mg y ClO4- en B3520/SSKd (Base de datos.xlsx, hoja con fechas) para validar.
