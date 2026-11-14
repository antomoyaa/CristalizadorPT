# Agua de dilución mínima y temperatura (notebook 06, B6–B8)

## Resultados (simulación con excursión de ClO4- +550 ppm a los 40 min y Mg2+ +35 % a los 140 min)
- **B6.** Con tope de agua de 36 t/h, las estrategias 1-3 (sin reciclo, reciclo fijo, AIT + AIC) superan los 40 ppm; la estrategia 4 (AIT + AIC + **supervisor de T2**) queda bajo la norma. El veredicto se calcula con el máximo simulado.
- **Agua.** Con el supervisor de T2 la calidad apenas depende del agua (pico 31,8 ppm con 36 t/h; 35,8 ppm con 20 t/h). El óptimo calculado con el criterio de `plant.barrido_agua_optimo` (pico ≤ 90 % de la norma y pérdida de KNO3 ≤ 1 t/h) es **20 t/h**, 46 % menos que los 36,8 t/h anteriores y 26 % menos que el nominal (27,14 t/h).
- **Temperatura.** Alimentación +10 °C satura el vacío de A3110 el 85 % del tiempo y pierde KNO3; agua de dilución a 35 °C (≈ 790 kW de enfriamiento) devuelve el margen de vacío; con prealimentación la T1 queda en 55,0 °C. Un cristalizador 2 caliente reduce el riesgo de precipitar KClO4 pero cuesta 1,82 t/h de KNO3 por °C.

## Límites de estos resultados (leer antes de usarlos)
1. Sobre 40 °C **no hay dato de solubilidad de KClO4 en la salmuera de SQM**: se usa la razón del agua pura (`KCLO4_SOL_AGUA_REF`, referencia tipo CRC). La eficacia del supervisor de T2 depende de ese dato.
2. Los modelos de coprecipitación y humedad son correlaciones empíricas (`params.EMPIRICAL_COEFFS`), no están calibrados con planta.
3. La pérdida de KNO3 por T usa el balance de la Entrega 2 (cristalizado por grado en A3210), no una isoterma de la salmuera.
4. El criterio del óptimo (90 % de la norma, 1 t/h) es una decisión de diseño; cámbiala en `barrido_agua_optimo(margen=..., tol_KNO3=...)`.
