# Bases de parámetros del repositorio (qué es cada una y qué NO se puede mezclar)

| Base | Dónde | KNO3 alimentación [kg/m3] | C*_ref (40 C) [kg/m3] | Pendiente solubilidad [1/K] | rho_liq [kg/m3] | Origen |
|---|---|---|---|---|---|---|
| `calibracion` (= `BASE` del notebook de sensibilidad y de `Estimacion_parametros_MSMPR_KNO3.ipynb`) | `params.BASE_CALIBRACION` | 220 | 180 | 0,010 | 1040 | Base **de trabajo** del notebook de Pareto; kG, kB, Kd calibrados solo valen sobre ella |
| `planta_sqm` | `params.BASE_PLANTA_SQM` | 525,2 (corriente 334 SSK-d) | 177,94 | 0,010 | 1518 | Alimentación: informe de parámetros SQM. **C*_ref y pendiente: SINTÉTICOS (placeholder)** |
| `literatura` | `params.BASE` | placeholder | Perry (curva g/g agua) | respaldada | 1040 | Valores con fuente donde existen |

Reglas:
1. `datos_cristalizador.xlsx` es un conjunto **sintético** (lo genera `Estimacion_parametros_MSMPR_KNO3.ipynb` con la base `calibracion` + ruido). Los parámetros que se "estiman" con él solo recuperan los valores usados para generarlo: no son una calibración de KNO3.
2. Los kG, kB, g, b, Kd de una base no son transferibles a otra (`BASIS_KEYS`); `params` lo verifica al aplicar una calibración.
3. Para dejar de usar placeholders: medir la curva de solubilidad real (protocolo `12.10.4`), reemplazar `Cstar_ref_kg_m3` y `solubility_slope_K_inv` de `planta_sqm`, y recalibrar la cinética con los protocolos de nucleación (`12.10.2`).
