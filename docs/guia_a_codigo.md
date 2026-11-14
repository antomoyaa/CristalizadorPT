# De la guía de ajuste (PDF 24-10-2026) al código

El PDF `Guia_ajuste_parametros_SQM.pdf` está en esta carpeta. Cada sección tiene su lugar en el repositorio:

| Guía | Contenido | Dónde vive |
|---|---|---|
| §1–2 | Vectores Y y P; MSMPR ideal (Gamma(4), d50 ∝ tau^((b−g)/(b+3g))) | `msmpr_core/model.py` (`msmpr`), `dpbe.py`; pruebas `test_gamma4_*`, `test_d50_exponent_*` |
| §3.1 | Ajuste de Pn con IC95, correlaciones, reparametrización y σ medido | Calibrador de referencia: `calibracion_do_crystal.ipynb` (equipo por equipo, dinámico, validación temporal). Chequeo estacionario rápido: `msmpr_core/fit.py`; `test_fit_recovers_*` |
| §3.2 | Sin datos: tornado y elasticidades | pendiente (usar los barridos del notebook de sensibilidad/Pareto existente sobre `msmpr`) |
| §4 | Recuperación por líquido y por sólido | `Y_liq`, `Y_sol`, `Y_cierre` en `model.py` y `plant.py`; `test_mass_balance_*` |
| §5 | Especificación wKNO3 mín. y ClO4 máx. | `config/specs.yaml` (PENDIENTE de SQM) |
| §6 | RSM Box-Behnken y modelo híbrido de ClO4 | pendiente (necesita datos; la plantilla de datos está en `data/`) |
| §7 | Frente de Pareto en t/h y ppm, con banda y rango de validez | pendiente de conectar con `msmpr`/RSM; hoy existe el notebook `MSMPR_KNO3_DPBE_sensibilidad_Pareto` |
| §8 | Prompt para Codex | `AGENTS.md` en la raíz |
| §9–10 | Plan de 2 semanas y lista de chequeo | la lista de chequeo se marca contra `pytest` y `params.placeholders()` |

Pendiente de la guía: faltan §3.2, §6 y §7 porque dependen del notebook de ejemplo `Ejemplo_ajuste_RSM_Pareto_SQM.ipynb` y de datos que aún no hay.
Cuando lleguen datos reales: llenar `data/experimentos_plantilla.csv` (con `origen`), correr `ajustar_cineticas` y recalibrar `Cstar_ref_kg_m3` con solubilidad de licor representativo (prioridad 1 del informe de parámetros).
