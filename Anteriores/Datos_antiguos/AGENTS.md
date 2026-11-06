# AGENTS.md — contexto permanente (guía de ajuste §8)

## Proyecto
Modelo de un cristalizador al vacío de KNO3 (NPT3, A3110/A3210, GEA Messo DTB) con reciclo de licor madre y
control PID. Código en `msmpr_core/`; notebooks numerados en `notebooks/`; datos en `data/`; parámetros con fuente en `parametros/`.

## Modelo (no cambiar sin preguntar)
- sigma = C/C* - 1;  G = kG exp(-Eg/R(1/T-1/Tref)) sigma^g [m/s];  B = kB sigma^b [#/(m3 s)]  (o forma centrada en sigma_ref).
- Cierre MSMPR: C_in - C = M_T = 6 rho_c kv B G^3 tau^4;  d50 = 3.672 G tau; (d90-d10)/d50 = 1.344.
- Unidades SI; concentraciones en kg KNO3 / m3 de pulpa; T en degC en la interfaz, K dentro.
- Recuperación de DOS formas que deben coincidir: Y_liq = 1 - C/C_in  y  Y_sol = M_T/C_in (`Y_cierre` ~ 1e-12).
- Vector de respuestas Y: d50, span, N_kg, Y, sigma, prod_tph.  Vector de parámetros P: tau_h, f_dil, T (operación) y Pn (cinética).

## Calibración (calibrador de referencia)
- Calibrar con `calibracion_do_crystal.ipynb` (A3110 y A3210 por separado). `msmpr_core/fit.py` es solo un chequeo estacionario rápido.
- Formato de datos canónico: el formato LARGO de ese notebook (hojas `mediciones` y `operacion`). Convertir Excel ancho con
  `python -m msmpr_core.datos` (exige declarar equipo y origen; no inventa operación).
- kG, kB y Kd calibrados solo valen sobre la base con que se calibraron. Usar `make_params("calibracion")` y `apply_calibration`.
  El notebook MSMPR_KNO3_DPBE_sensibilidad_Pareto NO se edita: su hash lo vigila la prueba "Original intacto".
- Nunca llevar a BASE valores calibrados con `origin = synthetic`.

## Reglas
1. No inventar datos ni parámetros. Todo valor tiene etiqueta en `params.PROVENANCE`: respaldado / interno / placeholder.
2. No cambiar unidades ni la base de concentración. No borrar código existente. No silenciar errores con try/except vacíos.
3. Todo dato sintético se rotula como tal (tablas, figuras, texto).
4. Un cambio del modelo = correr `pytest tests/` (balance, estacionario, ajuste de parámetros sintéticos, control).
5. Los tags de instrumentos vienen del P&ID 172164 hoja 10 (A3110). Lo no leído del plano se declara como supuesto.

## Criterios de aceptación mínimos
- Y entre 0.3 y 0.7 en el caso base calibrado (hoy NO se cumple: C*_ref y la alimentación son placeholders).
- Y_liq = Y_sol (|dif| < 1e-6).  El ajuste recupera parámetros sintéticos dentro de IC95.
- En lazo abierto el estacionario no deriva.

## Formato de respuesta de un agente
Primero un plan (5-8 pasos) y los supuestos; después archivos cambiados, salida de pytest y limitaciones.
