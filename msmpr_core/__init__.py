"""msmpr_core: núcleo compartido del modelo MSMPR/DPBE de KNO3 (planta NPT3, A3110/A3210).

Módulos
-------
params   : parámetros base con su procedencia (respaldado / interno / placeholder)
dpbe     : cinética, CSD estacionaria (analítica y en malla) y percentiles
plant    : modelo dinámico del cristalizador al vacío con reciclo y purga (única fuente de ecuaciones)
control  : PID industrial, cascada, sensor muestreado, identificación FOPDT y sintonía
model    : API msmpr(tau_h, f_dil, T, Pn) de la guía de ajuste (Y, d50, N, sigma, producción)
fit      : ajustar_cineticas(df) por mínimos cuadrados ponderados con IC95 y correlaciones
"""
__version__ = "0.1.0"
