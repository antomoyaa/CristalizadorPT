# Parámetros para DPBE-MSMPR evaporativo multicomponente de KNO₃

Fecha de corte: **15 de octubre de 2026**, America/Santiago. Alcance: documentación disponible del proyecto, carpeta de Drive indicada y búsquedas externas de originales. **La base resultante permite implementar cierres verificables y diseñar la calibración; no constituye una calibración ni una validación experimental de SQM.**

Se entregan [120 registros](parametros_DPBE.csv), incluidos **45 vacíos explícitos** y 6 casos no aplicables. Las demás filas comprenden mediciones publicadas o internas, ajustes, conversiones, definiciones geométricas y parámetros de simulación identificados como tales. Contar filas no equivale a contar parámetros independientes validados. La [bibliografía](bibliografia.bib) contiene 22 referencias con DOI verificado y 11 documentos o recursos sin DOI inventado.

## 1. Búsqueda, selección y límites de cobertura

Se inspeccionaron README, bases físicas y auditorías locales, código de construcción del notebook multicomponente, tabla previa de parámetros, corrientes industriales, protocolos y suplementos. El modelo original contiene correlaciones de **L-alanina** y valores ilustrativos; no se reutilizaron como datos de KNO₃. Los archivos preexistentes y Drive se usaron en modo de lectura.

La referencia de Drive contiene cuatro subcarpetas. La carpeta bibliográfica se listó completa: **136 archivos**; la carpeta adicional contiene **22 CSV**. Estos aportan **309 filas**, **239 candidatos únicos**, **224 DOI únicos** y 15 candidatos sin DOI. Se eliminaron duplicados por DOI normalizado y, en su ausencia, EID/título. Se corrigió la codificación `%23` de un DOI antiguo, sin alterar el CSV de referencia. Los **224 DOI fueron resueltos en Crossref**; también se verificaron referencias complementarias. Esto verifica identidad y metadatos, **no el contenido de las afirmaciones**. Los resúmenes y columnas Methods/Results de los CSV son pistas, no evidencia primaria.

La selección prioriza compatibilidad de compuesto, matriz, mecanismo, sobresaturación, temperatura, escala y calidad de medición. Se examinó evidencia reciente y se conservan trabajos antiguos cuando aportan mediciones directas difíciles de sustituir. Las citas de Crossref se archivan como información de priorización, sin equipararlas a calidad. No se adoptó ninguna fuente MDPI: había alternativas para los cierres utilizados.

Se consultaron originales completos de Rolfs, Kipp, Selzer —incluido su **suplemento oficial de ACS/Figshare**—, Bozorth, Hnedkovsky, Acevedo y Xu, además de documentos termodinámicos, referencias internas y la revisión de impurezas. Las tablas cuantitativas decisivas y las ecuaciones ambiguas se comprobaron visualmente. Los registros distinguen original completo, consulta dirigida, resumen editorial y acceso fallido. **Helt 1977, Graber 1999 y Kim 2023 siguen siendo pistas prioritarias con tablas originales no recuperadas**; no se inventaron sus constantes a partir de resúmenes. Tampoco se recuperó el original completo de Nielsen 1961 ni el apéndice A suplementario de Vetter.

La búsqueda tiene cobertura reproducible del corpus suministrado y expansión externa dirigida; las barreras de acceso impiden afirmar exhaustividad absoluta de toda la literatura mundial. El [registro de cribado](reproducibilidad/registro_cribado.csv) conserva una decisión para cada candidato; el [registro web](reproducibilidad/registro_busqueda_web.json) conserva consultas y respuestas, incluidos fallos. El [registro de fuentes](reproducibilidad/registro_fuentes.json) documenta el grado de acceso. No se evaluaron todos los 136 archivos de Drive con la misma profundidad.

**Categorías.** 1: sistema industrial multicomponente evaporativo, incluidos sus datos internos, aunque no estén validados. 2: KNO₃/agua, enfriamiento, caras o gotas. 3: sistemas parciales sin operación evaporativa equivalente, otros compuestos y propiedades auxiliares transversales. Así, un dato de KClO₄/agua es directo para esa sal, aunque sea categoría 3 respecto del objetivo principal. La categoría mide aplicabilidad, no calidad. El estado «calculado» también identifica transcripciones de escenarios numéricos y conversiones, cuya procedencia se explica; nunca implica medición nueva.

## 2. Referencia industrial y reconciliación

La corriente **334 SSK-d**, en la [transcripción local preservada](reproducibilidad/entradas/corrientes_usuario.csv), reporta 405.32 t/h de pulpa, 21.8376 t/h de sólidos, 383.4824 t/h de líquido, 70 °C, densidad líquida 1518 kg/m³ y Cp líquido 2248 J/(kg K). Las fracciones de sales del líquido incluyen agua 0.3994, KNO₃ equivalente 0.3460, NaNO₃ 0.0852, NaCl 0.0898 y KClO₄ equivalente 0.0051, junto con Mg, Ca, Li, sulfato y ácido bórico. Son una representación de composición; **no demuestran fases sólidas presentes**.

Con las fracciones iónicas del líquido, el diagnóstico estequiométrico calcula **I = 18.6654 mol/kg de agua** y un desequilibrio de carga de **0.000894%** mediante 100(Σz⁺m−Σ|z⁻|m)/(Σz⁺m+Σ|z⁻|m). Es un cálculo de composición analítica suponiendo disociación y H₃BO₃ neutro, no especiación validada. Esta salinidad descarta tratar el licor como solución ideal diluida. El cálculo también da Q líquido = 252.623 m³/h y volumen aditivo de pulpa = 262.550 m³/h, **0.401%** sobre los 261.5 m³/h tabulados; no se forzó una reconciliación alterando entradas.

Los documentos internos describen tanto 70→55→40 °C como 78→50→37 °C. Son escenarios distintos hasta confirmar campaña y etapas. También aparecen Cp de plantilla = 3500 J/(kg K), densidad de otras muestras ≈1460 kg/m³ y densidad de sólidos mezclados = 2200 kg/m³. No se promedian con las propiedades de la corriente ni se identifican con el cristal puro.

El contraste con el **XLSX original v6** identifica otra condición: `BM Detallado!X36` y `Desglose por corriente!M132` reportan **63.87574 °C**, y `BM Detallado!X31` reporta **260.66732 m³/h**. La densidad líquida original es1518.19038 kg/m³ (`X33`), coherente por redondeo con1518. La temperatura difiere6.12426 K de la transcripción, y esa diferencia no es redondeo. Se conservan ambas referencias; los cálculos de composición anteriores usan explícitamente la base local a70 °C. No se encontró en ese contraste una confirmación experimental del Cp2248. **Antes de calibrar se debe resolver cuál condición corresponde a cada campaña.** Las celdas leídas se archivan en `contraste_v6_columna_X.json` y `contraste_v6_filas.json`.

La base iónica de SSK-d implica 5.56379 g/L de **ClO₄⁻**, equivalentes a 7.75121 g/L de **KClO₄**, usando MClO₄/MKClO₄ = 0.7177966. Esta conversión no identifica precipitado. Los umbrales internos de calidad o concentración deben conservar campaña, base seca/húmeda y expresión analítica. El objetivo local de aproximadamente 40 ppm se interpreta provisionalmente como kg ClO₄ por kg de producto final seco.

## 3. Convenciones y ecuaciones implementables

Se recomienda SI: tiempo s, longitud m, energía J, masa kg y actividades adimensionales. Definir n(L,t) por volumen de **pulpa**, en m⁻⁴, y μk=∫Lᵏn dL. Para V constante y salida representativa Qp:

\[
\partial_t n+\partial_L(Gn)=-Q_p n/V+Q_f n_f/V+\mathcal A[n]+\mathcal F[n],
\qquad (Gn)_{L=0}=B.
\]

Con disolución, G<0 y la frontera L=0 es de salida: desaparecen partículas y su masa regresa al licor. No imponer nacimientos bajo saturación. n(0)=B/G solo sirve como densidad de frontera cuando G>0; el **flujo** es B. La entrada industrial contiene sólidos: n_f=0 requiere justificación experimental.

\[
M_T=\rho_s k_v\mu_3,\quad
M_s=VM_T,\quad
R_{cr}=3\rho_s k_vV\int L^2G(L)n(L)\,dL.
\]

Rcr es kg/s, con signo negativo para disolución. Incluir masa de núcleos si su tamaño de nacimiento es finito. En DPBE usar transferencia líquido-sólido consistente con el cambio discreto de μ₃ y registrar el flujo en Lmax. No calibrar parámetros físicos para compensar errores de frontera, tiempo o truncamiento.

Para KNO₃(s)⇌K⁺+NO₃⁻ y KClO₄(s)⇌K⁺+ClO₄⁻:

\[
a_i=\gamma_i m_i/m^\circ,\quad
S_a=(a_Ka_X/K_{sp})^{1/2},\quad \sigma_a=S_a-1.
\]

Ksp, actividades y estándar deben ser coherentes. ClO₄⁻ acuoso no tiene densidad cristalina, nucleación, factor de forma ni calor de cristalización propios. Si precipita KClO₄, se necesita una segunda población y consumo estequiométrico de K⁺ y ClO₄⁻. Examinar también NaCl y sales/hidratos de sulfato pertinentes; una sola población no describe cualquier conjunto de fases.

Un cierre empírico de referencia es:

\[
G=k_g(T,z)\max(\sigma,0)^g f_L,
\quad G_d=-k_d(T,z)\max(-\sigma,0)^d f_d(L),
\]
\[
B_p=k_{bp}\max(\sigma,0)^{b_p},\quad
B_s=k_{bs}\max(\sigma,0)^{b_s}(M_T/M_{T,ref})^j f_{ag}.
\]

Con σ y MT/MTref adimensionales, kg y kd tienen m/s y kbp/kbs m⁻³s⁻¹. Si se usa MT dimensional, kbs cambia a m⁻³s⁻¹(kg/m³)⁻ʲ. Una constante basada en μ₃ no equivale automáticamente a una basada en MT. Las tasas por volumen líquido deben multiplicarse por VL/V para usarlas por volumen de pulpa.

**No intercambiar sobresaturaciones.** Para q=kg de sal/kg de agua, Δq=q*σq, de modo que kσ=kΔ(q*)ᵍ. Esa conversión depende de T y composición. En un binario, w=q/(1+q), C=ρLw y m=q/M; sus cocientes de sobresaturación difieren cuando ρ y las bases cambian. Entre actividad y concentración se requiere la relación γ(T,z), no un factor universal. Tampoco es equivalente introducir 0.8 en una ley ajustada con σ=0.008.

Los balances de iones/inventarios incluyen alimentación, licor exterior de salida, evaporación de agua y transferencia a cada fase o piscina retenida con términos iguales y opuestos. Con nivel constante se necesita un cierre de volumen/densidad; Qp=Qf−E/ρagua es una simplificación, no una identidad para salmuera que cristaliza y cambia de densidad.

Para energía usar dH/dt=Ḣf−Ḣpulpa−Ḣvapor+UA(Tj−T)+Pagitación. Si H incluye entalpías de fase, no agregar otro calor de cristalización. En un balance sensible reducido, añadir −ṅcrΔHcr con ΔHcr<0 para cristalización exotérmica, y restar Eλ. Evitar mezclar Cp en kJ/(kg K) con U en W/(m² K).

## 4. Parámetros sustentados y alcance

| Bloque | Evidencia y selección | Uso permitido |
|---|---|---|
| Crecimiento KNO₃ | [Rolfs et al. 1997](https://doi.org/10.1016/S0022-0248(96)00430-7), ec.(2), tabla2 p.179: siete pares de k/g por cara; k=1.35×10⁻⁵ a 1.952×10⁻⁴ m/s, g=1.15–1.60. | Benchmark a saturación de 15 °C, σc=0.003–0.018 y licor a 0.16 m/s. No son siete estimaciones intercambiables de dL/dt. |
| Nucleación primaria | [Selzer et al. 2018](https://doi.org/10.1021/acs.cgd.7b01778), ecs.(3),(7),(8), tablas2–3: J₂=A₂Sa exp[−B₂/(ln Sa)²], A₂=(6.5±0.58)×10¹⁵ m⁻³s⁻¹, B₂=19.3±0.94. | Solo mecanismo2, gotas a 17.5 °C, Sa=3.4–3.8. No identifica kbp/bp industrial de ley potencial. |
| Disolución/secundaria | [Xu et al. 2025](https://doi.org/10.1016/j.partic.2024.10.018), tabla1 p.89 y ecs.(3)–(7). | Se conservan como parámetros prescritos de simulación con ambigüedades; no constituyen calibración experimental. |
| Equilibrio KClO₄ | [Bozorth 1923](https://doi.org/10.1021/ja01664a023), tablaII p.2655: solubilidad con KNO₃, NaNO₃ y NaCl añadidos a 25 °C. | Puntos primarios de validación de actividad y efecto de ion común; matriz parcial, sin evaporador. |
| Propiedades del sólido | Perry, tablas2-1 p.2-22 y2-151 p.2-161: densidades nominales 2110/2524 kg/m³ y CpKNO₃(T). [NIST WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7778747&Units=SI&Mask=2): CpKClO₄(T). | Propiedades de los sólidos puros; incertidumbre y expansión térmica no eliminadas. |
| Agua | [IAPWS SR6](https://iapws.org/technical-guidance/release/LiquidWater), tablas1–3 y8; [SR1](https://iapws.org/technical-guidance/release/Supp-sat), ecs.(1)–(3),(6)–(7). | Cp/rho de refrigeración y referencia de saturación de agua pura. |

**Crecimiento y forma.** Rolfs reporta errores de k y g y dispersión de velocidades ≈30%. [Kipp et al.](https://doi.org/10.1016/S0022-0248(96)00431-9) muestran efectos de aditivos e hidrodinámica; no proporcionan un Kd de ClO₄/KNO₃ utilizable aquí. Si L es distancia entre caras, podría relacionarse con velocidades de caras opuestas, pero un diámetro equivalente requiere volumen y evolución de hábito. kv=π/6 es exacto por definición para diámetro de esfera de igual volumen; kv=1 en un modelo de simulación no valida un hábito industrial.

**Nucleación.** Selzer distingue una mezcla de poblaciones con supervivencia P₀=αexp(−k₁t)+(1−α)exp(−k₂t); no corresponde sumar J₁+J₂ como dos mecanismos presentes simultáneamente en todas las gotas. Su suplemento confirma que una exponencial única describe mal los datos y ofrece Weibull como alternativa. El método de zona metaestable de [Acevedo y Nagy](https://doi.org/10.1016/j.jcrysgro.2014.02.024), §6.1, reconoce elevada variabilidad y finalidad de clasificación; no se tomó su constante como tasa numérica industrial. Las energías aparentes de DSC de [Stanford, McCulley y Vyazovkin](https://doi.org/10.1021/acs.jpcb.6b03860) tampoco equivalen a Ea de crecimiento DPBE.

**Equilibrio multicomponente.** Las curvas binarias de Perry son referencias, no cierres SQM. El volumen [IUPAC-NIST 89, parte2](https://doi.org/10.1063/1.4972807) trata nitrato de sodio e incluye evaluación del ternario NaNO₃–KNO₃–H₂O; no es una publicación dedicada al KNO₃ ni contiene toda la matriz industrial. [Marion](https://doi.org/10.1016/S0016-7037(02)00857-8) cubre bajas temperaturas hasta 298 K. [Carroll et al.](https://doi.org/10.1186/1467-4866-6-19) muestran deficiencias de parámetros Pitzer para mezclas con nitrato a 90–120 °C. Ninguno valida por sí solo el intervalo SQM.

La instantánea oficial **USGS pitzer.dat**, SHA256 `3640e62aee63a118f800b115b46a2760576e63e05e1792022315a28f75dbe9bb`, no contiene entradas NO₃ ni ClO₄. La disponibilidad del motor [PHREEQC/Pitzer](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-37.htm) no garantiza cobertura química. Se necesita un conjunto coherente de β⁽⁰⁾, β⁽¹⁾, Cφ, θ, ψ y dependencias térmicas, incluidos boratos/especies neutras cuando corresponda. No completar interacciones ausentes con cero sin evaluar su efecto.

**Información interna de perclorato.** Informaciones Varias §1.3 reporta una trayectoria de 74 a25 °C con inflexión cerca de40 °C. La concentración aumenta inicialmente mientras cambia el licor y después disminuye al precipitar. No se ajustó esa trayectoria como una solubilidad a composición constante. La hoja Curva KClO₄ se transcribe en CSV, pero carece de dominio validado de composición y presenta una etiqueta NO₃/NaNO₃ ambigua. Su ajuste interno no es validación independiente; límites operacionales de calidad tampoco son Ksp.

**Ostwald y Gibbs–Thomson.** Se utiliza como hipótesis γs-l,aparente=0.0227 J/m², inferida por Selzer mediante CNT, y Vm=M/ρ=4.79162×10⁻⁵ m³/mol. No es una medición interfacial independiente. Para esfera de diámetro L y sal1:1:

\[
\ln(IAP/K_{sp})=4\gamma V_m/(RTL),\quad
\ln S_{a,eq}=\ell_a/L,\quad \ell_a=2\gamma V_m/(RT).
\]

A 290.65 K resulta **ℓa=9.00189×10⁻¹⁰ m**; para el cociente IAP, la longitud es el doble. La incertidumbre de γ aporta 1.76%, sin incluir densidad, convención CNT ni transferencia a salmuera. A L=100 µm, Sa,eq−1≈9×10⁻⁶; a1 µm ≈9×10⁻⁴. Esto permite evaluar cuándo importa la curvatura frente a la sobresaturación medida, sin declarar insignificante todo ripening. Las energías cristal-vapor de Rolfs no sustituyen γ sólido-licor. El α adimensional de [Vetter et al.](https://doi.org/10.1021/cg4010714), ec.(3), usa número de unidades como coordenada: no equivale directamente a ℓ en metros ni justifica α=1.

## 5. Cálculos reproducibles y discrepancias

Los cálculos completos, bases y precisión numérica están en [reproducir.py](reproducibilidad/reproducir.py) y [calculos.json](reproducibilidad/resultados/calculos.json). La precisión del archivo no aumenta la precisión experimental.

- **Concentraciones:** q40,KNO₃=0.639 kg/kg agua; w40=0.639/1.639=0.38987. La pendiente local ln(0.855/0.458)/20=0.0312116 K⁻¹ corresponde a q, no a C en kg/m³. Rolfs da q40=0.64599, 1.09% mayor; son correlaciones distintas sin incertidumbre suficiente para resolver una discrepancia estadística. Para la matriz real se selecciona ninguna de las dos como equilibrio validado.
- **KClO₄ a25 °C:** Bozorth reporta 0.1491 mol/kg agua, equivalente a 0.0206577 kg/kg agua. Se prefiere ese punto directo a interpolar entre20 y30 °C para ese benchmark; no se reajusta toda la curva con un punto. Con KNO₃ añadido, el primer y segundo producto de actividad relativo se reproducen por s₀²/[s(s+mKNO₃)]. Para mKNO₃=0.6686, el cálculo da **0.472536**, frente a **0.490** tabulado. Se conservan ambos y se marca la discrepancia sin corregir arbitrariamente el original.
- **Calor cristalino:** CpKNO₃=(6.42+0.0530T_K)×4.184/0.1011032; a40 °C da **952.521 J/(kg K)**, incertidumbre publicada10%. Para KClO₄ se prefiere la dependencia Shomate NIST a la aproximación constante antigua de26.3 cal/(mol K), restringida a287–318 K; a40 °C Shomate da **851.492 J/(kg K)**. Ambas son compilaciones, no mediciones recientes del producto SQM.
- **Calor de cristalización:** invertir el calor de disolución a dilución infinita de Perry produce −36.1205 kJ/mol para KNO₃ y −54.1410 kJ/mol para KClO₄ a18 °C. El dato interno −7.73 kcal/mol se convierte a **−32.3423 kJ/mol**, redondeado en la hoja como −32.35. Los estados inicial/final difieren: no son estimaciones equivalentes del calor diferencial industrial. Faltan entalpías parciales y mezcla/dilución a composición real.
- **Agua:** SR6 a20 °C y0.1 MPa da **Cp=4184.189 J/(kg K)** y **ρ=998.20655 kg/m³**. SR1 más Clapeyron, λ=T(vv−vl)dp_sat/dT, da a40 °C **λ=2.406081 MJ/kg** y p_sat=7.38511 kPa. Condensación invierte el signo. El valor interno10.34 kcal/mol equivale a43.26256 kJ/mol, pero requiere T,p reales. La evaporación desde salmuera usa h_v−h̄agua,licor, que no se identifica únicamente con λ de agua pura.

Se comprobaron visualmente dos problemas adicionales de publicación. En **Xu**, ec.(4) omite L en D=kdS(1+α)^β, y la tabla1 rotula γ/α con unidades incompatibles con las sumas adimensionales; además, n y kb carecen de una normalización de volumen compatible directamente con nuestra DPBE. Se archivan valores sin reparar silenciosamente la ecuación. En el **SI de Selzer**, ec.(2) imprime τ=tmed(ln2)¹ᐟᵝ; derivando su ec.(1) resulta **tmed=τ(ln2)¹ᐟᵝ**, relación que sí reproduce la tablaS1. Se usa la derivación y se conserva la evidencia visual del error.

## 6. Incorporación de impurezas y transporte

La revisión de [Nyande et al. 2025](https://doi.org/10.1016/j.cherd.2025.08.002), §2, distingue inclusión en red, adsorción, oclusión de licor y retención exterior. Para cada piscina p usar dIp/dt=Jp,entrada−Jp,liberación−(Qp/V)Ip, retirando o restituyendo exactamente la misma cantidad en el licor. La oclusión incorpora también agua y sales acompañantes; el precipitado KClO₄ tiene su propio balance. La liberación durante disolución debe conservar impureza y masa.

Un Kd definido como wClO₄,cristal seco/wClO₄,licor es adimensional solo con esas bases. No equivale a una eficiencia de lavado ni a concentración de impureza en toda la pulpa. Un Kd aparente puede mezclar varios mecanismos y variar con G, tamaño, agitación y composición. **No se identificó Kd específico validado para ClO₄/KNO₃ en SQM.**

La pureza final requiere modelo de centrifugación, recuperación η(L), licor exterior retenido y lavado: ppm=10⁶ mClO₄,retenido/mproducto final seco. El denominador incluye residuos no volátiles depositados al secar el licor retenido. No deducir Kd a partir de un KPI ni igualar sal cristalina y residuo seco.

Para un cierre difusivo Jm=km(C−C*), con área kaL² y volumen kvL³:

\[
G=\frac{k_a}{3k_v\rho_s}k_m(C-C^*).
\]

Para diámetro esférico, el factor es2/ρs; para velocidad normal de cara es1/ρs. Por eso **km y kg no son el mismo parámetro**. Identificar viscosidad, difusión, Re/Sc y agitación antes de usar Sh(Re,Sc). Aglomeración/rotura necesitan kernels, eficiencia y distribución de fragmentos que conserven volumen; evitar contar la misma atrición como rotura y nucleación secundaria dos veces. Su omisión inicial es una hipótesis contrastable, no una tasa física igual a cero.

## 7. Calibración y validación necesarias

| Prioridad | Datos/experimento | Parámetros e identificación |
|---|---|---|
| 1 | Equilibrado desde subsaturación y sobresaturación a37,40,50,55,70 °C, con licor representativo; medir iones, agua, densidad, pH y sólidos por PXRD. Variar K,Na,Mg,Cl,NO₃ yClO₄ con electroneutralidad. | C*(T,z), actividades/aw, fases competidoras y umbral KClO₄. La lista de temperaturas es diseño propuesto, no rango ya validado. |
| 2 | Crecimiento y disolución sembrados, isoterma, baja generación de finos; CSD/imagen y masa frente a tiempo; varios tamaños y agitaciones. | kg,g,kd,d, dependencia L/T y transporte. Resolver hábito y base volumétrica antes del ajuste. |
| 3 | Réplicas de inducción sin semillas y ensayos con semillas variando σ,MT y agitación **independientemente**. | Primaria versus secundaria, k,b,j y función de agitación; distinguir finos por atrición. No identificar j si MT permanece constante. |
| 4 | Calorimetría de licor y sólidos, evaporación/condensado y balance térmico con T,p medidas. | Cp(T,z), calor diferencial, E,UA y efecto de actividad de agua. Separar sensible, mezcla y cambios de fase. |
| 5 | Comparar cristal sin lavar, lavado secuencial, licor y residuo seco; microscopía/PXRD y análisis ClO₄ específico; variar G, tamaño y wash ratio. | Kd de red, inclusiones, adsorción, precipitado y licor exterior. Fijar definición analítica de ppm. |

Validar fuera del conjunto de ajuste y por campaña: cierres de K,NO₃,ClO₄,agua/energía, producción, CSD, humedad y pureza. Reportar intervalos y correlaciones de parámetros; perfiles de verosimilitud o bootstrap cuando los datos lo permitan. Antes de ajustar, comprobar conservación DPBE, solución analítica MSMPR en el caso G constante, independencia de malla/paso y efectos de Lmax. Las nueve comprobaciones de este trabajo verifican cálculos documentales; **no validan un simulador ni desempeño industrial**.

## 8. Reproducción y preservación

Ver [instrucciones de reproducción](reproducibilidad/README.md). La corrida entregada ejecuta nueve comprobaciones: tres puntos oficiales de agua líquida, saturación normal, bases de concentración, factor iónico de Gibbs–Thomson, dos controles de Bozorth —incluida detección de discrepancia— y tablaS1 Weibull. Los scripts de cálculo y regeneración usan datos archivados, sin necesidad de acceso a Drive o a internet.

El inventario inicial abarca **8988 archivos preexistentes**, con SHA256, excluyendo `.venv-educational`, `.git` y `__pycache__`. La verificación final se conserva aparte. Todas las escrituras de la investigación se realizan dentro de la carpeta nueva. No se corrigieron notebooks, hojas, datos de referencia ni archivos originales.
