# Dinora — pipeline LC-MS/MS en R (xcms + RforMassSpectrometry)

Estado: **verificado y ejecutado**. El pipeline corre de punta a punta sobre datos
reales de LC-MS (4 archivos de ejemplo, 1 014 features, 46 s en 4 núcleos).

> No está listo para los datos del laboratorio todavía: falta que el doctor
> mande las muestras y confirme el tipo de experimento (ver
> [Qué necesitamos](#qué-necesitamos-del-doctor)).

---

## 1. Qué se pidió y qué se hizo

El doctor quiere **replicar una forma de analizar datos de un experimento LC-MS/MS**,
siguiendo un curso de R. Los enlaces compartidos fueron:

| Enlace | Qué es en realidad |
|---|---|
| `github.com/sneumann/xcms/tree/devel/R` | **Código fuente** del paquete `xcms` (66 archivos `.R`), no scripts de análisis. Los "comandos" que busca están en las **viñetas** y en [Metabonaut](https://rformassspectrometry.github.io/Metabonaut/). |
| PMC12750407 | Louail *et al.*, *xcms in Peak Form*, Anal. Chem. 2025. **Preprocesamiento**: detección de picos, alineamiento, agrupamiento, QC. |
| Metabolites 12(2):173 | Rainer *et al.*, *A Modular and Expandable Ecosystem for Metabolomics Data Annotation in R*, 2022. **Anotación**: `MetaboAnnotation`, `MetaboCoreUtils`, `CompoundDb`. |

Los dos artículos **no compiten**: uno es el preprocesamiento y el otro la
anotación posterior. Juntos forman la cadena completa, y es exactamente lo que
implementa este repositorio.

## 2. ¿Se puede hacer? Sí

Verificado ejecutando el stack completo, no solo leyendo documentación:

| Paso | Función | Resultado en la prueba |
|---|---|---|
| Importar | `readMsExperiment()` | 4 archivos, 5 112 espectros |
| Inspeccionar | `chromatogram()` | BPC + heatmap de similitud entre TIC |
| Detectar picos | `findChromPeaks(CentWaveParam())` | 1 385 picos |
| Refinar | `refineChromPeaks(MergeNeighboringPeaksParam())` | 1 292 picos |
| Alinear | `adjustRtime(ObiwarpParam())` | OK |
| Agrupar | `groupChromPeaks(PeakDensityParam())` | 1 014 features |
| Gap-filling | `fillChromPeaks(ChromPeakAreaParam())` | 588 valores recuperados |
| Cuantificar | `quantify()` | matriz 1 014 × 4 |
| QC | `filterFeatures(PercentMissingFilter())` | 872 features retenidas |
| Anotar MS1 | `matchMz(Mass2MzParam())` | 5/872 con candidato |
| Historial | `processHistory()` | 5 pasos, parámetros trazados |

Todo reproducible con `./run.sh`.

## 3. Estructura

```
config/config.yml        parámetros del pipeline (nada hardcodeado)
R/utils.R                lectura de config, tabla de muestras, helpers
R/01_preprocess.R        pasos 1-10: importar → matriz de features + QC
R/02_annotate.R          pasos 11-12: anotación MS1 y (opcional) MS2
env/environment.yml      versiones exactas del entorno
env/setup.sh             instala el entorno sin root (~2-3 GB)
run.sh                   corre todo
data/                    <- aquí van los archivos del doctor
results/                 salidas (CSV, RDS, PNG, historial)
tests/data/faahKO/       datos de ejemplo usados para validar
xcms_lcmsms_sintesis.md  notas de investigación y correcciones
```

## 4. Cómo se corre

```bash
bash env/setup.sh          # una sola vez
./run.sh                   # con config/config.yml
./run.sh config/otro.yml   # con otra configuración
```

Salidas en `results/`:

- `feature_table.csv` — m/z, RT e intensidad por muestra (la tabla que se entrega)
- `features_se.rds` — `SummarizedExperiment`, para estadística posterior
- `xcms_result.rds` — objeto `XcmsExperiment` completo (permite re-anotar sin reprocesar)
- `annotation_ms1.csv` — candidatos por feature
- `process_history.txt` / `.rds` — parámetros exactos de cada paso
- `plots/` — BPC, similitud entre TIC, alineamiento de RT

## 5. Cómo adaptarlo a los datos del laboratorio

1. Convertir los archivos crudos a **mzML** (ProteoWizard `msconvert`, en Windows).
2. Copiarlos a `data/` y crear un CSV de muestras:

   ```csv
   file,sample_name,sample_group,sample_type
   muestra01.mzML,M01,control,study
   muestra02.mzML,M02,control,study
   pool01.mzML,QC01,qc,qc
   ```

3. Ajustar `config/config.yml` → `input.format: "mzml"`, `input.directory: "data"`,
   `input.metadata: "data/muestra.csv"`, `input.pattern: "[.]mzML$"`.
4. **Calibrar** `peak_detection` y `grouping` con estándares internos antes de
   correr todo (ver abajo).
5. `./run.sh`.

### Calibración (esto no se puede saltar)

Los valores del `config.yml` son los de los datos de prueba, **no** los de este
laboratorio. La forma correcta de fijarlos, según Metabonaut:

1. Elegir 2-3 iones conocidos (estándares internos; si no hay, iones abundantes
   del matriz) repartidos a lo largo del gradiente.
2. Extraer sus EIC: `chromatogram(data, mz = ..., rt = ...)`.
3. Probar `CentWaveParam` sobre esos EIC y quedarse con el que detecte **un pico
   limpio** por ion, con el ancho esperado.
4. Repetir después de alinear y de agrupar: si un compuesto conocido se parte en
   dos features, `bw`/`binSize` están mal.

## 6. Notas de versión (importante para reproducibilidad)

El entorno instalado es **R 4.4.3 / xcms 4.4.0** (lo más nuevo empaquetado en
bioconda). El artículo de 2025 usó **R 4.5.1 / Bioconductor 3.22 / xcms 4.7.3**,
y Metabonaut (rama devel) ya pide **R ≥ 4.6**.

Diferencias encontradas al ejecutar (el código las maneja o las documenta):

| API | En 4.4.0 | En devel |
|---|---|---|
| `chromatogram()` | devuelve `MChromatograms` | devuelve `Chromatograms` |
| `filterFeatures()` sobre `SummarizedExperiment` | no | sí |
| `sampleData()` sobre `SummarizedExperiment` | no (usar `colData()`) | — |
| `RtMatcherParam` | no existe (usar `Mass2MzRtParam`) | sí |
| `matchMz(target = CompDb)` | no (extraer tabla primero) | sí |

Si se necesita clavar las versiones del artículo, hay que instalar R 4.5 y los
paquetes con `BiocManager::install(version = "3.22")` compilando desde fuente.

## 7. Limitaciones honestas

- **Tipo de experimento.** Todo esto asume *untargeted* con alta resolución
  (QTOF/Orbitrap). Si el experimento es *targeted* con triple cuadrupolo
  (MRM/SRM), xcms no es la herramienta: se usa Skyline. **Confirmar antes de
  invertir tiempo.**
- **MS2.** El paso de anotación MS2 está implementado pero **no se probó**: los
  datos de ejemplo no traen espectros MS2. Requiere una biblioteca (MGF/MSP) o
  estándares medidos en el mismo equipo.
- **Estadística.** El pipeline llega hasta la matriz de features y su QC. ANOVA,
  PCA/PLS-DA y figuras no están incluidos todavía: dependen del diseño
  experimental, que aún no se conoce.
- **DIA/SWATH.** Para datos DIA hay que cambiar el paso 3 por
  `findChromPeaksIsolationWindow()` + `reconstructChromPeakSpectra()`. No está
  implementado.

## 8. Qué necesitamos del doctor

Para dejar de adivinar:

1. **Los archivos crudos** (o 2-3 de muestra) y el formato/equipo de origen.
2. **Modo de adquisición**: DDA, DIA/SWATH o MSe.
3. **Diseño experimental**: grupos, réplicas, si hay QC agrupados, blancos,
   estándares internos.
4. **Qué entregable espera**: tabla de features, metabolitos identificados,
   comparación estadística entre grupos, figuras para publicación.
5. **La carpeta del curso** que le dieron: para comparar con este pipeline y no
   replicar algo distinto de lo que vio.

## 9. Referencias

- Louail P, *et al.* xcms in Peak Form. *Anal Chem* 2025;97:27639–27645.
  doi:10.1021/acs.analchem.5c04338
- Rainer J, *et al.* A Modular and Expandable Ecosystem for Metabolomics Data
  Annotation in R. *Metabolites* 2022;12(2):173. doi:10.3390/metabo12020173
- xcms (devel 4.11.5): https://github.com/sneumann/xcms
- Metabonaut, flujo completo: https://rformassspectrometry.github.io/Metabonaut/
- Tutoriales: https://jorainer.github.io/xcmsTutorials/
