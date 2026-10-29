#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## PASO 1-9: preprocesamiento de datos LC-MS(/MS) con xcms
##   1 importación  2 inspección  3 detección  4 refinamiento  5 alineamiento
##   6 correspondencia  7 gap-filling  8 matriz de features  9 QC de features
##
## Uso:  Rscript R/01_preprocess.R [config/config.yml]
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(xcms)
    library(MsExperiment)
    library(SummarizedExperiment)
    library(BiocParallel)
})

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_config(if (length(args)) args[1] else "config/config.yml")
outdir <- cfg$output$dir

## Paralelización: deja 1 núcleo libre para el sistema.
n_cores <- max(1, parallel::detectCores() - 1)
register(MulticoreParam(workers = n_cores))
cat("Núcleos usados:", n_cores, "\n")

## ---------------------------------------------------------------------------
step(1, "Importación de datos")
pd <- build_sample_table(cfg)
cat("Archivos:", nrow(pd), "\n")
print(pd[, c("sample_name", "sample_group")])

data <- readMsExperiment(spectraFiles = pd$file, sampleData = pd)
cat("Espectros totales:", length(spectra(data)), "\n")

## ---------------------------------------------------------------------------
step(2, "Inspección: BPC / TIC")
bpi <- chromatogram(data, aggregationFun = "max")
save_png(plot(bpi, main = "BPC (max) por muestra"),
         file.path(outdir, "plots", "01_bpc.png"))

## Heatmap de similitud entre TICs (detecta corridas problemáticas)
if (nrow(pd) >= 3) {
    tic <- chromatogram(data, aggregationFun = "sum")
    mat <- tryCatch({
        binned <- bin(tic, binSize = 2)
        do.call(cbind, lapply(seq_along(binned), function(i) intensity(binned[[i]])))
    }, error = function(e) NULL)

    if (!is.null(mat)) {
        sim <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))
        rownames(sim) <- colnames(sim) <- pd$sample_name
        ann <- data.frame(group = pd$sample_group, row.names = pd$sample_name)
        save_png(pheatmap::pheatmap(sim, annotation_col = ann,
                                    main = "Similitud entre TIC (Pearson)"),
                 file.path(outdir, "plots", "02_tic_similarity.png"))
        cat("Heatmap de similitud escrito.\n")
    }
}

## Recorte opcional del rango de retención
if (!is.na(cfg$rt$min) && !is.na(cfg$rt$max)) {
    data <- filterRt(data, rt = c(cfg$rt$min, cfg$rt$max))
    cat("Datos recortados a RT", cfg$rt$min, "-", cfg$rt$max, "s\n")
}

## ---------------------------------------------------------------------------
step(3, "Detección de picos cromatográficos (centWave, MS1)")
cwp <- CentWaveParam(
    ppm       = cfg$peak_detection$ppm,
    peakwidth = cfg$peak_detection$peakwidth,
    snthresh  = cfg$peak_detection$snthresh,
    noise     = cfg$peak_detection$noise,
    prefilter = cfg$peak_detection$prefilter,
    mzdiff    = cfg$peak_detection$mzdiff
)
data <- findChromPeaks(data, param = cwp, msLevel = 1L)
cat("Picos detectados:", nrow(chromPeaks(data)), "\n")

step(4, "Refinamiento de picos (unión de picos vecinos)")
data <- refineChromPeaks(data, param = MergeNeighboringPeaksParam(
    expandRt = cfg$refine$expandRt,
    expandMz = cfg$refine$expandMz
))
cat("Picos tras refinamiento:", nrow(chromPeaks(data)), "\n")

## ---------------------------------------------------------------------------
step(5, "Alineamiento de tiempos de retención")
if (identical(cfg$alignment$method, "peakgroups")) {
    ## Needs the correspondence first; see vignette "xcms.Rmd".
    stop("Método 'peakgroups' no implementado en este esqueleto.", call. = FALSE)
} else {
    data <- adjustRtime(data, param = ObiwarpParam(binSize = cfg$alignment$binSize))
}
save_png(plotAdjustedRtime(data),
         file.path(outdir, "plots", "03_rt_alignment.png"))

## ---------------------------------------------------------------------------
step(6, "Correspondencia entre muestras (agrupamiento de features)")
pdp <- PeakDensityParam(
    sampleGroups = sampleData(data)$sample_group,
    bw           = cfg$grouping$bw,
    minFraction  = cfg$grouping$minFraction,
    binSize      = cfg$grouping$binSize
)
data <- groupChromPeaks(data, param = pdp)
cat("Features:", nrow(featureDefinitions(data)), "\n")

step(7, "Gap-filling (recuperación de valores faltantes)")
data <- fillChromPeaks(data, param = ChromPeakAreaParam())
cat("Gap-filling terminado.\n")

## ---------------------------------------------------------------------------
step(8, "Matriz de features (SummarizedExperiment)")
res <- quantify(data, value = "into", method = "sum")
cat("Dimensiones:", nrow(res), "features x", ncol(res), "muestras\n")
cat("Valores faltantes:", sum(is.na(assay(res))), "\n")

## ---------------------------------------------------------------------------
step(9, "Control de calidad de features")
n_before <- nrow(res)
res_qc <- tryCatch(
    filterFeatures(res, filter = PercentMissingFilter(
        threshold = cfg$qc$max_percent_missing,
        f = colData(res)$sample_group
    )),
    error = function(e) {
        ## xcms < 4.5 solo acepta MsExperiment aquí -> filtrado manual
        message("filterFeatures() no soporta SummarizedExperiment en esta versión; ",
                "se aplica un filtro equivalente manual.")
        grp  <- colData(res)$sample_group
        miss <- sapply(split(seq_len(ncol(res)), grp), function(i)
            rowMeans(is.na(assay(res))[, i, drop = FALSE]) * 100)
        keep <- rowSums(miss > cfg$qc$max_percent_missing) == 0
        res[keep, ]
    }
)
cat("Features antes/después del QC:", n_before, "->", nrow(res_qc), "\n")

## ---------------------------------------------------------------------------
step(10, "Exportación")
ft <- as.data.frame(assay(res_qc))
fd <- as.data.frame(featureDefinitions(data))
fd <- fd[match(rownames(ft), rownames(fd)), c("mzmed", "rtmed")]
ft <- cbind(fd, ft)
write.csv(ft, file.path(outdir, "feature_table.csv"), row.names = TRUE)
saveRDS(data,    file.path(outdir, "xcms_result.rds"))
saveRDS(res_qc,  file.path(outdir, "features_se.rds"))

## Historial de procesamiento = trazabilidad de parámetros
ph <- processHistory(data)
writeLines(
    vapply(ph, function(p) paste(processType(p), "|", processDate(p)), character(1)),
    file.path(outdir, "process_history.txt")
)
saveRDS(ph, file.path(outdir, "process_history.rds"))

cat("\nListo. Resultados en:", normalizePath(outdir), "\n")
