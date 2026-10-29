#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## PASO 11-12: anotación de features
##   11 anotación MS1 (m/z + RT opcional) contra una tabla de masas exactas
##   12 (opcional) anotación MS2 contra una biblioteca de espectros (MGF/MSP)
##
## Uso:  Rscript R/02_annotate.R [config/config.yml]
## Requiere: haber corrido 01_preprocess.R
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(MetaboAnnotation)
    library(MetaboCoreUtils)
    library(SummarizedExperiment)
})

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg  <- read_config(if (length(args)) args[1] else "config/config.yml")
outdir <- cfg$output$dir

ann_cfg <- cfg$annotation
if (is.null(ann_cfg) || is.null(ann_cfg$database) || !nzchar(ann_cfg$database)) {
    stop("Falta la sección 'annotation: database:' en config.yml", call. = FALSE)
}
if (!file.exists(file.path(outdir, "features_se.rds")))
    stop("No existe results/features_se.rds. Corre primero R/01_preprocess.R",
         call. = FALSE)

res <- readRDS(file.path(outdir, "features_se.rds"))
fd  <- as.data.frame(rowData(res))
ft  <- as.data.frame(assay(res))

## Tabla de consulta: id + mz (+ rt si está disponible)
query <- data.frame(
    feature_id = rownames(ft),
    mz         = fd$mzmed
)
if (!is.null(fd$rtmed)) query$rt <- fd$rtmed

## ---------------------------------------------------------------------------
step(11, "Anotación MS1")
db <- read.csv(ann_cfg$database, stringsAsFactors = FALSE)
if (!"exactmass" %in% names(db))
    stop("La base de referencia necesita una columna 'exactmass' (masa neutra)",
         call. = FALSE)
cat("Base de referencia:", nrow(db), "compuestos |", ann_cfg$database, "\n")

## msLevel = 1 -> solo features con señal MS1; el aducto se aplica a la masa
param <- if (!is.null(db$rt)) {
    Mass2MzRtParam(adduct = ann_cfg$adduct, ppm = ann_cfg$ppm,
                   toleranceRt = ann_cfg$rt_tolerance)
} else {
    Mass2MzParam(adduct = ann_cfg$adduct, ppm = ann_cfg$ppm)
}

m <- matchMz(query, db, param)
md <- as.data.frame(matchedData(m))
n_hit <- sum(!is.na(md$target_name))
cat("Features con al menos un candidato:", n_hit, "de", nrow(query), "\n")

write.csv(md, file.path(outdir, "annotation_ms1.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
step(12, "Anotación MS2 (opcional)")
ms2_lib <- ann_cfg$ms2_library
if (is.null(ms2_lib) || !nzchar(ms2_lib) || !file.exists(ms2_lib)) {
    cat("Sin biblioteca MS2 configurada ('annotation: ms2_library:'). Se omite.\n")
} else {
    suppressPackageStartupMessages({
        library(Spectra)
        library(MsBackendMgf)
    })
    xdata <- readRDS(file.path(outdir, "xcms_result.rds"))
    qry   <- featureSpectra(xdata, return.type = "Spectra", msLevel = 2L)
    ref   <- Spectra(ms2_lib, source = MsBackendMgf())
    csp   <- CompareSpectraParam(ppm = ann_cfg$ppm, tolerance = 0.05,
                                 THRESHFUN = function(x) x >= ann_cfg$ms2_score)
    hits  <- matchSpectra(qry, ref, param = csp)
    write.csv(as.data.frame(matchedData(hits)),
              file.path(outdir, "annotation_ms2.csv"), row.names = FALSE)
    cat("Espectros MS2 consultados:", length(qry), "| coincidencias:",
        sum(!is.na(matchedData(hits)$target_name)), "\n")
}

cat("\nListo. Anotación en:", normalizePath(outdir), "\n")
