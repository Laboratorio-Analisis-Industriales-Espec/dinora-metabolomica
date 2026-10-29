## Utilidades compartidas por los scripts del pipeline.
suppressPackageStartupMessages({
    library(yaml)
})

## Lee la configuración y la devuelve como lista.
read_config <- function(path = "config/config.yml") {
    if (!file.exists(path))
        stop("No existe el archivo de configuración: ", path, call. = FALSE)
    cfg <- yaml::read_yaml(path)

    ## valores null de YAML -> NULL real
    cfg$rt$min <- if (is.null(cfg$rt$min)) NA_real_ else as.numeric(cfg$rt$min)
    cfg$rt$max <- if (is.null(cfg$rt$max)) NA_real_ else as.numeric(cfg$rt$max)

    outdir <- cfg$output$dir
    dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
    cfg
}

## Lista los archivos de entrada y construye la tabla de muestras.
## Si no hay metadata.csv, deduce el grupo del nombre de archivo/carpeta.
build_sample_table <- function(cfg) {
    files <- list.files(cfg$input$directory, pattern = cfg$input$pattern,
                        full.names = TRUE, recursive = TRUE)

    if (!is.null(cfg$input$metadata) && nzchar(cfg$input$metadata)) {
        pd <- read.csv(cfg$input$metadata, stringsAsFactors = FALSE)
        if (!all(c("file", "sample_name", "sample_group") %in% names(pd)))
            stop("La metadata debe tener las columnas: file, sample_name, sample_group",
                 call. = FALSE)
        pd$file <- ifelse(grepl("^(/|\\.)", pd$file), pd$file,
                          file.path(cfg$input$directory, pd$file))
        return(pd)
    }

    if (!length(files))
        stop("No se encontraron archivos en '", cfg$input$directory,
             "' con el patrón '", cfg$input$pattern, "'", call. = FALSE)

    ## Heurística: el nombre de la carpeta padre es el grupo; si no, "sample".
    group <- basename(dirname(files))
    group[group %in% c(".", cfg$input$directory)] <- "sample"

    data.frame(
        file         = files,
        sample_name  = tools::file_path_sans_ext(basename(files)),
        sample_group = group,
        sample_type  = "study",
        stringsAsFactors = FALSE
    )
}

## Imprime un encabezado de paso.
step <- function(n, txt) {
    cat(sprintf("\n[%s] %s\n", n, txt))
    utils::flush.console()
}

## Guarda un gráfico base-R en PNG.
save_png <- function(expr, file, width = 1600, height = 1200, res = 150) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off())
    force(expr)
}
