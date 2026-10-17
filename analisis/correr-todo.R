# Ejecutar desde la raiz del repositorio con:
# source("analisis/correr-todo.R")

archivos_requeridos <- c(
  file.path("R", "00-lectura.R"),
  file.path("R", "01-graficos.R"),
  file.path("R", "02-operadores.R"),
  file.path("R", "03-dependencia.R"),
  file.path("R", "04-metodos.R"),
  file.path("R", "05-evaluacion.R"),
  file.path("datos", "STUW9J-serie-1.csv"),
  file.path("datos", "STUW9J-serie-2.csv"),
  file.path("datos", "STUW9J-serie-3.csv")
)

faltantes <- archivos_requeridos[!file.exists(archivos_requeridos)]
if (length(faltantes) > 0L) {
  stop(
    "Ejecute este script desde la raiz del repositorio. No se encontraron: ",
    paste(faltantes, collapse = ", "),
    call. = FALSE
  )
}

paquetes_requeridos <- c("ggplot2", "patchwork", "tibble", "dplyr")
paquetes_faltantes <- paquetes_requeridos[
  !vapply(paquetes_requeridos, requireNamespace, logical(1), quietly = TRUE)
]
if (length(paquetes_faltantes) > 0L) {
  stop(
    "Faltan paquetes de R requeridos: ",
    paste(paquetes_faltantes, collapse = ", "),
    ". Instálelos antes de ejecutar el script.",
    call. = FALSE
  )
}

for (archivo in archivos_requeridos[seq_len(6L)]) {
  source(archivo)
}

dir.create("figs", showWarnings = FALSE, recursive = TRUE)

guardar_figura <- function(nombre, figura, ancho = 8, alto = 5) {
  ggplot2::ggsave(
    filename = file.path("figs", nombre),
    plot = figura,
    width = ancho,
    height = alto,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

series <- list(
  list(
    id = "serie-1",
    archivo = archivos_requeridos[7L],
    unidad = "casos",
    titulo = "Casos mensuales de enfermedad diarreica aguda notificados",
    rezago = 12L
  ),
  list(
    id = "serie-2",
    archivo = archivos_requeridos[8L],
    unidad = "afiliados",
    titulo = "Afiliados activos al programa de control de hipertension",
    rezago = 6L
  ),
  list(
    id = "serie-3",
    archivo = archivos_requeridos[9L],
    unidad = "consultas",
    titulo = "Consultas semanales por enfermedad respiratoria",
    rezago = 6L
  )
)

for (serie in series) {
  message("Procesando ", serie$id, "...")
  datos <- leer_serie(
    x = serie$archivo,
    unidad = serie$unidad,
    fuente = ""
  )

  guardar_figura(
    paste0(serie$id, "-serie-tiempo.png"),
    graficar_serie(datos, titulo = serie$titulo)
  )
  guardar_figura(
    paste0(serie$id, "-rezagos-", serie$rezago, ".png"),
    grafico_rezagos(datos, k = serie$rezago)
  )

  m <- min(24L, floor(length(datos$y) / 4))
  guardar_figura(
    paste0(serie$id, "-correlograma.png"),
    correlograma(datos, m = m),
    ancho = 10,
    alto = 7
  )

  if (serie$id == "serie-1") {
    guardar_figura(
      "serie-1-estacional-multiplicativo.png",
      grafico_estacional(datos, s = 12, tipo = "multiplicativo"),
      ancho = 10,
      alto = 8
    )
  }
}

message("Proceso terminado. Las figuras quedaron guardadas en figs/.")

writeLines(capture.output(sessionInfo()), "sesion-info.txt")
