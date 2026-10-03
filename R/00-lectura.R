leer_serie <- function(x, fuente, unidad) {
  stopifnot(is.character(fuente), is.character(unidad))
  
  # ---- Caso 1: objeto ts ---------------------------------------------------
  if (is.ts(x)) {
    frecuencia_val <- frequency(x)
    
    if (!frecuencia_val %in% c(1, 4, 12, 52, 365)) {
      stop("Solo se admiten objetos ts anuales, trimestrales, mensuales, semanales (52) o diarios (365).")
    }
    
    inicio <- start(x)
    anio_inicio <- inicio[1]
    periodo_inicio <- if (length(inicio) > 1L) inicio[2] else 1
    
    if (frecuencia_val == 365) {
      # Período = día del año
      if (periodo_inicio < 1 || periodo_inicio > 366) {
        stop("Período diario de inicio inválido.")
      }
      fecha_inicio <- as.Date(sprintf("%04d-01-01", anio_inicio)) + (periodo_inicio - 1)
      fechas <- seq.Date(fecha_inicio, by = "day", length.out = length(x))
      
    } else if (frecuencia_val == 52) {
      # Período = semana del año (semana 1 empieza el 1 de enero)
      if (periodo_inicio < 1 || periodo_inicio > 53) {
        stop("Período semanal de inicio inválido.")
      }
      fecha_inicio <- as.Date(sprintf("%04d-01-01", anio_inicio)) + 7 * (periodo_inicio - 1)
      fechas <- seq.Date(fecha_inicio, by = "week", length.out = length(x))
      
    } else if (frecuencia_val == 12) {
      if (periodo_inicio < 1 || periodo_inicio > 12) {
        stop("Período mensual de inicio inválido.")
      }
      fecha_inicio <- as.Date(sprintf("%04d-%02d-01", anio_inicio, periodo_inicio))
      fechas <- seq.Date(fecha_inicio, by = "month", length.out = length(x))
      
    } else if (frecuencia_val == 4) {
      if (periodo_inicio < 1 || periodo_inicio > 4) {
        stop("Período trimestral de inicio inválido.")
      }
      mes_inicio <- (periodo_inicio - 1) * 3 + 1
      fecha_inicio <- as.Date(sprintf("%04d-%02d-01", anio_inicio, mes_inicio))
      fechas <- seq.Date(fecha_inicio, by = "3 months", length.out = length(x))
      
    } else {
      if (periodo_inicio != 1) {
        stop("Período anual de inicio inválido.")
      }
      fecha_inicio <- as.Date(sprintf("%04d-01-01", anio_inicio))
      fechas <- seq.Date(fecha_inicio, by = "year", length.out = length(x))
    }
    
    datos <- tibble::tibble(
      t     = seq_len(length(x)),
      fecha = fechas,
      y     = as.numeric(x)
    )
    
    # ---- Caso 2: ruta a un CSV -----------------------------------------------
  } else {
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
      stop("x debe ser un objeto ts o una ruta válida a un archivo CSV.")
    }
    
    csv <- utils::read.csv(x, stringsAsFactors = FALSE)
    
    if (!all(c("fecha", "valor") %in% names(csv))) {
      stop("El CSV debe contener las columnas 'fecha' y 'valor'.")
    }
    
    if (nrow(csv) < 2L) {
      stop("El CSV debe contener al menos dos observaciones.")
    }
    
    fechas <- tryCatch(
      as.Date(csv$fecha),
      error = function(e) rep(as.Date(NA), nrow(csv))
    )
    valores <- suppressWarnings(as.numeric(csv$valor))
    
    if (anyNA(fechas)) {
      stop("La columna fecha contiene fechas inválidas.")
    }
    
    if (anyNA(valores)) {
      stop("La columna valor contiene valores no numéricos.")
    }
    
    dif_dias <- as.numeric(diff(fechas))
    
    if (any(dif_dias <= 0)) {
      stop("Las fechas no son crecientes.")
    }
    
    if (length(unique(dif_dias)) == 1L && dif_dias[1] %in% c(1, 7)) {
      # Diaria o semanal: se verifica en días
      frecuencia_val <- if (dif_dias[1] == 1) 365 else 52
    } else {
      # Mensual, trimestral o anual: se verifica en meses de calendario
      lt <- as.POSIXlt(fechas)
      meses <- lt$year * 12 + lt$mon
      dif_meses <- unique(diff(meses))
      
      if (length(dif_meses) != 1L || !dif_meses %in% c(1, 3, 12)) {
        stop("Las fechas no están equiespaciadas según una frecuencia admitida (diaria, semanal, mensual, trimestral o anual).")
      }
      
      frecuencia_val <- 12 / dif_meses  # 12, 4 o 1
    }
    
    datos <- tibble::tibble(
      t     = seq_len(nrow(csv)),
      fecha = fechas,
      y     = valores
    )
  }
  
  attr(datos, "frecuencia") <- frecuencia_val
  attr(datos, "fuente")     <- fuente
  attr(datos, "unidad")     <- unidad
  
  datos
}
