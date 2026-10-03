leer_serie <- function(x, fuente, unidad) {
  if (is.ts(x)) {
    frecuencia_val <- frequency(x)
    
    if (!frecuencia_val %in% c(1, 4, 12)) {
      stop("Solo se admiten objetos ts anuales, trimestrales o mensuales.")
    }
    
    inicio <- start(x)
    año_inicio <- inicio[1]
    periodo_inicio <- inicio[2]
    
    if (frecuencia_val == 12) {
      if (periodo_inicio < 1 || periodo_inicio > 12) {
        stop("Período mensual de inicio inválido.")
      }
      
      fecha_inicio <- as.Date(
        sprintf("%04d-%02d-01", año_inicio, periodo_inicio)
      )
      fechas_esperadas <- seq.Date(
        from = fecha_inicio,
        by = "month",
        length.out = length(x)
      )
    } else if (frecuencia_val == 4) {
      if (periodo_inicio < 1 || periodo_inicio > 4) {
        stop("Período trimestral de inicio inválido.")
      }
      
      mes_inicio <- (periodo_inicio - 1) * 3 + 1
      fecha_inicio <- as.Date(
        sprintf("%04d-%02d-01", año_inicio, mes_inicio)
      )
      fechas_esperadas <- seq.Date(
        from = fecha_inicio,
        by = "3 months",
        length.out = length(x)
      )
    } else {
      if (periodo_inicio != 1) {
        stop("Período anual de inicio inválido.")
      }
      
      fecha_inicio <- as.Date(sprintf("%04d-01-01", año_inicio))
      fechas_esperadas <- seq.Date(
        from = fecha_inicio,
        by = "year",
        length.out = length(x)
      )
    }
    
    datos <- tibble::tibble(
      t = seq_len(length(x)),
      fecha = fechas_esperadas,
      y = as.numeric(x)
    )
    
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
    
    fechas <- as.Date(csv$fecha)
    valores <- suppressWarnings(as.numeric(csv$valor))
    
    if (anyNA(fechas)) {
      stop("La columna fecha contiene fechas inválidas.")
    }
    
    if (anyNA(valores)) {
      stop("La columna valor contiene valores no numéricos.")
    }
    
    diferencias <- as.numeric(diff(fechas))
    
    if (any(diferencias <= 0)) {
      stop("Las fechas no son crecientes.")
    }
    
    if (length(unique(diferencias)) == 1L) {
      dias <- diferencias[1]
      
      frecuencia_val <- if (dias == 1) {
        365
      } else if (dias == 7) {
        52
      } else if (dias >= 28 && dias <= 31) {
        12
      } else if (dias >= 89 && dias <= 92) {
        4
      } else if (dias >= 365 && dias <= 366) {
        1
      } else {
        365 / dias
      }
    } else {
      stop("Las fechas no están equiespaciadas.")
    }
    
    datos <- tibble::tibble(
      t = seq_len(nrow(csv)),
      fecha = fechas,
      y = valores
    )
  }
  
  attr(datos, "frecuencia") <- frecuencia_val
  attr(datos, "fuente") <- fuente
  attr(datos, "unidad") <- unidad
  
  datos
}
