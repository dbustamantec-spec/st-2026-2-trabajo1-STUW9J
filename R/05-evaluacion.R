optimizar <- function(y, metodo, rejilla = NULL, s = NULL) {
  .validar_y(y, 2L)
  metodos <- c("mm", "dmm", "ses", "holt", "hw", "hw_mult")
  if (length(metodo) != 1L || is.na(metodo) || !metodo %in% metodos) {
    stop(sprintf("metodo debe ser uno de: %s.",
                 paste(metodos, collapse = ", ")), call. = FALSE)
  }
  if (metodo %in% c("hw", "hw_mult")) {
    if (is.null(s)) stop("Debe proporcionar s para Holt-Winters.",
                         call. = FALSE)
    s <- .validar_entero(s, "s", 2L)
  } else if (!is.null(s)) {
    stop("El argumento s solo se usa con hw o hw_mult.", call. = FALSE)
  }
  if (is.null(rejilla)) rejilla <- .rejilla_por_defecto(metodo)
  
  columnas <- switch(
    metodo,
    mm = "k", dmm = "k", ses = "alpha",
    holt = c("alpha", "beta"),
    hw = c("alpha", "beta", "gamma"),
    hw_mult = c("alpha", "beta", "gamma")
  )
  if (metodo %in% c("mm", "dmm", "ses")) {
    if (is.data.frame(rejilla)) {
      if (!columnas %in% names(rejilla)) {
        stop(sprintf("rejilla debe tener la columna %s.", columnas),
             call. = FALSE)
      }
      # CORRECCIÓN: Agregado drop = FALSE para que no se convierta en vector
      rejilla <- rejilla[, columnas, drop = FALSE] 
    } else {
      if (!is.numeric(rejilla) || !is.null(dim(rejilla)) ||
          length(rejilla) == 0L) {
        stop("rejilla debe ser un vector numerico no vacio.", call. = FALSE)
      }
      rejilla <- as.data.frame(setNames(list(rejilla), columnas))
    }
  } else if (!is.data.frame(rejilla) ||
             !all(columnas %in% names(rejilla)) || nrow(rejilla) == 0L) {
    stop(sprintf("rejilla debe ser un data.frame con: %s.",
                 paste(columnas, collapse = ", ")), call. = FALSE)
  }
  
  rejilla <- rejilla[, columnas, drop = FALSE]
  if (anyNA(rejilla) ||
      any(!vapply(rejilla, is.numeric, logical(1))) ||
      any(!vapply(rejilla, function(x) all(is.finite(x)), logical(1)))) {
    stop("Todos los parametros de rejilla deben ser numericos y finitos.",
         call. = FALSE)
  }
  if (anyDuplicated(rejilla)) {
    stop("rejilla no debe contener combinaciones duplicadas.", call. = FALSE)
  }
  
  if (metodo %in% c("mm", "dmm")) {
    minimo_k <- if (metodo == "mm") 2 else 3
    if (any(rejilla$k != round(rejilla$k)) ||
        any(rejilla$k < minimo_k)) {
      stop(sprintf("k debe ser entero y >= %d para %s.",
                   minimo_k, metodo), call. = FALSE)
    }
  } else {
    constantes <- unlist(rejilla, use.names = FALSE)
    if (any(constantes <= 0 | constantes >= 1)) {
      stop("Las constantes de suavizamiento deben estar entre 0 y 1.",
           call. = FALSE)
    }
  }
  
  # Construcción de modelos con manejo de errores (tryCatch)
  modelos <- lapply(seq_len(nrow(rejilla)), function(i) {
    tryCatch({
      .crear_metodo_rejilla(y, metodo, rejilla[i, , drop = FALSE], s)
    }, error = function(e) {
      list(yhat = rep(NA_real_, length(y)))
    })
  })
  
  # 1. Identificar cuáles modelos no fallaron
  modelos_validos <- vapply(modelos, function(m) any(is.finite(m$yhat)), logical(1))
  
  if (!any(modelos_validos)) {
    stop("Ninguna combinacion de la rejilla produjo un modelo valido.", call. = FALSE)
  }
  
  # 2. Calcular la máscara común SOLO con los modelos válidos
  mascara_comun <- Reduce(
    `&`,
    lapply(modelos[modelos_validos], function(modelo) {
      is.finite(modelo$yhat) & is.finite(y)
    })
  )
  
  if (!any(mascara_comun)) {
    stop("No hay origenes comunes con pronostico valido para comparar la rejilla.", call. = FALSE)
  }
  
  # 3. Calcular el MSE 
  mse <- vapply(seq_along(modelos), function(i) {
    if (!modelos_validos[i]) {
      return(Inf)
    }
    mean((y[mascara_comun] - modelos[[i]]$yhat[mascara_comun])^2)
  }, numeric(1))
  
  tabla <- cbind(rejilla, mse = mse)
  optimo <- tabla[which.min(tabla$mse), , drop = FALSE]
  
  list(
    rejilla = tabla,
    optimo = optimo,
    metodo = metodo,
    s = s,
    indices_mse = which(mascara_comun),
    numero_origenes = sum(mascara_comun)
  )
}
 


# origne
origen_movil <- function(y, metodo, h, n0, ...) {
  # n0: tamaño del tramo de estimación
  # h: horizonte de pronóstico a evaluar
  # metodo: función del método (ej. ajustar_holt, ajustar_ses)
  # ... : parámetros fijos elegidos previamente por optimizar()
  
  T_total <- length(y)
  origenes <- n0:(T_total - h)
  
  # Vectores para almacenar resultados
  reales <- numeric(length(origenes))
  predicciones <- numeric(length(origenes))
  errores <- numeric(length(origenes))
  
  for (i in seq_along(origenes)) {
    o <- origenes[i]
    
    # El método recibe estrictamente la información hasta el origen 'o'
    y_o <- y[1:o]
    
    # Se ajusta el modelo con la historia recortada y los parámetros fijos
    modelo <- metodo(y_o, ...)
    
    # Se pronostica hasta el horizonte 'h'
    pronosticos <- modelo$pronosticar(h)
    
    # Se extrae la predicción específica para el paso 'h'
    pred <- pronosticos[h]
    real <- y[o + h]
    
    predicciones[i] <- pred
    reales[i] <- real
    errores[i] <- real - pred
  }
  
  # Retorna una lista con la información completa de la evaluación
  list(
    origenes = origenes,
    reales = reales,
    predicciones = predicciones,
    errores = errores
  )
}

# medidas
medidas <- function(reales, predicciones) {
  errores <- reales - predicciones
  
  # RMSE: Root Mean Squared Error
  rmse <- sqrt(mean(errores^2, na.rm = TRUE))
  
  # MAD: Mean Absolute Deviation
  mad <- mean(abs(errores), na.rm = TRUE)
  
  # MAPE: Mean Absolute Percentage Error (expresado en porcentaje)
  mape <- mean(abs(errores / reales), na.rm = TRUE) * 100
  
  # Retorna un vector con los nombres para facilitar la creación de tablas
  c(RMSE = rmse, MAD = mad, MAPE = mape)
}
