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

errores_ganador <- function(y, metodo_ganador, n0) {
  eval <- origen_movil(y = y, metodo = metodo_ganador, h = 1, n0 = n0)
  data.frame(
    t = (n0 + 1):length(y),
    real = eval$reales,
    pronostico = eval$predicciones,
    error = eval$reales - eval$predicciones
  )
}

validar_ganador <- function(errores, titulo) {
  e <- errores$error
  n <- length(e)
  t_media <- t.test(e, mu = 0)
  
  m <- min(24L, floor(n / 4))
  r <- acf_muestral(e, m = m)
  Q <- ljung_box(r, T = n, m = m)
  banda <- 1.96 / sqrt(n)
  
  dw <- sum(diff(e)^2) / sum(e^2)
  centrados <- e - mean(e)
  momento_2 <- mean(centrados^2)
  asimetria <- mean(centrados^3) / momento_2^(3 / 2)
  curtosis <- mean(centrados^4) / momento_2^2
  jb <- n * (asimetria^2 / 6 + (curtosis - 3)^2 / 24)
  p_jb <- pchisq(jb, df = 2, lower.tail = FALSE)
  
  gl_t <- unname(t_media$parameter)
  critico_t <- qt(0.975, df = gl_t)
  critico_jb <- qchisq(0.95, df = 2)
  
  tabla_validacion <- data.frame(
    Prueba = c(
      "t de media cero",
      paste0("Ljung–Box ($Q$, $m=", m, ")"),
      "Durbin–Watson ($d$)",
      "Jarque–Bera ($JB$)"
    ),
    `Grados de libertad` = c(
      as.character(gl_t),
      as.character(Q$grados_libertad),
      "No aplica",
      "2"
    ),
    Estadístico = c(
      round(t_media$statistic, 4),
      round(Q$estadistico, 4),
      round(dw, 4),
      round(jb, 4)
    ),
    `Valor crítico (5%)` = c(
      paste0("±", round(critico_t, 4)),
      round(Q$valor_critico_5, 4),
      "Referencia: 2",
      round(critico_jb, 4)
    ),
    `Valor p` = c(
      round(t_media$p.value, 4),
      round(Q$p_valor, 4),
      "No aplica",
      round(p_jb, 4)
    ),
    Decisión = c(
      ifelse(t_media$p.value < 0.05, "Rechazar H₀", "No rechazar H₀"),
      ifelse(Q$p_valor < 0.05, "Rechazar H₀: autocorrelación conjunta",
             "No rechazar H₀"),
      ifelse(dw < 2, "Menor que 2: indicio de autocorrelación positiva",
             ifelse(dw > 2, "Mayor que 2: indicio de autocorrelación negativa",
                    "Igual a 2: sin indicio de autocorrelación de orden 1")),
      ifelse(p_jb < 0.05, "Rechazar H₀", "No rechazar H₀")
    )
  )
  
  list(
    tabla = tabla_validacion,
    errores = e,
    media = mean(e),
    desv_est = sd(e),
    acf = r,
    banda = banda,
    m = m,
    dw = dw,
    t_p = t_media$p.value,
    lb_p = Q$p_valor,
    jb_p = p_jb
  )
}



alpha_s1 <- opt_hw_mult_s1$optimo$alpha
beta_s1  <- opt_hw_mult_s1$optimo$beta
gamma_s1 <- opt_hw_mult_s1$optimo$gamma
alpha_s3 <- opt_se_s3$optimo$alpha 

ganador_s1 <- function(y) {
  ajustar_hw_mult(y, 
                  s = 12,
                  alpha = alpha_s1,
                  beta  = beta_s1,
                  gamma = gamma_s1)
}

ganador_s2 <- function(y) {
  ajustar_tendencia(y,
                    tipo = "cuadratica")
}


ganador_s3 <- function(y) {
  ajustar_ses(y, alpha = alpha_s3)
}

err_s1 <- errores_ganador(datos1$y, ganador_s1, n0 = n0_1)
val_s1 <- validar_ganador(err_s1, "Validación — Serie 1")

err_s2 <- errores_ganador(datos2$y, ganador_s2, n0 = n0_2)
val_s2 <- validar_ganador(err_s2, "Validación — Serie 2")

err_s3 <- errores_ganador(datos3$y, ganador_s3, n0 = n0_3)
val_s3 <- validar_ganador(err_s3, "Validación — Serie 3")


mostrar_tabla <- function(resultado, titulo) {
  cat("###", titulo, "\n\n")
  print(resultado$tabla, row.names = FALSE)
  cat("\n")
}
