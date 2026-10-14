# Metodos de pronostico del Trabajo 1.
# Cada ajuste devuelve yhat, pronosticar(h) y parametros.
# yhat[t] es un pronostico de origen t - 1; los calentamientos son NA.

.validar_y <- function(y, minimo = 1L) {
  if (!is.numeric(y) || !is.null(dim(y)) || length(y) < minimo ||
      anyNA(y) || any(!is.finite(y))) {
    stop(sprintf(
      "y debe ser un vector numerico finito sin NA y tener al menos %d observacion(es).",
      minimo
    ), call. = FALSE)
  }
  invisible(TRUE)
}

.validar_entero <- function(x, nombre, minimo = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != round(x) || x < minimo) {
    stop(sprintf("%s debe ser un entero mayor o igual a %d.", nombre, minimo),
         call. = FALSE)
  }
  as.integer(x)
}

.validar_h <- function(h) {
  .validar_entero(h, "h", 1L)
}

.validar_constante <- function(x, nombre) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop(sprintf("%s debe estar estrictamente entre 0 y 1.", nombre),
         call. = FALSE)
  }
  invisible(TRUE)
}

.resultado_metodo <- function(yhat, pronosticar, parametros) {
  list(yhat = as.numeric(yhat), pronosticar = pronosticar,
       parametros = parametros)
}

.diseno_tendencia <- function(t, tipo) {
  if (tipo == "cuadratica") {
    cbind(intercepto = 1, tendencia = t, cuadratica = t^2)
  } else {
    cbind(intercepto = 1, tendencia = t)
  }
}

.ajuste_ols <- function(X, y) {
  qr_x <- qr(X)
  if (qr_x$rank < ncol(X)) {
    stop("La matriz de regresores no tiene rango completo.", call. = FALSE)
  }
  beta <- as.numeric(qr.coef(qr_x, y))
  names(beta) <- colnames(X)
  ajustados <- as.numeric(X %*% beta)
  residuos <- as.numeric(y - ajustados)
  list(beta = beta, ajustados = ajustados, residuos = residuos)
}

.diagnosticos_ols <- function(X, y, ajuste) {
  n <- length(y)
  p <- ncol(X)
  if (n <= p) {
    stop("La regresion requiere mas observaciones que coeficientes.",
         call. = FALSE)
  }

  beta <- ajuste$beta
  e <- ajuste$residuos
  xtx_inv <- chol2inv(chol(crossprod(X)))
  sigma2 <- sum(e^2) / (n - p)
  ee_ols <- sqrt(pmax(diag(sigma2 * xtx_inv), 0))

  L <- min(floor(4 * (n / 100)^(2 / 9)), n - 1L)
  Xe <- X * e
  meat <- crossprod(Xe)
  if (L > 0L) {
    for (lag in seq_len(L)) {
      peso <- 1 - lag / (L + 1)
      gamma_l <- crossprod(
        Xe[(lag + 1L):n, , drop = FALSE],
        Xe[seq_len(n - lag), , drop = FALSE]
      )
      meat <- meat + peso * (gamma_l + t(gamma_l))
    }
  }
  vcov_hac <- xtx_inv %*% meat %*% xtx_inv
  ee_robusto <- sqrt(pmax(diag(vcov_hac), 0))
  t_val <- ifelse(ee_robusto > 0, beta / ee_robusto, NA_real_)
  p_val <- 2 * stats::pt(-abs(t_val), df = n - p)
  sst <- sum((y - mean(y))^2)
  r2 <- if (sst > 0) 1 - sum(e^2) / sst else NA_real_
  sse <- sum(e^2)
  dw <- if (sse > 0) sum(diff(e)^2) / sse else NA_real_

  list(
    tabla = data.frame(
      coeficiente = names(beta), estimacion = beta,
      ee_ols = ee_ols, ee_robusto = ee_robusto,
      t = t_val, p_valor = p_val, row.names = NULL
    ),
    sigma2 = sigma2,
    R2 = r2,
    dw = dw,
    L_hac = L
  )
}

# Media simple: media acumulada disponible antes de cada observacion.
ajustar_media <- function(y) {
  .validar_y(y)
  n <- length(y)
  yhat <- rep(NA_real_, n)
  if (n >= 2L) {
    acumuladas <- cumsum(y)
    yhat[2:n] <- acumuladas[seq_len(n - 1L)] / seq_len(n - 1L)
  }
  media_final <- mean(y)

  pronosticar <- function(h) {
    h <- .validar_h(h)
    rep(media_final, h)
  }

  .resultado_metodo(yhat, pronosticar, list(media = media_final))
}

# Ingenuo: el ultimo valor observado se repite.
ajustar_ingenuo <- function(y) {
  .validar_y(y)
  n <- length(y)
  yhat <- rep(NA_real_, n)
  if (n >= 2L) yhat[2:n] <- y[seq_len(n - 1L)]

  pronosticar <- function(h) {
    h <- .validar_h(h)
    rep(y[n], h)
  }

  .resultado_metodo(yhat, pronosticar, list(ultimo = y[n]))
}

# Ingenuo estacional: repite el ultimo ciclo observado.
ajustar_ingenuo_estacional <- function(y, s) {
  .validar_y(y)
  s <- .validar_entero(s, "s", 1L)
  n <- length(y)
  if (n <= s) stop("y debe contener mas de s observaciones.", call. = FALSE)

  yhat <- rep(NA_real_, n)
  yhat[(s + 1L):n] <- y[seq_len(n - s)]

  pronosticar <- function(h) {
    h <- .validar_h(h)
    ciclo <- tail(y, s)
    ciclo[((seq_len(h) - 1L) %% s) + 1L]
  }

  .resultado_metodo(yhat, pronosticar, list(s = s, ciclo_final = tail(y, s)))
}

# Media movil: pronostico como promedio de las k observaciones mas recientes.
ajustar_mm <- function(y, k) {
  .validar_y(y)
  k <- .validar_entero(k, "k", 2L)
  n <- length(y)
  if (n < k) stop("y debe tener al menos k observaciones.", call. = FALSE)

  yhat <- rep(NA_real_, n)
  if (n > k) {
    for (t in (k + 1L):n) {
      yhat[t] <- mean(y[(t - k):(t - 1L)])
    }
  }
  mm_final <- mean(tail(y, k))

  pronosticar <- function(h) {
    h <- .validar_h(h)
    rep(mm_final, h)
  }

  .resultado_metodo(yhat, pronosticar, list(k = k, mm_T = mm_final))
}

# Suavizamiento exponencial simple.
ajustar_ses <- function(y, alpha) {
  .validar_y(y)
  .validar_constante(alpha, "alpha")
  n <- length(y)
  yhat <- rep(NA_real_, n)
  nivel <- y[1L]

  if (n >= 2L) {
    for (t in 2:n) {
      yhat[t] <- nivel
      nivel <- alpha * y[t] + (1 - alpha) * nivel
    }
  }

  pronosticar <- function(h) {
    h <- .validar_h(h)
    rep(nivel, h)
  }

  .resultado_metodo(yhat, pronosticar,
                    list(alpha = alpha, nivel_final = nivel))
}

# Doble media movil, con nivel y pendiente de Brown.
ajustar_dmm <- function(y, k) {
  .validar_y(y)
  k <- .validar_entero(k, "k", 2L)
  n <- length(y)
  if (n < 2L * k) {
    stop("y debe tener al menos 2*k observaciones para ajustar la DMM.",
         call. = FALSE)
  }

  media_movil <- rep(NA_real_, n)
  for (t in k:n) {
    media_movil[t] <- mean(y[(t - k + 1L):t])
  }

  doble_media <- rep(NA_real_, n)
  inicio_doble <- 2L * k - 1L
  for (t in inicio_doble:n) {
    doble_media[t] <- mean(media_movil[(t - k + 1L):t])
  }

  nivel <- 2 * media_movil - doble_media
  pendiente <- (2 / (k - 1)) * (media_movil - doble_media)
  yhat <- rep(NA_real_, n)
  if (n >= inicio_doble) {
    for (t in inicio_doble:n) {
      if (t < n) yhat[t + 1L] <- nivel[t] + pendiente[t]
    }
  }

  nivel_final <- nivel[n]
  pendiente_final <- pendiente[n]
  pronosticar <- function(h) {
    h <- .validar_h(h)
    nivel_final + pendiente_final * seq_len(h)
  }

  .resultado_metodo(
    yhat, pronosticar,
    list(k = k, nivel = nivel, pendiente = pendiente)
  )
}

# Regresion de tendencia: los ajustes de un paso usan solo el prefijo anterior.
ajustar_tendencia <- function(y, tipo, corregir_sesgo = FALSE) {
  .validar_y(y)
  tipos <- c("lineal", "cuadratica", "exponencial")
  if (length(tipo) != 1L || is.na(tipo) || !tipo %in% tipos) {
    stop("tipo debe ser 'lineal', 'cuadratica' o 'exponencial'.",
         call. = FALSE)
  }
  if (!is.logical(corregir_sesgo) || length(corregir_sesgo) != 1L ||
      is.na(corregir_sesgo)) {
    stop("corregir_sesgo debe ser TRUE o FALSE.", call. = FALSE)
  }
  if (tipo == "exponencial" && any(y <= 0)) {
    stop("La tendencia exponencial requiere todos los valores de y > 0.",
         call. = FALSE)
  }

  n <- length(y)
  respuesta <- if (tipo == "exponencial") log(y) else y
  crear_x <- function(t) .diseno_tendencia(t, tipo)

  yhat <- rep(NA_real_, n)
  p <- ncol(crear_x(1))
  if (n >= p + 1L) {
    for (t in (p + 1L):n) {
      indices <- seq_len(t - 1L)
      ajuste_prev <- .ajuste_ols(crear_x(indices), respuesta[indices])
      pred <- sum(crear_x(t) * ajuste_prev$beta)
      if (tipo == "exponencial") {
        factor_sesgo <- 1
        if (corregir_sesgo && length(indices) > p) {
          factor_sesgo <- exp(sum(ajuste_prev$residuos^2) /
                                (2 * (length(indices) - p)))
        }
        pred <- exp(pred) * factor_sesgo
      }
      yhat[t] <- pred
    }
  }

  X <- crear_x(seq_len(n))
  ajuste_final <- .ajuste_ols(X, respuesta)
  diag_final <- .diagnosticos_ols(X, respuesta, ajuste_final)
  beta <- ajuste_final$beta
  sigma2 <- diag_final$sigma2
  pronosticar <- function(h) {
    h <- .validar_h(h)
    t_nuevo <- n + seq_len(h)
    pred <- as.numeric(crear_x(t_nuevo) %*% beta)
    if (tipo == "exponencial") {
      if (corregir_sesgo) pred <- pred + sigma2 / 2
      pred <- exp(pred)
    }
    pred
  }

  parametros <- c(
    list(tipo = tipo, coeficientes = beta),
    diag_final,
    list(corregir_sesgo = corregir_sesgo)
  )
  .resultado_metodo(yhat, pronosticar, parametros)
}

# Holt lineal con inicializacion L_1=y_1 y T_1=0.
ajustar_holt <- function(y, alpha, beta) {
  .validar_y(y, 2L)
  .validar_constante(alpha, "alpha")
  .validar_constante(beta, "beta")
  n <- length(y)
  yhat <- rep(NA_real_, n)
  nivel <- rep(NA_real_, n)
  tendencia <- rep(NA_real_, n)
  nivel[1L] <- y[1L]
  tendencia[1L] <- 0

  for (t in 2:n) {
    yhat[t] <- nivel[t - 1L] + tendencia[t - 1L]
    nivel[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
    tendencia[t] <- beta * (nivel[t] - nivel[t - 1L]) +
      (1 - beta) * tendencia[t - 1L]
  }

  nivel_final <- nivel[n]
  tendencia_final <- tendencia[n]
  pronosticar <- function(h) {
    h <- .validar_h(h)
    nivel_final + tendencia_final * seq_len(h)
  }

  .resultado_metodo(
    yhat, pronosticar,
    list(alpha = alpha, beta = beta, nivel = nivel,
         tendencia = tendencia)
  )
}

.diseno_tendencia_estacional <- function(t, s) {
  periodo <- ((t - 1L) %% s) + 1L
  X <- cbind(intercepto = 1, tendencia = t)
  if (s > 1L) {
    estacionales <- matrix(0, nrow = length(t), ncol = s - 1L)
    for (j in seq_len(s - 1L)) {
      estacionales[, j] <- ifelse(periodo == j, 1,
                                  ifelse(periodo == s, -1, 0))
    }
    colnames(estacionales) <- paste0("estacion_", seq_len(s - 1L))
    X <- cbind(X, estacionales)
  }
  X
}

# Regresion de tendencia con ficticias estacionales de suma cero.
ajustar_tendencia_estacional <- function(y, s) {
  .validar_y(y)
  s <- .validar_entero(s, "s", 2L)
  n <- length(y)
  p <- s + 1L
  if (n <= p) {
    stop("Se requieren mas de s+1 observaciones para la regresion estacional.",
         call. = FALSE)
  }

  yhat <- rep(NA_real_, n)
  for (t in (p + 1L):n) {
    indices <- seq_len(t - 1L)
    ajuste_prev <- .ajuste_ols(
      .diseno_tendencia_estacional(indices, s), y[indices]
    )
    yhat[t] <- sum(
      .diseno_tendencia_estacional(t, s) * ajuste_prev$beta
    )
  }

  X <- .diseno_tendencia_estacional(seq_len(n), s)
  ajuste <- .ajuste_ols(X, y)
  diagnosticos <- .diagnosticos_ols(X, y, ajuste)
  beta <- ajuste$beta
  pronosticar <- function(h) {
    h <- .validar_h(h)
    t_nuevo <- n + seq_len(h)
    as.numeric(.diseno_tendencia_estacional(t_nuevo, s) %*% beta)
  }

  parametros <- c(
    list(s = s, coeficientes = beta),
    diagnosticos
  )
  .resultado_metodo(yhat, pronosticar, parametros)
}

.validar_hw <- function(y, s, alpha, beta, gamma, multiplicativo = FALSE) {
  .validar_y(y)
  s <- .validar_entero(s, "s", 2L)
  .validar_constante(alpha, "alpha")
  .validar_constante(beta, "beta")
  .validar_constante(gamma, "gamma")
  if (length(y) < 2L * s) {
    stop("Holt-Winters requiere al menos dos ciclos completos (2*s observaciones).",
         call. = FALSE)
  }
  if (multiplicativo && any(y <= 0)) {
    stop("Holt-Winters multiplicativo requiere todos los valores de y > 0.",
         call. = FALSE)
  }
  s
}

.stop_hw_mult_inestable <- function() {
  stop(structure(
    list(
      message = paste(
        "La combinacion de alpha, beta y gamma hizo que el nivel",
        "o el indice estacional multiplicativo dejara de ser positivo."
      ),
      call = NULL
    ),
    class = c("hw_mult_inestable", "error", "condition")
  ))
}

.inicializar_holt_winters <- function(y, s, multiplicativo) {
  media_1 <- mean(y[seq_len(s)])
  media_2 <- mean(y[s + seq_len(s)])
  tendencia_2s <- (media_2 - media_1) / s

  if (multiplicativo) {
    est_1 <- y[seq_len(s)] / media_1
    est_2 <- y[s + seq_len(s)] / media_2
    est <- (est_1 + est_2) / 2
    est <- est / mean(est)
  } else {
    est_1 <- y[seq_len(s)] - media_1
    est_2 <- y[s + seq_len(s)] - media_2
    est <- (est_1 + est_2) / 2
    est <- est - mean(est)
  }

  list(nivel = media_2, tendencia = tendencia_2s,
       estacional = est)
}

# Holt-Winters aditivo. Los primeros 2*s datos inicializan el metodo;
# yhat empieza en 2*s+1 para no usar informacion futura en errores de ajuste.
ajustar_hw <- function(y, s, alpha, beta, gamma) {
  s <- .validar_hw(y, s, alpha, beta, gamma, FALSE)
  n <- length(y)
  inicio <- 2L * s
  inicial <- .inicializar_holt_winters(y, s, FALSE)
  nivel <- rep(NA_real_, n)
  tendencia <- rep(NA_real_, n)
  estacional <- rep(NA_real_, n)
  yhat <- rep(NA_real_, n)

  nivel[inicio] <- inicial$nivel
  tendencia[inicio] <- inicial$tendencia
  estacional[seq_len(inicio)] <- rep(inicial$estacional, 2L)

  if (n > inicio) {
    for (t in (inicio + 1L):n) {
      est_ant <- estacional[t - s]
      yhat[t] <- nivel[t - 1L] + tendencia[t - 1L] + est_ant
      nivel[t] <- alpha * (y[t] - est_ant) +
        (1 - alpha) * (nivel[t - 1L] + tendencia[t - 1L])
      tendencia[t] <- beta * (nivel[t] - nivel[t - 1L]) +
        (1 - beta) * tendencia[t - 1L]
      estacional[t] <- gamma * (y[t] - nivel[t]) +
        (1 - gamma) * est_ant
    }
  }

  nivel_final <- nivel[n]
  tendencia_final <- tendencia[n]
  pronosticar <- function(h) {
    h <- .validar_h(h)
    pasos <- seq_len(h)
    indices_estacionales <- n + ((pasos - 1L) %% s) + 1L - s
    nivel_final + tendencia_final * pasos + estacional[indices_estacionales]
  }

  .resultado_metodo(
    yhat, pronosticar,
    list(s = s, alpha = alpha, beta = beta, gamma = gamma,
         nivel = nivel, tendencia = tendencia, estacional = estacional,
         inicializacion = "dos ciclos; pronosticos evaluables desde 2*s+1")
  )
}

# Holt-Winters multiplicativo.
ajustar_hw_mult <- function(y, s, alpha, beta, gamma) {
  s <- .validar_hw(y, s, alpha, beta, gamma, TRUE)
  n <- length(y)
  inicio <- 2L * s
  inicial <- .inicializar_holt_winters(y, s, TRUE)
  nivel <- rep(NA_real_, n)
  tendencia <- rep(NA_real_, n)
  estacional <- rep(NA_real_, n)
  yhat <- rep(NA_real_, n)

  nivel[inicio] <- inicial$nivel
  tendencia[inicio] <- inicial$tendencia
  estacional[seq_len(inicio)] <- rep(inicial$estacional, 2L)

  if (n > inicio) {
    for (t in (inicio + 1L):n) {
      est_ant <- estacional[t - s]
      yhat[t] <- (nivel[t - 1L] + tendencia[t - 1L]) * est_ant
      nivel[t] <- alpha * (y[t] / est_ant) +
        (1 - alpha) * (nivel[t - 1L] + tendencia[t - 1L])
      tendencia[t] <- beta * (nivel[t] - nivel[t - 1L]) +
        (1 - beta) * tendencia[t - 1L]
      if (!is.finite(nivel[t]) || nivel[t] <= 0) {
        .stop_hw_mult_inestable()
      }
      estacional[t] <- gamma * (y[t] / nivel[t]) +
        (1 - gamma) * est_ant
      if (!is.finite(estacional[t]) || estacional[t] <= 0) {
        .stop_hw_mult_inestable()
      }
    }
  }

  nivel_final <- nivel[n]
  tendencia_final <- tendencia[n]
  pronosticar <- function(h) {
    h <- .validar_h(h)
    pasos <- seq_len(h)
    indices_estacionales <- n + ((pasos - 1L) %% s) + 1L - s
    (nivel_final + tendencia_final * pasos) *
      estacional[indices_estacionales]
  }

  .resultado_metodo(
    yhat, pronosticar,
    list(s = s, alpha = alpha, beta = beta, gamma = gamma,
         nivel = nivel, tendencia = tendencia, estacional = estacional,
         inicializacion = "dos ciclos; pronosticos evaluables desde 2*s+1")
  )
}

.rejilla_por_defecto <- function(metodo) {
  if (metodo == "mm") {
    return(data.frame(k = 2:12))
  }
  if (metodo == "dmm") {
    return(data.frame(k = 3:12))
  }
  if (metodo == "ses") {
    return(data.frame(alpha = seq(0.02, 0.98, by = 0.02)))
  }
  if (metodo == "holt") {
    valores <- seq(0.05, 0.95, by = 0.05)
    return(expand.grid(alpha = valores, beta = valores))
  }
  if (metodo %in% c("hw", "hw_mult")) {
    valores <- seq(0.1, 0.9, by = 0.1)
    return(expand.grid(alpha = valores, beta = valores, gamma = valores))
  }
  stop("metodo no tiene parametros optimizables en el catalogo.",
       call. = FALSE)
}

.crear_metodo_rejilla <- function(y, metodo, fila, s) {
  if (metodo == "mm") return(ajustar_mm(y, fila$k))
  if (metodo == "dmm") return(ajustar_dmm(y, fila$k))
  if (metodo == "ses") return(ajustar_ses(y, fila$alpha))
  if (metodo == "holt") return(ajustar_holt(y, fila$alpha, fila$beta))
  if (metodo == "hw") {
    return(ajustar_hw(y, s, fila$alpha, fila$beta, fila$gamma))
  }
  if (metodo == "hw_mult") {
    return(ajustar_hw_mult(y, s, fila$alpha, fila$beta, fila$gamma))
  }
  stop("metodo no implementado en optimizar().", call. = FALSE)
}
