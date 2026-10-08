# Serie
graficar_serie <- function(datos, titulo) {
  library(ggplot2)
  unidad <- attr(datos, "unidad")
  fuente <- attr(datos, "fuente")
  n_obs <- nrow(datos)
  
  g <- ggplot(datos, aes(x = fecha, y = y)) +
    geom_line(color = "black", linewidth = 0.5) +
    labs(title = titulo, x = "Fecha",y = sprintf("%s", unidad),
         caption = sprintf("Número de observaciones: %d", n_obs)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.caption = element_text(hjust = 0, color = "black", size = 9))
  
  return(g)
}


# Correlogramas
correlograma <- function(datos, m = NULL) {
  library(patchwork)
  library(ggplot2)
  y_vec <- if(is.data.frame(datos)) datos$y else as.numeric(datos)
  n <- length(y_vec)
  
  if(is.null(m)) {
    m <- min(24, floor(n/ 4))}
  
  media <- mean(y_vec)
  acfs <- numeric(m)
  denominador <- sum((y_vec - media)^2)
  
  ci <- 0.95
  climo <- qnorm((1 + ci)/2) / sqrt(n)
  banda_inf <- -climo
  banda_sup <- climo
  
  for (j in 1:m) {
    numerador <- 0
    for (i in (j+1):n) {
      numerador <- numerador + (y_vec[i] - media) * (y_vec[i-j] - media)
    }
    acfs[j] <- numerador / denominador
  }
  
  pacfs <- pacf(y_vec, lag.max = m, plot = FALSE)
  pacfs <- as.vector(pacfs$acf)
  df <- data.frame(Rezagos = 1:m, acfs = acfs, pacfs = pacfs)
  
  # Correlogramas
  g1 <- ggplot(df, aes(x = Rezagos, y = acfs)) +
    geom_segment(aes(xend = Rezagos, yend = 0), color = "black", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_hline(yintercept = c(banda_inf, banda_sup), color = "red", linetype = "dashed") + 
    labs(title ="Correlograma ACF", x = "Rezago", y = "Autocorrelación") +
    theme_minimal() + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  g2 <- ggplot(df, aes(x = Rezagos, y = pacfs)) +
    geom_segment(aes(xend = Rezagos, yend = 0), color = "black", linewidth = 1) + 
    geom_hline(yintercept = 0, color = "black") + 
    geom_hline(yintercept = c(banda_inf, banda_sup), color = "red", linetype = "dashed") + 
    labs(title ="Correlograma PACF", x = "Rezago", y = "Autocorrelación") +
    theme_minimal() + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  attr(acfs, "acfs") <- acfs
  return(g1 / g2)
  
}

# Rezagos
grafico_rezagos <- function(datos, k){
  library(ggplot2)
  yt <- datos$y[(k+1):(length(datos$y))]
  yt_k <- datos$y[1:(length(datos$y)-k)]
  df <- data.frame(yt = yt, yt_k = yt_k)
  
  ggplot(df, aes(x = yt_k, y = yt)) +
    geom_point(color = "black") +
    labs(title = paste("Gráfico de rezagos ( k =", k, ")"), 
         x = expression(y[t-k]),
         y = expression(y[t])) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}


# Estacionalidad 
media_movil_centrada <- function(y, s) {
  stopifnot(
    is.numeric(y),
    length(y) >= 2 * s,
    all(is.finite(y)),
    is.numeric(s),
    length(s) == 1L,
    is.finite(s),
    s >= 2,
    s == round(s)
  )
  
  n <- length(y)
  pesos <- if (s %% 2 == 1) {
    rep(1 / s, s)
  } else {
    c(0.5, rep(1, s - 1), 0.5) / s
  }
  
  h <- (length(pesos) - 1) / 2
  tendencia <- rep(NA_real_, n)
  
  for (i in seq.int(h + 1, n - h)) {
    tendencia[i] <- sum(pesos * y[(i - h):(i + h)])
  }
  
  tendencia
}


descomposicion_clasica <- function(datos, s,
                                   tipo = c("aditivo", "multiplicativo")) {
  tipo <- match.arg(tipo)
  
  stopifnot(
    is.data.frame(datos),
    all(c("fecha", "y") %in% names(datos)),
    inherits(datos$fecha, "Date"),
    is.numeric(datos$y),
    is.numeric(s),
    length(s) == 1L,
    is.finite(s),
    s %in% c(12, 52),
    nrow(datos) >= 2 * s,
    !anyNA(datos$fecha),
    all(is.finite(datos$y))
  )
  
  datos <- datos[order(datos$fecha), , drop = FALSE]
  
  if (any(diff(datos$fecha) <= 0)) {
    stop("Las fechas deben ser estrictamente crecientes.")
  }
  
  if (s == 52 && any(diff(as.numeric(datos$fecha)) != 7)) {
    stop("s = 52 exige observaciones semanales consecutivas.")
  }
  
  if (s == 12) {
    meses <- as.POSIXlt(datos$fecha)$year * 12 +
      as.POSIXlt(datos$fecha)$mon
    
    if (any(diff(meses) != 1)) {
      stop("s = 12 exige observaciones mensuales consecutivas.")
    }
  }
  
  if (tipo == "multiplicativo" && any(datos$y <= 0)) {
    stop("La descomposición multiplicativa exige valores positivos.")
  }
  
  n <- nrow(datos)
  posicion <- seq_len(n)
  indice <- ((posicion - 1) %% s) + 1
  ciclo <- ((posicion - 1) %/% s) + 1
  tendencia <- media_movil_centrada(datos$y, s)
  
  sin_tendencia <- if (tipo == "aditivo") {
    datos$y - tendencia
  } else {
    datos$y / tendencia
  }
  
  indices_crudos <- vapply(
    seq_len(s),
    function(j) mean(sin_tendencia[indice == j], na.rm = TRUE),
    numeric(1)
  )
  
  indices <- if (tipo == "aditivo") {
    indices_crudos - mean(indices_crudos)
  } else {
    indices_crudos / mean(indices_crudos)
  }
  
  estacional <- indices[indice]
  
  residuo <- if (tipo == "aditivo") {
    datos$y - tendencia - estacional
  } else {
    datos$y / (tendencia * estacional)
  }
  
  data.frame(
    fecha = datos$fecha,
    y = datos$y,
    indice = indice,
    ciclo = ciclo,
    tendencia = tendencia,
    indice_estacional = estacional,
    residuo = residuo
  )
}


grafico_estacional <- function(datos, s,
                               tipo = c("aditivo", "multiplicativo")) {
  library(ggplot2)
  library(patchwork)
  
  tipo <- match.arg(tipo)
  unidad <- attr(datos, "unidad")
  if (is.null(unidad)) unidad <- "valor"
  
  dec <- descomposicion_clasica(datos, s = s, tipo = tipo)
  
  etiquetas_indice <- if (s == 12) {
    meses <- c("ene", "feb", "mar", "abr", "may", "jun",
               "jul", "ago", "sep", "oct", "nov", "dic")
    mes_inicio <- as.POSIXlt(min(dec$fecha))$mon + 1
    meses[((mes_inicio - 1 + seq_len(s) - 1) %% 12) + 1]
  } else {
    as.character(seq_len(s))
  }
  
  g1 <- ggplot(
    dec,
    aes(x = indice, y = y, group = ciclo, colour = factor(ciclo))
  ) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    scale_x_continuous(breaks = seq_len(s), labels = etiquetas_indice) +
    labs(
      title = "Gráfico estacional",
      x = if (s == 12) "Mes del año" else "Semana del ciclo",
      y = unidad,
      colour = "Ciclo"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  g2 <- ggplot(dec, aes(x = fecha, y = tendencia)) +
    geom_line(color = "steelblue", linewidth = 0.7, na.rm = TRUE) +
    labs(title = "Tendencia estimada", x = "Fecha", y = unidad) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  indices <- unique(dec[, c("indice", "indice_estacional")])
  
  g3 <- ggplot(
    indices,
    aes(x = indice, y = indice_estacional, group = 1)
  ) +
    geom_line(color = "darkorange", linewidth = 0.8) +
    geom_point(color = "darkorange", size = 2) +
    scale_x_continuous(breaks = seq_len(s), labels = etiquetas_indice) +
    labs(
      title = "Índices estacionales normalizados",
      x = if (s == 12) "Mes del año" else "Semana del ciclo",
      y = if (tipo == "multiplicativo") "Índice (media = 1)" else "Índice (media = 0)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  g4 <- ggplot(dec, aes(x = fecha, y = residuo)) +
    geom_line(color = "black", linewidth = 0.7, na.rm = TRUE) +
    geom_hline(
      yintercept = if (tipo == "multiplicativo") 1 else 0,
      color = "red",
      linetype = "dashed"
    ) +
    labs(
      title = "Residuo",
      x = "Fecha",
      y = if (tipo == "multiplicativo") "y / (T * S)" else "y - T - S"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  g1 / (g2 | g3 | g4)
}
