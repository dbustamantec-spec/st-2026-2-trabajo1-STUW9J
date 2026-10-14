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
      rejilla <- rejilla[, columnas]
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
  
  modelos <- lapply(seq_len(nrow(rejilla)), function(i) {
    .crear_metodo_rejilla(y, metodo, rejilla[i, , drop = FALSE], s)
  })
  mascara_comun <- Reduce(
    `&`,
    lapply(modelos, function(modelo) {
      is.finite(modelo$yhat) & is.finite(y)
    })
  )
  if (!any(mascara_comun)) {
    stop("No hay origenes comunes con pronostico valido para comparar la rejilla.",
         call. = FALSE)
  }
  mse <- vapply(modelos, function(modelo) {
    mean((y[mascara_comun] - modelo$yhat[mascara_comun])^2)
  }, numeric(1))
  if (any(!is.finite(mse))) {
    stop("La rejilla produjo un MSE no finito; revise los datos y parametros.",
         call. = FALSE)
  }
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
