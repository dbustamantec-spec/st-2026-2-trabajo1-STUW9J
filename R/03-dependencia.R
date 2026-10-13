# pacf_dl
pacf_dl <- function(r) {
  m <- length(r)
  if (m == 0) return(numeric(0))
  
  phi <- matrix(0, nrow = m, ncol = m)
  phi[1, 1] <- r[1]
  
  if (m > 1) {
    for (k in 2:m) {
      suma_num <- sum(phi[k-1, 1:(k-1)] * r[(k-1):1])
      suma_den <- sum(phi[k-1, 1:(k-1)] * r[1:(k-1)])
      
      phi[k, k] <- (r[k] - suma_num) / (1 - suma_den)
      phi[k, 1:(k-1)] <- phi[k-1, 1:(k-1)] - phi[k, k] * phi[k-1, (k-1):1]
    }
  }
  
  return(diag(phi))
}


# acf_muestral
acf_muestral <- function(datos, m = NULL){
  y_vec <- if(is.data.frame(datos)) datos$y else as.numeric(datos)
  n <- length(y_vec)
  
  if(is.null(m)) {
    m <- min(24, floor(n/ 4))}
  
  media <- mean(y_vec)
  acfs <- numeric(m)
  denominador <- sum((y_vec - media)^2)
  
  for (j in 1:m){
    numerador <- 0
    for (i in (j+1):n) {
      numerador <- numerador + (y_vec[i] - media) * (y_vec[i-j] - media)
    }
    acfs[j] <- numerador / denominador}
  return(acfs)
}

# ljung-box
ljung_box <- function(r, T, m, p = 0) {
  T <- as.numeric(T)
  m <- as.numeric(m)
  p <- as.numeric(p)
  
  h <- 1:m
  r_m <- r[1:m]
  
  Q_m <- T * (T + 2) * sum((r_m^2) / (T - h))
  df <- m - p
  
  val_critico <- qchisq(0.95, df)
  p_val <- pchisq(Q_m, df, lower.tail = FALSE)
  
  return(list(
    estadistico = Q_m,
    grados_libertad = df,
    valor_critico_5 = val_critico,
    p_valor = p_val
  ))
}

# Varianza a largo plazo
varianza_largo_plazo <- function(y, L){
  if (anyNA(y)) stop("El vector 'y' contiene valores faltantes (NA).")
  
  n <- length(y)
  media <- mean(y)
  r <- numeric(L)
  
  denominador <- sum((y - media)^2)
  for (k in 1:L) {
    numerador <- sum((y[(k + 1):n] - media) * (y[1:(n - k)] - media))
    r[k] <- numerador / denominador
  }
  
  gamma0 <- denominador / n
  var_lp <- gamma0 * (1 + 2 * sum((1 - (1:L)/(L + 1)) * r))
  
  return(var_lp)
}