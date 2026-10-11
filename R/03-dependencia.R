# pacf_dl
pacf_dl <- function(r) {
  m <- length(r)
  if (m == 0) return(numeric(0))
  
  # Matriz para guardar los coeficientes phi_{k, j}
  phi <- matrix(0, nrow = m, ncol = m)
  
  # Inicialización k = 1
  phi[1, 1] <- r[1]
  
  if (m > 1) {
    for (k in 2:m) {
      # Calcular phi_{k, k} (el valor de la PACF en el rezago k)
      suma_num <- sum(phi[k-1, 1:(k-1)] * r[(k-1):1])
      suma_den <- sum(phi[k-1, 1:(k-1)] * r[1:(k-1)])
      
      phi[k, k] <- (r[k] - suma_num) / (1 - suma_den)
      
      # Actualizar los demás coeficientes phi_{k, j} de forma vectorizada
      phi[k, 1:(k-1)] <- phi[k-1, 1:(k-1)] - phi[k, k] * phi[k-1, (k-1):1]
    }
  }
  
  # La PACF corresponde a la diagonal principal de la matriz
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

