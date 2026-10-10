# diferenciar
diferenciar <- function(y, d, D, s){
  for (i in 1:D){
    n <- length(y)
    if (n <= s) stop()
    y <- y[(s+1):n] - y[1:(n-s)]
  }
  
  for (j in 1:d){
    n <- length(y)
    if (n <= 1) stop()
    y <- y[2:n] - y[1:(n-1)]
  }
  
  return(y)
}

# polinomio
polinomio_rezago <- function(coef){
  coefs <- c(1,coef)
  raices <- polyroot(coefs)
  modulos <- Mod(raices)
  df <- data.frame(Raiz = raices, Modulo = modulos)
  return(df)
}
