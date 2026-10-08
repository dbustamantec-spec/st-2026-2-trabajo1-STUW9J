source("R/00-lectura.R")
datos1 <- leer_serie(x = "datos/STUW9J-serie-1.csv", unidad = "casos", fuente = "")
datos2 <- leer_serie(x = "datos/STUW9J-serie-2.csv", unidad = "afiliados", fuente = "")
datos3 <- leer_serie(x = "datos/STUW9J-serie-3.csv", unidad = "consultas", fuente = "")
