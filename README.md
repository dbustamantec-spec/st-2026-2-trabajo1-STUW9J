https://github.com/dbustamantec-spec/st-2026-2-trabajo1-STUW9J.git — Código de equipo: STUW9J

# Trabajo 1 — Series de tiempo

## Integrantes y trabajo realizado

La mayoría de los puntos se hicieron de manera conjunta. En algunos puntos se dividió tareas  para avanzar, pero siempre se trabajó sobre el mismo punto y se compartió la revisión y las decisiones del análisis en cada momento. No se asignaron secciones exclusivas a una sola persona.

| Integrantes |
|---|---|
| Daniel Bustamante Cortés | 
| Andrés Felipe López Díaz | 
| Viviana Astrid León Oróstogue | 

## Contenido del repositorio

| Archivo o carpeta | Contenido |
|---|---|
| `datos/STUW9J-serie-1.csv` | Casos mensuales de enfermedad diarreica aguda. |
| `datos/STUW9J-serie-2.csv` | Afiliados semanales al programa de hipertensión. |
| `datos/STUW9J-serie-3.csv` | Consultas semanales por enfermedad respiratoria. |
| `R/00-lectura.R` | Lectura y validación de series temporales. |
| `R/01-graficos.R` | Gráficos de series, rezagos, estacionalidad y descomposición. |
| `R/02-operadores.R` | Diferenciación y polinomios de rezago. |
| `R/03-dependencia.R` | ACF, PACF, Ljung–Box y varianza de largo plazo. |
| `R/04-metodos.R` | Métodos de pronóstico implementados para el trabajo. |
| `R/05-evaluacion.R` | Optimización de parámetros, origen móvil y métricas de error. |
| `analisis/correr-todo.R` | Carga los módulos, procesa las series y genera las figuras de análisis. |
| `figs/` | Figuras generadas por `correr-todo.R` y por la nota para interesados. |
| `informe/informe.qmd` | Informe técnico del trabajo. |
| `informe/informe.pdf` | Informe técnico renderizado. |
| `informe/nota-interesados.qmd` | Nota de pronóstico dirigida a los interesados de las tres series. |
| `informe/nota-interesados.pdf` | Nota para interesados renderizada. |
| `sesion-info.txt` | Salida de `sessionInfo()` de R. |
| `st-2026-2-trabajo1-STUW9J.Rproj` | Proyecto de RStudio. |
| `exposicion/` | Carpeta para las diapositivas. |

## Cómo ejecutar

Ejecutar desde la raíz del repositorio, con R 4.3 o posterior y los paquetes `ggplot2`, `patchwork`, `tibble` y `dplyr` instalados:

```r
system.time(source("analisis/correr-todo.R"))
```

En una ejecución con R 4.5.3, el script tardó **11.70 segundos** de tiempo transcurrido (`elapsed`). El script genera las figuras de análisis en `figs/` y guarda `sesion-info.txt`. Para renderizar el informe:

```r
quarto render informe/informe.qmd
```

## Resumen de resultados

Los RMSE corresponden a las tablas de origen móvil del informe renderizado. Los referentes son ingenuo estacional para la serie 1 e ingenuo simple para las series 2 y 3. En las series 1 y 2, el método ganador es el mismo en \(h=1\) y en \(H\). Para la serie 3, el ganador entre los métodos comparados es suavizamiento exponencial simple en \(h=1\) y media simple en \(H=6\); el ingenuo se mantiene como referente y su RMSE a un paso es ligeramente menor que el de suavizamiento exponencial simple.

| Serie | Transformación elegida para el análisis | Método ganador y parámetros | RMSE del ganador / referente, \(h=1\) | RMSE del ganador / referente, \(H\) |
|---|---|---|---:|---:|
| 1 — Casos de EDA mensual | \((1-B^{12})\log(y_t)\) | Holt–Winters multiplicativo; \(s=12\), \(\alpha=0.10\), \(\beta=0.10\), \(\gamma=0.30\). | 12.1488 / 23.2386 | \(H=12\): 14.2661 / 21.7814 |
| 2 — Afiliados control de hipertensión semanal | \((1-B)y_t\) para el análisis de dependencia; los métodos de pronóstico se compararon en la escala original. | Tendencia cuadrática; sin parámetros de rejilla. Ajuste final en la serie original: intercepto \(\approx51.4133\), tendencia \(\approx0.7923\), término cuadrático \(\approx0.011954\). | 9.2492 / 13.4826 | \(H=6\): 9.7274 / 28.2725 |
| 3 — Consultas enfermedad respiratoria semanal | Sin transformación. | \(h=1\): suavizamiento exponencial simple, \(\alpha=0.98\). \(H=6\): media simple, sin parámetro ajustable. | 7.0520 / 7.0267 | \(H=6\): 9.0192 / 13.4173 |

## Declaración de uso de inteligencia artificial

Se utilizaron herramientas de inteligencia artificial para revisar, corregir y mejorar código, investigar errores de renderizado, mejorar la redacción y organizar documentación. Se revisaron e integraron los cambios del análisis, los resultados y el contenido final entregado. Se usaron Claude, Gemini y Github Copilot.
