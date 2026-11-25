# 00-lectura.R
# Lectura de las series y ficha verificada con código (sección 2a).

metadatos_series <- function() {
  data.frame( 
    serie      = 1:3,
    archivo    = sprintf("STMSYJ-serie-%d.csv", 1:3),
    nombre     = c("Material estéril acumulado en la escombrera norte",
                   "Precio de cierre semanal de la acción de la empresa",
                   "Costo directo de extracción por tonelada de mineral"),
    corto      = c("Estéril acumulado", "Precio de la acción", "Costo de extracción"),
    unidad     = c("millones de toneladas", "pesos por acción", "dólares por tonelada"),
    unidad_eje = c("Mt", "COP/acción", "USD/t"),
    frecuencia = c("trimestral", "semanal", "mensual"),
    paso       = c("3 months", "1 week", "1 month"),
    inicio     = c("2002-07-01", "2022-09-05", "2014-07-01"),
    fin        = c("2026-04-01", "2026-06-29", "2026-06-01"),
    T          = c(96L, 200L, 144L),
    stringsAsFactors = FALSE
  )
}

leer_serie <- function(i, dir_datos = "datos") {
  meta <- metadatos_series()[i, ]
  ruta <- file.path(dir_datos, meta$archivo)
  if (!file.exists(ruta)) stop("No existe el archivo ", ruta)
  x <- utils::read.csv(ruta, stringsAsFactors = FALSE)
  if (!identical(names(x), c("fecha", "valor")))
    stop("El archivo ", ruta, " no tiene las columnas fecha y valor")
  x$fecha <- as.Date(x$fecha)
  x$valor <- as.numeric(x$valor)
  attr(x, "meta") <- meta
  x
}


ficha_serie <- function(x) {
  meta <- attr(x, "meta")
  z <- x$valor
  f <- x$fecha
  esperadas <- seq(as.Date(meta$inicio), by = meta$paso, length.out = meta$T)
  chequeos <- c(
    T_coincide        = length(z) == meta$T,
    inicio_coincide   = f[1] == as.Date(meta$inicio),
    fin_coincide      = f[length(f)] == as.Date(meta$fin),
    sin_faltantes     = !anyNA(z) && !anyNA(f),
    fechas_regulares  = length(f) == length(esperadas) && all(f == esperadas),
    valores_positivos = all(z > 0)
  )
  if (!all(chequeos[1:5]))
    stop("La ficha de la serie ", meta$serie, " no coincide con lo declarado: ",
         paste(names(chequeos)[!chequeos[1:5]], collapse = ", "))
  list(
    tabla = data.frame(
      Serie = meta$serie, Nombre = meta$nombre, Unidad = meta$unidad,
      Frecuencia = meta$frecuencia, Inicio = format(f[1]),
      Fin = format(f[length(f)]), T = length(z),
      Minimo = min(z), Maximo = max(z), Media = mean(z),
      stringsAsFactors = FALSE
    ),
    chequeos = chequeos
  )
}
