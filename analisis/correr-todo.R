t_inicio <- Sys.time()

if (!file.exists("datos/STMSYJ-serie-1.csv"))
  stop("Ejecute este script desde la raíz del repositorio st-2026-2-trabajo2-STMSYJ")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
dir.create("figs", showWarnings = FALSE)

# serie 1 (trimestral, T = 96): 8 trimestres = 2 años  -> tiene 12 bloques

# serie 2 (semanal, T = 200):  10 semanas              -> tiene 20 bloques

# serie 3 (mensual, T = 144):  12 meses = 1 año        -> tiene 12 bloques

k_bloque <- c(8, 10, 12)
k_sensibilidad <- list(c(6, 8, 12), c(8, 10, 20), c(8, 12, 18))

resultados <- vector("list", 3)

for (i in 1:3) {
  message("== Serie ", i)
  x <- leer_serie(i)
  meta <- attr(x, "meta")
  z <- x$valor
  r <- list(datos = x, meta = meta)
  
  r$ficha <- ficha_serie(x)
  
  r$sd_nivel <- sd_nivel(z, k_bloque[i])
  r$sd_nivel_sens <- do.call(rbind, lapply(k_sensibilidad[[i]], function(k) {
    s <- sd_nivel(z, k)
    data.frame(k = k, bloques = s$n_bloques, pendiente = s$pendiente, lambda = s$lambda)
  }))
  if (all(z > 0)) {
    r$perfil <- perfil_lambda(z)
    r$lambda <- elegir_lambda(r$perfil)
    r$verif_perfil <- verificar_perfil(z, c(0, 0.5, 1))
    r$sens_modelo <- do.call(rbind, list(
      perfil_lambda_arima(z, c(0, 1, 0), deriva = TRUE),
      perfil_lambda_arima(z, c(1, 1, 0), deriva = TRUE),
      perfil_lambda_arima(z, c(2, 1, 0), deriva = TRUE),
      perfil_lambda_arima(z, c(1, 2, 0), deriva = FALSE),
      perfil_lambda_arima(z, c(0, 2, 1), deriva = FALSE)
    ))
  } else {
    stop("La serie ", i, " tiene valores no positivos; ver sección 3b")
  }
  r$g <- boxcox_tr(z, r$lambda$lambda)   # g(Z_t) que pasa a la sección 4
  
  pref <- sprintf("figs/serie%d-", i)
  guardar_png(paste0(pref, "tiempo.png"), graficar_serie(x))
  guardar_png(paste0(pref, "sd-nivel.png"), graficar_sd_nivel(r$sd_nivel, meta), 7, 5)
  guardar_png(paste0(pref, "perfil-lambda.png"), graficar_perfil(r$perfil, r$lambda, meta), 7, 5)
  
  
  resultados[[i]] <- r
}

# Resumen de la sección 3 en consola
print(do.call(rbind, lapply(resultados, function(r) data.frame(
  serie = r$meta$serie, lambda_hat = round(r$perfil$lambda_hat, 3),
  ic_inf = round(r$perfil$ic[1], 3), ic_sup = round(r$perfil$ic[2], 3),
  lambda_elegido = r$lambda$lambda, pendiente_sd = round(r$sd_nivel$pendiente, 3)))))

writeLines(utils::capture.output(utils::sessionInfo()), "sesion-info.txt")
message("Tiempo de ejecución: ", round(difftime(Sys.time(), t_inicio, units = "secs"), 1), " s")
