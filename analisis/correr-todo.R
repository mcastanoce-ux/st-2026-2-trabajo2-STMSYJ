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

tipos_df <- list(c("trend", "drift", "none"), c("trend", "drift", "none"),
                 c("trend", "drift", "none"))
etiquetas_g <- c("Z", "ln Z", "Z")
unidades_g <- c("Mt", "ln(COP)", "USD/t")

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
  
  g <- r$g
  r$tabla_dif <- tabla_diferencias(g, dmax = 3)
  r$df <- estrategia_df(g, tipos_df[[i]], criterio = "BIC")
  r$df_tabla <- tabla_df(r$df)
  r$df_sens <- sensibilidad_df(g)
  r$df_verif <- do.call(rbind, lapply(r$df, function(a) verificar_adf(diferenciar(g, a$d), a$tipo, a$p)))
  r$d <- r$df[[length(r$df)]]$d
  r$W <- diferenciar(g, r$d)
  r$K_sw <- elegir_K(r$W)
  r$sw <- prueba_media_sw(r$W, r$K_sw)
  r$sw_sens <- sensibilidad_sw(r$W)
  r$constante <- r$sw$rechaza
  guardar_png(paste0(pref, "diferencias.png"),
              graficar_diferencias(x$fecha, g, meta, etiquetas_g[i], unidades_g[i]), 8, 7)
  guardar_png(paste0(pref, "acf-diferencias.png"),
              graficar_acf_diferencias(g, meta, etiquetas_g[i]), 10, 3.8)
  
  r$m <- elegir_m(length(r$W))
  r$fm <- funciones_muestrales(r$W, r$m)
  r$eacf <- eacf_propia(r$W, ar.max = 7, ma.max = 13)
  r$eacf_verif <- verificar_eacf(r$W, ar.max = 7, ma.max = 13)
  r$eacf_vertices <- vertices_eacf(r$eacf$simbolo)
  etiqueta_w <- paste0(c("", "\u2207", "\u2207\u00b2")[r$d + 1], etiquetas_g[i])
  r$etiqueta_w <- etiqueta_w
  guardar_png(paste0(pref, "acf-pacf-W.png"), graficar_acf_pacf(r$W, r$m, meta, etiqueta_w), 10, 4.2)
  
  
  resultados[[i]] <- r
}

# Resumen de la sección 3 en consola
print(do.call(rbind, lapply(resultados, function(r) data.frame(
  serie = r$meta$serie, lambda_hat = round(r$perfil$lambda_hat, 3),
  ic_inf = round(r$perfil$ic[1], 3), ic_sup = round(r$perfil$ic[2], 3),
  lambda_elegido = r$lambda$lambda, pendiente_sd = round(r$sd_nivel$pendiente, 3)))))

writeLines(utils::capture.output(utils::sessionInfo()), "sesion-info.txt")
message("Tiempo de ejecución: ", round(difftime(Sys.time(), t_inicio, units = "secs"), 1), " s")
