boxcox_bc <- function(z, lambda) {
  if (any(z <= 0)) stop("Box y Cox exige valores positivos")
  if (abs(lambda) < 1e-10) log(z) else (z^lambda - 1) / lambda
}

boxcox_tr <- function(z, lambda) {
  if (any(z <= 0)) stop("Box y Cox exige valores positivos")
  if (abs(lambda - 1) < 1e-10) return(z)
  boxcox_bc(z, lambda)
}

#inversa de box cox
boxcox_inv <- function(g, lambda) {
  if (abs(lambda - 1) < 1e-10) return(g)
  if (abs(lambda) < 1e-10) return(exp(g))
  (lambda * g + 1)^(1 / lambda)
}

# Ajusta log(sd) = a + b log(media) por mínimos cuadrados; lambda = 1 - b.
sd_nivel <- function(z, k) {
  nb <- floor(length(z) / k)
  if (nb < 4) stop("Muy pocos bloques: use un k menor")
  bloque <- rep(seq_len(nb), each = k)
  zz <- z[seq_len(nb * k)]
  media <- as.numeric(tapply(zz, bloque, mean))
  desv  <- as.numeric(tapply(zz, bloque, function(u) stats::sd(diff(u))))
  ajuste <- stats::lm(log(desv) ~ log(media))
  b  <- unname(stats::coef(ajuste)[2])
  eb <- unname(summary(ajuste)$coefficients[2, 2])
  list(bloques = data.frame(bloque = seq_len(nb), media = media, sd_dif = desv),
       k = k, n_bloques = nb, descartadas = length(z) - nb * k,
       pendiente = b, ee_pendiente = eb, lambda = 1 - b,
       ajuste = ajuste)
}


# Modelo de trabajo: W_t = (1 - B) g_lambda(Z_t), t = 2..T, sigue un AR(1)

loglik_ar1_exacta <- function(w) {
  n <- length(w)
  w1 <- w[1]; wa <- w[-1]; wb <- w[-n]
  mu_gls <- function(phi) {
    ((1 - phi^2) * w1 + (1 - phi) * sum(wa - phi * wb)) /
      ((1 - phi^2) + (n - 1) * (1 - phi)^2)
  }
  lc <- function(phi) {
    mu <- mu_gls(phi)
    S <- (1 - phi^2) * (w1 - mu)^2 + sum(((wa - mu) - phi * (wb - mu))^2)
    -n / 2 * (log(2 * pi) + log(S / n) + 1) + 0.5 * log(1 - phi^2)
  }
  op <- stats::optimize(lc, interval = c(-0.999, 0.999), maximum = TRUE,
                        tol = 1e-10)
  phi <- op$maximum
  mu <- mu_gls(phi)
  S <- (1 - phi^2) * (w1 - mu)^2 + sum(((wa - mu) - phi * (wb - mu))^2)
  list(loglik = op$objective, phi = phi, mu = mu, sigma2 = S / n, n = n)
}

loglik_lambda <- function(z, lambda) {
  g <- boxcox_bc(z, lambda)
  loglik_ar1_exacta(diff(g))$loglik + (lambda - 1) * sum(log(z[-1]))
}


perfil_lambda <- function(z, rejilla = seq(-1, 2.5, by = 0.01), nivel = 0.95) {
  if (any(z <= 0))
    stop("La serie tiene valores no positivos: Box y Cox no se aplica tal cual")
  ll <- vapply(rejilla, function(l) loglik_lambda(z, l), numeric(1))
  j <- which.max(ll)
  if (j == 1 || j == length(rejilla)) warning("El máximo está en el borde de la rejilla")
  lo <- rejilla[max(1, j - 2)]; hi <- rejilla[min(length(rejilla), j + 2)]
  op <- stats::optimize(function(l) loglik_lambda(z, l), c(lo, hi),
                        maximum = TRUE, tol = 1e-8)
  lhat <- op$maximum; lmax <- op$objective
  corte <- lmax - stats::qchisq(nivel, df = 1) / 2
  f <- function(l) loglik_lambda(z, l) - corte
  inf <- if (f(min(rejilla)) < 0) stats::uniroot(f, c(min(rejilla), lhat))$root else NA
  sup <- if (f(max(rejilla)) < 0) stats::uniroot(f, c(lhat, max(rejilla)))$root else NA
  list(perfil = data.frame(lambda = rejilla, loglik = ll),
       lambda_hat = lhat, loglik_max = lmax, corte = corte,
       ic = c(inf = inf, sup = sup), nivel = nivel,
       ar1 = loglik_ar1_exacta(diff(boxcox_bc(z, lhat))))
}

elegir_lambda <- function(pl) {
  ic <- pl$ic
  dentro <- function(v) !is.na(ic[1]) && !is.na(ic[2]) && v >= ic[1] && v <= ic[2]
  if (dentro(1)) return(list(lambda = 1, razon = "1 está en el intervalo: no se transforma"))
  if (dentro(0)) return(list(lambda = 0, razon = "0 está en el intervalo y 1 no: logaritmo"))
  cand <- round(pl$lambda_hat * 2) / 2
  if (dentro(cand)) return(list(lambda = cand, razon = "múltiplo de 0.5 más cercano dentro del intervalo"))
  list(lambda = round(pl$lambda_hat, 2), razon = "ningún valor redondo está en el intervalo")
}


verificar_perfil <- function(z, lambdas = c(0, 0.5, 1)) {
  out <- lapply(lambdas, function(l) {
    g <- boxcox_bc(z, l)
    fit <- stats::arima(diff(g), order = c(1, 0, 0), include.mean = TRUE,
                        method = "ML",
                        optim.control = list(maxit = 1000, reltol = 1e-12))
    jac <- (l - 1) * sum(log(z[-1]))
    propio <- loglik_lambda(z, l)
    data.frame(lambda = l, propio = propio, arima_mas_jacobiano = fit$loglik + jac,
               diferencia = propio - (fit$loglik + jac))
  })
  do.call(rbind, out)
}

perfil_lambda_arima <- function(z, orden, deriva = TRUE,
                                rejilla = seq(-1, 2.5, by = 0.02), nivel = 0.95) {
  d <- orden[2]
  ll <- vapply(rejilla, function(l) {
    g <- boxcox_bc(z, l)
    fit <- try(suppressWarnings(stats::arima(diff(g, differences = d),
                                             order = c(orden[1], 0, orden[3]),
                                             include.mean = deriva, method = "ML")),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    fit$loglik + (l - 1) * sum(log(z[-seq_len(d)]))
  }, numeric(1))
  j <- which.max(ll)
  dentro <- rejilla[!is.na(ll) & ll >= max(ll, na.rm = TRUE) - stats::qchisq(nivel, 1) / 2]
  data.frame(modelo = sprintf("ARIMA(%d,%d,%d)%s", orden[1], orden[2], orden[3],
                              if (deriva) " con constante" else ""),
             lambda_hat = rejilla[j], ic_inf = min(dentro), ic_sup = max(dentro))
}


guardar_png <- function(archivo, expr, ancho = 8, alto = 4.5) {
  grDevices::png(archivo, width = ancho, height = alto, units = "in", res = 150)
  on.exit(grDevices::dev.off())
  force(expr)
  invisible(archivo)
}

graficar_serie <- function(x) {
  meta <- attr(x, "meta")
  op <- graphics::par(mar = c(6, 4.5, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(x$fecha, x$valor, type = "l", lwd = 1.4, col = "#1f4e79",
                 xaxt = "n", xlab = "Fecha",
                 ylab = sprintf("%s (%s)", meta$corto, meta$unidad_eje),
                 main = sprintf("Serie %d: %s", meta$serie, meta$nombre),
                 cex.main = 0.95)
  marcas <- pretty(x$fecha, n = 10)
  graphics::axis(1, at = marcas, labels = format(marcas, "%Y-%m"), cex.axis = 0.8)
  graphics::grid(nx = NA, ny = NULL, col = "grey85")
  graphics::mtext(sprintf("T = %d observaciones %ses, de %s a %s. Unidad: %s. Archivo %s.",
                          nrow(x), meta$frecuencia,
                          format(x$fecha[1]), format(x$fecha[nrow(x)]),
                          meta$unidad, meta$archivo),
                  side = 1, line = 4.5, cex = 0.75, adj = 0)
}

graficar_sd_nivel <- function(sn, meta) {
  b <- sn$bloques
  op <- graphics::par(mar = c(6, 4.5, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(b$media, b$sd_dif, log = "xy", pch = 19, col = "#1f4e79",
                 xlab = sprintf("Nivel medio del bloque (%s, escala log)", meta$unidad_eje),
                 ylab = sprintf("DE de la 1.ª diferencia (%s, escala log)", meta$unidad_eje),
                 main = sprintf("Serie %d: dispersión contra nivel, bloques de %d",
                                meta$serie, sn$k), cex.main = 0.95)
  xs <- exp(seq(log(min(b$media)), log(max(b$media)), length.out = 50))
  a <- stats::coef(sn$ajuste)
  graphics::lines(xs, exp(a[1] + a[2] * log(xs)), col = "#c0392b", lwd = 2)
  graphics::legend("topleft", bty = "n", cex = 0.85,
                   legend = sprintf("pendiente = %.3f (ee %.3f)  ->  lambda = 1 - pendiente = %.2f",
                                    sn$pendiente, sn$ee_pendiente, sn$lambda))
  graphics::mtext(sprintf("%d bloques consecutivos de %d observaciones; %d observaciones finales descartadas solo para este gráfico.",
                          sn$n_bloques, sn$k, sn$descartadas),
                  side = 1, line = 4.5, cex = 0.75, adj = 0)
}

graficar_perfil <- function(pl, eleccion, meta) {
  p <- pl$perfil
  vis <- p$loglik >= pl$loglik_max - 30
  op <- graphics::par(mar = c(6, 4.5, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(p$lambda[vis], p$loglik[vis], type = "l", lwd = 2, col = "#1f4e79",
                 xlab = expression(lambda), ylab = "Log-verosimilitud perfilada (con jacobiano)",
                 main = sprintf("Serie %d: perfil de lambda de Box y Cox", meta$serie),
                 cex.main = 0.95)
  graphics::abline(h = pl$corte, lty = 2, col = "grey40")
  graphics::abline(v = pl$ic, lty = 3, col = "grey40")
  graphics::abline(v = pl$lambda_hat, col = "#c0392b")
  graphics::abline(v = eleccion$lambda, col = "#27ae60", lwd = 2, lty = 4)
  graphics::legend("bottom", bty = "n", cex = 0.8, ncol = 2,
                   legend = c(sprintf("lambda_hat = %.3f", pl$lambda_hat),
                              sprintf("IC 95%%: [%.3f, %.3f]", pl$ic[1], pl$ic[2]),
                              sprintf("elegido = %g", eleccion$lambda),
                              "corte: l_max - qchisq(0.95,1)/2"),
                   col = c("#c0392b", "grey40", "#27ae60", "grey40"),
                   lty = c(1, 3, 4, 2), lwd = c(1, 1, 2, 1))
  graphics::mtext("Modelo de trabajo: AR(1) con media para (1-B) g(Z_t), MV exacta; jacobiano (lambda-1) sum ln Z_t, t = 2..T.",
                  side = 1, line = 4.5, cex = 0.75, adj = 0)
} 
