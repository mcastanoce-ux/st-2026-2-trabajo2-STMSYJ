# 02-integracion.R

diferenciar <- function(g, d) if (d == 0) g else diff(g, differences = d)

# r_k muestral (la misma definición de stats::acf: denominador sum (x - xbar)^2).
acf_k <- function(x, k) {
  x <- x - mean(x); n <- length(x)
  sum(x[1:(n - k)] * x[(1 + k):n]) / sum(x^2)
}

tabla_diferencias <- function(g, dmax = 3) {
  do.call(rbind, lapply(0:dmax, function(d) {
    y <- diferenciar(g, d)
    data.frame(d = d, n = length(y), varianza = stats::var(y),
               r1 = acf_k(y, 1), r2 = acf_k(y, 2))
  }))
}

graficar_acf_diferencias <- function(g, meta, etiqueta_g, lag.max = 24) {
  op <- graphics::par(mfrow = c(1, 3), mar = c(4.5, 4, 3.5, 1), oma = c(2, 0, 0, 0))
  on.exit(graphics::par(op))
  nombres <- c(etiqueta_g, paste0("∇", etiqueta_g), paste0("∇²", etiqueta_g))
  for (d in 0:2) {
    y <- diferenciar(g, d)
    a <- stats::acf(y, lag.max = lag.max, plot = FALSE)
    graphics::plot(a, main = sprintf("ACF de %s", nombres[d + 1]), ylim = c(-1, 1),
                   xlab = "Rezago", ylab = expression(r[k]), ci.col = "#c0392b")
  }
  graphics::mtext(sprintf("Serie %d. Bandas ±1,96/√n (ruido blanco). Rezagos en %s.",
                          meta$serie, sub("al$", "es", meta$frecuencia)),
                  side = 1, outer = TRUE, line = 0.5, cex = 0.75, adj = 0)
}

# Gráfico en el tiempo de g, ∇g y ∇²g
graficar_diferencias <- function(fecha, g, meta, etiqueta_g, unidad_g) {
  op <- graphics::par(mfrow = c(3, 1), mar = c(2.5, 4.5, 2.2, 1), oma = c(2.5, 0, 0, 0))
  on.exit(graphics::par(op))
  nombres <- c(etiqueta_g, paste0("\u2207", etiqueta_g), paste0("\u2207\u00b2", etiqueta_g))
  for (d in 0:2) {
    y <- diferenciar(g, d)
    f <- fecha[(d + 1):length(fecha)]
    graphics::plot(f, y, type = "l", col = "#1f4e79", xlab = "", ylab = unidad_g,
                   main = sprintf("Serie %d: %s", meta$serie, nombres[d + 1]), cex.main = 0.95)
    if (d > 0) graphics::abline(h = c(0, mean(y)), lty = c(1, 2), col = c("grey50", "#c0392b"))
  }
  graphics::mtext("Fecha. En las diferencias: línea continua en 0 y punteada en la media muestral.",
                  side = 1, outer = TRUE, line = 1, cex = 0.75, adj = 0)
}

# prueba de  Dickey y Fuller 
mco <- function(X, y) {
  q <- qr(X)
  b <- qr.coef(q, y)
  e <- y - X %*% b
  gl <- nrow(X) - ncol(X)
  s2 <- sum(e^2) / gl
  V <- s2 * chol2inv(qr.R(q))
  list(coef = drop(b), ee = sqrt(diag(V)), s2 = s2, sce = sum(e^2), n = nrow(X), gl = gl)
}

matriz_adf <- function(y, p, tipo, desde) {
  dy <- diff(y)                        # dy[s] = y[s+1] - y[s]
  idx <- (desde - 1):(length(dy))      # posiciones en dy de ∇y_t
  X <- cbind(y_lag = y[idx])           # y_{t-1} = y[s] si ∇y_t = dy[s]
  if (p > 0) for (j in 1:p) X <- cbind(X, dy[idx - j])
  if (p > 0) colnames(X)[-1] <- paste0("dy_lag", 1:p)
  if (tipo %in% c("drift", "trend")) X <- cbind(constante = 1, X)
  if (tipo == "trend") X <- cbind(X, tendencia = idx)
  list(X = X, y = dy[idx])
}

adf_propia <- function(y, tipo = c("none", "drift", "trend"), pmax = NULL,
                       criterio = c("BIC", "AIC")) {
  tipo <- match.arg(tipo); criterio <- match.arg(criterio)
  n <- length(y)
  if (is.null(pmax)) pmax <- floor(12 * (n / 100)^0.25)   # regla de Schwert (1989)
  desde_comun <- pmax + 2
  crit <- vapply(0:pmax, function(p) {
    m <- matriz_adf(y, p, tipo, desde_comun)
    f <- mco(m$X, m$y)
    N <- f$n; k <- ncol(m$X)
    pen <- if (criterio == "BIC") log(N) else 2
    N * log(f$sce / N) + pen * k
  }, numeric(1))
  p <- (0:pmax)[which.min(crit)]
  m <- matriz_adf(y, p, tipo, p + 2)
  f <- mco(m$X, m$y)
  j <- which(colnames(m$X) == "y_lag")
  list(tipo = tipo, p = p, pmax = pmax, criterio = criterio,
       criterios = data.frame(p = 0:pmax, valor = crit),
       gamma = unname(f$coef[j]), ee_gamma = unname(f$ee[j]),
       tau = unname(f$coef[j] / f$ee[j]), n_ef = f$n,
       coef = f$coef, ee = f$ee)
}

# Valores críticos de la tabla de Dickey y Fuller
criticos_df <- function(y, tipo, p) {
  u <- urca::ur.df(y, type = tipo, lags = p, selectlags = "Fixed")
  cv <- u@cval[1, ]
  list(cv = cv, tau_urca = unname(u@teststat[1]), modelo_urca = u)
}

valor_p_df <- function(tau, n_ef, tipo) {
  tr <- switch(tipo, none = "nc", drift = "c", trend = "ct")
  urca::punitroot(tau, N = n_ef, trend = tr, statistic = "t")
}

prueba_df <- function(y, tipo, pmax = NULL, criterio = "BIC") {
  a <- adf_propia(y, tipo, pmax, criterio)
  k <- criticos_df(y, tipo, a$p)
  a$cv <- k$cv
  a$cv5 <- unname(k$cv["5pct"])
  a$valor_p <- valor_p_df(a$tau, a$n_ef, tipo)
  a$rechaza <- a$tau < a$cv5
  a$tau_urca <- k$tau_urca
  a
}


verificar_adf <- function(y, tipo, p) {
  m <- matriz_adf(y, p, tipo, p + 2)   # misma muestra que ur.df(lags = p): t = p + 2..n
  f <- mco(m$X, m$y)
  j <- which(colnames(m$X) == "y_lag")
  tau_propio <- unname(f$coef[j] / f$ee[j])
  u <- urca::ur.df(y, type = tipo, lags = p, selectlags = "Fixed")
  cu <- u@testreg$coefficients
  fila <- if (tipo == "none") 1 else 2
  data.frame(tipo = tipo, rezagos = p,
             gamma_propio = unname(f$coef[j]), gamma_urca = cu[fila, 1],
             tau_propio = tau_propio, tau_urca = unname(u@teststat[1]),
             diferencia_tau = tau_propio - unname(u@teststat[1]))
}

estrategia_df <- function(g, tipos, pmax = NULL, criterio = "BIC") {
  res <- list()
  for (d in 0:2) {
    y <- diferenciar(g, d)
    r <- prueba_df(y, tipos[d + 1], pmax, criterio)
    r$d <- d
    res[[d + 1]] <- r
    if (r$rechaza) break
  }
  res
}

tabla_df <- function(res) {
  do.call(rbind, lapply(res, function(r) data.frame(
    d = r$d, especificacion = r$tipo, rezagos = r$p, n_ef = r$n_ef,
    gamma = r$gamma, ee = r$ee_gamma, tau = r$tau,
    cv1 = unname(r$cv["1pct"]), cv5 = r$cv5, cv10 = unname(r$cv["10pct"]),
    valor_p = r$valor_p, decision = ifelse(r$rechaza, "rechaza H0", "no rechaza H0"))))
}

sensibilidad_df <- function(g, dmax = 2) {
  out <- list()
  for (d in 0:dmax) for (tp in c("none", "drift", "trend")) for (cr in c("BIC", "AIC")) {
    r <- prueba_df(diferenciar(g, d), tp, criterio = cr)
    out[[length(out) + 1]] <- data.frame(d = d, especificacion = tp, criterio = cr,
                                         rezagos = r$p, tau = r$tau, cv5 = r$cv5,
                                         decision = ifelse(r$rechaza, "rechaza", "no rechaza"))
  }
  do.call(rbind, out)
}

texto_df <- function(r, nombre_y) {
  esp <- switch(r$tipo,
                none  = list(reg = "\\nabla y_t = \\gamma y_{t-1}", caso = "la especificación sin constante (caso 1, $\\tau$)"),
                drift = list(reg = "\\nabla y_t = \\beta_0 + \\gamma y_{t-1}", caso = "la especificación con constante (caso 2, $\\tau_\\mu$)"),
                trend = list(reg = "\\nabla y_t = \\beta_0 + \\beta_1 t + \\gamma y_{t-1}", caso = "la especificación con constante y tendencia (caso 4, $\\tau_\\tau$)"))
  rez <- if (r$p > 0) sprintf(" + \\sum_{j=1}^{%d}\\delta_j \\nabla y_{t-j}", r$p) else ""
  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",")
  fp_txt <- function(p) if (p < 1e-4) sprintf("$%s \\times 10^{%d}$",
                                              formatC(p / 10^floor(log10(p)), format = "f", digits = 1, decimal.mark = ","), floor(log10(p))) else fmt(p, 4)
  lect <- if (r$rechaza)
    sprintf("%s no tiene raíz unitaria: es estacionaria%s, y no necesita otra diferencia.", nombre_y,
            switch(r$tipo, none = " alrededor de cero", drift = " alrededor de una media", trend = " alrededor de una recta"))
  else sprintf("No hay evidencia contra la raíz unitaria en %s: se diferencia y se prueba de nuevo.", nombre_y)
  paste0(
    "Prueba sobre ", nombre_y, ", con ", esp$caso, ", con $p = ", r$p, "$ rezagos (", r$criterio,
    " sobre $p = 0,\\dots,", r$pmax, "$) y $n_{ef} = ", r$n_ef, "$:\n\n",
    "1. **Hipótesis.** $H_0: \\gamma = 0$ (raíz unitaria) frente a $H_1: \\gamma < 0$, en $", esp$reg, rez, " + e_t$.\n",
    "2. **Estadístico.** $\\hat\\tau = \\hat\\gamma / \\widehat{\\text{ee}}(\\hat\\gamma)$; bajo $H_0$ no sigue una $t$ sino la distribución de Dickey y Fuller para ", esp$caso,
    ". Valores críticos de Fuller (1976, tabla 8.5.2; Hamilton, 1994, tabla B.6) para el tamaño de muestra: 1 % = ", fmt(r$cv["1pct"], 2),
    ", 5 % = ", fmt(r$cv["5pct"], 2), ", 10 % = ", fmt(r$cv["10pct"], 2), ".\n",
    "3. **Región de rechazo al 5 %.** $\\hat\\tau < ", fmt(r$cv5, 2), "$.\n",
    "4. **Valor observado.** $\\hat\\gamma = ", fmt(r$gamma, 4), "$, ee $= ", fmt(r$ee_gamma, 4), "$, $\\hat\\tau = ", fmt(r$tau), "$; valor p aproximado de MacKinnon (1996) = ", fp_txt(r$valor_p), ".\n",
    "5. **Decisión.** ", if (r$rechaza) "Se rechaza $H_0$" else "No se rechaza $H_0$", " al 5 %.\n",
    "6. **Lectura.** ", lect, "\n")
}


#Prueba de la constante theta_0 con la regla S_W

prueba_media_sw <- function(w, K) {
  n <- length(w)
  c0 <- sum((w - mean(w))^2) / n
  r <- if (K > 0) vapply(1:K, function(k) acf_k(w, k), numeric(1)) else numeric(0)
  factor <- 1 + 2 * sum(r)
  sw <- sqrt(c0 / n * factor)
  tobs <- mean(w) / sw
  cv <- stats::qt(0.975, df = n - 1)
  list(n = n, media = mean(w), c0 = c0, K = K, r = r, factor = factor,
       sw = sw, ee_ingenuo = sqrt(c0 / n), t = tobs, cv = cv,
       valor_p = 2 * stats::pt(-abs(tobs), df = n - 1), rechaza = abs(tobs) > cv)
}

elegir_K <- function(w, kmax = 6) {
  n <- length(w)
  r <- vapply(1:kmax, function(k) acf_k(w, k), numeric(1))
  banda <- 1.96 * sqrt((1 + 2 * cumsum(c(0, r[-kmax]^2))) / n)
  fuera <- which(abs(r) > banda)
  if (length(fuera) == 0) 0L else max(fuera)
}

sensibilidad_sw <- function(w, Ks = 0:6) {
  do.call(rbind, lapply(Ks, function(K) {
    s <- prueba_media_sw(w, K)
    data.frame(K = K, factor = s$factor, S_W = s$sw, t = s$t, rechaza = s$rechaza)
  }))
}

texto_sw <- function(s, nombre_w, unidad_theta) {
  fmt <- function(x, d = 4) formatC(x, format = "f", digits = d, decimal.mark = ",")
  rs <- if (s$K > 0) paste0(" + 2(", sub("^\\+ ", "", paste0(ifelse(s$r < 0, "- ", "+ "), fmt(abs(s$r), 3), collapse = " ")), ")") else ""
  paste0(
    "1. **Hipótesis.** $H_0: \\mu_W = 0$ (es decir, $\\theta_0 = 0$) frente a $H_1: \\mu_W \\neq 0$, con $W_t = ", nombre_w, "$.\n",
    "2. **Estadístico.** $t = \\bar W / S_{\\bar W}$, con $S_{\\bar W} = \\sqrt{\\hat\\gamma_0/n\\,(1 + 2r_1 + \\dots + 2r_K)}$, $K = ", s$K,
    "$; bajo $H_0$, aproximadamente $t_{n-1} = t_{", s$n - 1, "}$ (≈ normal estándar).\n",
    "3. **Región de rechazo al 5 %.** $|t| > t_{0,975;", s$n - 1, "} = ", fmt(s$cv, 3), "$.\n",
    "4. **Valor observado.** $\\bar W = ", fmt(s$media, 5), "$, $\\hat\\gamma_0 = ", fmt(s$c0, 6), "$, factor $= 1", rs, " = ", fmt(s$factor, 3),
    "$, $S_{\\bar W} = ", fmt(s$sw, 5), "$, $t = ", fmt(s$t, 3), "$, valor p $= ", fmt(s$valor_p, 4), "$.\n",
    "5. **Decisión.** ", if (s$rechaza) "Se rechaza $H_0$" else "No se rechaza $H_0$", " al 5 %.\n",
    "6. **Lectura.** ", if (s$rechaza)
      paste0("La media de $W_t$ es distinta de cero: el modelo lleva constante, y $\\bar W$ estima ", unidad_theta, ".")
    else paste0("No hay evidencia de que la media de $W_t$, ", unidad_theta, ", sea distinta de cero: el modelo no lleva constante."), "\n")
}