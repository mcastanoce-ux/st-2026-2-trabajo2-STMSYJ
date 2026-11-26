# 03-identificacion.R

# Prueba de  Box y Jenkins
elegir_m <- function(n, tope = 24) min(floor(n / 4), tope)

funciones_muestrales <- function(w, m) {
  n <- length(w)
  r <- stats::acf(w, lag.max = m, plot = FALSE)$acf[-1]
  phi <- stats::pacf(w, lag.max = m, plot = FALSE)$acf[, 1, 1]
  data.frame(k = 1:m, r = r, phi_kk = phi,
             banda_rb = 1.96 / sqrt(n),
             banda_bartlett = 1.96 * sqrt((1 + 2 * cumsum(c(0, r[-m]^2))) / n),
             fuera_acf_rb = abs(r) > 1.96 / sqrt(n),
             fuera_acf_bartlett = abs(r) > 1.96 * sqrt((1 + 2 * cumsum(c(0, r[-m]^2))) / n),
             fuera_pacf = abs(phi) > 1.96 / sqrt(n))
}

graficar_acf_pacf <- function(w, m, meta, etiqueta_w) {
  fm <- funciones_muestrales(w, m)
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4, 3.5, 1), oma = c(2.2, 0, 0, 0))
  on.exit(graphics::par(op))
  yl <- range(c(-0.5, 0.5, fm$r, fm$phi_kk, fm$banda_bartlett, -fm$banda_bartlett))
  graphics::plot(fm$k, fm$r, type = "h", lwd = 2, ylim = yl, xlab = "Rezago k",
                 ylab = expression(r[k]), main = sprintf("ACF de %s", etiqueta_w))
  graphics::abline(h = 0)
  graphics::abline(h = c(-1, 1) * fm$banda_rb[1], lty = 2, col = "#c0392b")
  graphics::lines(fm$k, fm$banda_bartlett, lty = 3, col = "#1f4e79", lwd = 1.5)
  graphics::lines(fm$k, -fm$banda_bartlett, lty = 3, col = "#1f4e79", lwd = 1.5)
  graphics::legend("topright", bty = "n", cex = 0.75, lty = c(2, 3), col = c("#c0392b", "#1f4e79"),
                   legend = c("±1,96/√n (ruido blanco)", "Bartlett (corte de MA)"))
  graphics::plot(fm$k, fm$phi_kk, type = "h", lwd = 2, ylim = yl, xlab = "Rezago k",
                 ylab = expression(hat(phi)[kk]), main = sprintf("PACF de %s", etiqueta_w))
  graphics::abline(h = 0)
  graphics::abline(h = c(-1, 1) * fm$banda_rb[1], lty = 2, col = "#c0392b")
  graphics::legend("topright", bty = "n", cex = 0.75, lty = 2, col = "#c0392b",
                   legend = "±1,96/√n (corte de AR)")
  graphics::mtext(sprintf("Serie %d: n = %d, m = %d rezagos.", meta$serie, length(w), m),
                  side = 1, outer = TRUE, line = 0.8, cex = 0.75, adj = 0)
}


ar_mco_sin_media <- function(z, k) {
  n <- length(z)
  X <- sapply(1:k, function(i) z[(k + 1 - i):(n - i)])
  if (k == 1) X <- matrix(X, ncol = 1)
  qr.coef(qr(X), z[(k + 1):n])
}

eacf_propia <- function(w, ar.max = 7, ma.max = 13) {
  z <- w - mean(w)
  n <- length(z)
  kmax <- ar.max + ma.max + 1           # órdenes AR necesarios para la recursión
  # iteración 0: lista con el vector de coeficientes de cada AR(k), k = 1..kmax
  fi <- lapply(1:kmax, function(k) ar_mco_sin_media(z, k))
  val <- matrix(NA_real_, ar.max + 1, ma.max + 1,
                dimnames = list(AR = 0:ar.max, MA = 0:ma.max))
  for (q in 0:ma.max) {
    j <- q + 1
    # una iteración de la recursión: AR(k) en j a partir de AR(k) y AR(k+1) en j-1
    fi <- lapply(seq_len(length(fi) - 1), function(k) {
      a <- fi[[k]]; b <- fi[[k + 1]]
      prev <- c(-1, a)                  # phi_{k,0} = -1, phi_{k,1..k}
      b[1:k] - prev[1:k] * b[k + 1] / a[k]
    })
    val[1, q + 1] <- acf_k(z, j)        # k = 0: ACF muestral de z
    for (k in 1:ar.max) {
      phi <- fi[[k]]
      t <- (k + 1):n
      res <- z[t] - vapply(t, function(s) sum(phi * z[s - (1:k)]), numeric(1))
      val[k + 1, q + 1] <- acf_k(res, j)
    }
  }
  umbral <- outer(0:ar.max, 0:ma.max, function(k, q) 2 / sqrt(n - k - q - 1))
  simbolo <- ifelse(abs(val) > umbral, "x", "o")
  dimnames(simbolo) <- dimnames(val)
  list(eacf = val, simbolo = simbolo, n = n)
}

verificar_eacf <- function(w, ar.max = 7, ma.max = 13) {
  propia <- eacf_propia(w, ar.max, ma.max)
  tsa <- NULL
  utils::capture.output(tsa <- TSA::eacf(w, ar.max = ar.max, ma.max = ma.max))
  list(max_dif_valores = max(abs(propia$eacf - tsa$eacf)),
       simbolos_iguales = all(propia$simbolo == tsa$symbol),
       celdas_distintas = sum(propia$simbolo != tsa$symbol))
}

vertices_eacf <- function(simbolo, max_suma = 5, prof = 2, ancho = 3) {
  K <- nrow(simbolo) - 1; Q <- ncol(simbolo) - 1
  out <- NULL
  for (k in 0:K) for (q in 0:Q) {
    if (k + q > max_suma) next
    ok <- TRUE
    for (i in 0:min(prof, K - k)) for (j in i:min(i + ancho, Q - q)) {
      if (simbolo[k + i + 1, q + j + 1] != "o") ok <- FALSE
    }
    if (ok) out <- rbind(out, data.frame(p = k, q = q, suma = k + q))
  }
  if (!is.null(out)) out <- out[order(out$suma, out$p), ]
  out
}