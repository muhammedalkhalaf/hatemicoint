# Hand computations of the Phillips statistics, Hatemi-J (2008) eq. (3) to (6)

ns <- asNamespace("hatemicoint")

set.seed(1)
n <- 100
x <- cumsum(rnorm(n))
y <- 1 + 2 * x + cumsum(rnorm(n))
D1 <- as.numeric(seq_len(n) > 30)
D2 <- as.numeric(seq_len(n) > 60)
u <- lm.fit(cbind(1, D1, D2, x, D1 * x, D2 * x), y)$residuals

S <- sum(u[-n]^2)                       # sum_{t=1}^{n-1} u_t^2
cr <- sum(u[-n] * u[-1])                # sum_{t=1}^{n-1} u_t u_{t+1}
rho <- cr / S
v <- u[-1] - rho * u[-n]                # residuals of u_t on u_{t-1}
N <- length(v)
g <- function(e, j) sum(e[(j + 1):length(e)] * e[1:(length(e) - j)]) / n
bart <- function(z) ifelse(z < 1, 1 - z, 0)
qs <- function(z) { w <- 6 * pi * z / 5; 3 / w^2 * (sin(w) / w - cos(w)) }
hand <- function(e, w, jmax, recolour = 1) {
  s2 <- recolour * (g(e, 0) + 2 * sum(sapply(seq_len(jmax), function(j) w(j) * g(e, j))))
  lam <- (s2 - g(v, 0)) / 2
  rs <- (cr - n * lam) / S
  list(za = n * (rs - 1), zt = (rs - 1) / sqrt(s2 / S), lrvar = s2)
}

test_that("iid Zt and Za match eq. (5)-(6) and the DF t-ratio", {
  pk <- ns$.compute_pp(u, NULL, "iid", FALSE)
  expect_equal(pk$za, n * (rho - 1), tolerance = 1e-10)
  expect_equal(pk$zt, (rho - 1) / sqrt(g(v, 0) / S), tolerance = 1e-10)
  # DF t-ratio uses the (n - 2) divisor for the residual variance, eq. (4) uses n
  df_t <- coef(summary(lm(diff(u) ~ u[-n] - 1)))[1, 3]
  expect_equal(pk$zt * sqrt((n - 2) / n), df_t, tolerance = 1e-10)
  expect_equal(pk$zt, -3.3077536731, tolerance = 1e-8)
  expect_true(is.na(pk$bw))
})

test_that("Bartlett kernel without prewhitening matches the hand computation", {
  h <- hand(v, function(j) bart(j / 4), 4)
  pk <- ns$.compute_pp(u, 4, "bartlett", FALSE)
  expect_equal(pk$za, h$za, tolerance = 1e-10)
  expect_equal(pk$zt, h$zt, tolerance = 1e-10)
  expect_equal(pk$lrvar, h$lrvar, tolerance = 1e-10)
  expect_equal(pk$bw, 4)
})

test_that("QS kernel (untruncated, Andrews weights) matches the hand computation", {
  h <- hand(v, function(j) qs(j / 4), N - 1)
  pk <- ns$.compute_pp(u, 4, "qs", FALSE)
  expect_equal(pk$za, h$za, tolerance = 1e-10)
  expect_equal(pk$zt, h$zt, tolerance = 1e-10)
  # weight at j = B has z = 6 pi / 5 (the 1/5 factor of Andrews 1991)
  z <- 6 * pi / 5
  expect_equal(ns$.kernel_weight(1, "qs"), 3 / z^2 * (sin(z) / z - cos(z)))
  expect_equal(ns$.kernel_weight(0, "qs"), 1)
  expect_equal(ns$.kernel_weight(c(0.25, 1, 2), "bartlett"), c(0.75, 0, 0))
})

test_that("prewhitened QS with the Andrews bandwidth matches the hand computation", {
  a1 <- sum(v[-1] * v[-N]) / sum(v[-N]^2)
  a1 <- max(min(a1, 0.97), -0.97)
  e <- v[-1] - a1 * v[-N]
  Te <- length(e)
  phi <- sum(e[-1] * e[-Te]) / sum(e[-Te]^2)
  phi <- max(min(phi, 0.97), -0.97)
  alpha2 <- 4 * phi^2 / (1 - phi)^4
  B <- 1.3221 * (alpha2 * Te)^(1 / 5)
  h <- hand(e, function(j) qs(j / B), Te - 1, 1 / (1 - a1)^2)
  pk <- ns$.compute_pp(u, NULL, "qs", TRUE)
  expect_equal(pk$bw, B, tolerance = 1e-10)
  expect_equal(pk$za, h$za, tolerance = 1e-10)
  expect_equal(pk$zt, h$zt, tolerance = 1e-10)
  # Bartlett with prewhitening and a user bandwidth
  h2 <- hand(e, function(j) bart(j / 4), 4, 1 / (1 - a1)^2)
  pk2 <- ns$.compute_pp(u, 4, "bartlett", TRUE)
  expect_equal(pk2$za, h2$za, tolerance = 1e-10)
  expect_equal(pk2$zt, h2$zt, tolerance = 1e-10)
})

test_that("the one-sided correction is used (not twice the one-sided sum)", {
  pk <- ns$.compute_pp(u, 4, "bartlett", FALSE)
  lam <- sum(sapply(1:4, function(j) bart(j / 4) * g(v, j)))
  expect_equal(pk$rho_star, (cr - n * lam) / S, tolerance = 1e-10)
})

test_that("degenerate residuals give NA, never a substitute value", {
  pk <- ns$.compute_pp(rep(0, 50), NULL, "qs", TRUE)
  expect_true(is.na(pk$zt) && is.na(pk$za))
  expect_true(is.character(pk$reason))
})

test_that("a zero Andrews bandwidth gives finite statistics (no NaN weights)", {
  # residual series with an AR(1) coefficient of exactly zero after prewhitening
  set.seed(5)
  e <- rnorm(60)
  expect_equal(ns$.lrvar_kernel(e, 0, "qs", 60), sum(e^2) / 60)
  expect_equal(ns$.andrews_bandwidth_qs(c(1, 0, -1, 0, 1, 0, -1, 0)), 0)
  pk <- ns$.compute_pp(cumsum(rnorm(80)), 0, "qs", FALSE)
  expect_true(is.finite(pk$zt) && is.finite(pk$za))
})
