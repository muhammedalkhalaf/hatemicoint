# Hand computations of the ADF statistic and lag selection on a common sample

ns <- asNamespace("hatemicoint")

set.seed(1)
n <- 100
x <- cumsum(rnorm(n))
y <- 1 + 2 * x + cumsum(rnorm(n))
D1 <- as.numeric(seq_len(n) > 30)
D2 <- as.numeric(seq_len(n) > 60)
u <- lm.fit(cbind(1, D1, D2, x, D1 * x, D2 * x), y)$residuals
du <- diff(u)
ul <- u[-n]
adf_hand <- function(p, start) {
  idx <- start:(n - 1)
  X <- ul[idx]
  if (p > 0) for (i in 1:p) X <- cbind(X, du[idx - i])
  coef(summary(lm(du[idx] ~ X - 1)))[1, 3]
}
ic_hand <- function(p, maxlags, pen) {
  idx <- (maxlags + 1):(n - 1)
  X <- ul[idx]
  if (p > 0) for (i in 1:p) X <- cbind(X, du[idx - i])
  f <- lm(du[idx] ~ X - 1)
  T <- length(idx)
  log(sum(resid(f)^2) / T) + pen(T) * (p + 1) / T
}

test_that("fixed lag ADF equals the OLS t-ratio", {
  expect_equal(ns$.compute_adf(u, 0, "fixed")$stat, adf_hand(0, 1), tolerance = 1e-10)
  expect_equal(ns$.compute_adf(u, 3, "fixed")$stat, adf_hand(3, 4), tolerance = 1e-10)
  expect_equal(ns$.compute_adf(u, 3, "fixed")$lag, 3)
})

test_that("AIC and SIC choose the lag on the common sample", {
  ml <- 4
  p_sic <- which.min(sapply(0:ml, ic_hand, maxlags = ml, pen = log)) - 1
  p_aic <- which.min(sapply(0:ml, ic_hand, maxlags = ml, pen = function(T) 2)) - 1
  r_sic <- ns$.compute_adf(u, ml, "sic")
  r_aic <- ns$.compute_adf(u, ml, "aic")
  expect_equal(r_sic$lag, p_sic)
  expect_equal(r_aic$lag, p_aic)
  expect_equal(r_sic$stat, adf_hand(p_sic, p_sic + 1), tolerance = 1e-10)
  expect_equal(r_aic$stat, adf_hand(p_aic, p_aic + 1), tolerance = 1e-10)
})

test_that("tstat rule is general-to-specific at |t| > 1.645 on the common sample", {
  # Build residuals with strong dependence in differences so that a lag is kept
  set.seed(7)
  e <- as.numeric(stats::filter(rnorm(150), c(0.6, -0.3), method = "recursive"))
  z <- cumsum(e)
  ml <- 4
  nz <- length(z)
  dz <- diff(z); zl <- z[-nz]
  tsel <- 0
  for (p in ml:1) {
    idx <- (ml + 1):(nz - 1)
    X <- zl[idx]
    for (i in 1:p) X <- cbind(X, dz[idx - i])
    tt <- coef(summary(lm(dz[idx] ~ X - 1)))[p + 1, 3]
    if (abs(tt) > 1.645) { tsel <- p; break }
  }
  expect_true(tsel > 0)
  r <- ns$.compute_adf(z, ml, "tstat")
  expect_equal(r$lag, tsel)
  idx <- (tsel + 1):(nz - 1)
  X <- zl[idx]
  for (i in 1:tsel) X <- cbind(X, dz[idx - i])
  expect_equal(r$stat, coef(summary(lm(dz[idx] ~ X - 1)))[1, 3], tolerance = 1e-10)
})

test_that("singular ADF regression returns NA with a reason", {
  r <- ns$.compute_adf(rep(1, 40), 2, "fixed")
  expect_true(is.na(r$stat))
  expect_true(is.character(r$reason))
})
