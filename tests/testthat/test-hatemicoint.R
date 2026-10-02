ns <- asNamespace("hatemicoint")

set.seed(3)
n <- 60
x <- cumsum(rnorm(n))
y <- 1 + 2 * x + rnorm(n)

test_that("the infimum equals the minimum over the grid of the per-pair statistics", {
  r <- suppressWarnings(hatemicoint(y, x, maxlags = 1, lag_selection = "fixed",
                                    kernel = "bartlett", bwl = 2, prewhite = FALSE))
  g <- ns$.break_grid(n, 0.15)
  stats <- t(apply(g, 1, function(tb) {
    d1 <- as.numeric(seq_len(n) > tb[1]); d2 <- as.numeric(seq_len(n) > tb[2])
    u <- lm.fit(cbind(1, d1, d2, x, d1 * x, d2 * x), y)$residuals
    pp <- ns$.compute_pp(u, 2, "bartlett", FALSE)
    c(ns$.compute_adf(u, 1, "fixed")$stat, pp$zt, pp$za)
  }))
  expect_equal(r$adf_min, min(stats[, 1]), tolerance = 1e-10)
  expect_equal(r$zt_min, min(stats[, 2]), tolerance = 1e-10)
  expect_equal(r$za_min, min(stats[, 3]), tolerance = 1e-10)
  i <- which.min(stats[, 2])
  expect_equal(c(r$tb1_zt, r$tb2_zt), unname(g[i, ]))
  expect_equal(r$lag_adf, 1)
  expect_equal(r$bwl_zt, 2)
  expect_equal(r$n_skipped, 0)
})

test_that("defaults and arguments are recorded", {
  r <- suppressWarnings(hatemicoint(y, x))
  expect_equal(r$maxlags, floor(4 * (n / 100)^(1 / 4)))
  expect_equal(r$kernel, "qs")
  expect_true(r$prewhite)
  expect_true(is.na(r$bwl))
  expect_true(is.finite(r$bwl_zt))
  expect_equal(r$cv_adfzt, c(-6.503, -6.015, -5.653))
  expect_equal(r$cv_za, c(-90.794, -76.003, -52.232))
  expect_output(print(r), "Bandwidth: automatic")
  expect_output(s <- summary(r), "Prewhitening")
  expect_type(s, "list")
})

test_that("input checks", {
  expect_error(hatemicoint(y, x[-1]), "must equal")
  expect_error(hatemicoint(y, cbind(x, x, x, x, x)), "Maximum of 4")
  expect_error(hatemicoint(y, x, trimming = 0.6), "trimming")
  expect_error(hatemicoint(y, x, bwl = -1), "bwl")
  expect_error(hatemicoint(y, x, maxlags = -1), "maxlags")
  expect_warning(hatemicoint(y[1:40], x[1:40]), "small")
})

test_that("a constant regressor is skipped with a warning and NA statistics", {
  expect_warning(r <- hatemicoint(y, rep(1, n), kernel = "iid"), "skipped")
  expect_true(is.na(r$adf_min) && is.na(r$zt_min) && is.na(r$za_min))
  expect_true(r$n_skipped > 0)
  expect_output(print(r), "NA")
})
