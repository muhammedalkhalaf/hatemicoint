ns <- asNamespace("hatemicoint")

test_that("break grid respects the trimming bounds for every n = 20 to 2000", {
  bad <- integer(0)
  for (n in 20:2000) {
    g <- ns$.break_grid(n, 0.15)
    ok <- nrow(g) > 0 &&
      all(g[, "tb1"] >= 0.15 * n) &&
      all(g[, "tb1"] <= 0.70 * n) &&
      all(g[, "tb2"] <= 0.85 * n) &&
      all(g[, "tb2"] - g[, "tb1"] >= 0.15 * n)
    if (!ok) bad <- c(bad, n)
  }
  expect_length(bad, 0)
  g <- ns$.break_grid(100, 0.15)
  expect_equal(range(g[, "tb1"]), c(15, 70))
  expect_equal(range(g[, "tb2"]), c(30, 85))
  expect_equal(nrow(g), sum(sapply(15:70, function(t1) max(0, 85 - (t1 + 15) + 1))))
})

test_that("the cases n = 24, 50, 57 no longer leave the bounds", {
  for (n in c(24, 50, 57)) {
    g <- ns$.break_grid(n, 0.15)
    expect_true(all(g[, "tb2"] <= 0.85 * n))
    expect_true(all(g[, "tb2"] - g[, "tb1"] >= 0.15 * n))
  }
})
