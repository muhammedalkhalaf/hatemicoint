# hatemicoint 1.1.0

This release corrects the Phillips statistics and changes several defaults.
Results differ from version 1.0.1 for every data set.

* Zt statistic (R/utils.R): the statistic was multiplied by an extra factor
  sqrt(n - 1). It is now Zt = (rho* - 1) / sqrt(sigma2 / sum_{t=1}^{n-1} u_t^2),
  the form of Gregory and Hansen (1996, p. 105) and Phillips (1987); equation
  (6) of Hatemi-J (2008) as printed omits the square root in the denominator. With `kernel = "iid"`
  Zt now equals the Dickey-Fuller t-ratio up to the residual-variance divisor
  (n instead of n - 2); the old value was about ten times too large in
  absolute value for n = 100.
* Serial-correlation correction: the correction subtracted from the
  cross-product term was twice the one-sided sum sum_{j>=1} w(j/B) gamma(j)
  of equation (3). It is now the one-sided sum, lambda = (sigma2 - gamma(0))/2,
  multiplied by n to undo the 1/n in the autocovariances of equation (4), as
  in Phillips (1987). Sums use the normalisation of the paper: sum over
  t = 1, ..., n - 1 of u_t^2, autocovariances divided by n, and
  Za = n (rho* - 1).
* Quadratic spectral kernel: the weight used z = 6 pi j / B instead of
  z = 6 pi (j/B) / 5 (Andrews 1991), and the sum was truncated at the
  bandwidth. The QS sum now runs over all lags with the correct weight. The
  Bartlett weight changed from 1 - j/(B + 1) (Newey-West form, used in
  1.0.1) to 1 - j/B, the kernel k(x) = 1 - |x| evaluated at x = j/B, matching
  the paper's w(j/B) notation and the way the QS kernel is evaluated. This is
  the package's choice, not a correction from the paper; results with
  `kernel = "bartlett"` and a given `bwl` therefore differ from 1.0.1. With
  integer B only lags 1, ..., B - 1 enter, so `bwl = 1` gives no correction.
* New defaults for the Phillips statistics: `kernel = "qs"` with AR(1)
  prewhitening (new argument `prewhite = TRUE`) and, when `bwl = NULL`, the
  Andrews (1991) automatic bandwidth for the QS kernel, following footnote 3
  of the paper and Andrews and Monahan (1992). The AR(1) coefficient is capped
  at 0.97 in absolute value and the kernel estimate is recoloured by
  1/(1 - a1)^2. The bandwidth is chosen for every break pair; the values at
  the selected break pairs are returned as `bwl_zt` and `bwl_za`. Versions
  before 1.1.0 used `kernel = "iid"` (no correction) by default. For
  `kernel = "bartlett"` with `bwl = NULL` the old rule
  `round(4 (n/100)^(2/9))` is kept.
* ADF lag selection: candidate lags are now compared on a common sample
  (the first `maxlags` usable observations are dropped for every candidate)
  for the "tstat", "aic" and "sic" rules; the final statistic for the chosen
  lag is computed on the largest sample available for that lag. The lag is
  still re-selected for every break pair. New rule `lag_selection = "fixed"`
  uses `maxlags` lags. The default `maxlags` is now
  `floor(4 (n/100)^(1/4))` (Schwert rule; the paper does not give a lag rule)
  instead of 8. The selected lag at the ADF* break pair is returned as
  `lag_adf`.
* Break grid: the bounds were rounded separately, so for some n (for example
  24, 50, 57) pairs outside the trimming region were used. tb1 now runs from
  ceiling(0.15 n) to floor(0.70 n) and tb2 from tb1 + ceiling(0.15 n) to
  floor(0.85 n), so that tb1 >= 0.15 n, tb2 <= 0.85 n and tb2 - tb1 >= 0.15 n
  for every n (the test suite checks every n from 20 to 2000).
* A singular regressor matrix in the ADF regression was replaced by a huge
  diagonal (giving a statistic near zero) and the long-run variance was
  floored at 1e-10. Both now give NA for that break pair, the pair is skipped
  in the infimum and a warning reports the number of skipped pairs
  (`n_skipped` in the result). Missing values in `y` or `x` are an error.
* `print` and `summary` show the prewhitening, the bandwidth and `NA` break
  dates when no statistic could be computed.
* Added a testthat suite pinning every statistic to hand computations and
  checking the break grid bounds.

# hatemicoint 1.0.0

## Initial CRAN Release

* Implements the Hatemi-J (2008) cointegration test with two unknown regime shifts
* Three test statistics: ADF*, Zt*, and Za*
* Lag selection methods: t-statistic, AIC, SIC
* Kernel options: IID, Bartlett, Quadratic Spectral
* Critical values for k = 1 to 4 regressors
* Print and summary methods for test results

## References

* Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts
  with an application to financial market integration. *Empirical Economics*,
  35, 497-505. DOI: 10.1007/s00181-007-0175-9
