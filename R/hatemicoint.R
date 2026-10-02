#' Hatemi-J Cointegration Test with Two Unknown Regime Shifts
#'
#' @description
#' Performs the Hatemi-J (2008) cointegration test which allows for two unknown
#' structural breaks (regime shifts) in the cointegrating relationship. The test
#' searches over all admissible break date pairs and returns the minimum
#' test statistics along with the endogenously determined break dates.
#'
#' @details
#' The Hatemi-J (2008) test extends the Gregory and Hansen (1996) cointegration
#' test by allowing for two structural breaks instead of one. The test is based
#' on the residuals from the cointegrating regression with regime shift dummies:
#'
#' \deqn{y_t = \alpha_0 + \alpha_1 D_{1t} + \alpha_2 D_{2t} + \beta_0' x_t +
#'       \beta_1' D_{1t} x_t + \beta_2' D_{2t} x_t + u_t}
#'
#' where \eqn{D_{1t} = 1} if \eqn{t > tb_1} and \eqn{D_{2t} = 1} if
#' \eqn{t > tb_2}.
#'
#' Three test statistics are computed for every break pair and the smallest
#' value over the grid is reported (equations 7 to 9 of the paper):
#' \itemize{
#'   \item \strong{ADF*}: t-ratio of \eqn{\delta} in
#'     \eqn{\Delta \hat u_t = \delta \hat u_{t-1} + \sum_{i=1}^p \phi_i
#'     \Delta \hat u_{t-i} + e_t} (no intercept).
#'   \item \strong{Zt*}: Phillips \eqn{Z_t} statistic, equation (6) (see below
#'     for the square root).
#'   \item \strong{Za*}: Phillips \eqn{Z_\alpha} statistic, equation (5).
#' }
#'
#' \strong{Phillips statistics (equations 3 to 6).} With \eqn{\hat\rho} the OLS
#' coefficient of \eqn{\hat u_{t-1}} on \eqn{\hat u_t} (no intercept),
#' \eqn{v_t = \hat u_t - \hat\rho \hat u_{t-1}}, autocovariances
#' \eqn{\hat\gamma(j) = n^{-1} \sum_t v_t v_{t-j}} (equation 4), and
#' \eqn{S = \sum_{t=1}^{n-1} \hat u_t^2}, the long-run variance is
#' \eqn{\hat\sigma^2 = \hat\gamma(0) + 2 \sum_{j \ge 1} w(j/B) \hat\gamma(j)},
#' the one-sided correction is \eqn{\hat\lambda = (\hat\sigma^2 - \hat\gamma(0))/2},
#' \eqn{\hat\rho^* = (\sum_{t=1}^{n-1} \hat u_t \hat u_{t+1} - n \hat\lambda)/S},
#' \eqn{Z_\alpha = n(\hat\rho^* - 1)} and
#' \eqn{Z_t = (\hat\rho^* - 1)/\sqrt{\hat\sigma^2 / S}}. The factor \eqn{n} in
#' front of \eqn{\hat\lambda} compensates the \eqn{1/n} in equation (4), as in
#' Phillips (1987); equation (3) of the paper prints the difference without it.
#' The \eqn{Z_t} formula is the one of Gregory and Hansen (1996, p. 105) and
#' Phillips (1987); equation (6) of the paper as printed omits the square root
#' in the denominator (that form would grow with \eqn{n}).
#' All sums use the normalisation of the paper (\eqn{n}, not \eqn{n - 1}).
#'
#' \strong{Long-run variance.} The paper (footnote 3) estimates
#' \eqn{\hat\sigma^2} with a prewhitened quadratic spectral (QS) kernel with
#' first-order autoregressive prewhitening and the automatic bandwidth of
#' Andrews (1991), following Andrews and Monahan (1992). This is the default
#' since version 1.1.0 (\code{kernel = "qs"}, \code{prewhite = TRUE},
#' \code{bwl = NULL}) and is implemented as follows for every break pair:
#' an AR(1) without intercept is fitted to \eqn{v_t}, its coefficient
#' \eqn{a_1} is capped at 0.97 in absolute value, the AR(1) residuals
#' \eqn{e_t} are used to compute the Andrews bandwidth
#' \eqn{B = 1.3221 (\hat\alpha(2) T)^{1/5}} with
#' \eqn{\hat\alpha(2) = 4 \hat\phi^2 / (1 - \hat\phi)^4}, where \eqn{\hat\phi}
#' is the AR(1) coefficient of \eqn{e_t}, also capped at 0.97 in absolute
#' value (a package choice, not part of Andrews 1991), and \eqn{T} is the
#' length of \eqn{e_t}; if the bandwidth is zero (\eqn{\hat\phi = 0}) the
#' kernel estimate reduces to \eqn{\hat\gamma(0)} of \eqn{e_t}. The QS
#' kernel \eqn{k(x) = 3/z^2 (\sin z / z - \cos z)}, \eqn{z = 6\pi x/5}, is
#' applied to all lags \eqn{j = 1, \ldots, T - 1} of \eqn{e_t} (the QS kernel
#' is not truncated), and the result is recoloured by dividing by
#' \eqn{(1 - a_1)^2}. With \code{kernel = "bartlett"} the weights are
#' \eqn{1 - j/B} for \eqn{j < B} and zero beyond, that is the kernel
#' \eqn{k(x) = 1 - |x|} evaluated at \eqn{x = j/B} as in the paper's
#' \eqn{w(j/B)} notation; versions before 1.1.0 used the Newey-West form
#' \eqn{1 - j/(B + 1)}. This is the package's choice, not a correction from
#' the paper. With integer \eqn{B} only lags \eqn{1, \ldots, B - 1} enter, so
#' \code{bwl = 1} gives no correction. If \code{bwl} is \code{NULL} the
#' Bartlett bandwidth is \code{round(4 (n/100)^(2/9))} (the rule used by
#' earlier versions). With \code{kernel = "iid"} no serial-correlation
#' correction is applied and \eqn{Z_t} equals the Dickey-Fuller t-ratio up to
#' the residual-variance divisor. Versions before 1.1.0 used
#' \code{kernel = "iid"} by default.
#'
#' \strong{ADF lag selection.} The paper does not state a lag rule. The
#' package compares the candidate lags \eqn{0, \ldots, maxlags} on a common
#' sample (the first \code{maxlags} usable observations are dropped for every
#' candidate) and computes the final statistic for the chosen lag on the
#' largest sample available for that lag. \code{"tstat"} is the
#' general-to-specific rule (largest lag whose last coefficient has
#' \eqn{|t| > 1.645}), \code{"aic"} and \code{"sic"} minimise
#' \eqn{\log(RSS/T) + c (p + 1)/T} with \eqn{c = 2} and \eqn{c = \log T}, and
#' \code{"fixed"} uses \code{maxlags} lags. The lag is re-selected for every
#' break pair, as in Gregory and Hansen (1996) style infimum tests. The
#' default \code{maxlags = floor(4 (n/100)^(1/4))} is the Schwert rule and is
#' the package's choice, not the paper's.
#'
#' \strong{Break grid.} With trimming \eqn{\tau} (default 0.15), \eqn{tb_1}
#' runs from \eqn{\lceil \tau n \rceil} to \eqn{\lfloor (1 - 2\tau) n \rfloor}
#' and \eqn{tb_2} from \eqn{tb_1 + \lceil \tau n \rceil} to
#' \eqn{\lfloor (1 - \tau) n \rfloor}, so that \eqn{tb_1 \ge \tau n},
#' \eqn{tb_2 \le (1 - \tau) n} and \eqn{tb_2 - tb_1 \ge \tau n} for every
#' \eqn{n}. Break pairs for which a statistic cannot be computed (singular
#' regressor matrix or non-positive long-run variance) are skipped with a
#' warning and never replaced by a substitute value.
#'
#' The null hypothesis is no cointegration. Rejection occurs when the test
#' statistic is smaller (more negative) than the critical value. The critical
#' values of Table 1 are asymptotic (response-surface intercepts). In a
#' Monte Carlo with the data generating process of Section 4 of the paper
#' (n = 100, m = 1, 400 replications, package defaults) the empirical 5
#' percent null quantiles were about -6.6 for ADF*, -6.5 for Zt* and -62 for
#' Za*, against the asymptotic values -6.015, -6.015 and -76.003, so at
#' n = 100 ADF* and Zt* reject more often than the nominal level and Za*
#' much less often. With the asymptotic Table 1 critical values the package
#' therefore does not reproduce Table 2 of the paper at n = 100: a 200
#' replication run with the defaults gave rejection frequencies of 0.175
#' (ADF*), 0.140 (Zt*) and 0.000 (Za*) under the null against 0.032, 0.069
#' and 0.069 in Table 2, and the ADF* size depends strongly on the lag rule
#' (0.005 with a fixed lag of 4, 0.105 with lag 0, 0.175 with the default
#' "tstat" rule). The package does not adjust for this.
#'
#' @param y Numeric vector. The dependent variable (must be I(1)).
#' @param x Numeric matrix or vector. The independent variable(s) (must be I(1)).
#'   Maximum of 4 regressors allowed (k <= 4).
#' @param maxlags Integer. Maximum number of lags for the ADF regression (the
#'   lag itself when \code{lag_selection = "fixed"}). Default \code{NULL}
#'   means \code{floor(4 * (n/100)^(1/4))}.
#' @param lag_selection Character. Lag selection rule: \code{"tstat"}
#'   (default), \code{"aic"}, \code{"sic"}, or \code{"fixed"}.
#' @param kernel Character. Kernel for the long-run variance in the Phillips
#'   statistics: \code{"qs"} (quadratic spectral, default), \code{"bartlett"},
#'   or \code{"iid"} (no correction).
#' @param bwl Numeric. Bandwidth for the kernel. If \code{NULL} (default), the
#'   Andrews (1991) automatic bandwidth is used for \code{"qs"} and
#'   \code{round(4 * (n/100)^(2/9))} for \code{"bartlett"}. Ignored for
#'   \code{"iid"}.
#' @param prewhite Logical. AR(1) prewhitening of the kernel estimator
#'   (Andrews and Monahan 1992). Default \code{TRUE}. Ignored for
#'   \code{"iid"}.
#' @param trimming Numeric. Trimming parameter for the break point search. Must
#'   be between 0 and 0.5 (exclusive). Default is 0.15.
#'
#' @return An object of class \code{"hatemicoint"} containing:
#' \describe{
#'   \item{adf_min}{Minimum ADF* test statistic}
#'   \item{tb1_adf}{First break location (observation number) for ADF*}
#'   \item{tb2_adf}{Second break location (observation number) for ADF*}
#'   \item{lag_adf}{ADF lag used at the ADF* break pair}
#'   \item{zt_min}{Minimum Zt* test statistic}
#'   \item{tb1_zt}{First break location for Zt*}
#'   \item{tb2_zt}{Second break location for Zt*}
#'   \item{bwl_zt}{Bandwidth used at the Zt* break pair (\code{NA} for "iid")}
#'   \item{za_min}{Minimum Za* test statistic}
#'   \item{tb1_za}{First break location for Za*}
#'   \item{tb2_za}{Second break location for Za*}
#'   \item{bwl_za}{Bandwidth used at the Za* break pair (\code{NA} for "iid")}
#'   \item{cv_adfzt}{Critical values for ADF* and Zt* tests (1\%, 5\%, 10\%)}
#'   \item{cv_za}{Critical values for Za* test (1\%, 5\%, 10\%)}
#'   \item{nobs}{Number of observations}
#'   \item{k}{Number of regressors}
#'   \item{maxlags}{Maximum lags used}
#'   \item{lag_selection}{Lag selection method used}
#'   \item{kernel}{Kernel used}
#'   \item{bwl}{Bandwidth argument (\code{NA} when chosen automatically)}
#'   \item{prewhite}{Whether prewhitening was used}
#'   \item{trimming}{Trimming parameter}
#'   \item{n_skipped}{Number of break pairs skipped because a statistic could
#'     not be computed}
#' }
#'
#' @references
#' Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts
#' with an application to financial market integration. \emph{Empirical Economics},
#' 35, 497-505. \doi{10.1007/s00181-007-0175-9}
#'
#' Gregory, A.W. and Hansen, B.E. (1996). Residual-based tests for cointegration
#' in models with regime shifts. \emph{Journal of Econometrics}, 70(1), 99-126.
#' \doi{10.1016/0304-4076(69)41685-7}
#'
#' Andrews, D.W.K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' Andrews, D.W.K. and Monahan, J.C. (1992). An improved heteroskedasticity and
#' autocorrelation consistent covariance matrix estimator. \emph{Econometrica},
#' 60(4), 953-966. \doi{10.2307/2951574}
#'
#' Phillips, P.C.B. (1987). Time series regression with a unit root.
#' \emph{Econometrica}, 55(2), 277-301. \doi{10.2307/1913237}
#'
#' @examples
#' \donttest{
#' # Generate example data with structural breaks
#' set.seed(123)
#' n <- 200
#' x <- cumsum(rnorm(n))
#'
#' # Create cointegrated series with two breaks
#' y <- numeric(n)
#' y[1:70] <- 1 + 0.8 * x[1:70] + rnorm(70, sd = 0.5)
#' y[71:140] <- 3 + 1.2 * x[71:140] + rnorm(70, sd = 0.5)
#' y[141:200] <- 2 + 0.6 * x[141:200] + rnorm(60, sd = 0.5)
#'
#' # Run the test
#' result <- hatemicoint(y, x)
#' print(result)
#' summary(result)
#' }
#'
#' @export
#' @importFrom stats lm.fit
hatemicoint <- function(y, x,
                        maxlags = NULL,
                        lag_selection = c("tstat", "aic", "sic", "fixed"),
                        kernel = c("qs", "bartlett", "iid"),
                        bwl = NULL,
                        prewhite = TRUE,
                        trimming = 0.15) {

  # Input validation
  lag_selection <- match.arg(lag_selection)
  kernel <- match.arg(kernel)

  # Convert to matrix if vector
  y <- as.numeric(y)
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  } else {
    x <- as.matrix(x)
  }

  n <- length(y)
  k <- ncol(x)

  # Validation checks
  if (length(y) != nrow(x)) {
    stop("Length of 'y' must equal number of rows in 'x'")
  }

  if (k > 4) {
    stop("Maximum of 4 independent variables allowed (k <= 4)")
  }

  if (anyNA(y) || anyNA(x)) {
    stop("'y' and 'x' must not contain missing values")
  }

  if (is.null(maxlags)) {
    maxlags <- floor(4 * (n / 100)^(1 / 4))
  }
  if (length(maxlags) != 1 || maxlags < 0 || maxlags != round(maxlags)) {
    stop("'maxlags' must be a non-negative integer")
  }
  maxlags <- as.integer(maxlags)

  if (trimming <= 0 || trimming >= 0.5) {
    stop("'trimming' must be between 0 and 0.5 (exclusive)")
  }

  if (!is.null(bwl) && (length(bwl) != 1 || !is.finite(bwl) || bwl <= 0)) {
    stop("'bwl' must be a positive number")
  }
  prewhite <- isTRUE(prewhite)

  if (n < 50) {
    warning("Sample size is small (< 50). Results may be unreliable.")
  }

  # Break point grid
  grid <- .break_grid(n, trimming)
  if (nrow(grid) == 0) {
    stop("No admissible break pairs for this sample size and trimming")
  }

  # Initialize minimum statistics
  adf_min <- Inf
  zt_min <- Inf
  za_min <- Inf
  tb1_adf <- tb2_adf <- NA_integer_
  tb1_zt <- tb2_zt <- NA_integer_
  tb1_za <- tb2_za <- NA_integer_
  lag_adf <- NA_integer_
  bwl_zt <- bwl_za <- NA_real_
  skipped <- character(0)

  # Grid search over break points
  for (r in seq_len(nrow(grid))) {
    tb1 <- unname(grid[r, 1])
    tb2 <- unname(grid[r, 2])

    # Create dummy variables for regime shifts
    du1 <- as.numeric(seq_len(n) > tb1)
    du2 <- as.numeric(seq_len(n) > tb2)

    # Model: y = a0 + a1*DU1 + a2*DU2 + b0'*x + b1'*DU1*x + b2'*DU2*x + u
    X_reg <- cbind(1, du1, du2)
    for (j in seq_len(k)) {
      X_reg <- cbind(X_reg, x[, j], du1 * x[, j], du2 * x[, j])
    }

    fit <- stats::lm.fit(X_reg, y)
    if (fit$rank < ncol(X_reg)) {
      skipped <- c(skipped, sprintf("(%d, %d): singular cointegrating regression",
                                    tb1, tb2))
      next
    }
    resid <- fit$residuals

    # Test statistics on residuals
    adf <- .compute_adf(resid, maxlags, lag_selection)
    pp <- .compute_pp(resid, bwl, kernel, prewhite)
    if (!is.null(adf$reason)) {
      skipped <- c(skipped, sprintf("(%d, %d): %s", tb1, tb2, adf$reason))
    }
    if (!is.null(pp$reason)) {
      skipped <- c(skipped, sprintf("(%d, %d): %s", tb1, tb2, pp$reason))
    }

    if (is.finite(adf$stat) && adf$stat < adf_min) {
      adf_min <- adf$stat
      tb1_adf <- tb1
      tb2_adf <- tb2
      lag_adf <- adf$lag
    }

    if (is.finite(pp$zt) && pp$zt < zt_min) {
      zt_min <- pp$zt
      tb1_zt <- tb1
      tb2_zt <- tb2
      bwl_zt <- pp$bw
    }

    if (is.finite(pp$za) && pp$za < za_min) {
      za_min <- pp$za
      tb1_za <- tb1
      tb2_za <- tb2
      bwl_za <- pp$bw
    }
  }

  if (length(skipped) > 0) {
    warning(sprintf(paste0("%d break pair(s) skipped because a statistic could ",
                           "not be computed, e.g. %s"),
                    length(skipped), skipped[1]))
  }
  if (!is.finite(adf_min)) adf_min <- NA_real_
  if (!is.finite(zt_min)) zt_min <- NA_real_
  if (!is.finite(za_min)) za_min <- NA_real_

  # Get critical values
  cv <- .get_critical_values(k)

  # Create result object
  result <- list(
    adf_min = adf_min,
    tb1_adf = tb1_adf,
    tb2_adf = tb2_adf,
    lag_adf = lag_adf,
    zt_min = zt_min,
    tb1_zt = tb1_zt,
    tb2_zt = tb2_zt,
    bwl_zt = bwl_zt,
    za_min = za_min,
    tb1_za = tb1_za,
    tb2_za = tb2_za,
    bwl_za = bwl_za,
    cv_adfzt = cv$cv_adfzt,
    cv_za = cv$cv_za,
    nobs = n,
    k = k,
    maxlags = maxlags,
    lag_selection = lag_selection,
    kernel = kernel,
    bwl = if (is.null(bwl)) NA_real_ else bwl,
    prewhite = prewhite,
    trimming = trimming,
    n_skipped = length(skipped),
    call = match.call()
  )

  class(result) <- "hatemicoint"
  return(result)
}


#' Print Method for hatemicoint Objects
#'
#' @param x An object of class \code{"hatemicoint"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.hatemicoint <- function(x, ...) {
  cat("\n")
  cat("Hatemi-J Cointegration Test with Two Unknown Regime Shifts\n")
  cat("===========================================================\n")
  cat("Reference: Hatemi-J (2008), Empirical Economics, 35, 497-505.\n")
  cat("           DOI: 10.1007/s00181-007-0175-9\n")
  cat("\n")
  cat(sprintf("Observations: %d | Regressors (k): %d | Trimming: %.2f\n",
              x$nobs, x$k, x$trimming))
  cat(sprintf("Lag selection: %s (max lags %d) | Kernel: %s",
              x$lag_selection, x$maxlags, x$kernel))
  if (x$kernel != "iid") {
    cat(sprintf(" | Prewhitening: %s | Bandwidth: %s",
                if (isTRUE(x$prewhite)) "AR(1)" else "none",
                .fmt_bw(x$bwl, x$bwl_zt, x$bwl_za)))
  }
  cat("\n\n")
  
  cat("Test Results:\n")
  cat("-------------\n")
  cat(sprintf("%-8s %12s %12s %12s\n", "Test", "Statistic", "Break 1", "Break 2"))
  cat(sprintf("%-8s %12.4f %12s %12s\n", "ADF*", x$adf_min,
              .fmt_int(x$tb1_adf), .fmt_int(x$tb2_adf)))
  cat(sprintf("%-8s %12.4f %12s %12s\n", "Zt*", x$zt_min,
              .fmt_int(x$tb1_zt), .fmt_int(x$tb2_zt)))
  cat(sprintf("%-8s %12.4f %12s %12s\n", "Za*", x$za_min,
              .fmt_int(x$tb1_za), .fmt_int(x$tb2_za)))
  cat("\n")
  
  cat("Critical Values (ADF* and Zt*):\n")
  cat(sprintf("  1%%: %.3f | 5%%: %.3f | 10%%: %.3f\n",
              x$cv_adfzt[1], x$cv_adfzt[2], x$cv_adfzt[3]))
  cat("\n")
  
  cat("Critical Values (Za*):\n")
  cat(sprintf("  1%%: %.3f | 5%%: %.3f | 10%%: %.3f\n",
              x$cv_za[1], x$cv_za[2], x$cv_za[3]))
  cat("\n")
  
  cat("H0: No cointegration\n")
  cat("Ha: Cointegration with two regime shifts\n")
  cat("Note: Reject H0 if test statistic < critical value\n")
  cat("\n")
  
  invisible(x)
}


#' Summary Method for hatemicoint Objects
#'
#' @param object An object of class \code{"hatemicoint"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a summary list with inference at various significance levels.
#'
#' @export
summary.hatemicoint <- function(object, ...) {
  cat("\n")
  cat("================================================================\n")
  cat("        HATEMI-J COINTEGRATION TEST WITH TWO REGIME SHIFTS\n")
  cat("================================================================\n")
  cat("\n")
  cat("Reference:\n")
  cat("  Hatemi-J, A. (2008). Tests for cointegration with two unknown\n")
  cat("  regime shifts with an application to financial market integration.\n")
  cat("  Empirical Economics, 35, 497-505.\n")
  cat("  DOI: 10.1007/s00181-007-0175-9\n")
  cat("\n")
  
  cat("Specifications:\n")
  cat("---------------\n")
  cat(sprintf("  Number of observations:    %d\n", object$nobs))
  cat(sprintf("  Number of regressors (k):  %d\n", object$k))
  cat(sprintf("  Trimming parameter:        %.2f\n", object$trimming))
  cat(sprintf("  Maximum ADF lags:          %d\n", object$maxlags))
  cat(sprintf("  Lag selection method:      %s\n", toupper(object$lag_selection)))
  cat(sprintf("  Kernel:                    %s\n", object$kernel))
  if (object$kernel != "iid") {
    cat(sprintf("  Prewhitening:              %s\n",
                if (isTRUE(object$prewhite)) "AR(1)" else "none"))
    cat(sprintf("  Bandwidth:                 %s\n",
                .fmt_bw(object$bwl, object$bwl_zt, object$bwl_za)))
  }
  if (isTRUE(object$n_skipped > 0)) {
    cat(sprintf("  Break pairs skipped:       %d\n", object$n_skipped))
  }
  cat("\n")
  
  # Determine significance
  .sig_stars <- function(stat, cv) {
    if (is.na(stat)) return("")
    if (stat < cv[1]) return("***")
    if (stat < cv[2]) return("**")
    if (stat < cv[3]) return("*")
    return("")
  }
  
  cat("Test Results:\n")
  cat("-------------\n")
  cat(sprintf("%-8s %12s %10s %10s %8s\n",
              "Test", "Statistic", "Break 1", "Break 2", "Signif."))
  cat("--------------------------------------------------------\n")
  
  adf_sig <- .sig_stars(object$adf_min, object$cv_adfzt)
  zt_sig <- .sig_stars(object$zt_min, object$cv_adfzt)
  za_sig <- .sig_stars(object$za_min, object$cv_za)
  
  cat(sprintf("%-8s %12.4f %10s %10s %8s\n",
              "ADF*", object$adf_min, .fmt_int(object$tb1_adf),
              .fmt_int(object$tb2_adf), adf_sig))
  cat(sprintf("%-8s %12.4f %10s %10s %8s\n",
              "Zt*", object$zt_min, .fmt_int(object$tb1_zt),
              .fmt_int(object$tb2_zt), zt_sig))
  cat(sprintf("%-8s %12.4f %10s %10s %8s\n",
              "Za*", object$za_min, .fmt_int(object$tb1_za),
              .fmt_int(object$tb2_za), za_sig))
  cat("\n")
  cat("Significance codes: '***' 1%, '**' 5%, '*' 10%\n")
  cat("\n")
  
  cat("Critical Values:\n")
  cat("----------------\n")
  cat(sprintf("%-12s %12s %12s %12s\n", "Test", "1%", "5%", "10%"))
  cat("----------------------------------------------------\n")
  cat(sprintf("%-12s %12.3f %12.3f %12.3f\n",
              "ADF*/Zt*", object$cv_adfzt[1], object$cv_adfzt[2], object$cv_adfzt[3]))
  cat(sprintf("%-12s %12.3f %12.3f %12.3f\n",
              "Za*", object$cv_za[1], object$cv_za[2], object$cv_za[3]))
  cat("\n")
  
  cat("Interpretation:\n")
  cat("---------------\n")
  cat("H0: No cointegration (null hypothesis)\n")
  cat("Ha: Cointegration exists with two structural breaks\n")
  cat("\n")
  cat("Decision rule: Reject H0 if test statistic < critical value\n")
  cat("             (more negative values provide stronger evidence against H0)\n")
  cat("\n")
  
  # Return summary invisibly
  summary_list <- list(
    adf_reject_1pct = object$adf_min < object$cv_adfzt[1],
    adf_reject_5pct = object$adf_min < object$cv_adfzt[2],
    adf_reject_10pct = object$adf_min < object$cv_adfzt[3],
    zt_reject_1pct = object$zt_min < object$cv_adfzt[1],
    zt_reject_5pct = object$zt_min < object$cv_adfzt[2],
    zt_reject_10pct = object$zt_min < object$cv_adfzt[3],
    za_reject_1pct = object$za_min < object$cv_za[1],
    za_reject_5pct = object$za_min < object$cv_za[2],
    za_reject_10pct = object$za_min < object$cv_za[3]
  )
  
  invisible(summary_list)
}


# Formatting helpers for print and summary methods
.fmt_int <- function(v) if (is.na(v)) "NA" else sprintf("%d", as.integer(v))

.fmt_bw <- function(bwl, bw_zt, bw_za) {
  if (!is.na(bwl)) return(format(bwl))
  if (is.na(bw_zt) && is.na(bw_za)) return("automatic")
  sprintf("automatic (Zt*: %.2f, Za*: %.2f)", bw_zt, bw_za)
}
