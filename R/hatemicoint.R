#' Hatemi-J Cointegration Test with Two Unknown Regime Shifts
#'
#' @description
#' Performs the Hatemi-J (2008) cointegration test which allows for two unknown
#' structural breaks (regime shifts) in the cointegrating relationship. The test
#' searches over all possible break date combinations and returns the minimum
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
#' where \eqn{D_{1t}} and \eqn{D_{2t}} are dummy variables for the two regime
#' shifts.
#'
#' Three test statistics are computed:
#' \itemize{
#'   \item \strong{ADF*}: Modified Augmented Dickey-Fuller test on residuals
#'   \item \strong{Zt*}: Phillips-Perron Z_t test on residuals
#'   \item \strong{Za*}: Phillips-Perron Z_alpha test on residuals
#' }
#'
#' The null hypothesis is no cointegration. Rejection occurs when the test
#' statistic is smaller (more negative) than the critical value.
#'
#' @param y Numeric vector. The dependent variable (must be I(1)).
#' @param x Numeric matrix or vector. The independent variable(s) (must be I(1)).
#'   Maximum of 4 regressors allowed (k <= 4).
#' @param maxlags Integer. Maximum number of lags for ADF test. Default is 8.
#' @param lag_selection Character. Lag selection criterion: \code{"tstat"} (default),
#'   \code{"aic"}, or \code{"sic"}.
#' @param kernel Character. Kernel for long-run variance estimation in PP tests:
#'   \code{"iid"} (default), \code{"bartlett"}, or \code{"qs"} (quadratic spectral).
#' @param bwl Integer. Bandwidth for kernel estimation. If NULL (default), computed
#'   as \code{round(4 * (n/100)^(2/9))}.
#' @param trimming Numeric. Trimming parameter for break point search. Must be
#'   between 0 and 0.5 (exclusive). Default is 0.15.
#'
#' @return An object of class \code{"hatemicoint"} containing:
#' \describe{
#'   \item{adf_min}{Minimum ADF* test statistic}
#'   \item{tb1_adf}{First break location (observation number) for ADF*}
#'   \item{tb2_adf}{Second break location (observation number) for ADF*}
#'   \item{zt_min}{Minimum Zt* test statistic}
#'   \item{tb1_zt}{First break location for Zt*}
#'   \item{tb2_zt}{Second break location for Zt*}
#'   \item{za_min}{Minimum Za* test statistic}
#'   \item{tb1_za}{First break location for Za*}
#'   \item{tb2_za}{Second break location for Za*}
#'   \item{cv_adfzt}{Critical values for ADF* and Zt* tests (1\%, 5\%, 10\%)}
#'   \item{cv_za}{Critical values for Za* test (1\%, 5\%, 10\%)}
#'   \item{nobs}{Number of observations}
#'   \item{k}{Number of regressors}
#'   \item{maxlags}{Maximum lags used}
#'   \item{lag_selection}{Lag selection method used}
#'   \item{kernel}{Kernel used}
#'   \item{bwl}{Bandwidth used}
#'   \item{trimming}{Trimming parameter}
#' }
#'
#' @references
#' Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts
#' with an application to financial market integration. \emph{Empirical Economics},
#' 35, 497-505. \doi{10.1007/s00181-007-0175-9}
#'
#' Gregory, A.W., & Hansen, B.E. (1996). Residual-based tests for cointegration
#' in models with regime shifts. \emph{Journal of Econometrics}, 70(1), 99-126.
#' \doi{10.1016/0304-4076(69)41685-7}
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
#' @importFrom stats lm residuals coef var lm.fit complete.cases
hatemicoint <- function(y, x,
                        maxlags = 8,
                        lag_selection = c("tstat", "aic", "sic"),
                        kernel = c("iid", "bartlett", "qs"),
                        bwl = NULL,
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
  
  if (maxlags < 0) {
    stop("'maxlags' must be non-negative")
  }
  
  if (trimming <= 0 || trimming >= 0.5) {
    stop("'trimming' must be between 0 and 0.5 (exclusive)")
  }
  
  if (n < 50) {
    warning("Sample size is small (< 50). Results may be unreliable.")
  }
  
  # Compute bandwidth if not provided
  if (is.null(bwl)) {
    bwl <- round(4 * (n / 100)^(2/9))
  } else if (bwl <= 0) {
    stop("'bwl' must be positive")
  }
  
  # Compute break point search bounds
  T1 <- round(trimming * n)
  T2 <- round((1 - 2 * trimming) * n)
  T3 <- round((1 - trimming) * n)
  
  # Initialize minimum statistics
  adf_min <- Inf
  zt_min <- Inf
  za_min <- Inf
  tb1_adf <- tb2_adf <- 0
  tb1_zt <- tb2_zt <- 0
  tb1_za <- tb2_za <- 0
  
  # Grid search over break points
  for (tb1 in T1:T2) {
    for (tb2 in (tb1 + T1):T3) {
      
      # Create dummy variables for regime shifts
      du1 <- as.numeric(seq_len(n) > tb1)
      du2 <- as.numeric(seq_len(n) > tb2)
      
      # Build regressor matrix with interactions
      # Model: y = a0 + a1*DU1 + a2*DU2 + b0'*x + b1'*DU1*x + b2'*DU2*x + u
      X_reg <- cbind(1, du1, du2)
      
      for (j in seq_len(k)) {
        X_reg <- cbind(X_reg, x[, j], du1 * x[, j], du2 * x[, j])
      }
      
      # Estimate cointegrating regression
      fit <- lm.fit(X_reg, y)
      resid <- fit$residuals
      
      # Compute test statistics on residuals
      adf_stat <- .compute_adf(resid, maxlags, lag_selection)
      pp_stats <- .compute_pp(resid, bwl, kernel)
      zt_stat <- pp_stats$zt
      za_stat <- pp_stats$za
      
      # Update minimum statistics
      if (adf_stat < adf_min) {
        adf_min <- adf_stat
        tb1_adf <- tb1
        tb2_adf <- tb2
      }
      
      if (zt_stat < zt_min) {
        zt_min <- zt_stat
        tb1_zt <- tb1
        tb2_zt <- tb2
      }
      
      if (za_stat < za_min) {
        za_min <- za_stat
        tb1_za <- tb1
        tb2_za <- tb2
      }
    }
  }
  
  # Get critical values
  cv <- .get_critical_values(k)
  
  # Create result object
  result <- list(
    adf_min = adf_min,
    tb1_adf = tb1_adf,
    tb2_adf = tb2_adf,
    zt_min = zt_min,
    tb1_zt = tb1_zt,
    tb2_zt = tb2_zt,
    za_min = za_min,
    tb1_za = tb1_za,
    tb2_za = tb2_za,
    cv_adfzt = cv$cv_adfzt,
    cv_za = cv$cv_za,
    nobs = n,
    k = k,
    maxlags = maxlags,
    lag_selection = lag_selection,
    kernel = kernel,
    bwl = bwl,
    trimming = trimming,
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
  cat(sprintf("Lag selection: %s | Kernel: %s",
              x$lag_selection, x$kernel))
  if (x$kernel != "iid") {
    cat(sprintf(" | Bandwidth: %d", x$bwl))
  }
  cat("\n\n")
  
  cat("Test Results:\n")
  cat("-------------\n")
  cat(sprintf("%-8s %12s %12s %12s\n", "Test", "Statistic", "Break 1", "Break 2"))
  cat(sprintf("%-8s %12.4f %12d %12d\n", "ADF*", x$adf_min, x$tb1_adf, x$tb2_adf))
  cat(sprintf("%-8s %12.4f %12d %12d\n", "Zt*", x$zt_min, x$tb1_zt, x$tb2_zt))
  cat(sprintf("%-8s %12.4f %12d %12d\n", "Za*", x$za_min, x$tb1_za, x$tb2_za))
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
    cat(sprintf("  Bandwidth:                 %d\n", object$bwl))
  }
  cat("\n")
  
  # Determine significance
  .sig_stars <- function(stat, cv) {
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
  
  cat(sprintf("%-8s %12.4f %10d %10d %8s\n",
              "ADF*", object$adf_min, object$tb1_adf, object$tb2_adf, adf_sig))
  cat(sprintf("%-8s %12.4f %10d %10d %8s\n",
              "Zt*", object$zt_min, object$tb1_zt, object$tb2_zt, zt_sig))
  cat(sprintf("%-8s %12.4f %10d %10d %8s\n",
              "Za*", object$za_min, object$tb1_za, object$tb2_za, za_sig))
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
