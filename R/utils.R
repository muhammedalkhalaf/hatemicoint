#' Compute ADF Test Statistic on Residuals
#'
#' @description
#' Internal function to compute the Augmented Dickey-Fuller test statistic
#' on regression residuals, without intercept or trend. The regression is
#' \deqn{\Delta u_t = \delta u_{t-1} + \sum_{i=1}^{p} \phi_i \Delta u_{t-i} + e_t}
#' and the statistic is the t-ratio of \eqn{\delta} (OLS standard error with
#' the degrees-of-freedom divisor).
#'
#' Lag selection: all candidate lags \code{0, ..., maxlags} are compared on the
#' common sample \code{t = maxlags + 2, ..., n} (the first \code{maxlags} usable
#' rows are dropped for every candidate). Once the lag \code{p} is chosen, the
#' final statistic is computed on the largest sample available for that lag,
#' \code{t = p + 2, ..., n}. With \code{ic = "fixed"} the lag is \code{maxlags}.
#'
#' @param resid Numeric vector of residuals.
#' @param maxlags Maximum number of lags to consider (or the fixed lag).
#' @param ic Lag selection rule: "tstat", "aic", "sic" or "fixed".
#'
#' @return A list with \code{stat} (the ADF t-statistic, \code{NA} if the
#'   regressor matrix is singular) and \code{lag} (the lag used).
#'
#' @keywords internal
#' @noRd
.compute_adf <- function(resid, maxlags, ic) {
  n <- length(resid)
  delta_resid <- diff(resid)          # Delta u_t, t = 2..n  (length n - 1)
  resid_lag <- resid[-n]              # u_{t-1},   t = 2..n  (length n - 1)
  T_eff <- n - 1

  # Regressor matrix for lag p on rows t = start..n (row index in 1..T_eff)
  .design <- function(p, start) {
    idx <- start:T_eff
    X <- matrix(resid_lag[idx], ncol = 1)
    if (p > 0) {
      for (i in 1:p) {
        X <- cbind(X, delta_resid[idx - i])
      }
    }
    list(y = delta_resid[idx], X = X)
  }

  .fit <- function(p, start) {
    d <- .design(p, start)
    if (length(d$y) < ncol(d$X) + 2) return(NULL)
    fit <- stats::lm.fit(d$X, d$y)
    if (fit$rank < ncol(d$X)) return(NULL)
    df_resid <- length(d$y) - ncol(d$X)
    sigma2 <- sum(fit$residuals^2) / df_resid
    XtX_inv <- tryCatch(solve(crossprod(d$X)), error = function(e) NULL)
    if (is.null(XtX_inv)) return(NULL)
    se <- sqrt(sigma2 * diag(XtX_inv))
    list(coef = fit$coefficients, se = se,
         rss = sum(fit$residuals^2), nobs = length(d$y))
  }

  # Common sample for all candidates: rows maxlags + 1 .. T_eff
  common_start <- maxlags + 1
  if (common_start > T_eff - 2) {
    return(list(stat = NA_real_, lag = NA_integer_,
                reason = "too few observations for the ADF regression"))
  }

  if (ic == "fixed") {
    lag_optimal <- maxlags
  } else if (ic == "tstat") {
    # General-to-specific: start at maxlags, keep the first lag whose last
    # coefficient is significant at the 10 percent level (|t| > 1.645)
    lag_optimal <- 0
    if (maxlags > 0) {
      for (p in maxlags:1) {
        f <- .fit(p, common_start)
        if (is.null(f)) next
        t_stat <- abs(f$coef[p + 1] / f$se[p + 1])
        if (is.finite(t_stat) && t_stat > 1.645) {
          lag_optimal <- p
          break
        }
      }
    }
  } else {
    best_ic <- Inf
    lag_optimal <- 0
    for (p in 0:maxlags) {
      f <- .fit(p, common_start)
      if (is.null(f)) next
      pen <- if (ic == "aic") 2 else log(f$nobs)
      ic_val <- log(f$rss / f$nobs) + pen * (1 + p) / f$nobs
      if (is.finite(ic_val) && ic_val < best_ic) {
        best_ic <- ic_val
        lag_optimal <- p
      }
    }
  }

  # Final ADF regression on the largest sample for the chosen lag
  f <- .fit(lag_optimal, lag_optimal + 1)
  if (is.null(f)) {
    return(list(stat = NA_real_, lag = lag_optimal,
                reason = "singular regressor matrix in the ADF regression"))
  }
  list(stat = unname(f$coef[1] / f$se[1]), lag = lag_optimal, reason = NULL)
}


#' Kernel weights
#'
#' @description
#' Bartlett kernel \code{1 - x} for \code{x < 1} (zero beyond the bandwidth)
#' and the quadratic spectral (QS) kernel of Andrews (1991),
#' \code{k(x) = 3 / z^2 * (sin(z) / z - cos(z))} with \code{z = 6 pi x / 5}.
#'
#' @param x Numeric vector of lag / bandwidth ratios.
#' @param kernel "bartlett" or "qs".
#' @keywords internal
#' @noRd
.kernel_weight <- function(x, kernel) {
  if (kernel == "bartlett") {
    w <- 1 - x
    w[x >= 1] <- 0
    return(w)
  }
  z <- 6 * pi * x / 5
  w <- 3 / z^2 * (sin(z) / z - cos(z))
  w[abs(z) < 1e-10] <- 1
  w
}


#' Andrews (1991) automatic bandwidth for the QS kernel
#'
#' @description
#' Fits an AR(1) (no intercept) to \code{e}, obtains \code{rho} and computes
#' \code{alpha2 = 4 rho^2 / (1 - rho)^4} and \code{S = 1.3221 (alpha2 T)^(1/5)},
#' where \code{T} is the length of \code{e}. The AR(1) coefficient is capped at
#' 0.97 in absolute value.
#'
#' @param e Numeric vector.
#' @return The bandwidth (a non-negative number, not rounded; zero when the
#'   AR(1) coefficient is exactly zero).
#' @keywords internal
#' @noRd
.andrews_bandwidth_qs <- function(e) {
  T <- length(e)
  rho <- sum(e[-1] * e[-T]) / sum(e[-T]^2)
  if (!is.finite(rho)) rho <- 0
  # cap at 0.97 in absolute value (package choice, as for the prewhitening)
  rho <- max(min(rho, 0.97), -0.97)
  alpha2 <- 4 * rho^2 / (1 - rho)^4
  1.3221 * (alpha2 * T)^(1 / 5)
}


#' Long-run variance of a series with a kernel
#'
#' @description
#' Returns \code{gamma(0) + 2 sum_j w(j / B) gamma(j)} with
#' \code{gamma(j) = (1 / norm) sum_t e_t e_{t-j}}. The Bartlett sum stops at
#' \code{j = B}; the QS sum runs over all \code{j = 1, ..., length(e) - 1}
#' because the QS kernel is not truncated.
#'
#' @keywords internal
#' @noRd
.lrvar_kernel <- function(e, bw, kernel, norm) {
  N <- length(e)
  gamma0 <- sum(e^2) / norm
  if (kernel == "iid" || N < 2) return(gamma0)
  # A zero bandwidth (Andrews rule with an AR(1) coefficient of exactly zero)
  # gives no weight to any lag: the estimate is gamma(0)
  if (!is.finite(bw) || bw <= 0) return(gamma0)
  jmax <- if (kernel == "bartlett") min(floor(bw), N - 1) else N - 1
  if (jmax < 1) return(gamma0)
  j <- seq_len(jmax)
  # stats::acf divides the cross products by N; rescale to the 1 / norm form
  ac <- stats::acf(e, lag.max = jmax, type = "covariance", demean = FALSE,
                   plot = FALSE)$acf[-1]
  gammaj <- ac * N / norm
  gamma0 + 2 * sum(.kernel_weight(j / bw, kernel) * gammaj)
}


#' Compute Phillips-Perron Test Statistics on Residuals
#'
#' @description
#' Internal function to compute the Phillips-Perron \eqn{Z_t} and \eqn{Z_\alpha}
#' statistics on regression residuals following Hatemi-J (2008), equations
#' (3) to (6), with the normalisation of the paper: sums over
#' \eqn{t = 1, \ldots, n - 1} of \eqn{u_t^2}, autocovariances divided by
#' \eqn{n}, and \eqn{Z_\alpha = n (\rho^* - 1)}.
#'
#' Steps:
#' \enumerate{
#'   \item \eqn{\hat\rho} from the regression of \eqn{u_t} on \eqn{u_{t-1}}
#'     (no intercept); \eqn{v_t = u_t - \hat\rho u_{t-1}}, \eqn{t = 2, \ldots, n}.
#'   \item \eqn{\hat\gamma(j) = n^{-1} \sum_t v_t v_{t-j}} (eq. 4).
#'   \item Long-run variance \eqn{\hat\sigma^2 = \hat\gamma(0) + 2 \sum_j w(j/B)
#'     \hat\gamma(j)}. With \code{prewhite = TRUE} the kernel is applied to the
#'     residuals of an AR(1) fitted to \eqn{v_t} (coefficient capped at 0.97 in
#'     absolute value) and the result is recoloured by dividing by
#'     \eqn{(1 - a_1)^2}.
#'   \item One-sided correction \eqn{\hat\lambda = (\hat\sigma^2 - \hat\gamma(0)) / 2},
#'     which equals \eqn{\sum_{j \ge 1} w(j/B) \hat\gamma(j)} without prewhitening.
#'   \item \eqn{\hat\rho^* = (\sum_{t=1}^{n-1} u_t u_{t+1} - n \hat\lambda) /
#'     \sum_{t=1}^{n-1} u_t^2}. The factor \eqn{n} undoes the \eqn{1/n} in
#'     \eqn{\hat\gamma(j)}, as in Phillips (1987); without it the correction
#'     would be of order \eqn{1/n^2} relative to the leading term.
#'   \item \eqn{Z_\alpha = n(\hat\rho^* - 1)} and
#'     \eqn{Z_t = (\hat\rho^* - 1) / \sqrt{\hat\sigma^2 / \sum_{t=1}^{n-1} u_t^2}}.
#' }
#' With \code{kernel = "iid"} there is no correction and \eqn{Z_t} equals the
#' Dickey-Fuller t-ratio up to the residual-variance divisor (\eqn{n} here,
#' \eqn{n - 2} in the OLS t-ratio).
#'
#' @param resid Numeric vector of residuals.
#' @param bwl Bandwidth. \code{NULL} means the automatic rule for "qs".
#' @param kernel Kernel type: "iid", "bartlett", or "qs".
#' @param prewhite Logical. AR(1) prewhitening (ignored for "iid").
#'
#' @return A list with components \code{zt}, \code{za}, \code{bw} (bandwidth
#'   used, \code{NA} for "iid"), \code{lrvar}, \code{rho_star} and
#'   \code{reason} (\code{NULL} or a string when the statistics are \code{NA}).
#'
#' @keywords internal
#' @noRd
.compute_pp <- function(resid, bwl, kernel, prewhite = TRUE) {
  u <- as.numeric(resid)
  n <- length(u)
  bad <- function(reason) list(zt = NA_real_, za = NA_real_, bw = NA_real_,
                               lrvar = NA_real_, rho_star = NA_real_,
                               reason = reason)

  S <- sum(u[-n]^2)                       # sum_{t=1}^{n-1} u_t^2
  if (!is.finite(S) || S <= 0) return(bad("zero residual sum of squares"))
  cross <- sum(u[-n] * u[-1])             # sum_{t=1}^{n-1} u_t u_{t+1}
  rho_hat <- cross / S
  v <- u[-1] - rho_hat * u[-n]            # t = 2..n, length n - 1
  gamma0 <- sum(v^2) / n                  # eq. (4), j = 0

  bw <- NA_real_
  if (kernel == "iid") {
    lrvar <- gamma0
  } else {
    if (isTRUE(prewhite)) {
      N <- length(v)
      a1 <- sum(v[-1] * v[-N]) / sum(v[-N]^2)
      if (!is.finite(a1)) a1 <- 0
      a1 <- max(min(a1, 0.97), -0.97)
      e <- v[-1] - a1 * v[-N]
      recolour <- 1 / (1 - a1)^2
    } else {
      e <- v
      recolour <- 1
    }
    if (is.null(bwl)) {
      if (kernel == "qs") {
        bw <- .andrews_bandwidth_qs(e)
      } else {
        bw <- round(4 * (n / 100)^(2 / 9))
      }
    } else {
      bw <- bwl
    }
    lrvar <- recolour * .lrvar_kernel(e, bw, kernel, norm = n)
  }

  if (!is.finite(lrvar) || lrvar <= 0) {
    return(bad("non-positive long-run variance estimate"))
  }

  lambda <- (lrvar - gamma0) / 2           # one-sided sum of weighted gamma(j)
  rho_star <- (cross - n * lambda) / S     # eq. (3)
  za <- n * (rho_star - 1)                 # eq. (5)
  zt <- (rho_star - 1) / sqrt(lrvar / S)   # eq. (6)

  list(zt = unname(zt), za = unname(za), bw = bw, lrvar = lrvar,
       rho_star = unname(rho_star), reason = NULL)
}


#' Break-point grid
#'
#' @description
#' Returns the admissible pairs (tb1, tb2) of break observations for a sample
#' of size n and trimming tau: tb1 from \code{ceiling(tau n)} to
#' \code{floor((1 - 2 tau) n)}, tb2 from \code{tb1 + ceiling(tau n)} to
#' \code{floor((1 - tau) n)}. This guarantees tb1 >= tau n, tb2 <= (1 - tau) n
#' and tb2 - tb1 >= tau n for every n.
#'
#' @keywords internal
#' @noRd
.break_grid <- function(n, trimming) {
  lo <- ceiling(trimming * n)
  hi1 <- floor((1 - 2 * trimming) * n)
  hi2 <- floor((1 - trimming) * n)
  empty <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("tb1", "tb2")))
  if (hi1 < lo) return(empty)
  tb1 <- lo:hi1
  counts <- pmax(0L, hi2 - (tb1 + lo) + 1L)
  if (sum(counts) == 0) return(empty)
  tb1 <- tb1[counts > 0]
  counts <- counts[counts > 0]
  out <- cbind(tb1 = rep(tb1, counts),
               tb2 = sequence(counts, from = tb1 + lo))
  storage.mode(out) <- "integer"
  out
}


#' Get Critical Values for Hatemi-J Test
#'
#' @description
#' Internal function to retrieve critical values based on the number of
#' regressors. Critical values are from Hatemi-J (2008), Table 1.
#'
#' @param k Number of regressors (1 to 4).
#'
#' @return A list with \code{cv_adfzt} and \code{cv_za} vectors containing
#'   critical values at 1\%, 5\%, and 10\% significance levels.
#'
#' @references
#' Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts
#' with an application to financial market integration. \emph{Empirical Economics},
#' 35, 497-505. \doi{10.1007/s00181-007-0175-9}
#'
#' @keywords internal
#' @noRd
.get_critical_values <- function(k) {
  # Critical values from Hatemi-J (2008), Table 1
  # Columns: 1%, 5%, 10%

  cv_adfzt_table <- matrix(c(
    -6.503, -6.015, -5.653,  # k = 1
    -6.928, -6.458, -6.224,  # k = 2
    -7.833, -7.352, -7.118,  # k = 3
    -8.353, -7.903, -7.705   # k = 4
  ), nrow = 4, ncol = 3, byrow = TRUE)

  cv_za_table <- matrix(c(
    -90.794, -76.003, -52.232,   # k = 1
    -99.458, -83.644, -76.806,   # k = 2
    -118.577, -104.860, -97.749, # k = 3
    -140.135, -123.870, -116.169 # k = 4
  ), nrow = 4, ncol = 3, byrow = TRUE)

  return(list(
    cv_adfzt = cv_adfzt_table[k, ],
    cv_za = cv_za_table[k, ]
  ))
}
