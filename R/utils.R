#' Compute ADF Test Statistic on Residuals
#'
#' @description
#' Internal function to compute the Augmented Dickey-Fuller test statistic
#' on regression residuals, without intercept or trend.
#'
#' @param resid Numeric vector of residuals.
#' @param maxlags Maximum number of lags to consider.
#' @param ic Information criterion for lag selection: "tstat", "aic", or "sic".
#'
#' @return The ADF t-statistic.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats lm.fit complete.cases
.compute_adf <- function(resid, maxlags, ic) {
  n <- length(resid)
  
  # Create first difference
  delta_resid <- diff(resid)
  # Lagged level (for time 2:n, the lagged value is resid[1:(n-1)])
  resid_lag <- resid[-n]
  
  # Now both delta_resid and resid_lag have length n-1
  T_eff <- n - 1
  
  # Function to create lagged differences matrix
  .create_lag_matrix <- function(p) {
    if (p == 0) {
      return(matrix(resid_lag, ncol = 1))
    }
    
    # Start with lagged level
    X <- matrix(NA, nrow = T_eff, ncol = 1 + p)
    X[, 1] <- resid_lag
    
    # Add lagged differences
    for (i in 1:p) {
      lag_diff <- c(rep(NA, i), delta_resid[1:(T_eff - i)])
      X[, i + 1] <- lag_diff
    }
    return(X)
  }
  
  # Determine optimal lag using specified criterion
  if (ic == "tstat") {
    # Sequential t-stat approach: start from maxlags and reduce
    lag_optimal <- 0
    
    for (p in maxlags:0) {
      X <- .create_lag_matrix(p)
      
      # Find valid observations (no NAs)
      valid <- complete.cases(X, delta_resid)
      
      if (sum(valid) < (p + 3)) next
      
      y_reg <- delta_resid[valid]
      X_reg <- X[valid, , drop = FALSE]
      
      # Fit regression without intercept
      fit <- stats::lm.fit(X_reg, y_reg)
      
      if (p > 0) {
        # Check t-stat of last lag
        df_resid <- length(y_reg) - ncol(X_reg)
        if (df_resid <= 0) next
        
        sigma2 <- sum(fit$residuals^2) / df_resid
        XtX_inv <- tryCatch(solve(crossprod(X_reg)), error = function(e) NULL)
        
        if (!is.null(XtX_inv)) {
          se <- sqrt(sigma2 * diag(XtX_inv))
          t_stat <- abs(fit$coefficients[p + 1] / se[p + 1])
          
          if (!is.na(t_stat) && is.finite(t_stat) && t_stat > 1.645) {
            lag_optimal <- p
            break
          }
        }
      } else {
        lag_optimal <- 0
        break
      }
    }
  } else if (ic == "aic") {
    # AIC selection
    best_ic <- Inf
    lag_optimal <- 0
    
    for (p in 0:maxlags) {
      X <- .create_lag_matrix(p)
      valid <- complete.cases(X, delta_resid)
      
      if (sum(valid) < (p + 3)) next
      
      y_reg <- delta_resid[valid]
      X_reg <- X[valid, , drop = FALSE]
      
      fit <- stats::lm.fit(X_reg, y_reg)
      rss <- sum(fit$residuals^2)
      n_valid <- length(y_reg)
      
      aic_val <- log(rss / n_valid) + 2 * (1 + p) / n_valid
      
      if (is.finite(aic_val) && aic_val < best_ic) {
        best_ic <- aic_val
        lag_optimal <- p
      }
    }
  } else {
    # SIC (BIC) selection
    best_ic <- Inf
    lag_optimal <- 0
    
    for (p in 0:maxlags) {
      X <- .create_lag_matrix(p)
      valid <- complete.cases(X, delta_resid)
      
      if (sum(valid) < (p + 3)) next
      
      y_reg <- delta_resid[valid]
      X_reg <- X[valid, , drop = FALSE]
      
      fit <- stats::lm.fit(X_reg, y_reg)
      rss <- sum(fit$residuals^2)
      n_valid <- length(y_reg)
      
      sic_val <- log(rss / n_valid) + (1 + p) * log(n_valid) / n_valid
      
      if (is.finite(sic_val) && sic_val < best_ic) {
        best_ic <- sic_val
        lag_optimal <- p
      }
    }
  }
  
  # Final ADF regression with optimal lag
  X <- .create_lag_matrix(lag_optimal)
  valid <- complete.cases(X, delta_resid)
  
  y_reg <- delta_resid[valid]
  X_reg <- X[valid, , drop = FALSE]
  
  fit <- stats::lm.fit(X_reg, y_reg)
  
  # Compute t-statistic for rho_hat (first coefficient)
  df_resid <- length(y_reg) - ncol(X_reg)
  sigma2 <- sum(fit$residuals^2) / df_resid
  XtX_inv <- tryCatch(solve(crossprod(X_reg)), error = function(e) {
    diag(1e10, ncol(X_reg))
  })
  se <- sqrt(sigma2 * XtX_inv[1, 1])
  
  adf_stat <- fit$coefficients[1] / se
  
  return(adf_stat)
}


#' Compute Phillips-Perron Test Statistics on Residuals
#'
#' @description
#' Internal function to compute the Phillips-Perron Zt and Za test statistics
#' on regression residuals.
#'
#' @param resid Numeric vector of residuals.
#' @param bwl Bandwidth for kernel estimation.
#' @param kernel Kernel type: "iid", "bartlett", or "qs".
#'
#' @return A list with components \code{zt} and \code{za}.
#'
#' @keywords internal
#' @noRd
#' @importFrom stats lm.fit
.compute_pp <- function(resid, bwl, kernel) {
  n <- length(resid)
  
  # Create first difference and lagged level
  delta_resid <- diff(resid)
  resid_lag <- resid[-n]
  T_eff <- n - 1
  
  # Fit regression: delta_u = rho * u_{t-1} + e (no constant)
  X_reg <- matrix(resid_lag, ncol = 1)
  fit <- stats::lm.fit(X_reg, delta_resid)
  rho_hat <- fit$coefficients[1]
  u_hat <- fit$residuals
  
  # Compute gamma_0 (variance of innovations)
  gamma0 <- sum(u_hat^2) / T_eff
  
  # Compute long-run variance
  lrvar <- gamma0
  
  if (kernel != "iid") {
    for (j in 1:min(bwl, T_eff - 1)) {
      # Compute autocovariance at lag j
      gammaj <- sum(u_hat[(j + 1):T_eff] * u_hat[1:(T_eff - j)]) / T_eff
      
      # Kernel weight
      if (kernel == "bartlett") {
        weight <- 1 - j / (bwl + 1)
      } else {
        # Quadratic spectral kernel
        x <- 6 * pi * j / bwl
        if (abs(x) < 1e-10) {
          weight <- 1
        } else {
          weight <- 3 / (x^2) * (sin(x) / x - cos(x))
        }
      }
      
      lrvar <- lrvar + 2 * weight * gammaj
    }
  }
  
  # Ensure lrvar is positive
  lrvar <- max(lrvar, 1e-10)
  
  # Compute sum of squared residuals (levels)
  sum_resid2 <- sum((resid - mean(resid))^2)
  
  # Lambda squared
  lambda2 <- lrvar - gamma0
  
  # Phillips-Perron adjustments
  rho_star <- rho_hat - lambda2 / sum_resid2
  
  # Za statistic
  za_stat <- T_eff * rho_star
  
  # Zt statistic
  zt_stat <- rho_star * sqrt(T_eff / lrvar * sum_resid2)
  
  return(list(zt = zt_stat, za = za_stat))
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
