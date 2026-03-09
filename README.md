# hatemicoint

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/hatemicoint)](https://CRAN.R-project.org/package=hatemicoint)
<!-- badges: end -->

## Overview

**hatemicoint** implements the Hatemi-J (2008) cointegration test which allows for **two unknown structural breaks** (regime shifts) in the cointegrating relationship.

Standard cointegration tests may fail to detect cointegration when the relationship between variables has changed over time due to structural breaks. The Hatemi-J test addresses this by:

1. Searching over all possible pairs of break dates
2. Finding the break dates that minimize the test statistics
3. Providing critical values that account for the data-driven break selection

## Installation

```r
# Install from CRAN (when available)
install.packages("hatemicoint")

# Or install the development version from GitHub
# install.packages("devtools")
```

## Usage

```r
library(hatemicoint)

# Generate example data with structural breaks
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))  # I(1) regressor

# Create cointegrated series with two breaks
y <- numeric(n)
y[1:70] <- 1 + 0.8 * x[1:70] + rnorm(70, sd = 0.5)
y[71:140] <- 3 + 1.2 * x[71:140] + rnorm(70, sd = 0.5)
y[141:200] <- 2 + 0.6 * x[141:200] + rnorm(60, sd = 0.5)

# Run the Hatemi-J test
result <- hatemicoint(y, x)
print(result)
summary(result)
```

## Test Statistics

The package computes three test statistics:

| Test | Description |
|------|-------------|
| **ADF*** | Augmented Dickey-Fuller test on residuals |
| **Zt*** | Phillips-Perron Z_t test |
| **Za*** | Phillips-Perron Z_α test |

All statistics test the null hypothesis of **no cointegration** against the alternative of **cointegration with two structural breaks**.

## Options

```r
hatemicoint(y, x,
  maxlags = 8,              # Maximum lags for ADF
  lag_selection = "tstat",  # Lag selection: "tstat", "aic", or "sic"
  kernel = "iid",           # Kernel: "iid", "bartlett", or "qs"
  bwl = NULL,               # Bandwidth (auto-computed if NULL)
  trimming = 0.15           # Trimming for break search
)
```

## Critical Values

Critical values from Hatemi-J (2008, Table 1):

### ADF* and Zt* Tests

| k | 1% | 5% | 10% |
|---|----|----|-----|
| 1 | -6.503 | -6.015 | -5.653 |
| 2 | -6.928 | -6.458 | -6.224 |
| 3 | -7.833 | -7.352 | -7.118 |
| 4 | -8.353 | -7.903 | -7.705 |

### Za* Test

| k | 1% | 5% | 10% |
|---|----|----|-----|
| 1 | -90.794 | -76.003 | -52.232 |
| 2 | -99.458 | -83.644 | -76.806 |
| 3 | -118.577 | -104.860 | -97.749 |
| 4 | -140.135 | -123.870 | -116.169 |

**Decision rule:** Reject H₀ (no cointegration) if test statistic < critical value.

## References

Hatemi-J, A. (2008). Tests for cointegration with two unknown regime shifts with an application to financial market integration. *Empirical Economics*, 35, 497-505. [DOI: 10.1007/s00181-007-0175-9](https://doi.org/10.1007/s00181-007-0175-9)

Gregory, A.W., & Hansen, B.E. (1996). Residual-based tests for cointegration in models with regime shifts. *Journal of Econometrics*, 70(1), 99-126. [DOI: 10.1016/0304-4076(69)41685-7](https://doi.org/10.1016/0304-4076(69)41685-7)

## Author

## License
MIT
