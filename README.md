
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package mcst: Overview

<!-- badges: start -->
<!-- badges: end -->

**mcst** is an R package that delivers end-to-end workflows for
portfolio risk and return simulations, risk analysis, and stress
testing. Key features include: - Modeling of endogenous and exogenous
shocks - Multiple return distributions (Normal, t, Skewed Generalized
T) - Option for Copula-based dependence structures (Gaussian and
t-copula) - Choice of forward-looking variance estimation for covariance
matrix (Sample, EWMA, GARCH)

## Installation

You can install the development version of mcst from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("calebong/mcst")
```

## Example

This is a basic example which shows you how to solve a common problem:

## Quick Start:

A sample portfolio of 3 equal-weighted stocks as stocks (endogenous
factors) and 2 market-wide (exogenous) factors - SPY and Crude Oil.

    #> # A tibble: 6 × 4
    #>   week          NVDA      JPM     XOM
    #>   <date>       <dbl>    <dbl>   <dbl>
    #> 1 2016-01-03 -0.101  -0.102   -0.0418
    #> 2 2016-01-10 -0.0850 -0.0319   0.0387
    #> 3 2016-01-17  0.0494 -0.00158 -0.0130
    #> 4 2016-01-24  0.0295  0.0448   0.0167
    #> 5 2016-01-31 -0.0976 -0.0294   0.0286
    #> 6 2016-02-07 -0.0265 -0.00450  0.0210
    #> # A tibble: 6 × 3
    #>   week         Crude      SPY
    #>   <date>       <dbl>    <dbl>
    #> 1 2016-01-03 -0.105  -0.0586 
    #> 2 2016-01-10 -0.113  -0.0214 
    #> 3 2016-01-17  0       0      
    #> 4 2016-01-24  0       0      
    #> 5 2016-01-31 -0.0812 -0.0298 
    #> 6 2016-02-07 -0.0469 -0.00702

Basic simulation with equal portfolio weights

    #>   ℹ️ historical_returns converted from tibble to data frame
    #>   ℹ️ Date column detected: 'week' (536 dates, 2016-01-03 to 2026-04-05)
    #>   ℹ️ Data frequency inferred: Weekly (52 periods/year)
    #> 
    #> ========================================================================
    #> DATA PREPROCESSING SUMMARY
    #> ========================================================================
    #>   Assets         : 3
    #>   Asset names    : NVDA, JPM, XOM
    #>   Observations   : 536 (from 536 raw rows)
    #>   Distribution   : rmvnorm
    #>   Copula         : NO
    #>   Cov estimation : SAMPLE
    #> 
    #>   Per-asset empirical summary:
    #>     NVDA: mean=0.0121  sd=0.0637  skew=0.362  kurt=4.498  min=-0.2005  max=0.3019
    #>     JPM: mean=0.0040  sd=0.0374  skew=-0.004  kurt=7.418  min=-0.1964  max=0.2226
    #>     XOM: mean=0.0030  sd=0.0384  skew=-0.111  kurt=6.673  min=-0.2007  max=0.1674
    #> 
    #>   📌 Reference asset defaulted to: NVDA
    #> 
    #> ========================================================================
    #> PORTFOLIO RISK SIMULATION (WEEKLY RETURNS)
    #> ========================================================================
    #>   Assets                : 3 — NVDA, JPM, XOM
    #>   Reference asset       : NVDA (index 1)
    #>   Observations          : 536
    #>   Distribution          : rmvnorm
    #>   Copula                : NONE
    #>   Covariance estimation : SAMPLE
    #>   Calibrate params      : NO
    #>   Propagation           : YES
    #> 
    #>   Weights: Equal weighted (1/N)
    #>   ✅ Long-only portfolio
    #> 
    #> ========================================================================
    #> WEIGHTS SUMMARY
    #> ========================================================================
    #>   Portfolios            : 1
    #>   Simulations each      : 10,000
    #>   Total simulations     : 10,000
    #>   Portfolio type        : Long-Only
    #> 
    #>   Portfolio weights:
    #>     NVDA                 : +0.3333 (+33.33%)  ███████
    #>     JPM                  : +0.3333 (+33.33%)  ███████
    #>     XOM                  : +0.3333 (+33.33%)  ███████
    #>     SUM                  : +1.0000
    #> 
    #> ========================================================================
    #> FITTING DISTRIBUTIONS
    #> ========================================================================
    #> 
    #> Path: NO COPULA
    #> 
    #> Distribution: Multivariate Normal
    #>   Means : 0.0121 (1.21%), 0.0040 (0.40%), 0.0030 (0.30%)
    #>   Vols  : 0.0637 (6.37%), 0.0374 (3.74%), 0.0384 (3.84%)
    #> 
    #> Historical correlation matrix:
    #>        NVDA    JPM    XOM
    #> NVDA 1.0000 0.3552 0.1281
    #> JPM  0.3552 1.0000 0.4864
    #> XOM  0.1281 0.4864 1.0000
    #> 
    #>   Final Σ condition number : 6.53
    #>   Final Σ min eigenvalue   : 6.83e-04
    #> 
    #> ========================================================================
    #> SIMULATING PORTFOLIO RETURNS
    #> ========================================================================
    #> 
    #>   Simulating 10,000 scenarios across 1 weight set(s)...
    #> 
    #> ========================================================================
    #> SIMULATION SUMMARY
    #> ========================================================================
    #>   Distribution   : rmvnorm
    #>   Weight sets    : 1
    #>   Sims per set   : 10,000
    #>   Valid pre-obs  : 10,000 / 10,000
    #> 
    #>   Sanity check:
    #>     Portfolio μ          Target:   0.6367%  Simulated:   0.6163%
    #>     Portfolio σ          Target:   3.4386%  Simulated:   3.4434%
    #>     ✅ Simulated vol consistent with target
    #> 
    #> ========================================================================
    #> PRE-STRESS PORTFOLIO STATISTICS (Weekly Returns)
    #> ========================================================================
    #> 
    #>   Sample size: 10,000
    #> 
    #>   RETURN (Weekly)
    #>   Expected Return (mean)         : +0.6163%
    #>   Expected Return (median)       : +0.6038%
    #> 
    #>   RETURN (Annualised)
    #>   Annual Return (geometric)      : +37.6447%
    #>   Annual Volatility              : 24.8309%
    #>   Sharpe Ratio (annual)          : 1.2907
    #> 
    #>   RISK (Weekly)
    #>   Volatility                     : 3.4434%
    #>   Skewness                       : 0.0178
    #>   Excess Kurtosis                : -0.0510
    #> 
    #>   VALUE-AT-RISK (Weekly, negative = loss)
    #>   90% VaR (10th pct)             : -3.8537%
    #>   95% VaR (5th pct)              : -5.0433%
    #>   99% VaR (1st pct)              : -7.2927%
    #> 
    #>   EXPECTED SHORTFALL (Weekly)
    #>   90% CVaR                       : -5.4002%
    #>   95% CVaR                       : -6.3897%
    #>   99% CVaR                       : -8.3681%
    #> 
    #>   TAIL DIAGNOSTICS
    #>   % Positive months              : 56.70%
    #>   % Negative months              : 43.30%
    #>   Worst single return            : -11.4769%
    #>   Best single return             : +14.4975%
    #> 
    #> ========================================================================
    #> GENERATING PORTFOLIO DIAGNOSTIC PLOTS
    #> ========================================================================
    #>   ✅ Plot 1: Asset distributions
    #>   ✅ Plot 2: Portfolio distribution
    #>   ✅ Plot 3: QQ plot
    #>   ✅ Plot 4: Historical correlation
    #>   ✅ Plot 6: Covariance matrix
    #>   ℹ️ Plot 7: GARCH time-series skipped (cov_estimation != 'garch')
    #> 
    #> ========================================================================
    #> GENERATING PAIRWISE COPULA DIAGNOSTIC PLOTS
    #>   Reference asset: NVDA
    #> ========================================================================
    #>   ℹ️ Pairwise plots skipped — copula = FALSE
    #>   ℹ️  Tail dependency plot skipped — copula = FALSE or fitting failed
    #> 
    #> ========================================================================
    #> COMPILING RESULTS
    #> ========================================================================
    #> 
    #> ========================================================================
    #> SIMULATION COMPLETE
    #> ========================================================================
    #>   Pre-stress samples  : 10,000
    #>   Diagnostic plots    : 5

<img src="man/figures/README-s1 portfolio distribution-1.png" width="100%" />

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure plot-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
