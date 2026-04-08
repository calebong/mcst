
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

The illustration uses a sample portfolio of 3 equal-weighted stocks
(endogenous factors) and 2 market-wide, exogenous factors - SPY and
Crude Oil - as shocks.

``` r
head(returns_data)
#> # A tibble: 6 × 4
#>   week          NVDA      JPM     XOM
#>   <date>       <dbl>    <dbl>   <dbl>
#> 1 2016-01-03 -0.101  -0.102   -0.0418
#> 2 2016-01-10 -0.0850 -0.0319   0.0387
#> 3 2016-01-17  0.0494 -0.00158 -0.0130
#> 4 2016-01-24  0.0295  0.0448   0.0167
#> 5 2016-01-31 -0.0976 -0.0294   0.0286
#> 6 2016-02-07 -0.0265 -0.00450  0.0210
head(exo_data)
#> # A tibble: 6 × 3
#>   week         Crude      SPY
#>   <date>       <dbl>    <dbl>
#> 1 2016-01-03 -0.105  -0.0586 
#> 2 2016-01-10 -0.113  -0.0214 
#> 3 2016-01-17  0       0      
#> 4 2016-01-24  0       0      
#> 5 2016-01-31 -0.0812 -0.0298 
#> 6 2016-02-07 -0.0469 -0.00702

tail(returns_data)
#> # A tibble: 6 × 4
#>   week           NVDA     JPM      XOM
#>   <date>        <dbl>   <dbl>    <dbl>
#> 1 2026-03-01  0.00356 -0.0360 -0.00846
#> 2 2026-03-08  0.0137  -0.0209  0.0325 
#> 3 2026-03-15 -0.0419   0.0110  0.0227 
#> 4 2026-03-22 -0.0300  -0.0130  0.0709 
#> 5 2026-03-29  0.0589   0.0416 -0.0602 
#> 6 2026-04-05  0.00400  0.0147  0.0200
tail(exo_data)
#> # A tibble: 6 × 3
#>   week          Crude      SPY
#>   <date>        <dbl>    <dbl>
#> 1 2026-03-01  0.356   -0.0198 
#> 2 2026-03-08  0.0859  -0.0150 
#> 3 2026-03-15 -0.00395 -0.0180 
#> 4 2026-03-22  0.0134  -0.0223 
#> 5 2026-03-29  0.119    0.0343 
#> 6 2026-04-05  0.0126   0.00517
```

Illustration of a stagflationary-like stress scenario: energy costs
surge while equity markets sell off simultaneously.

``` r
library(mcst)

exogenous_shock            <- list(SPY = -0.05, Crude = 0.20) # stagflationary-like shock; -5% to SPY, +20% to crude oil
exogenous_volatility_shock <- list(SPY = 1.1,     Crude = 1.5) # SPY and Crude volatility increases by a factor of 1.1x and 1.5x respectively 

# Basic simulation with equal weights
result_stag <- portfolio_risk_simulation(
  historical_returns          = returns_data,
  sim_returns_dist            = "rmvt",
  cov_estimation              = "garch",
  calibrate_params            = TRUE,
  calibration_method          = "full",
  copula                      = TRUE,
  copula_type                 = "t",
  copula_df                   = 4,
  n_sim_returns               = NULL,
  propagation                 = TRUE,
  exogenous_returns           = exo_data,
  exogenous_shock             = exogenous_shock,
  exogenous_volatility_shock  = exogenous_volatility_shock,
  df                          = 4,
  seed                        = 123
)
```

<figure>
<img src="man/figures/result_stag_portfolio_distribution_comparison.png"
alt="“”" />
<figcaption aria-hidden="true">“”</figcaption>
</figure>

<figure>
<img src="man/figures/result_stag_exogenous_propagation.png" alt="“”" />
<figcaption aria-hidden="true">“”</figcaption>
</figure>

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
