
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

    #> Warning in var.mpl(copula, u): the covariance matrix of the parameter estimates
    #> is computed as if 'df.fixed = TRUE' with df = 8
    #> Warning in data.frame(Date = plot_dates, Return = returns_numeric[, i] * : row
    #> names were found from a short variable and have been discarded
    #> Warning in data.frame(Date = plot_dates, Return = returns_numeric[, i] * : row
    #> names were found from a short variable and have been discarded
    #> Warning in data.frame(Date = plot_dates, Return = returns_numeric[, i] * : row
    #> names were found from a short variable and have been discarded

## Solution 3: Use a Separate Chunk to Save and Display

<figure>
<img
src="/Users/Caleb/Documents/Working%20Directory/mcst/man/figures/result_stag_portfolio_distribution_comparison.png"
alt="“”" />
<figcaption aria-hidden="true">“”</figcaption>
</figure>

<figure>
<img
src="/Users/Caleb/Documents/Working%20Directory/mcst/man/figures/result_stag_exogenous_propagation.png"
alt="“”" />
<figcaption aria-hidden="true">“”</figcaption>
</figure>

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
