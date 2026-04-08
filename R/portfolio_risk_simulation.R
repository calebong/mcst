#' Check if a matrix is positive definite
#'
#' Internal helper function that check for positive definiteness of a matrix.
#' Used internally by \code{fix_covariance()}.
#' @param x A matrix.
#' @return A Boolean TRUE or FALSE.
#' @noRd
is.positive.definite <- function(x) {
  if(any(is.na(x)) || any(is.infinite(x))) return(FALSE)
  if(!is.matrix(x)) return(FALSE)
  if(nrow(x) != ncol(x)) return(FALSE)
  eigen_values <- tryCatch({
    eigen(x, symmetric = TRUE, only.values = TRUE)$values
  }, error = function(e) {
    return(rep(-1, nrow(x)))
  })
  return(min(eigen_values) > 1e-8)
}
#'
#' Fix covariance matrix to be positive definite
#'
#' This function checks if a covariance matrix is positive definite and
#' repairs it if necessary using the nearest positive definite matrix.
#' Used internally by \code{simulate_portfolio_returns()}.
#' @param cov_mat A numeric matrix. The covariance matrix to be checked and
#'   potentially repaired.
#' @return A positive definite covariance matrix.
#' @noRd
fix_covariance <- function(cov_mat) {
  if(!is.matrix(cov_mat)) cov_mat <- as.matrix(cov_mat)
  if(any(is.na(cov_mat)) || any(is.infinite(cov_mat))) {
    diag_vals <- diag(cov_mat)
    diag_vals[is.na(diag_vals) | is.infinite(diag_vals)] <- 0.01
    cov_mat <- diag(diag_vals)
  }
  if(!is.positive.definite(cov_mat)) {
    cov_mat <- as.matrix(nearPD(cov_mat, corr = FALSE)$mat)
  }
  return(cov_mat)
}
#'
#' Calibrate optimal degrees of freedom for multivariate t-distribution
#'
#' Internal helper that finds the optimal degrees of freedom parameter for a
#' t-distribution fitted to return data using maximum likelihood estimation.
#'
#' @param returns A numeric vector of asset returns
#' @param df_range Numeric vector. Range of degrees of freedom to search over.
#'   Default is seq(2.5, 10, by = 0.5).
#'
#' @return Integer. Optimal degrees of freedom value.
#'
#' @details
#' The function standardizes returns and evaluates the log-likelihood across
#' a grid of df values. The df that maximizes the log-likelihood is returned.
#'
#' @noRd
calibrate_df_mvt <- function(returns, df_range = seq(2.5, 10, by = 0.5)) {
  returns_std <- (returns - mean(returns)) / sd(returns)
  best_df     <- 4
  best_loglik <- -Inf
  for(df in df_range) {
    tryCatch({
      loglik <- sum(dt(returns_std, df = df, log = TRUE) - log(sd(returns)))
      if(loglik > best_loglik) {
        best_loglik <- loglik
        best_df     <- df
      }
    }, error = function(e) {})
  }
  return(best_df)
}
#' Calibrate parameters for Skewed Generalized T distribution
#'
#' Internal helper that estimates parameters for the Skewed Generalized T (SGT)
#' distribution using one of several calibration methods.
#'
#' @param returns A numeric vector of asset returns
#' @param method Character string specifying calibration method:
#'   \itemize{
#'     \item \code{"mle"}: Maximum likelihood estimation (full optimization)
#'     \item \code{"tail"}: Heuristic based on tail ratios (default)
#'     \item \code{"empirical"}: Based on empirical skewness
#'     \item \code{"default"}: Fallback to conservative parameters
#'   }
#'
#' @return A list containing:
#'   \item{lambda}{Skewness parameter (between -0.9 and 0.9)}
#'   \item{p}{Left tail parameter (degrees of freedom for left tail)}
#'   \item{q}{Right tail parameter (degrees of freedom for right tail)}
#'   \item{converged}{Logical indicating successful convergence (for MLE method)}
#'
#' @details
#' The SGT distribution is a flexible distribution that can capture skewness and
#' heavy tails. The \code{"tail"} method uses heuristics based on the ratio of
#' 99\% VaR to standard deviation to select parameters that match observed tail
#' behavior. This is much faster than MLE and often sufficient for risk modeling.
#'
#' @importFrom moments skewness
#' @importFrom sgt dsgt psgt rsgt qsgt
#' @noRd
calibrate_sgt <- function(returns, method = "tail") {

  if(method == "mle") {
    init_lambda <- skewness(returns) / 3
    init_lambda <- max(-0.8, min(0.8, init_lambda))

    neg_log_lik <- function(params, data) {
      mu     <- params[1]
      sigma  <- params[2]
      lambda <- params[3]
      p      <- params[4]
      q      <- params[5]
      if(sigma <= 0 || abs(lambda) >= 0.99 || p <= 2 || q <= 2) return(1e10)
      ll <- sum(dsgt(data, mu = mu, sigma = sigma,
                     lambda = lambda, p = p, q = q, log = TRUE))
      return(-ll)
    }

    tryCatch({
      opt_result <- optim(
        c(mean(returns), sd(returns), init_lambda, 4, 4),
        neg_log_lik, data = returns, method = "L-BFGS-B",
        lower = c(-0.5, 0.01, -0.9, 2.1, 2.1),
        upper = c(0.5,  1,     0.9, 10,  10),
        control = list(maxit = 500)
      )
      return(list(lambda = opt_result$par[3],
                  p = opt_result$par[4],
                  q = opt_result$par[5],
                  converged = TRUE))
    }, error = function(e) {
      # Fall back to tail method if MLE fails
      return(calibrate_sgt(returns, method = "tail"))
    })

  } else if(method == "tail") {
    var_95     <- quantile(returns, 0.05)
    var_99     <- quantile(returns, 0.01)
    tail_ratio <- abs(var_99) / sd(returns)

    # Heuristic mapping from tail ratio to SGT parameters
    if(tail_ratio > 3.5)      { lambda <- -0.7; p <- 2.5; q <- 6   }
    else if(tail_ratio > 3.0) { lambda <- -0.6; p <- 2.8; q <- 5.5 }
    else if(tail_ratio > 2.5) { lambda <- -0.5; p <- 3.2; q <- 5   }
    else if(tail_ratio > 2.0) { lambda <- -0.4; p <- 3.5; q <- 4.5 }
    else if(tail_ratio > 1.5) { lambda <- -0.2; p <- 4;   q <- 4   }
    else                      { lambda <-  0;   p <- 5;   q <- 5   }

    return(list(lambda = lambda, p = p, q = q, converged = TRUE))

  } else if(method == "empirical") {
    emp_skew <- skewness(returns)
    lambda   <- max(-0.8, min(0.8, emp_skew / 2.5))
    return(list(lambda = lambda, p = 4, q = 4, converged = TRUE))

  } else {
    # Default conservative parameters
    return(list(lambda = -0.3, p = 4, q = 4, converged = TRUE))
  }
}
#'
#' Calculate GJR-GARCH forward volatility
#'
#' Internal helper that estimates conditional volatility using a GJR-GARCH(1,1)
#' model for each asset in the return series. This captures leverage effects
#' where negative returns have a larger impact on volatility than positive returns.
#'
#' @param returns A numeric matrix or data frame of asset returns with
#'   dimensions n_obs × n_assets. Columns should be named with asset identifiers.
#' @param forecast_horizon Integer. Number of steps ahead to forecast volatility.
#'   Default is 1 (one-period-ahead forecast).
#'
#' @return A list containing four elements:
#'   \item{forecast_cov}{n_assets × n_assets diagonal covariance matrix with
#'     forecast variances on the diagonal}
#'   \item{cond_vol_mat}{n_obs × n_assets matrix of in-sample conditional
#'     volatilities (sigma_t) for each asset}
#'   \item{uncond_vol}{Named vector of unconditional (long-run) volatilities
#'     for each asset}
#'   \item{forecast_vol}{Named vector of 1-step-ahead volatility forecasts}
#'
#' @details
#' The function fits a separate GJR-GARCH(1,1) model to each asset using the
#' \code{rugarch} package with a skewed Student-t distribution. The GJR-GARCH
#' model incorporates leverage effects through an asymmetry parameter (gamma):
#'
#' \deqn{\sigma_t^2 = \omega + (\alpha + \gamma \cdot I_{t-1}) \varepsilon_{t-1}^2 + \beta \sigma_{t-1}^2}
#'
#' where \eqn{I_{t-1} = 1} if \eqn{\varepsilon_{t-1} < 0} and 0 otherwise.
#'
#' The function includes robust fallback logic for unconditional volatility
#' estimation: if \code{rugarch::uncvariance()} fails or returns implausible
#' values, it falls back to the empirical mean of in-sample conditional variances.
#'
#' @importFrom rugarch ugarchspec ugarchfit ugarchforecast sigma uncvariance
#' @importFrom stats sd
#' @noRd
calculate_garch_vol <- function(returns, forecast_horizon = 1) {
  if(!requireNamespace("rugarch", quietly = TRUE)) {
    stop("Package 'rugarch' is required for GARCH volatility estimation. ",
         "Please install it with: install.packages('rugarch')")
  }

  n_assets    <- ncol(returns)
  n_obs       <- nrow(returns)
  asset_names <- colnames(returns)

  garch_spec <- ugarchspec(
    variance.model    = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model        = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "sstd"
  )

  var_forecasts <- numeric(n_assets)
  cond_vol_mat  <- matrix(NA_real_, n_obs, n_assets,
                          dimnames = list(NULL, asset_names))
  uncond_vol    <- numeric(n_assets)
  forecast_vol  <- numeric(n_assets)
  names(uncond_vol) <- names(forecast_vol) <- asset_names

  for(i in seq_len(n_assets)) {
    fit   <- ugarchfit(spec = garch_spec, data = returns[, i], solver = "hybrid")
    fcast <- ugarchforecast(fit, n.ahead = forecast_horizon)

    fcast_sigma_1         <- sigma(fcast)[1, 1]
    var_forecasts[i]      <- fcast_sigma_1^2
    forecast_vol[i]       <- fcast_sigma_1

    # In-sample conditional volatilities (one per observation)
    cond_vol_mat[, i]     <- as.numeric(sigma(fit))

    # Unconditional (long-run) volatility.
    #
    # We use rugarch's own uncvariance() to extract unconditional variance.
    #
    # uncvariance(fit) accounts for the correct kappa = E[z² I(z<0)] for
    # the fitted innovation distribution (sstd here), not the hardcoded 0.5.
    #
    # If uncvariance() fails or returns a non-finite/implausible value, we
    # fall back to the empirical mean conditional variance — the time-average
    # of sigma²(t) from the in-sample fit.  This is always finite and
    # sensibly bounded between the minimum and maximum conditional variance.
    uncond_var_i <- tryCatch({
      uv <- rugarch::uncvariance(fit)
      if(is.finite(uv) && uv > 0 && sqrt(uv) < 5 * sd(returns[, i]))
        uv
      else
        mean(as.numeric(sigma(fit))^2)   # empirical fallback
    }, error = function(e) mean(as.numeric(sigma(fit))^2))

    uncond_vol[i] <- sqrt(pmax(uncond_var_i, 0))
  }

  forecast_cov <- matrix(0, n_assets, n_assets)
  diag(forecast_cov) <- var_forecasts
  rownames(forecast_cov) <- colnames(forecast_cov) <- asset_names
  forecast_cov <- fix_covariance(forecast_cov)

  return(list(
    forecast_cov = forecast_cov,
    cond_vol_mat = cond_vol_mat,
    uncond_vol   = uncond_vol,
    forecast_vol = forecast_vol
  ))
}
#' Calculate EWMA covariance matrix
#'
#' Internal helper that estimates covariance using the Exponentially Weighted
#' Moving Average (EWMA) approach, as used in RiskMetrics.
#'
#' @param returns A numeric matrix or data frame of asset returns with
#'   dimensions n_obs × n_assets. Columns should be named with asset identifiers.
#' @param lambda Numeric. Decay factor for the EWMA model. Default is 0.94,
#'   which is the RiskMetrics recommended value for daily data. Lower values
#'   give more weight to recent observations.
#'
#' @return A positive definite covariance matrix of dimensions n_assets × n_assets.
#'
#' @details
#' The EWMA covariance estimator updates the covariance matrix recursively:
#'
#' \deqn{\Sigma_t = \lambda \Sigma_{t-1} + (1-\lambda) \mathbf{r}_t \mathbf{r}_t^\top}
#'
#' where \eqn{\mathbf{r}_t} is the vector of returns at time t. The recursion
#' is initialized with the sample covariance matrix.
#'
#' The EWMA approach gives exponentially declining weights to past observations,
#' making it more responsive to recent market conditions than the equally-weighted
#' sample covariance. The half-life of the weights is approximately \eqn{\log(0.5)/\log(\lambda)}.
#'
#' The resulting covariance matrix is passed through \code{fix_covariance()} to
#' ensure positive definiteness.
#'
#' @seealso \code{\link{fix_covariance}} for positive definiteness repair
#' @noRd
calculate_ewma_covariance <- function(returns, lambda = 0.94) {
  n_assets <- ncol(returns)
  n_obs    <- nrow(returns)
  ewma_cov <- cov(returns)

  for(t in 2:n_obs) {
    r_t <- if(n_assets == 1) matrix(as.matrix(returns[t, ]), ncol = 1) else
      matrix(as.matrix(returns[t, ]), nrow = 1)
    outer_product <- t(r_t) %*% r_t
    ewma_cov <- lambda * ewma_cov + (1 - lambda) * outer_product
  }

  ewma_cov <- fix_covariance(ewma_cov)
  return(ewma_cov)
}
#' Monte Carlo Portfolio Risk Simulation with Stress Testing
#'
#' @description
#' An end-to-end function for multi-asset portfolio risk modelling.
#'
#' @param historical_returns A data frame, tibble, or matrix of asset returns.
#'   One column per asset, one row per observation period. An optional date
#'   column (class `Date`, `POSIXct`, or character coercible to `Date`) is
#'   used to infer return frequency automatically.
#' @param weights One of: `NULL` (equal weights, default); a named numeric
#'   vector summing to 1; a named list with asset names or 1-based column
#'   indices as keys; or a data frame / matrix where each row is a distinct
#'   weight allocation.
#' @param sim_returns_dist Character. Marginal return distribution.
#'   One of `"rmvnorm"` (Gaussian), `"rmvt"` (Student-t), or `"rsgt"`
#'   (Skewed Generalised-t). Default `"rmvnorm"`.
#' @param cov_estimation Character. Forward-looking covariance method.
#'   One of `"sample"`, `"ewma"`, or `"garch"`. Default `"sample"`.
#' @param cov_lambda Numeric in (0, 1). EWMA decay parameter. Only used
#'   when `cov_estimation = "ewma"`. Defaults to `0.97` if `NULL`.
#' @param calibrate_params Logical. If `TRUE`, fits distribution shape
#'   parameters (degrees of freedom for `rmvt`; lambda, p, q for `rsgt`)
#'   to each asset's historical returns. Default `FALSE`.
#' @param calibration_method Character. Moment-matching strategy when
#'   `calibrate_params = TRUE`. One of `"tail"` (match excess kurtosis only)
#'   or `"full"` (match all four moments equally). Default `"tail"`.
#' @param asset_dfs Numeric vector. Per-asset degrees of freedom for `rmvt`
#'   when `calibrate_params = FALSE`. Length must equal number of assets.
#' @param asset_lambdas Numeric vector. Per-asset skewness parameters for
#'   `rsgt`. Values in (-1, 1).
#' @param asset_p Numeric vector. Per-asset peakedness parameters for `rsgt`.
#' @param asset_q Numeric vector. Per-asset tail-decay parameters for `rsgt`.
#' @param copula Logical. If `TRUE`, fits a copula to the PIT residuals and
#'   uses it to model dependence separately from the marginals. Default `FALSE`.
#' @param copula_type Character. Copula family. One of `"gaussian"` or `"t"`.
#'   Only used when `copula = TRUE`. Default `"gaussian"`.
#' @param copula_df Numeric. Starting degrees of freedom for the t-copula.
#'   Only used when `copula = TRUE` and `copula_type = "t"`. Default `4`.
#' @param n_sim_returns Integer. Number of simulation draws per weight set.
#'   If `NULL` (default), uses `nrow(historical_returns)`.
#' @param propagation Logical. If `TRUE` (default), return shocks are
#'   propagated to unshocked assets via the Schur complement conditional
#'   mean formula. If `FALSE`, each shocked asset is modified in isolation.
#' @param asset_shock Named list or named numeric vector. Return shock
#'   (decimal) for each endogenous asset. Unspecified assets default to 0.
#'   Accepts names or 1-based column indices as keys.
#' @param volatility_shock Named list or named numeric vector. Volatility
#'   multiplier for each endogenous asset. Unspecified assets default to 1
#'   (no change). Accepts names or 1-based column indices as keys.
#' @param exogenous_returns A data frame of returns for non-portfolio
#'   instruments (e.g. SPY, Crude Oil). Must include a date column aligned
#'   with `historical_returns`. These instruments carry zero portfolio weight
#'   and are used as shock transmission channels only.
#' @param exogenous_shock Named list or named numeric vector. Return shock
#'   (decimal) for each exogenous instrument. Propagated to all portfolio
#'   asset means via `Sigma_pe %*% solve(Sigma_ee) %*% delta_e`.
#' @param exogenous_volatility_shock Named list or named numeric vector.
#'   Volatility multiplier for each exogenous instrument. Widens portfolio
#'   covariance via the Law of Total Variance addon.
#' @param df Numeric. Global degrees of freedom for `rmvt` when
#'   `calibrate_params = FALSE`. Default `4`.
#' @param seed Integer. Random seed for reproducibility. Default `123`.
#' @param reference_asset Character or integer. The asset shown on the
#'   y-axis of the pairwise latent-space diagnostic scatter plots.
#' @param n_pairwise_plots Integer. Maximum number of asset pairs to plot
#'   in the pairwise diagnostics panel. Default `12`.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{`portfolio_stats`}{14-row data frame: Metric, Pre_Stress, Post_Stress.}
#'   \item{`asset_draws`}{`n_sim x n_assets` data frame of pre-stress per-asset draws.}
#'   \item{`asset_draws_stress`}{`n_sim x n_assets` data frame of post-stress draws; `NULL` if not stressed.}
#'   \item{`asset_empirical_summary`}{Per-asset empirical mean, SD, skewness, kurtosis, min, max.}
#'   \item{`exogenous_asset_empirical_summary`}{Empirical summary for exogenous instruments; `NULL` if none.}
#'   \item{`exogenous_historical_correlation`}{Portfolio x exogenous Pearson correlation data frame.}
#'   \item{`ks_test`}{KS uniformity test results: Asset, KS_stat, p_value, Flag.}
#'   \item{`stress_impact_summary`}{7-row data frame: Metric, Pre_Stress, Post_Stress, Change.}
#'   \item{`stressed_covariance`}{n x n stressed covariance matrix after shock application.}
#'   \item{`stressed_means`}{Named vector of stressed conditional means per asset.}
#'   \item{`final_covariance`}{n x n final forward-looking covariance matrix.}
#'   \item{`copula_correlation`}{Fitted copula correlation matrix; `NULL` if `copula = FALSE`.}
#'   \item{`tail_dependence`}{Tail dependence matrix; `NULL` if Gaussian copula or `copula = FALSE`.}
#'   \item{`calibrated_params`}{Calibrated shape parameters per asset; `NULL` if not calibrated.}
#'   \item{`plots`}{List with `portfolio_diagnostics` containing named ggplot/patchwork objects.}
#'   \item{`stressed_plots`}{Stress-scenario visualisations; `NULL` if no shocks specified.}
#' }
#'
#' @examples
#' \dontrun{
#' # Minimal call — Gaussian, equal weights, no stress
#' result <- portfolio_risk_simulation(returns_df)
#'
#' # Fat-tailed with t-copula and EWMA covariance
#' result <- portfolio_risk_simulation(
#'   returns_df,
#'   sim_returns_dist   = "rmvt",
#'   cov_estimation     = "ewma",
#'   cov_lambda         = 0.97,
#'   calibrate_params   = TRUE,
#'   copula             = TRUE,
#'   copula_type        = "t"
#' )
#'
#' # Stress: SPY falls 5%, vol doubles
#' result <- portfolio_risk_simulation(
#'   returns_df,
#'   exogenous_returns          = exo_df,
#'   exogenous_shock            = list(SPY = -0.05),
#'   exogenous_volatility_shock = list(SPY = 2)
#' )
#' }
#'
#' @references
#' Glosten, L.R., Jagannathan, R., & Runkle, D.E. (1993).
#' On the relation between the expected value and the volatility of the
#' nominal excess return on stocks. \emph{Journal of Finance}, 48(5), 1779--1801.
#'
#' Ledoit, O., & Wolf, M. (2004). A well-conditioned estimator for
#' large-dimensional covariance matrices.
#' \emph{Journal of Multivariate Analysis}, 88(2), 365--411.
#'
#' Nelsen, R.B. (2006). \emph{An Introduction to Copulas} (2nd ed.). Springer.
#'
#' Sklar, A. (1959). Fonctions de repartition a n dimensions et leurs marges.
#' \emph{Publications de l'Institut Statistique de l'Universite de Paris}, 8, 229--231.
#'
#' Theodossiou, P. (1998). Financial data and the skewed generalised t distribution.
#' \emph{Management Science}, 44(12), 1650--1661.
#'
#' @seealso
#' \itemize{
#'   \item \code{vignette("portfolio_risk_simulation_vignette", package = "mcst")}
#'   \item \code{vignette("portfolio_risk_simulation_scenarios", package = "mcst")}
#'   \item \code{\link{portfolio_risk_simulation}} for main function details
#' }
#'
#' @importFrom dplyr select filter mutate group_by summarise ungroup arrange
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_bar
#'   geom_line geom_point geom_vline geom_hline labs theme theme_minimal
#'   element_text scale_fill_manual scale_colour_manual coord_cartesian
#'   after_stat position_dodge
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom Matrix nearPD
#' @importFrom MASS ginv mvrnorm
#' @importFrom copula normalCopula tCopula fitCopula dCopula pobs getTheta
#' @importFrom rugarch ugarchspec ugarchfit ugarchforecast sigma uncvariance
#' @importFrom mvtnorm rmvnorm rmvt
#' @importFrom moments skewness kurtosis
#'
#' @export
portfolio_risk_simulation <- function(
    historical_returns,
    weights           = NULL,
    sim_returns_dist  = "rmvnorm",
    cov_estimation    = "sample",
    cov_lambda        = 0.97,
    calibrate_params  = FALSE,
    calibration_method = "tail",
    asset_dfs         = NULL,
    asset_lambdas     = NULL,
    asset_p           = NULL,
    asset_q           = NULL,
    copula            = FALSE,
    copula_type       = "gaussian",
    copula_df         = 4,
    n_sim_returns     = NULL,
    propagation       = TRUE,
    asset_shock                = NULL,
    volatility_shock           = NULL,
    exogenous_returns          = NULL,
    exogenous_shock            = NULL,
    exogenous_volatility_shock = NULL,
    df                         = 4,
    seed                       = 123,
    reference_asset            = NULL,
    n_pairwise_plots           = 12
) {

  set.seed(seed)

  # ===========================================================================
  # 0. ARGUMENT NORMALISATION — list → named numeric vector / matrix
  # ===========================================================================
  #
  # Five arguments accept flexible list input in addition to their existing
  # numeric-vector / data-frame forms:
  #
  #   weights, asset_shock, volatility_shock  → resolved against portfolio asset names
  #   exogenous_shock, exogenous_volatility_shock → resolved against exogenous asset names
  #
  # Accepted list formats:
  #
  #   Format 1 — names  : list(AXP = 0.3, TSM = 0.3, NVDA = 0.4)
  #     Names must match asset names from the returns / exogenous data frame.
  #     Unspecified assets default to default_val (0 for shocks, 1 for vol multipliers).
  #     Order is inferred automatically.
  #
  #   Format 2 — indices : list(`1` = 0.3, `2` = 0.3, `3` = 0.4)
  #     Numeric-string keys (after date column removal) reference the 1-based
  #     column position in the returns / exogenous data frame.
  #     Unspecified positions default to default_val.
  #
  # resolve_arg() is called BEFORE any validation, so the downstream code
  # still receives a plain numeric vector and is completely unaffected.
  # ──────────────────────────────────────────────────────────────────────────

  resolve_arg <- function(arg, all_names, arg_name, default_val, full_length = TRUE) {
    # Returns a named numeric vector of length(all_names), or NULL if arg is NULL.
    # full_length = TRUE  → always return all n assets (pad missing with default_val)
    # full_length = FALSE → return only the entries that were specified (for exog shocks
    #                       where partial specification is fine)

    if(is.null(arg)) return(NULL)

    # Already a plain numeric vector or matrix — pass through untouched
    if(!is.list(arg)) return(arg)

    n <- length(all_names)
    keys <- names(arg)

    if(is.null(keys))
      stop(sprintf(
        "`%s`: list must have named elements.\n  Use either asset names (e.g. list(AXP=0.3)) or\n  1-based column indices (e.g. list(`1`=0.3)).",
        arg_name))

    vals <- unlist(arg, use.names = FALSE)
    if(!is.numeric(vals))
      stop(sprintf("`%s`: list values must be numeric.", arg_name))

    # Detect format: are all keys integer-like strings?
    key_nums <- suppressWarnings(as.integer(keys))
    is_index_format <- !any(is.na(key_nums))

    if(is_index_format) {
      # ── Format 2: integer indices ────────────────────────────────────────────
      bad_idx <- key_nums[key_nums < 1 | key_nums > n]
      if(length(bad_idx) > 0)
        stop(sprintf(
          "`%s`: index %s out of range [1, %d] (available assets: %s).",
          arg_name, paste(bad_idx, collapse=", "), n, paste(all_names, collapse=", ")))

      dup_idx <- key_nums[duplicated(key_nums)]
      if(length(dup_idx) > 0)
        stop(sprintf("`%s`: duplicate indices: %s.", arg_name, paste(unique(dup_idx), collapse=", ")))

      if(full_length) {
        out <- setNames(rep(default_val, n), all_names)
        out[key_nums] <- vals
      } else {
        # Partial — only return specified assets (named)
        out <- setNames(vals, all_names[key_nums])
      }

    } else {
      # ── Format 1: name keys ──────────────────────────────────────────────────
      bad_names <- setdiff(keys, all_names)
      if(length(bad_names) > 0)
        stop(sprintf(
          "`%s`: name(s) not found: %s.\n  Available: %s.",
          arg_name, paste(bad_names, collapse=", "), paste(all_names, collapse=", ")))

      dup_names <- keys[duplicated(keys)]
      if(length(dup_names) > 0)
        stop(sprintf("`%s`: duplicate names: %s.", arg_name, paste(unique(dup_names), collapse=", ")))

      if(full_length) {
        out <- setNames(rep(default_val, n), all_names)
        out[keys] <- vals
      } else {
        out <- setNames(vals, keys)
      }
    }

    cat(sprintf(
      "  ℹ️ `%s` resolved from list (%s format): %s\n",
      arg_name,
      if(is_index_format) "index" else "name",
      paste(sprintf("%s=%g", names(out), out), collapse=", ")))

    out
  }

  # Portfolio-argument names are resolved after preprocessing (when asset_names
  # is known). We store the raw inputs here and resolve them in Step 2.
  .weights_raw    <- weights
  .shock_raw      <- asset_shock
  .vol_shock_raw  <- volatility_shock
  .exog_shock_raw <- exogenous_shock
  .exog_vs_raw    <- exogenous_volatility_shock

  # ── Argument normalisation: cov_lambda ───────────────────────────────────────
  # cov_lambda is only consumed by EWMA. If the user passes NULL (or an invalid
  # value) the EWMA loop silently produces a 0×0 matrix and crashes.
  # Default: 0.97 for weekly/monthly (RiskMetrics), 0.94 for daily.
  if(cov_estimation == "ewma") {
    if(is.null(cov_lambda) || !is.numeric(cov_lambda) ||
       length(cov_lambda) != 1 || !is.finite(cov_lambda) ||
       cov_lambda <= 0 || cov_lambda >= 1) {
      cov_lambda <- 0.97
      cat(sprintf("  ℹ️ cov_lambda NULL or invalid — defaulting to %.2f for EWMA\n", cov_lambda))
    }
  }

  # ── Argument normalisation: calibration_method ───────────────────────────────
  # Only "tail" and "full" are valid. Anything else (e.g. "mle") is unsupported;
  # default to "tail" with a warning so the function continues cleanly.
  if(calibrate_params && !is.null(calibration_method)) {
    valid_methods <- c("tail", "full")
    if(!calibration_method %in% valid_methods) {
      warning(sprintf(
        paste0("calibration_method = '%s' is not recognised (valid: %s). ",
               "Defaulting to 'tail'."),
        calibration_method, paste(valid_methods, collapse = ", ")))
      calibration_method <- "tail"
    }
  } else if(calibrate_params && is.null(calibration_method)) {
    calibration_method <- "tail"
    cat("  ℹ️ calibration_method = NULL — defaulting to 'tail'\n")
  }

  # ===========================================================================
  # 1. DATA PREPROCESSING
  # ===========================================================================

  # 1a. Input type validation — accept tibble, data.frame, or matrix
  if(inherits(historical_returns, "tbl") || inherits(historical_returns, "tbl_df")) {
    historical_returns <- as.data.frame(historical_returns)
    cat("  ℹ️ historical_returns converted from tibble to data frame\n")
  }
  if(is.matrix(historical_returns)) {
    historical_returns <- as.data.frame(historical_returns)
    cat("  ℹ️ historical_returns converted from matrix to data frame\n")
  }
  if(!is.data.frame(historical_returns))
    stop("historical_returns must be a data frame, tibble, or matrix")

  # 1b. Detect and extract date column
  # The date column is the first column of class Date, POSIXct/lt, or character
  # that can be coerced to Date. All remaining columns are treated as asset returns.
  detect_date_col <- function(df) {
    for(nm in colnames(df)) {
      col <- df[[nm]]
      if(inherits(col, c("Date", "POSIXct", "POSIXlt"))) return(nm)
      if(is.character(col) || is.factor(col)) {
        trial <- suppressWarnings(as.Date(as.character(col)))
        if(sum(!is.na(trial)) / length(trial) > 0.8) return(nm)
      }
    }
    return(NULL)
  }

  date_col_name <- detect_date_col(historical_returns)
  dates_raw     <- NULL

  if(!is.null(date_col_name)) {
    col <- historical_returns[[date_col_name]]
    dates_raw <- if(inherits(col, "Date")) col
                 else if(inherits(col, c("POSIXct","POSIXlt"))) as.Date(col)
                 else suppressWarnings(as.Date(as.character(col)))
    dates_raw <- dates_raw[!is.na(dates_raw)]
    cat(sprintf("  ℹ️ Date column detected: '%s' (%d dates, %s to %s)\n",
                date_col_name, length(dates_raw),
                format(min(dates_raw)), format(max(dates_raw))))
  } else {
    cat("  ℹ️ No date column detected — frequency defaulting to 'monthly'\n")
  }

  # 1b2. Infer data frequency from median gap between consecutive dates
  infer_frequency <- function(dates) {
    if(is.null(dates) || length(dates) < 2) return(list(freq = "monthly", label = "Monthly", periods_per_year = 12L))
    gaps <- as.numeric(diff(sort(dates)))
    med_gap <- median(gaps, na.rm = TRUE)
    if(med_gap <= 5)        return(list(freq = "daily",   label = "Daily",   periods_per_year = 252L))
    else if(med_gap <= 10)  return(list(freq = "weekly",  label = "Weekly",  periods_per_year = 52L))
    else                    return(list(freq = "monthly", label = "Monthly", periods_per_year = 12L))
  }

  freq_info        <- infer_frequency(dates_raw)
  data_freq        <- freq_info$freq          # "daily" | "weekly" | "monthly"
  data_freq_label  <- freq_info$label         # "Daily" | "Weekly" | "Monthly"
  periods_per_year <- freq_info$periods_per_year  # 252 | 52 | 12

  cat(sprintf("  ℹ️ Data frequency inferred: %s (%d periods/year)\n",
              data_freq_label, periods_per_year))

  # 1b3. Extract numeric return columns (excluding date column)
  returns_numeric <- historical_returns %>% select_if(is.numeric)
  if(!is.null(date_col_name) && date_col_name %in% colnames(returns_numeric))
    returns_numeric <- returns_numeric[, colnames(returns_numeric) != date_col_name, drop = FALSE]
  if(ncol(returns_numeric) == 0)
    stop("No numeric return columns found in historical_returns.")

  # 1c. Handle missing values
  n_rows_raw <- nrow(returns_numeric)
  n_missing  <- sum(!complete.cases(returns_numeric))

  if(n_missing > 0) {
    miss_pct <- round(100 * n_missing / n_rows_raw, 2)
    cat(sprintf("\n  ⚠️ Missing values: %d rows (%.2f%%) contain NA\n",
                n_missing, miss_pct))
    col_miss <- colSums(is.na(returns_numeric))
    col_miss <- col_miss[col_miss > 0]
    for(nm in names(col_miss))
      cat(sprintf("    %s: %d NAs (%.1f%%)\n",
                  nm, col_miss[nm], 100 * col_miss[nm] / n_rows_raw))
    if(miss_pct > 20)
      warning(sprintf("%.1f%% of rows dropped due to NA.", miss_pct))
  }

  returns_numeric <- na.omit(returns_numeric) %>% data.frame()

  # 1d. Constant column check
  col_sds       <- apply(returns_numeric, 2, sd)
  constant_cols <- names(col_sds[col_sds < 1e-10])
  if(length(constant_cols) > 0)
    stop(sprintf("Near-zero variance assets: %s",
                 paste(constant_cols, collapse = ", ")))

  # 1e. Infinite value check
  if(sum(!is.finite(as.matrix(returns_numeric))) > 0)
    stop("Infinite values found after NA removal.")

  # 1f. Core dimensions
  asset_names <- colnames(returns_numeric)
  n_assets    <- ncol(returns_numeric)
  n_obs       <- nrow(returns_numeric)

  if(n_obs < 30)
    warning(sprintf("Only %d observations — covariance estimates may be unreliable.", n_obs))
  if(n_obs < n_assets * 5)
    warning(sprintf("n_obs (%d) < 5 × n_assets (%d) — covariance may be ill-conditioned.",
                    n_obs, n_assets * 5))

  cat("\n========================================================================\n")
  cat("DATA PREPROCESSING SUMMARY\n")
  cat("========================================================================\n")
  cat(sprintf("  Assets         : %d\n", n_assets))
  cat(sprintf("  Asset names    : %s\n", paste(asset_names, collapse = ", ")))
  cat(sprintf("  Observations   : %d (from %d raw rows)\n", n_obs, n_rows_raw))
  cat(sprintf("  Distribution   : %s\n", sim_returns_dist))
  cat(sprintf("  Copula         : %s%s\n",
              ifelse(copula, "YES", "NO"),
              if(copula) sprintf(" (%s)", toupper(copula_type)) else ""))
  cat(sprintf("  Cov estimation : %s\n", toupper(cov_estimation)))

  cat("\n  Per-asset empirical summary:\n")
  asset_empirical_summary <- data.frame(
    Asset    = asset_names,
    Mean     = NA_real_, SD       = NA_real_, Skewness = NA_real_,
    Kurtosis = NA_real_, Min      = NA_real_, Max      = NA_real_,
    stringsAsFactors = FALSE
  )
  for(i in seq_len(n_assets)) {
    x <- returns_numeric[, i]
    asset_empirical_summary[i, "Mean"]     <- mean(x)
    asset_empirical_summary[i, "SD"]       <- sd(x)
    asset_empirical_summary[i, "Skewness"] <- skewness(x)
    asset_empirical_summary[i, "Kurtosis"] <- kurtosis(x)
    asset_empirical_summary[i, "Min"]      <- min(x)
    asset_empirical_summary[i, "Max"]      <- max(x)
    cat(sprintf("    %s: mean=%.4f  sd=%.4f  skew=%.3f  kurt=%.3f  min=%.4f  max=%.4f\n",
                asset_names[i], mean(x), sd(x), skewness(x), kurtosis(x),
                min(x), max(x)))
  }


  # ===========================================================================
  # 1.7 EXOGENOUS INSTRUMENTS PREPROCESSING
  #
  # exogenous_returns holds returns for assets NOT in the portfolio (e.g. SPY).
  # These instruments propagate shocks onto portfolio assets via:
  #
  #   μ_p | X_e = x_e  =  μ_p + Σ_pe Σ_ee⁻¹ (x_e − μ_e)
  #   Σ_p | X_e         =  Σ_pp − Σ_pe Σ_ee⁻¹ Σ_ep
  #
  # Exogenous assets carry zero portfolio weight and are NEVER simulated.
  # ===========================================================================

  has_exogenous <- !is.null(exogenous_returns)
  exog_numeric  <- NULL
  exog_names    <- character(0)
  n_exog        <- 0L
  exog_means    <- NULL
  Sigma_pe_hist <- NULL   # Cov(X_port, X_exog)
  Sigma_ee_hist <- NULL   # Cov(X_exog, X_exog)

  if(has_exogenous) {

    # 1.7a. Coerce to data frame
    if(inherits(exogenous_returns, c("tbl", "tbl_df")))
      exogenous_returns <- as.data.frame(exogenous_returns)
    if(is.matrix(exogenous_returns))
      exogenous_returns <- as.data.frame(exogenous_returns)
    if(!is.data.frame(exogenous_returns))
      stop("exogenous_returns must be a data frame, tibble, or matrix")

    # 1.7b. Strip date column from exogenous df
    exog_date_col <- detect_date_col(exogenous_returns)
    exog_num      <- exogenous_returns %>% select_if(is.numeric)
    if(!is.null(exog_date_col) && exog_date_col %in% colnames(exog_num))
      exog_num <- exog_num[, colnames(exog_num) != exog_date_col, drop = FALSE]
    exog_num <- na.omit(exog_num) %>% data.frame()

    # 1.7c. No overlap with portfolio assets
    overlap <- intersect(colnames(exog_num), asset_names)
    if(length(overlap) > 0)
      stop(sprintf(
        paste0("exogenous_returns shares assets with historical_returns: %s\n",
               "  Each asset must appear in exactly one of the two data frames."),
        paste(overlap, collapse = ", ")))

    # 1.7d. Row count must match
    if(nrow(exog_num) != n_obs)
      stop(sprintf(
        paste0("exogenous_returns has %d rows but portfolio returns has %d.\n",
               "  Align both to the same dates before calling."),
        nrow(exog_num), n_obs))

    # 1.7e. Validity checks
    exog_sds <- apply(exog_num, 2, sd)
    const_e  <- names(exog_sds[exog_sds < 1e-10])
    if(length(const_e) > 0)
      stop(sprintf("Near-zero variance exogenous assets: %s",
                   paste(const_e, collapse = ", ")))
    if(sum(!is.finite(as.matrix(exog_num))) > 0)
      stop("Infinite values found in exogenous_returns after NA removal.")

    exog_numeric <- exog_num
    exog_names   <- colnames(exog_numeric)
    n_exog       <- ncol(exog_numeric)
    exog_means   <- colMeans(exog_numeric)

    cat("\n========================================================================\n")
    cat("EXOGENOUS INSTRUMENTS\n")
    cat("========================================================================\n")
    cat(sprintf("  Count        : %d\n", n_exog))
    cat(sprintf("  Names        : %s\n", paste(exog_names, collapse = ", ")))
    cat(sprintf("  Rows matched : %d (aligned with portfolio observations)\n", n_obs))
    cat("\n  Exogenous empirical summary:\n")
    exogenous_asset_empirical_summary <- data.frame(
      Asset    = exog_names,
      Mean     = NA_real_, SD       = NA_real_,
      Skewness = NA_real_, Kurtosis = NA_real_,
      stringsAsFactors = FALSE
    )
    for(i in seq_len(n_exog)) {
      x <- exog_numeric[, i]
      exogenous_asset_empirical_summary[i, "Mean"]     <- mean(x)
      exogenous_asset_empirical_summary[i, "SD"]       <- sd(x)
      exogenous_asset_empirical_summary[i, "Skewness"] <- skewness(x)
      exogenous_asset_empirical_summary[i, "Kurtosis"] <- kurtosis(x)
      cat(sprintf("    %s: mean=%.4f  sd=%.4f  skew=%.3f  kurt=%.3f\n",
                  exog_names[i], mean(x), sd(x), skewness(x), kurtosis(x)))
    }

    # 1.7f. Joint historical covariance: Σ_aug = cov([X_port, X_exog])
    #   Sample cov is always used for the cross-block Σ_pe regardless of
    #   cov_estimation, since GARCH/EWMA apply only to marginal volatilities
    #   of portfolio assets. The diagonal Σ_pp is later replaced by the
    #   forward-looking covariance from the chosen cov_estimation method.
    joint_returns  <- cbind(returns_numeric, exog_numeric)
    joint_cov_full <- cov(joint_returns)

    p_idx <- seq_len(n_assets)
    e_idx <- seq(n_assets + 1L, n_assets + n_exog)

    Sigma_pe_hist <- joint_cov_full[p_idx, e_idx, drop = FALSE]
    Sigma_ee_hist <- joint_cov_full[e_idx, e_idx, drop = FALSE]
    rownames(Sigma_pe_hist) <- asset_names
    colnames(Sigma_pe_hist) <- exog_names
    rownames(Sigma_ee_hist) <- colnames(Sigma_ee_hist) <- exog_names

    if(!is.positive.definite(Sigma_ee_hist)) {
      cat("  ⚠️ Σ_ee (exogenous cov) not PD — applying nearPD\n")
      Sigma_ee_hist <- as.matrix(nearPD(Sigma_ee_hist, corr = FALSE)$mat)
      rownames(Sigma_ee_hist) <- colnames(Sigma_ee_hist) <- exog_names
    }

    # Report portfolio × exogenous correlations
    joint_cor_full <- cor(joint_returns)
    pe_cor         <- joint_cor_full[p_idx, e_idx, drop = FALSE]
    cat("\n  Portfolio x Exogenous historical correlations:\n")
    cat(sprintf("  %-20s", ""))
    for(en in exog_names) cat(sprintf("  %8s", en))
    cat("\n")
    for(pi in seq_len(n_assets)) {
      cat(sprintf("  %-20s", asset_names[pi]))
      for(ei in seq_len(n_exog)) cat(sprintf("  %8.3f", pe_cor[pi, ei]))
      cat("\n")
    }
    # Store as data frame: Asset | exog1 | exog2 | ...
    exogenous_historical_correlation <- as.data.frame(pe_cor, stringsAsFactors = FALSE)
    exogenous_historical_correlation <- cbind(
      data.frame(Asset = asset_names, stringsAsFactors = FALSE),
      exogenous_historical_correlation
    )
    rownames(exogenous_historical_correlation) <- NULL
  } else {
    exogenous_asset_empirical_summary    <- NULL
    exogenous_historical_correlation     <- NULL
  }

  # ===========================================================================
  # 1.5 CALIBRATE DISTRIBUTION PARAMETERS
  # ===========================================================================

  calibrated_dfs     <- NULL
  calibrated_lambdas <- NULL
  calibrated_ps      <- NULL
  calibrated_qs      <- NULL

  if(calibrate_params && sim_returns_dist == "rmvnorm") {
    cat("\n  ℹ️ calibrate_params ignored: nothing to calibrate for rmvnorm\n")
    calibrate_params <- FALSE
  }

  if(calibrate_params) {

    cat("\n========================================================================\n")
    cat("CALIBRATING DISTRIBUTION PARAMETERS\n")
    cat("========================================================================\n")
    cat(sprintf("  Distribution      : %s\n", sim_returns_dist))
    cat(sprintf("  Calibration method: %s\n", toupper(calibration_method)))

    if(sim_returns_dist == "rmvt") {

      calibrated_dfs        <- numeric(n_assets)
      names(calibrated_dfs) <- asset_names
      cat("\nCalibrating degrees of freedom:\n")

      for(i in seq_len(n_assets)) {
        calibrated_dfs[i] <- tryCatch({
          df_i <- calibrate_df_mvt(returns_numeric[, i])
          if(!is.finite(df_i) || df_i <= 2) { cat(sprintf("  ⚠️ %s: df invalid — using 4\n", asset_names[i])); 4 }
          else df_i
        }, error = function(e) {
          cat(sprintf("  ⚠️ %s: calibration failed — using 4\n", asset_names[i])); 4
        })
        cat(sprintf("  %s: df = %.4f\n", asset_names[i], calibrated_dfs[i]))
      }

      df_joint <- mean(calibrated_dfs)
      df       <- df_joint
      if(df_joint <= 2) stop(sprintf("Mean calibrated df = %.2f <= 2", df_joint))
      cat(sprintf("\n  Joint df (mean): %.4f\n", df_joint))

    } else if(sim_returns_dist == "rsgt") {

      calibrated_lambdas        <- numeric(n_assets)
      calibrated_ps             <- numeric(n_assets)
      calibrated_qs             <- numeric(n_assets)
      names(calibrated_lambdas) <- asset_names
      names(calibrated_ps)      <- asset_names
      names(calibrated_qs)      <- asset_names

      cat("\nCalibrating SGT parameters:\n")

      for(i in seq_len(n_assets)) {
        result <- tryCatch(
          calibrate_sgt(returns_numeric[, i], method = calibration_method),
          error = function(e) {
            cat(sprintf("  ⚠️ %s: failed — using defaults\n", asset_names[i]))
            list(lambda = -0.3, p = 4, q = 4, converged = FALSE)
          }
        )

        lam_ok <- is.finite(result$lambda) && abs(result$lambda) < 1
        p_ok   <- is.finite(result$p)      && result$p > 0
        q_ok   <- is.finite(result$q)      && result$q > 1 / result$p

        if(!lam_ok || !p_ok || !q_ok) {
          cat(sprintf("  ⚠️ %s: invalid params — using defaults\n", asset_names[i]))
          result <- list(lambda = -0.3, p = 4, q = 4)
        }

        calibrated_lambdas[i] <- result$lambda
        calibrated_ps[i]      <- result$p
        calibrated_qs[i]      <- result$q

        if(result$q <= 2 / result$p)
          warning(sprintf("Asset %s: q <= 2/p — SGT variance may be infinite", asset_names[i]))

        converged_str <- if(!is.null(result$converged) && !result$converged) " [not converged]" else ""
        cat(sprintf("  %s: λ=%.4f, p=%.4f, q=%.4f%s\n",
                    asset_names[i], calibrated_lambdas[i],
                    calibrated_ps[i], calibrated_qs[i], converged_str))
      }
    }

    cat("\n========================================================================\n")
    cat("CALIBRATION SUMMARY\n")
    cat("========================================================================\n")
    if(sim_returns_dist == "rmvt") {
      print(data.frame(Asset = asset_names, Optimal_DF = round(calibrated_dfs, 4),
                       Valid = calibrated_dfs > 2), row.names = FALSE)
      cat(sprintf("\n  Joint df for simulation: %.4f\n", df_joint))
    } else if(sim_returns_dist == "rsgt") {
      print(data.frame(Asset      = asset_names,
                       Lambda     = round(calibrated_lambdas, 4),
                       p          = round(calibrated_ps, 4),
                       q          = round(calibrated_qs, 4),
                       Finite_Var = calibrated_qs > 2 / calibrated_ps),
            row.names = FALSE)
    }
  }

  # Resolve final distribution parameters
  if(sim_returns_dist == "rmvt") {
    if(calibrate_params && !is.null(calibrated_dfs)) {
      cat(sprintf("\n  rmvt: using calibrated df_joint = %.4f\n", df_joint))
    } else if(!is.null(asset_dfs)) {
      if(length(asset_dfs) != n_assets)
        stop(sprintf("asset_dfs length %d != n_assets %d", length(asset_dfs), n_assets))
      if(any(asset_dfs <= 2)) stop("All asset_dfs must be > 2")
      df_joint <- mean(asset_dfs)
      df       <- df_joint
      cat(sprintf("\n  rmvt: user-supplied asset_dfs, df_joint = %.4f\n", df_joint))
    } else {
      if(df <= 2) stop(sprintf("df = %.2f must be > 2", df))
      df_joint <- df
      cat(sprintf("\n  rmvt: fixed df = %.4f\n", df_joint))
    }
  }

  if(sim_returns_dist == "rsgt") {
    if(calibrate_params && !is.null(calibrated_lambdas)) {
      asset_lambdas <- calibrated_lambdas
      asset_p       <- calibrated_ps
      asset_q       <- calibrated_qs
      cat("\n  rsgt: using calibrated SGT parameters\n")
    } else if(!is.null(asset_lambdas)) {
      if(length(asset_lambdas) != n_assets)
        stop(sprintf("asset_lambdas length %d != n_assets %d", length(asset_lambdas), n_assets))
      if(any(abs(asset_lambdas) >= 1)) stop("|lambda| must be < 1")
      asset_p <- if(!is.null(asset_p)) asset_p else rep(4, n_assets)
      asset_q <- if(!is.null(asset_q)) asset_q else rep(4, n_assets)
      if(length(asset_p) != n_assets || length(asset_q) != n_assets)
        stop("asset_p and asset_q must have length == n_assets")
      if(any(asset_p <= 0)) stop("All asset_p must be > 0")
      if(any(asset_q <= 1 / asset_p)) stop("All asset_q must be > 1/p")
      cat("\n  rsgt: using user-supplied SGT parameters\n")
    } else {
      asset_lambdas <- rep(-0.3, n_assets)
      asset_p       <- rep(4,    n_assets)
      asset_q       <- rep(4,    n_assets)
      cat("\n  rsgt: using conservative defaults (λ=-0.3, p=4, q=4)\n")
    }
    names(asset_lambdas) <- asset_names
    names(asset_p)       <- asset_names
    names(asset_q)       <- asset_names
  }

  # ===========================================================================
  # REFERENCE ASSET SELECTION
  # ===========================================================================

  # Resolve list arguments now (asset_names is known). This is the primary
  # resolution point — the Step 2a-pre block re-runs resolve_arg idempotently.
  asset_shock      <- resolve_arg(.shock_raw,     asset_names, "asset_shock",
                                  default_val = 0, full_length = TRUE)
  volatility_shock <- resolve_arg(.vol_shock_raw, asset_names, "volatility_shock",
                                  default_val = 1, full_length = TRUE)
  weights          <- resolve_arg(.weights_raw,   asset_names, "weights",
                                  default_val = 0, full_length = TRUE)

  if(!is.null(asset_shock)) {
    if(!is.numeric(asset_shock) || length(asset_shock) != n_assets)
      stop(sprintf(
        "asset_shock must be numeric of length %d (one per portfolio asset). Got length %d.",
        n_assets, length(asset_shock)))
    names(asset_shock) <- asset_names
  }

  if(is.null(reference_asset)) {
    if(!is.null(asset_shock) && sum(asset_shock != 0) == 1) {
      reference_asset <- asset_names[which(asset_shock != 0)[1]]
      cat(sprintf("\n  🔹 Reference asset auto-set to shocked asset: %s\n", reference_asset))
    } else if(!is.null(asset_shock) && sum(asset_shock != 0) > 1) {
      reference_asset <- asset_names[1]
      cat(sprintf("\n  ⚠️ Multiple shocks — reference asset defaulted to: %s\n", reference_asset))
    } else {
      reference_asset <- asset_names[1]
      cat(sprintf("\n  📌 Reference asset defaulted to: %s\n", reference_asset))
    }
  } else {
    if(!is.character(reference_asset) || length(reference_asset) != 1)
      stop("reference_asset must be a single character string")
    if(!(reference_asset %in% asset_names))
      stop(sprintf("reference_asset '%s' not found in: %s",
                   reference_asset, paste(asset_names, collapse = ", ")))
    cat(sprintf("\n  🎯 Reference asset: %s\n", reference_asset))
  }

  ref_idx <- which(asset_names == reference_asset)

  cat("\n========================================================================\n")
  cat(sprintf("PORTFOLIO RISK SIMULATION (%s RETURNS)\n", toupper(data_freq_label)))
  cat("========================================================================\n")
  cat(sprintf("  Assets                : %d — %s\n", n_assets, paste(asset_names, collapse = ", ")))
  cat(sprintf("  Reference asset       : %s (index %d)\n", reference_asset, ref_idx))
  cat(sprintf("  Observations          : %d\n", n_obs))
  cat(sprintf("  Distribution          : %s\n", sim_returns_dist))
  cat(sprintf("  Copula                : %s\n", ifelse(copula, toupper(copula_type), "NONE")))
  cat(sprintf("  Covariance estimation : %s\n", toupper(cov_estimation)))
  cat(sprintf("  Calibrate params      : %s\n", ifelse(calibrate_params, "YES", "NO")))
  cat(sprintf("  Propagation           : %s\n", ifelse(propagation, "YES", "NO")))
  if(!is.null(asset_shock) && any(asset_shock != 0))
    cat(sprintf("  Shock applied to      : %s\n",
                paste(asset_names[asset_shock != 0], collapse = ", ")))

  # ===========================================================================
  # 2. WEIGHTS PROCESSING
  # ===========================================================================

  # 2a-pre. Resolve exogenous list arguments (portfolio args resolved earlier)
  # ── Exogenous arguments (partial lists allowed — unspecified assets unshocked)
  exogenous_shock            <- if(has_exogenous && length(exog_names) > 0)
    resolve_arg(.exog_shock_raw, exog_names, "exogenous_shock",
                default_val = 0, full_length = FALSE) else .exog_shock_raw
  exogenous_volatility_shock <- if(has_exogenous && length(exog_names) > 0)
    resolve_arg(.exog_vs_raw,    exog_names, "exogenous_volatility_shock",
                default_val = 1, full_length = FALSE) else .exog_vs_raw

  # 2a. Build weights_matrix
  if(is.null(weights)) {
    weights_matrix           <- matrix(rep(1 / n_assets, n_assets), nrow = 1)
    rownames(weights_matrix) <- "Equal_Weight"
    colnames(weights_matrix) <- asset_names
    n_weight_sims            <- 1
    cat("\n  Weights: Equal weighted (1/N)\n")

  } else if(is.numeric(weights) && is.vector(weights) && !is.data.frame(weights)) {
    if(length(weights) != n_assets)
      stop(sprintf("weights vector length %d != n_assets %d", length(weights), n_assets))
    if(!is.null(names(weights))) {
      missing_names <- setdiff(asset_names, names(weights))
      extra_names   <- setdiff(names(weights), asset_names)
      if(length(missing_names) > 0)
        stop(sprintf("weights missing: %s", paste(missing_names, collapse = ", ")))
      if(length(extra_names) > 0)
        stop(sprintf("weights unknown: %s", paste(extra_names, collapse = ", ")))
      weights <- weights[asset_names]
    }
    weights_matrix           <- matrix(weights, nrow = 1)
    colnames(weights_matrix) <- asset_names
    rownames(weights_matrix) <- "Custom_Weight"
    n_weight_sims            <- 1
    cat("\n  Weights: Custom static (single portfolio)\n")

  } else if(is.data.frame(weights) || is.matrix(weights)) {
    weights_matrix <- as.matrix(weights)
    missing_cols   <- setdiff(asset_names, colnames(weights_matrix))
    extra_cols     <- setdiff(colnames(weights_matrix), asset_names)
    if(length(missing_cols) > 0)
      stop(sprintf("weights missing columns: %s", paste(missing_cols, collapse = ", ")))
    if(length(extra_cols) > 0)
      cat(sprintf("  ℹ️ Dropping unrecognised weight columns: %s\n",
                  paste(extra_cols, collapse = ", ")))
    weights_matrix <- weights_matrix[, asset_names, drop = FALSE]
    n_weight_sims  <- nrow(weights_matrix)
    if(n_weight_sims == 0) stop("weights data frame has 0 rows")
    cat(sprintf("\n  Weights: Data frame (%d portfolios)\n", n_weight_sims))

  } else {
    stop(sprintf("weights must be NULL, numeric vector, or data frame. Got: %s",
                 paste(class(weights), collapse = "/")))
  }

  # 2b. Numeric integrity
  if(any(!is.finite(weights_matrix))) {
    bad_rows <- which(apply(weights_matrix, 1, function(r) any(!is.finite(r))))
    stop(sprintf("Non-finite weights in rows: %s", paste(bad_rows, collapse = ", ")))
  }

  row_sums    <- rowSums(weights_matrix)
  severe_rows <- which(abs(row_sums - 1) > 0.05)
  minor_rows  <- which(abs(row_sums - 1) > 1e-6 & abs(row_sums - 1) <= 0.05)

  if(length(severe_rows) > 0)
    stop(sprintf("%d weight rows deviate >5%% from sum-to-1: rows %s\n  Sums: %s",
                 length(severe_rows), paste(severe_rows, collapse = ", "),
                 paste(round(row_sums[severe_rows], 4), collapse = ", ")))

  if(length(minor_rows) > 0) {
    warning(sprintf("%d weight rows renormalised (max dev: %.2e)",
                    length(minor_rows), max(abs(row_sums[minor_rows] - 1))))
    weights_matrix[minor_rows, ] <- weights_matrix[minor_rows, , drop = FALSE] /
                                     row_sums[minor_rows]
  }

  # 2c. Portfolio type
  has_shorts   <- any(weights_matrix < 0)
  has_leverage <- any(rowSums(weights_matrix) > 1 + 1e-6)
  is_long_only <- !has_shorts && !has_leverage
  if(has_shorts)   cat("  ℹ️ Short positions detected — long-short portfolio\n")
  if(has_leverage) cat("  ℹ️ Gross exposure > 1 — leveraged portfolio\n")
  if(is_long_only) cat("  ✅ Long-only portfolio\n")

  # 2d. Resolve n_sim_returns
  if(!is.null(n_sim_returns)) {
    if(!is.numeric(n_sim_returns) || length(n_sim_returns) != 1 ||
       n_sim_returns != floor(n_sim_returns) || n_sim_returns < 1)
      stop("n_sim_returns must be a single positive integer")
    n_sim_returns <- as.integer(n_sim_returns)
  } else {
    n_sim_returns <- if(n_weight_sims > 1) as.integer(n_weight_sims) else 10000L
  }

  total_sims <- n_weight_sims * n_sim_returns
  if(total_sims > 1e8)
    warning(sprintf("Total simulations = %s — may require significant memory.",
                    format(total_sims, big.mark = ",")))

  # 2e. Summary
  cat("\n========================================================================\n")
  cat("WEIGHTS SUMMARY\n")
  cat("========================================================================\n")
  cat(sprintf("  Portfolios            : %d\n", n_weight_sims))
  cat(sprintf("  Simulations each      : %s\n", format(n_sim_returns, big.mark = ",")))
  cat(sprintf("  Total simulations     : %s\n", format(total_sims, big.mark = ",")))
  cat(sprintf("  Portfolio type        : %s\n",
              if(has_shorts) "Long-Short" else if(has_leverage) "Leveraged" else "Long-Only"))

  if(n_weight_sims == 1) {
    cat("\n  Portfolio weights:\n")
    for(i in seq_len(n_assets)) {
      bar_len <- round(abs(weights_matrix[1, i]) * 20)
      bar     <- paste(rep(if(weights_matrix[1, i] >= 0) "█" else "░", bar_len), collapse = "")
      cat(sprintf("    %-20s : %+7.4f (%+6.2f%%)  %s\n",
                  asset_names[i], weights_matrix[1, i],
                  weights_matrix[1, i] * 100, bar))
    }
    cat(sprintf("    %-20s : %+7.4f\n", "SUM", sum(weights_matrix[1, ])))
  } else {
    cat("\n  Weight distribution across portfolios:\n")
    cat(sprintf("    %-20s %8s %8s %8s %8s\n", "Asset", "Min", "Mean", "Max", "SD"))
    for(i in seq_len(n_assets)) {
      w_col <- weights_matrix[, i]
      cat(sprintf("    %-20s %8.4f %8.4f %8.4f %8.4f\n",
                  asset_names[i], min(w_col), mean(w_col), max(w_col), sd(w_col)))
    }
  }

  # ===========================================================================
  # 3. FIT DISTRIBUTIONS TO HISTORICAL RETURNS
  # ===========================================================================

  cat("\n========================================================================\n")
  cat("FITTING DISTRIBUTIONS\n")
  cat("========================================================================\n")

  historical_cor <- cor(returns_numeric)
  historical_cov <- cov(returns_numeric)

  if(!is.positive.definite(historical_cov)) {
    cat("  ⚠️ Historical covariance not PD — applying nearPD\n")
    historical_cov <- as.matrix(nearPD(historical_cov, corr = FALSE)$mat)
  }
  if(!is.positive.definite(historical_cor)) {
    cat("  ⚠️ Historical correlation not PD — applying nearPD\n")
    historical_cor <- as.matrix(nearPD(historical_cor, corr = TRUE)$mat)
  }
  rownames(historical_cov) <- colnames(historical_cov) <- asset_names
  rownames(historical_cor) <- colnames(historical_cor) <- asset_names

  # Internal helper: moment-match SGT location/scale
  fit_sgt_moments <- function(target_mean, target_sd, lam, p, q,
                               n_mc = 100000, seed_offset = 0) {
    tryCatch({
      # Save and restore global RNG state so that set.seed() here does not
      # corrupt the random stream for subsequent copula fitting calls.
      old_seed <- if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
                    .GlobalEnv$.Random.seed else NULL
      on.exit({
        if(!is.null(old_seed)) .GlobalEnv$.Random.seed <- old_seed
        else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      }, add = TRUE)
      set.seed(seed + seed_offset)
      sgt_std  <- rsgt(n_mc, mu = 0, sigma = 1, lambda = lam, p = p, q = q)
      raw_mean <- mean(sgt_std)
      raw_sd   <- sd(sgt_std)
      if(!is.finite(raw_mean) || !is.finite(raw_sd) || raw_sd <= 0)
        stop("Non-finite moments from SGT draw")
      sigma_param <- target_sd / raw_sd
      mu_param    <- target_mean - sigma_param * raw_mean
      list(mu = mu_param, sigma = sigma_param, success = TRUE)
    }, error = function(e) {
      cat(sprintf("    ⚠️ SGT moment matching failed: %s\n", e$message))
      list(mu = target_mean, sigma = target_sd, success = FALSE)
    })
  }

  # garch_result is initialised here (function scope) so it is always accessible
  # in the plotting section, regardless of the copula or cov_estimation path taken.
  garch_result <- NULL   # populated in Step 4 when cov_estimation == "garch"

  # =========================================================================
  # PATH A: NO COPULA
  # =========================================================================
  if(!copula) {

    cat("\nPath: NO COPULA\n")
    base_mu    <- colMeans(returns_numeric)
    base_sigma <- historical_cov

    # --- A1: Multivariate Normal ---
    if(sim_returns_dist == "rmvnorm") {

      cat("\nDistribution: Multivariate Normal\n")
      fitted_params <- lapply(seq_len(n_assets), function(i)
        list(mean = base_mu[i], sd = sqrt(base_sigma[i, i]), dist = "normal"))
      names(fitted_params) <- asset_names

      cat(sprintf("  Means : %s\n",
                  paste(sprintf("%.4f (%.2f%%)", base_mu, base_mu * 100), collapse = ", ")))
      cat(sprintf("  Vols  : %s\n",
                  paste(sprintf("%.4f (%.2f%%)", sqrt(diag(base_sigma)),
                                sqrt(diag(base_sigma)) * 100), collapse = ", ")))
      cat("\nHistorical correlation matrix:\n"); print(round(historical_cor, 4))

    # --- A2: Multivariate Student-t ---
    } else if(sim_returns_dist == "rmvt") {

      cat("\nDistribution: Multivariate Student-t\n")

      df_vec <- if(calibrate_params && !is.null(calibrated_dfs)) {
        calibrated_dfs
      } else if(!is.null(asset_dfs)) {
        asset_dfs
      } else {
        rep(df, n_assets)
      }

      df_joint <- mean(df_vec)
      if(df_joint <= 2) stop(sprintf("df_joint = %.2f must be > 2", df_joint))
      cat(sprintf("  Joint df: %.4f\n", df_joint))

      delta_vec    <- colMeans(returns_numeric)
      scale_matrix <- historical_cov * (df_joint - 2) / df_joint
      scale_matrix <- (scale_matrix + t(scale_matrix)) / 2
      if(!is.positive.definite(scale_matrix)) {
        cat("  ⚠️ Scale matrix not PD — applying nearPD\n")
        scale_matrix <- as.matrix(nearPD(scale_matrix, corr = FALSE)$mat)
      }
      rownames(scale_matrix) <- colnames(scale_matrix) <- asset_names

      fitted_params <- lapply(seq_len(n_assets), function(i) {
        marginal_sd <- sqrt(df_joint / (df_joint - 2) * scale_matrix[i, i])
        list(mean = delta_vec[i], sd = marginal_sd, df = df_vec[i], dist = "t")
      })
      names(fitted_params) <- asset_names

      for(i in seq_len(n_assets))
        cat(sprintf("  %s: t(μ=%.4f, σ=%.4f, df=%.2f)\n",
                    asset_names[i], fitted_params[[i]]$mean,
                    fitted_params[[i]]$sd, fitted_params[[i]]$df))

      cat("\nHistorical correlation matrix:\n"); print(round(historical_cor, 4))

      base_mu    <- delta_vec
      base_sigma <- scale_matrix
      df         <- df_joint
      rmvt_delta <- delta_vec
      rmvt_scale <- scale_matrix
      rmvt_df    <- df_joint

    # --- A3: SGT marginals + Cholesky ---
    } else if(sim_returns_dist == "rsgt") {

      cat("\nDistribution: SGT marginals + Cholesky correlation\n")
      cat("Note: No closed-form multivariate SGT — Gaussian copula approximation\n\n")

      lam_vec <- if(calibrate_params && !is.null(calibrated_lambdas)) calibrated_lambdas
                 else if(!is.null(asset_lambdas)) asset_lambdas
                 else rep(-0.3, n_assets)
      p_vec   <- if(calibrate_params && !is.null(calibrated_ps)) calibrated_ps
                 else if(!is.null(asset_p)) asset_p
                 else rep(4, n_assets)
      q_vec   <- if(calibrate_params && !is.null(calibrated_qs)) calibrated_qs
                 else if(!is.null(asset_q)) asset_q
                 else rep(4, n_assets)

      fitted_params <- vector("list", n_assets)
      names(fitted_params) <- asset_names

      for(i in seq_len(n_assets)) {
        returns_i   <- returns_numeric[, i]
        lam_i <- lam_vec[i]; p_i <- p_vec[i]; q_i <- q_vec[i]

        if(abs(lam_i) >= 1) stop(sprintf("Asset %s: |lambda| must be < 1", asset_names[i]))
        if(p_i <= 0)         stop(sprintf("Asset %s: p must be > 0", asset_names[i]))
        if(q_i <= 1/p_i)     stop(sprintf("Asset %s: q must be > 1/p", asset_names[i]))
        if(q_i <= 2/p_i)     warning(sprintf("Asset %s: q <= 2/p — variance may be infinite", asset_names[i]))

        sgt_fit <- fit_sgt_moments(mean(returns_i), sd(returns_i), lam_i, p_i, q_i, seed_offset = i)

        fitted_params[[i]] <- list(
          mean = sgt_fit$mu, sd = sgt_fit$sigma,
          lambda = lam_i, p = p_i, q = q_i, dist = "sgt",
          moment_match_success = sgt_fit$success
        )
        cat(sprintf("  %s: SGT(μ=%.4f, σ=%.4f, λ=%.3f, p=%.2f, q=%.2f)%s\n",
                    asset_names[i], sgt_fit$mu, sgt_fit$sigma, lam_i, p_i, q_i,
                    if(!sgt_fit$success) " [fallback]" else ""))
      }

      chol_cor <- tryCatch(chol(historical_cor), error = function(e) {
        cat("  ⚠️ Cholesky failed — using nearPD\n")
        chol(as.matrix(nearPD(historical_cor, corr = TRUE)$mat))
      })

      cat("\nHistorical correlation matrix:\n"); print(round(historical_cor, 4))
      base_mu <- sapply(fitted_params, function(x) x$mean)
      names(base_mu) <- asset_names

    } else {
      stop("sim_returns_dist must be 'rmvnorm', 'rmvt', or 'rsgt'")
    }

    copula_success   <- FALSE
    copula_cor       <- NULL
    tail_dependence  <- NULL
    latent_Z         <- NULL
    ks_results_df    <- NULL   # populated in Step 6.5 when copula = TRUE
    copula_df_fitted <- NA

  # =========================================================================
  # PATH B: COPULA
  # =========================================================================
  } else {

    cat("\nPath: COPULA — marginals + copula separately\n")
    cat(sprintf("Copula type: %s\n", toupper(copula_type)))

    base_mu    <- colMeans(returns_numeric)
    base_sigma <- historical_cov

    fitted_params <- vector("list", n_assets)
    names(fitted_params) <- asset_names

    for(i in seq_len(n_assets)) {
      returns_i <- returns_numeric[, i]

      # B1: Normal
      if(sim_returns_dist == "rmvnorm") {
        fit_norm <- tryCatch(fitdistr(returns_i, "normal"), error = function(e) NULL)
        mu_i    <- if(!is.null(fit_norm)) fit_norm$estimate["mean"] else mean(returns_i)
        sigma_i <- if(!is.null(fit_norm)) fit_norm$estimate["sd"]   else sd(returns_i)
        fitted_params[[i]] <- list(mean = mu_i, sd = sigma_i, dist = "normal")
        cat(sprintf("  %s: N(μ=%.4f, σ=%.4f)\n", asset_names[i], mu_i, sigma_i))

      # B2: Student-t
      } else if(sim_returns_dist == "rmvt") {
        df_i <- if(calibrate_params && !is.null(calibrated_dfs)) calibrated_dfs[i]
                else if(!is.null(asset_dfs)) asset_dfs[i]
                else df
        if(df_i <= 2) stop(sprintf("Asset %s: df = %.2f must be > 2", asset_names[i], df_i))

        mu_i          <- mean(returns_i)
        sigma_i       <- sd(returns_i) * sqrt((df_i - 2) / df_i)
        marginal_sd_i <- sqrt(df_i / (df_i - 2)) * sigma_i

        fitted_params[[i]] <- list(
          mean = mu_i, sd = sigma_i, marginal_sd = marginal_sd_i, df = df_i, dist = "t")
        cat(sprintf("  %s: t(μ=%.4f, σ_scale=%.4f, σ_marginal=%.4f, df=%.2f)\n",
                    asset_names[i], mu_i, sigma_i, marginal_sd_i, df_i))

      # B3: SGT
      } else if(sim_returns_dist == "rsgt") {
        lam_i <- if(calibrate_params && !is.null(calibrated_lambdas)) calibrated_lambdas[i]
                 else if(!is.null(asset_lambdas)) asset_lambdas[i]
                 else {
                   sgt_calib <- tryCatch(
                     calibrate_sgt(returns_i, method = calibration_method),
                     error = function(e) list(lambda = -0.3, p = 4, q = 4))
                   sgt_calib$lambda
                 }
        p_i <- if(calibrate_params && !is.null(calibrated_ps)) calibrated_ps[i]
               else if(!is.null(asset_p)) asset_p[i] else 4
        q_i <- if(calibrate_params && !is.null(calibrated_qs)) calibrated_qs[i]
               else if(!is.null(asset_q)) asset_q[i] else 4

        if(abs(lam_i) >= 1) stop(sprintf("Asset %s: |lambda| must be < 1", asset_names[i]))
        if(p_i <= 0)         stop(sprintf("Asset %s: p must be > 0", asset_names[i]))
        if(q_i <= 1/p_i)     stop(sprintf("Asset %s: q must be > 1/p", asset_names[i]))
        if(q_i <= 2/p_i)     warning(sprintf("Asset %s: q <= 2/p — variance may be infinite", asset_names[i]))

        sgt_fit <- fit_sgt_moments(mean(returns_i), sd(returns_i), lam_i, p_i, q_i, seed_offset = i)

        fitted_params[[i]] <- list(
          mean = sgt_fit$mu, sd = sgt_fit$sigma, lambda = lam_i, p = p_i, q = q_i,
          dist = "sgt", moment_match_success = sgt_fit$success)
        cat(sprintf("  %s: SGT(μ=%.4f, σ=%.4f, λ=%.3f, p=%.2f, q=%.2f)%s\n",
                    asset_names[i], sgt_fit$mu, sgt_fit$sigma, lam_i, p_i, q_i,
                    if(!sgt_fit$success) " [fallback]" else ""))

      } else {
        stop("sim_returns_dist must be 'rmvnorm', 'rmvt', or 'rsgt'")
      }
    }

    cat("\n========================================================================\n")
    cat("MARGINAL FITTING SUMMARY\n")
    cat("========================================================================\n")
    cat(sprintf("  Distribution : %s\n", sim_returns_dist))
    cat(sprintf("  Assets fitted: %d\n", n_assets))
    cat(sprintf("  Fitted means : %s\n",
                paste(sprintf("%.4f", sapply(fitted_params, `[[`, "mean")), collapse = ", ")))
    cat(sprintf("  Fitted SDs   : %s\n",
                paste(sprintf("%.4f", sapply(fitted_params, `[[`, "sd")), collapse = ", ")))

    # =========================================================================
    # STEP 4: FORWARD-LOOKING VOLATILITIES
    # =========================================================================

    if(cov_estimation == "sample") {
      forward_cov <- cov(returns_numeric)
      cat("\n📊 Using SAMPLE covariance for forward volatilities\n")
    } else if(cov_estimation == "ewma") {
      forward_cov <- calculate_ewma_covariance(returns_numeric, lambda = cov_lambda)
      cat(sprintf("\n📊 Using EWMA covariance (lambda = %.2f)\n", cov_lambda))
    } else if(cov_estimation == "garch") {
      cat("\n📊 Estimating GARCH forward volatilities...\n")
      garch_result <- calculate_garch_vol(returns_numeric)
      forward_cov  <- garch_result$forecast_cov
      cat("GARCH estimation complete\n")
    } else {
      stop("cov_estimation must be 'sample', 'ewma', or 'garch'")
    }

    forward_vols <- sqrt(diag(forward_cov))
    names(forward_vols) <- asset_names

    cat("\n📈 Forward-looking volatilities:\n")
    for(i in seq_len(n_assets))
      cat(sprintf("  %s: %.2f%% (historical: %.2f%%)\n",
                  asset_names[i], forward_vols[i] * 100,
                  fitted_params[[i]]$sd * 100))

    # =========================================================================
    # STEP 5: FIT COPULA
    # =========================================================================

    cat("\n📈 Fitting copula for dependence structure...\n")

    # PIT using fitted marginals
    u_fitted <- matrix(NA, nrow = n_obs, ncol = n_assets)
    for(i in seq_len(n_assets)) {
      tryCatch({
        fp <- fitted_params[[i]]
        u_fitted[, i] <- if(sim_returns_dist == "rmvnorm") {
          pnorm(returns_numeric[, i], mean = fp$mean, sd = fp$sd)
        } else if(sim_returns_dist == "rmvt") {
          pt((returns_numeric[, i] - fp$mean) / fp$sd, df = fp$df)
        } else {
          psgt(returns_numeric[, i], mu = fp$mean, sigma = fp$sd,
               lambda = fp$lambda, p = fp$p, q = fp$q)
        }
        u_fitted[, i] <- pmax(pmin(u_fitted[, i], 1 - 1e-7), 1e-7)
      }, error = function(e) {
        cat(sprintf("  ⚠️ PIT failed for %s — using empirical ranks\n", asset_names[i]))
        u_fitted[, i] <<- pobs(returns_numeric[, i, drop = FALSE])[, 1]
      })
    }

    cat(sprintf("  PIT method: fitted %s marginals\n", toupper(sim_returns_dist)))
    cat(sprintf("  Uniform range: [%.4f, %.4f]\n", min(u_fitted), max(u_fitted)))

    # -------------------------------------------------------------------------
    # 5b. Determine feasible copula structure based on n_obs / n_params ratio.
    #
    # fitCopula with method="mpl" (max pseudo-likelihood) is a numerical
    # optimisation. When observations per parameter is too low, the likelihood
    # surface is nearly flat and the optimiser diverges or fails.
    #
    # Parameter counts per copula structure:
    #   "un"   (unstructured) : n*(n-1)/2 correlation params  — fully flexible
    #   "toep" (Toeplitz)     : n-1 params                    — lag-structured
    #   "ex"   (exchangeable) : 1 param                       — single common ρ
    #
    # Tiered strategy based on obs/params ratio:
    #   ratio >= 10 → "un"   + mpl   (fully flexible, well-identified)
    #   ratio >= 5  → "un"   + itau  (rank-based inversion, closed-form)
    #   ratio >= 3  → "toep" + itau  (structured + rank-based)
    #   ratio <  3  → "ex"   + itau  (most parsimonious)
    #
    # itau (inversion of Kendall's tau) is closed-form and NEVER fails.
    # -------------------------------------------------------------------------
    n_rho_un   <- n_assets * (n_assets - 1L) / 2L
    n_params_t <- n_rho_un + 1L          # +1 for df
    n_params_g <- n_rho_un

    obs_per_param <- n_obs / if(copula_type == "t") n_params_t else n_params_g

    if(obs_per_param >= 10) {
      fit_dispstr <- "un";   fit_method <- "mpl";  fit_label <- "unstructured + MPL"
    } else if(obs_per_param >= 5) {
      fit_dispstr <- "un";   fit_method <- "itau"; fit_label <- "unstructured + iTau"
    } else if(obs_per_param >= 3) {
      fit_dispstr <- "toep"; fit_method <- "itau"; fit_label <- "Toeplitz + iTau"
    } else {
      fit_dispstr <- "ex";   fit_method <- "itau"; fit_label <- "exchangeable + iTau"
    }

    cat(sprintf("  n_obs=%d | n_params=%d | ratio=%.1f → strategy: %s\n",
                n_obs,
                if(copula_type == "t") n_params_t else n_params_g,
                obs_per_param, fit_label))

    # -------------------------------------------------------------------------
    # 5c. Warm start: convert Kendall's τ → copula ρ analytically.
    #   t-copula:        ρ = sin(π/2 × τ)        (Lindskog et al. 2003)
    #   Gaussian copula: ρ = 2 sin(π/6 × ρ_s)   (van der Waerden / Pearson)
    # These always lie in (-1,1) and satisfy model assumptions — far better
    # than the default zero-start used by fitCopula.
    # -------------------------------------------------------------------------
    kendall_tau     <- cor(u_fitted, method = "kendall")
    spearman_cor_u  <- cor(u_fitted, method = "spearman")

    warmstart_cor <- if(copula_type == "t") sin(pi / 2 * kendall_tau) else
                       2 * sin(pi / 6 * spearman_cor_u)
    diag(warmstart_cor) <- 1
    warmstart_cor <- (warmstart_cor + t(warmstart_cor)) / 2
    if(!is.positive.definite(warmstart_cor))
      warmstart_cor <- as.matrix(nearPD(warmstart_cor, corr = TRUE)$mat)

    warmstart_params <- pmax(pmin(warmstart_cor[upper.tri(warmstart_cor)], 0.98), -0.98)
    cat(sprintf("  Warm-start ρ: [%.3f, %.3f] (from %s)\n",
                min(warmstart_params), max(warmstart_params),
                if(copula_type == "t") "Kendall's τ" else "Spearman's ρ"))

    # -------------------------------------------------------------------------
    # 5d. Extract correlation matrix from a fitted copula object.
    # Works for "un", "ex", and "toep" dispstr values.
    # -------------------------------------------------------------------------
    extract_cor_matrix <- function(params_est, dispstr, n) {
      cm <- diag(n)
      if(dispstr == "un") {
        k <- 1
        for(ii in 1:(n - 1))
          for(jj in (ii + 1):n) {
            cm[ii, jj] <- cm[jj, ii] <- params_est[k]; k <- k + 1
          }
      } else if(dispstr == "ex") {
        cm[cm == 0] <- params_est[1]
      } else if(dispstr == "toep") {
        for(ii in 1:(n - 1))
          for(jj in (ii + 1):n) {
            lag <- jj - ii
            if(lag <= length(params_est))
              cm[ii, jj] <- cm[jj, ii] <- params_est[lag]
          }
      }
      cm
    }

    copula_result <- tryCatch({

      if(copula_type == "gaussian") {

        gauss_cop <- normalCopula(param = warmstart_params,
                                   dim = n_assets, dispstr = fit_dispstr)
        fit <- if(fit_method == "mpl") {
          fitCopula(gauss_cop, u_fitted, method = "mpl",
                     start = warmstart_params,
                     optim.method = "L-BFGS-B",
                     optim.control = list(maxit = 1000L, factr = 1e7))
        } else {
          fitCopula(gauss_cop, u_fitted, method = "itau")
        }

        cor_matrix <- extract_cor_matrix(fit@estimate, fit_dispstr, n_assets)
        if(!is.positive.definite(cor_matrix)) {
          cat("  ⚠️ Gaussian copula correlation not PD — applying nearPD\n")
          cor_matrix <- as.matrix(nearPD(cor_matrix, corr = TRUE)$mat)
        }

        ll_val <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
        cat(sprintf("  Gaussian copula fitted [%s]. Log-lik: %.2f\n", fit_label,
                    if(is.finite(ll_val)) ll_val else NA))

        list(success = TRUE, cor_matrix = cor_matrix, fit = fit, type = "gaussian",
             dispstr = fit_dispstr, method_used = fit_method,
             df = NULL, tail_dep_mat = NULL,
             loglik = ll_val,
             aic = tryCatch(AIC(fit), error = function(e) NA_real_),
             bic = tryCatch(BIC(fit), error = function(e) NA_real_))

      } else if(copula_type == "t") {

        df_init        <- 8
        warmstart_full <- c(warmstart_params, df_init)

        t_cop_obj <- tCopula(param = warmstart_params, dim = n_assets,
                              dispstr = fit_dispstr, df = df_init, df.fixed = FALSE)

        # ── Fit correlation structure ─────────────────────────────────────────
        # For method="itau": the copula package estimates the CORRELATION params
        # via rank-based Kendall tau inversion (always succeeds), then attempts
        # to estimate df via profile MLE conditional on those correlations.
        # The profile MLE for df can silently return NA with large n_assets or
        # near-1 pairwise correlations (ill-conditioned density evaluation).
        #
        # Strategy: fit correlations first (itau always works), then recover df:
        #   1. Try df from getTheta() — use it if finite and in (2, 200)
        #   2. Otherwise: 1-D profile MLE over df ∈ [2.1, 50] with rho fixed
        #   3. Final fallback: df_init (= 8)
        fit <- if(fit_method == "mpl") {
          fitCopula(t_cop_obj, u_fitted, method = "mpl",
                     start = warmstart_full,
                     optim.method = "L-BFGS-B",
                     optim.control = list(maxit = 1000L, factr = 1e7))
        } else {
          fitCopula(t_cop_obj, u_fitted, method = "itau")
        }

        theta    <- getTheta(fit@copula)
        n_rho_fit <- switch(fit_dispstr,
                            "un"   = n_assets * (n_assets - 1L) / 2L,
                            "ex"   = 1L,
                            "toep" = n_assets - 1L,
                            n_assets * (n_assets - 1L) / 2L)
        cor_params <- theta[seq_len(n_rho_fit)]
        df_est_raw <- theta[n_rho_fit + 1L]

        cor_matrix <- extract_cor_matrix(cor_params, fit_dispstr, n_assets)
        if(!is.positive.definite(cor_matrix)) {
          cat("  ⚠️ t-copula correlation not PD — applying nearPD\n")
          cor_matrix <- as.matrix(nearPD(cor_matrix, corr = TRUE)$mat)
        }

        # ── Recover df when itau returns NA ──────────────────────────────────
        # itau fits rho analytically but df via profile MLE may fail → NA.
        # When that happens, do a 1-D grid search over df with rho fixed.
        df_est <- if(is.finite(df_est_raw) && df_est_raw > 2 && df_est_raw < 500) {
          df_est_raw
        } else {
          if(fit_method == "itau") {
            cat("  ⚠️ df = NA from itau — estimating df via 1-D profile MLE\n")
            # Fix rho at the itau estimate; search over df only
            df_search <- tryCatch({
              neg_ll_df <- function(df_val) {
                if(df_val <= 2) return(1e12)
                t_fixed <- tCopula(param = cor_params, dim = n_assets,
                                    dispstr = fit_dispstr, df = df_val,
                                    df.fixed = TRUE)
                tryCatch(-sum(dCopula(u_fitted, t_fixed, log = TRUE)),
                         error = function(e) 1e12)
              }
              df_grid  <- c(3, 4, 5, 6, 8, 10, 15, 20, 30)
              ll_grid  <- sapply(df_grid, neg_ll_df)
              best_df  <- df_grid[which.min(ll_grid)]
              # Refine with Brent's method in [2.1, 50]
              opt <- optimize(neg_ll_df, interval = c(2.1, 50), tol = 0.1)
              if(is.finite(opt$minimum) && opt$minimum > 2) opt$minimum else best_df
            }, error = function(e) {
              cat(sprintf("  ⚠️ 1-D df MLE failed (%s) — using df_init=%.0f\n",
                          e$message, df_init))
              df_init
            })
            cat(sprintf("  ℹ️ df estimated via profile MLE: %.2f\n", df_search))
            df_search
          } else {
            cat(sprintf("  ⚠️ df = NA from mpl — using df_init=%.0f\n", df_init))
            df_init
          }
        }

        ll_val <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
        cat(sprintf("  t-copula fitted [%s]. df=%.2f | Log-lik: %s\n",
                    fit_label, df_est,
                    if(is.finite(ll_val)) sprintf("%.2f", ll_val) else "NA"))
        if(fit_dispstr != "un")
          cat(sprintf("  ℹ️ Structured copula (%s): correlations are constrained\n",
                      fit_dispstr))
        if(is.finite(df_est) && df_est > 30)
          cat("  ℹ️ df > 30: t-copula near-Gaussian, tail dependence minimal\n")

        # ── Tail dependence (only valid when df is finite) ────────────────────
        tail_dep_mat <- matrix(NA_real_, n_assets, n_assets)
        diag(tail_dep_mat) <- 1
        if(is.finite(df_est) && df_est > 2) {
          for(ii in 1:(n_assets - 1))
            for(jj in (ii + 1):n_assets) {
              rho_ij <- cor_matrix[ii, jj]
              if(is.finite(rho_ij) && abs(rho_ij) < 1) {
                td <- 2 * pt(-sqrt((df_est + 1) * (1 - rho_ij) / (1 + rho_ij)),
                             df = df_est + 1)
                tail_dep_mat[ii, jj] <- tail_dep_mat[jj, ii] <- td
              }
            }
        }
        colnames(tail_dep_mat) <- rownames(tail_dep_mat) <- asset_names

        cat("\n  Tail Dependence Coefficients:\n")
        for(ii in 1:(n_assets - 1))
          for(jj in (ii + 1):n_assets)
            cat(sprintf("    %s \u2014 %s: \u03c1=%.3f, \u03bb=%.3f\n",
                        asset_names[ii], asset_names[jj],
                        cor_matrix[ii, jj], tail_dep_mat[ii, jj]))

        list(success = TRUE, cor_matrix = cor_matrix, fit = fit, type = "t",
             dispstr = fit_dispstr, method_used = fit_method,
             df = df_est, tail_dep_mat = tail_dep_mat,
             loglik = ll_val,
             aic = tryCatch(AIC(fit), error = function(e) NA_real_),
             bic = tryCatch(BIC(fit), error = function(e) NA_real_))

      } else {
        stop("copula_type must be 'gaussian' or 't'")
      }

    }, error = function(e) {
      # If the tiered strategy still fails, fall back to Kendall's τ inversion.
      # This is closed-form and always produces a valid correlation matrix.
      cat(sprintf("  ⚠️ Copula fitting failed: %s\n  Falling back to Kendall's \u03c4 inversion\n",
                  e$message))
      fallback_cor <- if(copula_type == "t") sin(pi / 2 * kendall_tau) else
                        2 * sin(pi / 6 * spearman_cor_u)
      diag(fallback_cor) <- 1
      fallback_cor <- (fallback_cor + t(fallback_cor)) / 2
      if(!is.positive.definite(fallback_cor))
        fallback_cor <- as.matrix(nearPD(fallback_cor, corr = TRUE)$mat)
      colnames(fallback_cor) <- rownames(fallback_cor) <- asset_names
      cat("  ℹ️ Kendall \u03c4-based correlation used as dependence structure\n")
      list(success = FALSE, cor_matrix = fallback_cor, fit = NULL,
           type = "fallback_kendall", dispstr = "un", method_used = "kendall_inversion",
           df = if(copula_type == "t") 8 else NULL, tail_dep_mat = NULL,
           loglik = NA_real_, aic = NA_real_, bic = NA_real_)
    })

    copula_success <- copula_result$success
    copula_cor     <- copula_result$cor_matrix
    colnames(copula_cor) <- rownames(copula_cor) <- asset_names

    if(copula_type == "t" && copula_success && !is.null(copula_result$df)) {
      copula_df_fitted <- copula_result$df
      if(abs(copula_df_fitted - copula_df) > 0.5)
        cat(sprintf("\n  ℹ️ Fitted copula df (%.2f) differs from input copula_df (%d).\n",
                    copula_df_fitted, copula_df))
    } else {
      copula_df_fitted <- copula_df
    }

    tail_dependence <- if(copula_type == "t" && copula_success) copula_result$tail_dep_mat else NULL

    if(copula_success) {
      cat(sprintf("\n  Copula fit: Log-lik=%.2f  AIC=%.2f  BIC=%.2f\n",
                  copula_result$loglik, copula_result$aic, copula_result$bic))
    }

    # =========================================================================
    # STEP 6: BUILD FINAL COVARIANCE MATRIX Σ = D × ρ × D
    # =========================================================================

    # 6a. Validate inputs
    if(any(!is.finite(forward_vols)) || any(forward_vols <= 0)) {
      cat("  ⚠️ Non-positive forward vols detected — replacing with historical\n")
      historical_vols <- sqrt(diag(historical_cov))
      forward_vols[!is.finite(forward_vols) | forward_vols <= 0] <-
        historical_vols[!is.finite(forward_vols) | forward_vols <= 0]
    }
    if(!is.positive.definite(copula_cor)) {
      cat("  ⚠️ copula_cor not PD — applying nearPD\n")
      copula_cor <- as.matrix(nearPD(copula_cor, corr = TRUE)$mat)
      colnames(copula_cor) <- rownames(copula_cor) <- asset_names
    }

    # 6b. Construct
    D_forward  <- diag(forward_vols, nrow = n_assets, ncol = n_assets)
    base_sigma <- D_forward %*% copula_cor %*% D_forward
    base_sigma <- (base_sigma + t(base_sigma)) / 2
    rownames(base_sigma) <- colnames(base_sigma) <- asset_names

    if(!is.positive.definite(base_sigma)) {
      cat("  ⚠️ Final Σ not PD — applying nearPD\n")
      base_sigma <- as.matrix(nearPD(base_sigma, corr = FALSE, ensureSymmetry = TRUE)$mat)
      rownames(base_sigma) <- colnames(base_sigma) <- asset_names
    }

    implied_cor   <- cov2cor(base_sigma)
    cor_deviation <- max(abs(implied_cor - copula_cor))
    if(cor_deviation > 0.01)
      cat(sprintf("  ⚠️ Max copula/implied correlation deviation: %.4f\n", cor_deviation))

    # 6c. Forward-looking means
    if(cov_estimation == "ewma") {
      ewma_weights <- (1 - cov_lambda) * cov_lambda ^ (rev(seq_len(n_obs) - 1))
      ewma_weights <- ewma_weights / sum(ewma_weights)
      base_mu      <- apply(returns_numeric, 2, function(x) sum(ewma_weights * x))
      cat(sprintf("\n📈 EWMA expected returns (lambda=%.2f)\n", cov_lambda))
    } else {
      base_mu <- sapply(fitted_params, function(x) x$mean)
      cat(sprintf("\n📈 Using fitted marginal means (%s path)\n", cov_estimation))
    }
    names(base_mu) <- asset_names

    # 6d. rmvt scale matrix for copula path
    if(sim_returns_dist == "rmvt") {
      df_for_scale <- if(exists("df_joint") && is.finite(df_joint)) df_joint else copula_df_fitted
      rmvt_scale   <- base_sigma * (df_for_scale - 2) / df_for_scale
      rmvt_scale   <- fix_covariance(rmvt_scale)
      rownames(rmvt_scale) <- colnames(rmvt_scale) <- asset_names
      cat(sprintf("\n📐 rmvt scale matrix (df=%.2f)\n", df_for_scale))
    }

    # 6e. Diagnostics
    cat("\n========================================================================\n")
    cat("STEP 6 DIAGNOSTICS\n")
    cat("========================================================================\n")
    hist_vols <- sqrt(diag(historical_cov))
    cat("\nForward vs historical volatilities:\n")
    for(i in seq_len(n_assets)) {
      dir <- if(forward_vols[i] > hist_vols[i]) "↑" else if(forward_vols[i] < hist_vols[i]) "↓" else "="
      cat(sprintf("  %s: %.4f%% %s (historical: %.4f%%)\n",
                  asset_names[i], forward_vols[i]*100, dir, hist_vols[i]*100))
    }
    cat("\nHistorical correlation:\n");    print(round(historical_cor, 4))
    cat("\nCopula correlation:\n");        print(round(copula_cor, 4))
    cat("\nDifference (Copula − Historical):\n")
    diff_mat <- copula_cor - historical_cor
    print(round(diff_mat, 4))
    cat(sprintf("  Max |diff|: %.4f\n", max(abs(diff_mat[upper.tri(diff_mat)]))))
    cat("\nFinal Σ:\n"); print(round(base_sigma, 6))
    cat(sprintf("\n  Condition number : %.2f\n", kappa(base_sigma)))
    cat(sprintf("  Min eigenvalue   : %.2e\n",
                min(eigen(base_sigma, only.values = TRUE)$values)))

    # =========================================================================
    # STEP 6.5: LATENT VARIABLES FOR DIAGNOSTIC PLOTS
    # =========================================================================

    if(is.null(fitted_params) || length(fitted_params) == 0)
      stop("fitted_params empty entering Step 6.5")

    u_data      <- matrix(NA_real_, nrow = n_obs, ncol = n_assets,
                          dimnames = list(NULL, asset_names))
    pit_warnings <- character(0)

    for(i in seq_len(n_assets)) {
      x_i <- returns_numeric[, i]; fp <- fitted_params[[i]]; asset_i <- asset_names[i]

      u_i <- tryCatch({
        if(sim_returns_dist == "rmvnorm") {
          if(is.null(fp$mean) || is.null(fp$sd) || fp$sd <= 0)
            stop("Invalid normal params")
          pnorm(x_i, mean = fp$mean, sd = fp$sd)

        } else if(sim_returns_dist == "rmvt") {
          if(is.null(fp$mean) || is.null(fp$sd) || is.null(fp$df) || fp$sd <= 0 || fp$df <= 0)
            stop("Invalid t params")
          pt((x_i - fp$mean) / fp$sd, df = fp$df)

        } else if(sim_returns_dist == "rsgt") {
          if(is.null(fp$mean) || is.null(fp$sd) || is.null(fp$lambda) ||
             is.null(fp$p) || is.null(fp$q) || fp$sd <= 0 || abs(fp$lambda) >= 1)
            stop("Invalid SGT params")
          psgt(x_i, mu = fp$mean, sigma = fp$sd, lambda = fp$lambda, p = fp$p, q = fp$q)
        }
      }, error = function(e) {
        pit_warnings <<- c(pit_warnings,
                           sprintf("%s: PIT failed (%s) — empirical ranks", asset_i, e$message))
        as.vector(pobs(matrix(x_i, ncol = 1)))
      })

      u_i <- pmax(pmin(u_i, 1 - 1e-7), 1e-7)
      if(any(!is.finite(u_i))) { u_i[!is.finite(u_i)] <- 0.5 }
      u_data[, i] <- u_i
    }

    if(length(pit_warnings) > 0) {
      cat("\n  ⚠️ PIT warnings:\n")
      for(w in pit_warnings) cat(sprintf("    %s\n", w))
    }

    # KS uniformity test
    # Ties in the PIT scores (from identical raw returns or boundary values of
    # psgt) cause ks.test() to emit "ties should not be present" warnings.
    # The test is still valid but the warning is noisy. We break ties by adding
    # a tiny jitter (1e-7 scale) before the test only — u_data itself is unchanged.
    cat("\n  Uniformity diagnostic (KS test):\n")
    ks_results <- data.frame(Asset = asset_names, KS_stat = NA_real_,
                              p_value = NA_real_, Flag = "", stringsAsFactors = FALSE)
    for(i in seq_len(n_assets)) {
      u_ks <- u_data[, i]
      if(any(duplicated(u_ks))) {
        set.seed(seed + i)   # reproducible jitter
        u_ks <- u_ks + runif(length(u_ks), -1e-7, 1e-7)
        u_ks <- pmax(pmin(u_ks, 1 - 1e-10), 1e-10)
      }
      ks_i <- tryCatch(suppressWarnings(ks.test(u_ks, "punif", 0, 1)),
                       error = function(e) NULL)
      if(!is.null(ks_i)) {
        ks_results$KS_stat[i] <- round(ks_i$statistic, 4)
        ks_results$p_value[i] <- round(ks_i$p.value,   4)
        ks_results$Flag[i]    <- if(ks_i$p.value < 0.05) "⚠️ REJECT" else "✅ OK"
      }
    }
    print(ks_results, row.names = FALSE)
    ks_results_df <- ks_results   # exposed as result$ks_test

    # Latent Z
    latent_Z <- apply(u_data, 2, qnorm)
    colnames(latent_Z) <- asset_names
    n_inf_Z <- sum(!is.finite(latent_Z))
    if(n_inf_Z > 0) {
      warning(sprintf("%d non-finite values in latent_Z", n_inf_Z))
      for(i in seq_len(n_assets)) {
        bad_idx <- !is.finite(latent_Z[, i])
        if(any(bad_idx)) latent_Z[bad_idx, i] <- median(latent_Z[!bad_idx, i], na.rm = TRUE)
      }
    }

    # Tail dependence
    if(copula_type == "gaussian" && copula_success) {
      cat("\n  ℹ️ Gaussian copula: asymptotic tail dependence = 0\n")
    } else if(copula_type == "t" && copula_success) {
      df_used <- if(!is.null(copula_result$df) && is.finite(copula_result$df))
                   copula_result$df else copula_df

      if(!is.null(copula_result$tail_dep_mat)) {
        tail_dependence <- copula_result$tail_dep_mat
        cat(sprintf("\n  Tail dependence retrieved from Step 5 (df=%.4f)\n", df_used))
      } else {
        tail_dependence <- matrix(NA_real_, n_assets, n_assets,
                                   dimnames = list(asset_names, asset_names))
        diag(tail_dependence) <- 1
        for(ii in seq_len(n_assets - 1))
          for(jj in seq(ii + 1, n_assets)) {
            rho_ij <- copula_cor[ii, jj]
            if(!is.na(rho_ij) && is.finite(rho_ij) && abs(rho_ij) < 1) {
              td <- 2 * pt(-sqrt((df_used + 1) * (1 - rho_ij) / (1 + rho_ij)), df = df_used + 1)
              tail_dependence[ii, jj] <- tail_dependence[jj, ii] <- td
            }
          }
      }

      cat(sprintf("\n  %-25s %8s %12s\n", "Pair", "ρ_copula", "λ"))
      cat(sprintf("  %s\n", strrep("-", 48)))
      for(ii in seq_len(n_assets - 1))
        for(jj in seq(ii + 1, n_assets)) {
          td_val  <- tail_dependence[ii, jj]
          rho_val <- copula_cor[ii, jj]
          cat(sprintf("  %-25s %8.4f %12.4f%s\n",
                      paste(asset_names[ii], "—", asset_names[jj]),
                      rho_val, td_val,
                      if(!is.na(td_val) && td_val > 0.3) " ⚠️ HIGH" else ""))
        }
    }

    # Latent Z vs copula correlation check
    latent_cor <- cor(latent_Z)
    if(!is.null(copula_cor)) {
      max_dev <- max(abs(latent_cor - copula_cor))
      cat(sprintf("\n  Max |latent_cor − copula_cor|: %.4f%s\n",
                  max_dev,
                  if(max_dev > 0.05) " ⚠️" else " ✅"))
    }

  }  # end PATH B / copula else block

  # Final base_sigma validation (unconditional)
  base_sigma <- (base_sigma + t(base_sigma)) / 2
  rownames(base_sigma) <- colnames(base_sigma) <- asset_names
  if(!is.positive.definite(base_sigma)) {
    cat("\n  ⚠️ base_sigma not PD before stress testing — applying nearPD\n")
    base_sigma <- as.matrix(nearPD(base_sigma, corr = FALSE, ensureSymmetry = TRUE)$mat)
    rownames(base_sigma) <- colnames(base_sigma) <- asset_names
  }
  cat(sprintf("\n  Final Σ condition number : %.2f\n", kappa(base_sigma)))
  cat(sprintf("  Final Σ min eigenvalue   : %.2e\n",
              min(eigen(base_sigma, only.values = TRUE)$values)))

  # ===========================================================================
  # 7. STRESS TESTING
  # ===========================================================================

  mu_stressed    <- base_mu
  sigma_stressed <- base_sigma
  is_stressed    <- FALSE

  # Three shock scenarios supported:
  #   1. Endogenous only  : asset_shock / volatility_shock target portfolio assets
  #   2. Exogenous only   : exogenous_shock / exogenous_volatility_shock target
  #                         assets in exogenous_returns (not in portfolio)
  #   3. Mixed            : both sets of shocks applied simultaneously
  #
  # All scenarios use the Schur complement conditional propagation formula.
  # The "shocked block" s is the union of endogenous and exogenous shocked assets.
  # The "non-shocked block" n is the remaining portfolio assets.

  has_endo_shock <- (!is.null(asset_shock)      && any(asset_shock != 0)) ||
                   (!is.null(volatility_shock)  && any(volatility_shock != 1))
  has_exog_shock <- !is.null(exogenous_shock) || !is.null(exogenous_volatility_shock)

  if(has_endo_shock || has_exog_shock) {

    is_stressed <- TRUE
    cat("\n========================================================================\n")
    cat("STRESS TESTING\n")
    cat("========================================================================\n")

    # 7b. Validate endogenous shocks
    if(!is.null(asset_shock)) {
      if(!is.numeric(asset_shock) || length(asset_shock) != n_assets)
        stop(sprintf("asset_shock must be numeric length %d (one per portfolio asset)", n_assets))
      if(any(!is.finite(asset_shock))) stop("asset_shock contains non-finite values")
      names(asset_shock) <- asset_names
      large_s <- which(abs(asset_shock) > 0.5)
      if(length(large_s) > 0)
        warning(sprintf("asset_shock > 50%% for: %s — verify decimal units",
                        paste(asset_names[large_s], collapse = ", ")))
    }
    if(!is.null(volatility_shock)) {
      if(!is.numeric(volatility_shock) || length(volatility_shock) != n_assets)
        stop(sprintf("volatility_shock must be numeric length %d", n_assets))
      if(any(!is.finite(volatility_shock))) stop("volatility_shock contains non-finite values")
      if(any(volatility_shock < 0))  stop("volatility_shock must be non-negative")
      zero_vs <- which(volatility_shock == 0)
      if(length(zero_vs) > 0)
        stop(sprintf(paste0(
          "volatility_shock = 0 for: %s.\n",
          "  A multiplier of 0 collapses asset variance to zero, making the\n",
          "  covariance matrix singular and the Schur complement undefined.\n",
          "  Use a small positive value (e.g. 0.01) if you want near-zero vol,\n",
          "  or set to 1 to leave that asset's volatility unchanged."),
          paste(asset_names[zero_vs], collapse = ", ")))
      names(volatility_shock) <- asset_names
    }

    # 7b2. Validate exogenous shocks
    if(!is.null(exogenous_shock)) {
      if(!has_exogenous)
        stop("exogenous_shock supplied but exogenous_returns = NULL.")
      if(!is.numeric(exogenous_shock) || is.null(names(exogenous_shock)))
        stop("exogenous_shock must be a NAMED numeric vector")
      bad_e <- setdiff(names(exogenous_shock), exog_names)
      if(length(bad_e) > 0)
        stop(sprintf("exogenous_shock names not in exogenous_returns: %s",
                     paste(bad_e, collapse = ", ")))
      if(any(!is.finite(exogenous_shock))) stop("exogenous_shock contains non-finite values")
    }
    if(!is.null(exogenous_volatility_shock)) {
      if(!has_exogenous)
        stop("exogenous_volatility_shock supplied but exogenous_returns = NULL.")
      if(!is.numeric(exogenous_volatility_shock) || is.null(names(exogenous_volatility_shock)))
        stop("exogenous_volatility_shock must be a NAMED numeric vector")
      bad_ev <- setdiff(names(exogenous_volatility_shock), exog_names)
      if(length(bad_ev) > 0)
        stop(sprintf("exogenous_volatility_shock names not in exogenous_returns: %s",
                     paste(bad_ev, collapse = ", ")))
      if(any(!is.finite(exogenous_volatility_shock)))
        stop("exogenous_volatility_shock contains non-finite values")
      if(any(exogenous_volatility_shock < 0))
        stop("exogenous_volatility_shock must be non-negative")
      zero_evs <- names(exogenous_volatility_shock)[exogenous_volatility_shock == 0]
      if(length(zero_evs) > 0)
        stop(sprintf(paste0(
          "exogenous_volatility_shock = 0 for: %s.\n",
          "  A multiplier of 0 collapses exogenous variance to zero, making\n",
          "  the cross-covariance block Sigma_ee singular and non-invertible.\n",
          "  Use a small positive value or omit that asset from the shock."),
          paste(zero_evs, collapse = ", ")))
    }

    # 7c. Identify shocked index sets
    shocked_assets     <- if(!is.null(asset_shock)) which(asset_shock != 0) else integer(0)
    non_shocked_assets <- setdiff(seq_len(n_assets), shocked_assets)
    vol_shocked_assets <- if(!is.null(volatility_shock)) which(volatility_shock != 1) else integer(0)
    n_shocked          <- length(shocked_assets)
    n_vol_shocked      <- length(vol_shocked_assets)

    exog_shocked_idx     <- if(!is.null(exogenous_shock))
                              match(names(exogenous_shock), exog_names) else integer(0)
    exog_vol_shocked_idx <- if(!is.null(exogenous_volatility_shock))
                              match(names(exogenous_volatility_shock), exog_names) else integer(0)
    n_exog_shocked       <- length(exog_shocked_idx)

    cat(sprintf("  Portfolio return-shocked : %d — %s\n", n_shocked,
                if(n_shocked > 0) paste(asset_names[shocked_assets], collapse=", ") else "none"))
    cat(sprintf("  Portfolio vol-shocked    : %d — %s\n", n_vol_shocked,
                if(n_vol_shocked > 0) paste(asset_names[vol_shocked_assets], collapse=", ") else "none"))
    if(has_exogenous) {
      cat(sprintf("  Exogenous return-shocked : %d — %s\n", n_exog_shocked,
                  if(n_exog_shocked > 0) paste(exog_names[exog_shocked_idx], collapse=", ") else "none"))
      cat(sprintf("  Exogenous vol-shocked    : %d — %s\n", length(exog_vol_shocked_idx),
                  if(length(exog_vol_shocked_idx) > 0)
                    paste(exog_names[exog_vol_shocked_idx], collapse=", ") else "none"))
    }
    cat(sprintf("  Propagation              : %s\n", ifelse(propagation, "ENABLED", "DISABLED")))

    # ─────────────────────────────────────────────────────────────────────────
    # 7d. Portfolio volatility shock: Σ_pp_s = K_p × Σ_pp × K_p
    # ─────────────────────────────────────────────────────────────────────────
    if(!is.null(volatility_shock) && n_vol_shocked > 0) {
      cat("\n--- PORTFOLIO VOLATILITY SHOCK ---\n")
      cat(sprintf("  %-20s %10s %12s %12s\n", "Asset", "Multiplier", "\u03c3_base", "\u03c3_stressed"))
      cat(sprintf("  %s\n", strrep("-", 58)))
      for(i in seq_len(n_assets)) {
        k_i <- volatility_shock[i]; sb <- sqrt(sigma_stressed[i, i])
        if(k_i != 1)
          cat(sprintf("  %-20s %10.4fx %11.4f%% %11.4f%%\n",
                      asset_names[i], k_i, sb*100, sb*k_i*100))
        else
          cat(sprintf("  %-20s %10.4fx %11.4f%%  [unchanged]\n",
                      asset_names[i], k_i, sb*100))
      }
      sigma_stressed <- sigma_stressed * outer(volatility_shock, volatility_shock)
      cat("  \u2705 Portfolio covariance scaled\n")
    }

    # ─────────────────────────────────────────────────────────────────────────
    # 7e. Exogenous volatility shock: scale Σ_pe and Σ_ee
    #   Σ_pe_s[i,j] = Σ_pe[i,j] × k_e[j]
    #   Σ_ee_s[i,j] = Σ_ee[i,j] × k_e[i] × k_e[j]
    # ─────────────────────────────────────────────────────────────────────────
    Sigma_pe_use <- Sigma_pe_hist
    Sigma_ee_use <- Sigma_ee_hist

    if(!is.null(exogenous_volatility_shock) && length(exog_vol_shocked_idx) > 0 && has_exogenous) {
      cat("\n--- EXOGENOUS VOLATILITY SHOCK ---\n")
      k_e_vec <- rep(1, n_exog); names(k_e_vec) <- exog_names
      k_e_vec[names(exogenous_volatility_shock)] <- exogenous_volatility_shock
      cat(sprintf("  %-20s %10s\n", "Exogenous Asset", "Multiplier"))
      for(en in exog_names)
        cat(sprintf("  %-20s %10.4fx%s\n", en, k_e_vec[en],
                    if(k_e_vec[en]!=1) "" else "  [unchanged]"))
      Sigma_pe_use <- sweep(Sigma_pe_hist, 2, k_e_vec, `*`)
      Sigma_ee_use <- Sigma_ee_hist * outer(k_e_vec, k_e_vec)
      if(!is.positive.definite(Sigma_ee_use)) {
        cat("  \u26a0\ufe0f Stressed \u03a3_ee not PD — applying nearPD\n")
        Sigma_ee_use <- as.matrix(nearPD(Sigma_ee_use, corr = FALSE)$mat)
      }
      cat("  \u2705 Exogenous covariance blocks scaled\n")
    }

    # ─────────────────────────────────────────────────────────────────────────
    # 7f. Return shock propagation via Schur complement
    #
    # Augmented system: X_aug = [X_port, X_exog]
    #   Σ_aug = [ Σ_pp   Σ_pe ]
    #           [ Σ_ep   Σ_ee ]
    #
    # Shocked block s  = endogenous shocked ∪ exogenous shocked (in aug indices)
    # Non-shocked n    = portfolio assets NOT in s
    #
    # Conditional mean:  μ_n|s = μ_n + Σ_ns Σ_ss⁻¹ δ_s
    # Conditional cov :  Σ_n|s = Σ_nn − Σ_ns Σ_ss⁻¹ Σ_sn
    # ─────────────────────────────────────────────────────────────────────────
    if(propagation) {

      any_endo     <- n_shocked > 0
      any_exog_ret <- n_exog_shocked > 0 && has_exogenous               # exog RETURN shock
      any_exog_vol <- length(exog_vol_shocked_idx) > 0 && has_exogenous # exog VOL shock
      any_exog     <- any_exog_ret || any_exog_vol   # either kind

      # ── Schur mean propagation: only when there are RETURN shocks ─────────────
      if(any_endo || any_exog_ret) {

        cat("\n--- RETURN SHOCK & PROPAGATION ---\n")

        # Build augmented covariance
        Sigma_pe_blk <- if(!is.null(Sigma_pe_use) && n_exog > 0) Sigma_pe_use else
                          matrix(0, n_assets, max(n_exog, 1))
        Sigma_ee_blk <- if(!is.null(Sigma_ee_use) && n_exog > 0) Sigma_ee_use else
                          diag(max(n_exog, 1))

        if(n_exog > 0) {
          Sig_aug <- rbind(cbind(sigma_stressed, Sigma_pe_blk),
                           cbind(t(Sigma_pe_blk), Sigma_ee_blk))
          rownames(Sig_aug) <- colnames(Sig_aug) <- c(asset_names, exog_names)
          mu_aug  <- c(base_mu, exog_means)
        } else {
          Sig_aug <- sigma_stressed
          rownames(Sig_aug) <- colnames(Sig_aug) <- asset_names
          mu_aug  <- base_mu
        }

        # Shocked block: non-empty because any_endo || any_exog_ret is TRUE here
        s_aug_endo <- shocked_assets
        s_aug_exog <- if(any_exog_ret) n_assets + exog_shocked_idx else integer(0)
        s_aug      <- c(s_aug_endo, s_aug_exog)
        n_aug_idx  <- setdiff(seq_len(n_assets), shocked_assets)

        delta_endo <- if(any_endo)     asset_shock[shocked_assets] else numeric(0)
        delta_exog <- if(any_exog_ret) exogenous_shock             else numeric(0)
        delta_s    <- c(delta_endo, delta_exog)

        if(any_endo) {
          cat("\nEndogenous (portfolio) return shocks:\n")
          for(i in shocked_assets)
            cat(sprintf("  %s: %+.4f (%+.2f%%)\n",
                        asset_names[i], asset_shock[i], asset_shock[i]*100))
        }
        if(any_exog_ret) {
          cat("\nExogenous return shocks (zero portfolio weight \u2014 not simulated):\n")
          for(en in names(exogenous_shock))
            cat(sprintf("  %s: %+.4f (%+.2f%%)\n",
                        en, exogenous_shock[en], exogenous_shock[en]*100))
        }

        Sigma_ss <- Sig_aug[s_aug,     s_aug,     drop = FALSE]
        Sigma_ns <- Sig_aug[n_aug_idx, s_aug,     drop = FALSE]
        Sigma_sn <- Sig_aug[s_aug,     n_aug_idx, drop = FALSE]
        Sigma_nn <- Sig_aug[n_aug_idx, n_aug_idx, drop = FALSE]

        Sigma_ss_inv <- if(qr(Sigma_ss)$rank < nrow(Sigma_ss)) {
          cat("  \u26a0\ufe0f Shocked block \u03a3_ss rank-deficient \u2014 using ginv\n")
          MASS::ginv(Sigma_ss)
        } else {
          tryCatch(solve(Sigma_ss), error = function(e) {
            cat(sprintf("  \u26a0\ufe0f solve(\u03a3_ss) failed: %s \u2014 using ginv\n", e$message))
            MASS::ginv(Sigma_ss)
          })
        }

        mu_stressed[shocked_assets] <- base_mu[shocked_assets] + delta_endo
        if(length(n_aug_idx) > 0)
          mu_stressed[n_aug_idx] <- mu_aug[n_aug_idx] +
                                    as.vector(Sigma_ns %*% Sigma_ss_inv %*% delta_s)

        cat("\nConditional mean propagation:\n")
        cat(sprintf("  %-22s  %10s %10s %10s  %s\n",
                    "Asset", "Base \u03bc", "Stressed \u03bc", "\u0394\u03bc", "Type"))
        cat(sprintf("  %s\n", strrep("-", 68)))
        for(i in shocked_assets)
          cat(sprintf("  %-22s  %10.5f %10.5f %10.5f  [ENDO SHOCKED]\n",
                      asset_names[i], base_mu[i], mu_stressed[i],
                      mu_stressed[i] - base_mu[i]))
        for(i in n_aug_idx)
          cat(sprintf("  %-22s  %10.5f %10.5f %10.5f  [propagated]\n",
                      asset_names[i], base_mu[i], mu_stressed[i],
                      mu_stressed[i] - base_mu[i]))
        if(any_exog_ret) {
          cat("  Exogenous (not in portfolio, shown for reference):\n")
          for(en in names(exogenous_shock))
            cat(sprintf("  %-22s  %10.5f %10.5f %10.5f  [EXOG SHOCKED]\n",
                        en, exog_means[en], exog_means[en] + exogenous_shock[en],
                        exogenous_shock[en]))
        }

      } else if(!any_exog_vol) {
        cat("  \u2139\ufe0f No return shocks \u2014 mean propagation skipped\n")
      }

      # ── LTV covariance addon: independent of return shocks ───────────────────────
      # Fires whenever any exog vol shock is set, regardless of whether there is
      # also a return shock. Mean propagation (Schur) and this addon are orthogonal:
      #   Return shock   \u2192 conditional mean shifts; covariance unchanged
      #   Exog vol shock \u2192 marginal portfolio covariance widens via LTV; mean unchanged
      #
      # Formula:  \u03a3_pp_stressed = \u03a3_pp + \u03a3_e [ (k_e\u00b2 \u2212 1) \u00d7 \u03a3_pe[.,e] \u03a3_ee[e,e]\u207b\u00b9 \u03a3_ep[e,.] ]
      # ──────────────────────────────────────────────────────────────────────────────
      if(any_exog_vol && has_exogenous &&
         !is.null(Sigma_pe_hist) && !is.null(Sigma_ee_hist)) {

        k_vec_all <- setNames(rep(1, n_exog), exog_names)
        for(en in names(exogenous_volatility_shock))
          if(en %in% exog_names) k_vec_all[en] <- exogenous_volatility_shock[en]

        exog_vol_addon <- matrix(0, n_assets, n_assets)
        addon_applied  <- FALSE

        for(ei in seq_len(n_exog)) {
          k_e    <- k_vec_all[ei]
          factor <- k_e^2 - 1
          if(factor > 1e-12) {
            addon_applied <- TRUE
            sig_pe_col    <- Sigma_pe_hist[, ei, drop = FALSE]
            var_ee        <- Sigma_ee_hist[ei, ei]
            if(var_ee > 1e-12)
              exog_vol_addon <- exog_vol_addon +
                factor * (sig_pe_col %*% t(sig_pe_col)) / var_ee
          }
        }

        if(addon_applied) {
          sigma_stressed <- sigma_stressed + exog_vol_addon
          sigma_stressed <- (sigma_stressed + t(sigma_stressed)) / 2

          cat("\n--- EXOGENOUS VOLATILITY ADDON (LTV) ---\n")
          for(ei in seq_len(n_exog)) {
            k_e <- k_vec_all[ei]
            if(k_e != 1)
              cat(sprintf("  %s: k_e=%.2f \u2192 (k_e\u00b2\u22121)=%.4f \u00d7 outer product added\n",
                          exog_names[ei], k_e, k_e^2 - 1))
          }
          cat(sprintf("  Asset marginal vols: pre=%s | post=%s\n",
              paste(sprintf("%.3f%%", sqrt(diag(base_sigma)) * 100), collapse=", "),
              paste(sprintf("%.3f%%", sqrt(pmax(diag(sigma_stressed), 0)) * 100),
                    collapse=", ")))
        }
      }

    } else {
      # Propagation disabled
      if(n_shocked > 0) {
        cat("\nPropagation DISABLED — portfolio means shifted only\n")
        mu_stressed[shocked_assets] <- base_mu[shocked_assets] + asset_shock[shocked_assets]
      }
      if(n_exog_shocked > 0)
        cat("  \u2139\ufe0f Propagation DISABLED: exogenous shocks do not affect portfolio\n")
    }

    # 7g. Final validation
    sigma_stressed <- (sigma_stressed + t(sigma_stressed)) / 2
    rownames(sigma_stressed) <- colnames(sigma_stressed) <- asset_names
    if(!is.positive.definite(sigma_stressed)) {
      cat("  \u26a0\ufe0f sigma_stressed not PD — applying nearPD\n")
      sigma_stressed <- as.matrix(nearPD(sigma_stressed, corr = FALSE, ensureSymmetry = TRUE)$mat)
      rownames(sigma_stressed) <- colnames(sigma_stressed) <- asset_names
    }
    names(mu_stressed) <- asset_names

    # 7h. Summary
    cat("\n========================================================================\n")
    cat("STRESSED PARAMETERS SUMMARY\n")
    cat("========================================================================\n")
    cat(sprintf("\n  %-22s %12s %12s %12s\n", "Asset", "\u03bc_base", "\u03bc_stressed", "\u0394\u03bc"))
    cat(sprintf("  %s\n", strrep("-", 62)))
    for(i in seq_len(n_assets)) {
      is_endo_s <- i %in% shocked_assets
      is_prop   <- !is_endo_s && abs(mu_stressed[i] - base_mu[i]) > 1e-10
      flag <- if(is_endo_s) "[ENDO SHOCKED]" else if(is_prop) "[propagated]" else ""
      cat(sprintf("  %-22s %12.5f %12.5f %12.5f  %s\n",
                  asset_names[i], base_mu[i], mu_stressed[i],
                  mu_stressed[i] - base_mu[i], flag))
    }
    cat(sprintf("\n  %-22s %12s %12s %12s\n", "Asset", "\u03c3_base", "\u03c3_stressed", "\u0394\u03c3"))
    cat(sprintf("  %s\n", strrep("-", 62)))
    for(i in seq_len(n_assets)) {
      sb   <- sqrt(base_sigma[i, i]); ss <- sqrt(sigma_stressed[i, i])
      flag <- if(i %in% vol_shocked_assets) "[VOL SHOCKED]"
              else if(abs(ss - sb) > 1e-8) "[propagated]" else ""
      cat(sprintf("  %-22s %12.4f %12.4f %12.4f  %s\n",
                  asset_names[i], sb*100, ss*100, (ss-sb)*100, flag))
    }
    if(n_exog_shocked > 0 || length(exog_vol_shocked_idx) > 0) {
      cat("\n  Exogenous shock applied (zero weight — not simulated):\n")
      all_exog_t <- unique(c(exog_names[exog_shocked_idx], exog_names[exog_vol_shocked_idx]))
      for(en in all_exog_t) {
        rd  <- if(!is.null(exogenous_shock) && en %in% names(exogenous_shock))
                 exogenous_shock[en] else 0
        vm  <- if(!is.null(exogenous_volatility_shock) &&
                  en %in% names(exogenous_volatility_shock))
                 exogenous_volatility_shock[en] else 1
        cat(sprintf("    %-20s  \u0394ret=%+.4f  vol_mult=%.4fx\n", en, rd, vm))
      }
    }
    cat(sprintf("\n  sigma_stressed condition number : %.2f\n", kappa(sigma_stressed)))
    cat(sprintf("  sigma_stressed min eigenvalue   : %.2e\n",
                min(eigen(sigma_stressed, only.values = TRUE)$values)))

  }  # end stress testing block

  # ===========================================================================
  # 8. SIMULATE RETURNS
  # ===========================================================================

  cat("\n========================================================================\n")
  cat("SIMULATING PORTFOLIO RETURNS\n")
  cat("========================================================================\n")

  # 8a. Resolve and validate parameters
  use_mu_pre    <- base_mu
  use_sigma_pre <- (base_sigma + t(base_sigma)) / 2
  if(!is.positive.definite(use_sigma_pre)) {
    cat("  ⚠️ base_sigma not PD — applying nearPD\n")
    use_sigma_pre <- as.matrix(nearPD(use_sigma_pre, corr = FALSE, ensureSymmetry = TRUE)$mat)
  }
  rownames(use_sigma_pre) <- colnames(use_sigma_pre) <- asset_names

  if(is_stressed) {
    use_mu_post    <- mu_stressed
    use_sigma_post <- (sigma_stressed + t(sigma_stressed)) / 2
    if(!is.positive.definite(use_sigma_post)) {
      cat("  ⚠️ sigma_stressed not PD — applying nearPD\n")
      use_sigma_post <- as.matrix(nearPD(use_sigma_post, corr = FALSE, ensureSymmetry = TRUE)$mat)
    }
    rownames(use_sigma_post) <- colnames(use_sigma_post) <- asset_names
  } else {
    use_mu_post    <- use_mu_pre
    use_sigma_post <- use_sigma_pre
  }

  # rmvt scale matrices
  if(sim_returns_dist == "rmvt") {
    df_sim <- if(exists("rmvt_df") && is.finite(rmvt_df)) rmvt_df
              else if(exists("df_joint") && is.finite(df_joint)) df_joint
              else df
    if(df_sim <= 2) stop(sprintf("df_sim = %.4f must be > 2", df_sim))

    scale_pre <- if(exists("rmvt_scale") && !is.null(rmvt_scale) && !copula) {
      cat(sprintf("  rmvt: using pre-computed scale matrix (df=%.4f)\n", df_sim))
      rmvt_scale
    } else {
      cat(sprintf("  rmvt: deflating covariance to scale matrix (df=%.4f)\n", df_sim))
      use_sigma_pre * (df_sim - 2) / df_sim
    }
    scale_pre <- (scale_pre + t(scale_pre)) / 2
    if(!is.positive.definite(scale_pre))
      scale_pre <- as.matrix(nearPD(scale_pre, corr = FALSE, ensureSymmetry = TRUE)$mat)

    scale_post <- use_sigma_post * (df_sim - 2) / df_sim
    scale_post <- (scale_post + t(scale_post)) / 2
    if(!is.positive.definite(scale_post))
      scale_post <- as.matrix(nearPD(scale_post, corr = FALSE, ensureSymmetry = TRUE)$mat)

    cat(sprintf("  rmvt df: %.4f\n", df_sim))
  }

  # rsgt copula df
  if(sim_returns_dist == "rsgt") {
    copula_df_sim <- if(exists("copula_df_fitted") && is.finite(copula_df_fitted))
                       copula_df_fitted else copula_df
    cat(sprintf("  rsgt copula df: %.4f\n", copula_df_sim))
  }

  # 8b. Packages
  if(sim_returns_dist == "rsgt") {
    if(!requireNamespace("sgt", quietly = TRUE)) install.packages("sgt")
    library(sgt)
  }
  if(sim_returns_dist == "rmvt") {
    if(!requireNamespace("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")
    library(mvtnorm)
  }

  # 8c. Pre-allocate
  portfolio_returns_pre  <- matrix(NA_real_, nrow = n_weight_sims, ncol = n_sim_returns)
  portfolio_returns_post <- matrix(NA_real_, nrow = n_weight_sims, ncol = n_sim_returns)

  # 8d. Simulation helpers
  simulate_mvnorm <- function(n_sim, mu, Sigma)
    mvrnorm(n_sim, mu = mu, Sigma = Sigma)

  simulate_mvt <- function(n_sim, delta, sigma, df)
    rmvt(n_sim, delta = delta, sigma = sigma, df = df)

  simulate_rsgt <- function(n_sim, mu_vec, sigma_mat, ref_sigma_mat, fp_list,
                             lam_vec, p_vec, q_vec,
                             copula_type_sim, copula_df_sim, n_assets) {
    # sigma_mat     : current (possibly stressed) covariance matrix
    # ref_sigma_mat : reference (pre-stress) covariance matrix — used to compute
    #                 vol scaling ratios relative to the fitted SGT parameters
    #
    # WHY vol scaling is needed:
    #   cov2cor(sigma_mat) discards marginal volatility information.
    #   Without explicit scaling, all SGT marginals would use fp$sd regardless
    #   of any volatility shock applied to sigma_mat. The vol ratio corrects this:
    #     sigma_j* = fp$sd × sqrt(sigma_mat[j,j]) / sqrt(ref_sigma_mat[j,j])
    #   This scales the SGT width proportionally to the stressed marginal vol
    #   while preserving the fitted shape parameters (lambda, p, q).

    cor_mat <- cov2cor(sigma_mat)
    cor_mat <- (cor_mat + t(cor_mat)) / 2
    if(!is.positive.definite(cor_mat))
      cor_mat <- as.matrix(nearPD(cor_mat, corr = TRUE)$mat)

    Z <- mvrnorm(n_sim, mu = rep(0, n_assets), Sigma = cor_mat)
    U <- if(copula_type_sim == "t") pt(Z, df = copula_df_sim) else pnorm(Z)

    ref_vols     <- sqrt(pmax(diag(ref_sigma_mat), 1e-12))
    current_vols <- sqrt(pmax(diag(sigma_mat),     1e-12))

    sim_mat <- matrix(NA_real_, nrow = n_sim, ncol = n_assets)
    for(j in seq_len(n_assets)) {
      vol_ratio_j    <- current_vols[j] / ref_vols[j]
      sigma_j_scaled <- fp_list[[j]]$sd * vol_ratio_j
      mean_shift     <- mu_vec[j] - fp_list[[j]]$mean
      sim_mat[, j]   <- qsgt(U[, j],
                              mu     = fp_list[[j]]$mean,
                              sigma  = sigma_j_scaled,
                              lambda = lam_vec[j],
                              p      = p_vec[j],
                              q      = q_vec[j]) + mean_shift
    }
    return(sim_mat)
  }

  # 8e. SGT shape parameters
  if(sim_returns_dist == "rsgt") {
    lam_sim <- sapply(seq_len(n_assets), function(j)
      if(!is.null(fitted_params[[j]]$lambda)) fitted_params[[j]]$lambda
      else if(!is.null(asset_lambdas)) asset_lambdas[j] else -0.3)
    p_sim <- sapply(seq_len(n_assets), function(j)
      if(!is.null(fitted_params[[j]]$p)) fitted_params[[j]]$p
      else if(!is.null(asset_p)) asset_p[j] else 4)
    q_sim <- sapply(seq_len(n_assets), function(j)
      if(!is.null(fitted_params[[j]]$q)) fitted_params[[j]]$q
      else if(!is.null(asset_q)) asset_q[j] else 4)
    cat("  SGT shape parameters:\n")
    for(j in seq_len(n_assets))
      cat(sprintf("    %s: λ=%.4f, p=%.4f, q=%.4f\n",
                  asset_names[j], lam_sim[j], p_sim[j], q_sim[j]))
  }

  # 8f. Main simulation loop
  cat(sprintf("\n  Simulating %s scenarios across %d weight set(s)...\n",
              format(n_sim_returns, big.mark = ","), n_weight_sims))
  sim_errors <- 0L
  asset_draws_pre  <- NULL   # n_sim_returns x n_assets (first weight set)
  asset_draws_post <- NULL   # n_sim_returns x n_assets (first weight set, stressed)

  for(i in seq_len(n_weight_sims)) {

    if(n_weight_sims > 1 && (i %% max(1L, n_weight_sims %/% 10L) == 0L || i == n_weight_sims))
      cat(sprintf("  Progress: %d / %d (%.0f%%)\n", i, n_weight_sims, 100*i/n_weight_sims))

    # Pre-stress
    sim_returns_pre <- tryCatch({
      if(sim_returns_dist == "rmvnorm") {
        simulate_mvnorm(n_sim_returns, mu = use_mu_pre, Sigma = use_sigma_pre)
      } else if(sim_returns_dist == "rmvt") {
        simulate_mvt(n_sim_returns, delta = use_mu_pre, sigma = scale_pre, df = df_sim)
      } else {
        simulate_rsgt(n_sim_returns, use_mu_pre, use_sigma_pre,
                      use_sigma_pre,
                      fitted_params, lam_sim, p_sim, q_sim,
                      copula_type, copula_df_sim, n_assets)
      }
    }, error = function(e) {
      sim_errors <<- sim_errors + 1L
      if(sim_errors <= 5L)
        cat(sprintf("  ⚠️ Pre-stress error (set %d): %s\n", i, e$message))
      NULL
    })

    if(is.null(sim_returns_pre)) next
    if(any(!is.finite(sim_returns_pre)))
      sim_returns_pre[!is.finite(sim_returns_pre)] <- NA_real_
    if(i == 1L) {
      asset_draws_pre <- as.data.frame(sim_returns_pre)
      colnames(asset_draws_pre) <- asset_names
    }

    # Post-stress
    sim_returns_post <- NULL
    if(is_stressed) {
      sim_returns_post <- tryCatch({
        if(sim_returns_dist == "rmvnorm") {
          simulate_mvnorm(n_sim_returns, mu = use_mu_post, Sigma = use_sigma_post)
        } else if(sim_returns_dist == "rmvt") {
          simulate_mvt(n_sim_returns, delta = use_mu_post, sigma = scale_post, df = df_sim)
        } else {
          # ref_sigma_mat = use_sigma_pre: vol ratio captures the stress vol change
          simulate_rsgt(n_sim_returns, use_mu_post, use_sigma_post,
                        use_sigma_pre,
                        fitted_params, lam_sim, p_sim, q_sim,
                        copula_type, copula_df_sim, n_assets)
        }
      }, error = function(e) {
        sim_errors <<- sim_errors + 1L
        if(sim_errors <= 5L)
          cat(sprintf("  ⚠️ Post-stress error (set %d): %s\n", i, e$message))
        NULL
      })
      if(!is.null(sim_returns_post) && any(!is.finite(sim_returns_post)))
        sim_returns_post[!is.finite(sim_returns_post)] <- NA_real_
      if(i == 1L && !is.null(sim_returns_post)) {
        asset_draws_post <- as.data.frame(sim_returns_post)
        colnames(asset_draws_post) <- asset_names
      }
    }

    w_current <- weights_matrix[i, ]
    portfolio_returns_pre[i, ] <- as.vector(sim_returns_pre %*% w_current)
    if(is_stressed && !is.null(sim_returns_post))
      portfolio_returns_post[i, ] <- as.vector(sim_returns_post %*% w_current)

  }  # end weight loop

  if(sim_errors > 0)
    warning(sprintf("%d simulation error(s) — affected rows excluded.", sim_errors))

  # 8g. Flatten and validate
  portfolio_returns_pre_flat <- as.vector(portfolio_returns_pre)
  portfolio_returns_pre_flat <- portfolio_returns_pre_flat[is.finite(portfolio_returns_pre_flat)]
  if(length(portfolio_returns_pre_flat) == 0)
    stop("No valid pre-stress portfolio returns generated.")

  n_dropped_pre <- n_weight_sims * n_sim_returns - length(portfolio_returns_pre_flat)
  if(n_dropped_pre > 0)
    warning(sprintf("%d pre-stress values excluded (%.2f%%)", n_dropped_pre,
                    100 * n_dropped_pre / (n_weight_sims * n_sim_returns)))

  if(is_stressed) {
    portfolio_returns_post_flat <- as.vector(portfolio_returns_post)
    portfolio_returns_post_flat <- portfolio_returns_post_flat[is.finite(portfolio_returns_post_flat)]
    if(length(portfolio_returns_post_flat) == 0)
      stop("No valid post-stress portfolio returns generated.")
  }

  # 8h. Simulation summary
  cat("\n========================================================================\n")
  cat("SIMULATION SUMMARY\n")
  cat("========================================================================\n")
  cat(sprintf("  Distribution   : %s\n", sim_returns_dist))
  if(sim_returns_dist == "rmvt")
    cat(sprintf("  Simulation df  : %.4f\n", df_sim))
  if(sim_returns_dist == "rsgt")
    cat(sprintf("  Copula df      : %.4f [%s]\n", copula_df_sim, toupper(copula_type)))
  cat(sprintf("  Weight sets    : %d\n", n_weight_sims))
  cat(sprintf("  Sims per set   : %s\n", format(n_sim_returns, big.mark = ",")))
  cat(sprintf("  Valid pre-obs  : %s / %s\n",
              format(length(portfolio_returns_pre_flat), big.mark = ","),
              format(n_weight_sims * n_sim_returns, big.mark = ",")))
  if(is_stressed)
    cat(sprintf("  Valid post-obs : %s / %s\n",
                format(length(portfolio_returns_post_flat), big.mark = ","),
                format(n_weight_sims * n_sim_returns, big.mark = ",")))
  if(sim_errors > 0)
    cat(sprintf("  ⚠️ Errors      : %d weight set(s) failed\n", sim_errors))

  # Sanity check (single portfolio)
  if(n_weight_sims == 1) {
    w          <- weights_matrix[1, ]
    target_vol <- sqrt(as.numeric(t(w) %*% use_sigma_pre %*% w)) * 100
    sim_vol    <- sd(portfolio_returns_pre_flat) * 100
    target_ret <- sum(w * use_mu_pre) * 100
    sim_ret    <- mean(portfolio_returns_pre_flat) * 100
    cat(sprintf("\n  Sanity check:\n"))
    cat(sprintf("    %-20s  Target: %8.4f%%  Simulated: %8.4f%%\n",
                "Portfolio μ", target_ret, sim_ret))
    cat(sprintf("    %-20s  Target: %8.4f%%  Simulated: %8.4f%%\n",
                "Portfolio σ", target_vol, sim_vol))
    if(abs(sim_vol / target_vol - 1) > 0.1)
      cat(sprintf("    ⚠️ Vol ratio = %.2f — check scale matrix\n", sim_vol / target_vol))
    else
      cat("    ✅ Simulated vol consistent with target\n")
  }

  # ===========================================================================
  # 9. CALCULATE STATISTICS
  # ===========================================================================

  # 9a. Risk statistics helper
  compute_risk_stats <- function(returns_flat, label, rf_monthly = 0) {
    n <- length(returns_flat)
    if(n == 0) stop(sprintf("No returns for '%s'", label))
    if(any(!is.finite(returns_flat)))
      stop(sprintf("%d non-finite values for '%s'", sum(!is.finite(returns_flat)), label))

    ret_mean   <- mean(returns_flat)
    ret_median <- median(returns_flat)
    ret_vol    <- sd(returns_flat)
    ret_skew   <- skewness(returns_flat)

    # Excess kurtosis: compute explicitly to avoid convention ambiguity
    xbar            <- mean(returns_flat)
    m4              <- mean((returns_flat - xbar)^4)
    m2              <- mean((returns_flat - xbar)^2)
    excess_kurt_explicit <- (m4 / m2^2) - 3

    kurt_pkg <- tryCatch(kurtosis(returns_flat), error = function(e) NA)
    ret_excess_kurt <- if(!is.na(kurt_pkg)) {
      if(abs(kurt_pkg - excess_kurt_explicit) < 0.5) kurt_pkg else kurt_pkg - 3
    } else excess_kurt_explicit

    var_90 <- quantile(returns_flat, 0.10, names = FALSE)
    var_95 <- quantile(returns_flat, 0.05, names = FALSE)
    var_99 <- quantile(returns_flat, 0.01, names = FALSE)

    compute_cvar <- function(rets, thr, min_obs = 10) {
      tail_s <- rets[rets < thr]
      if(length(tail_s) >= min_obs) mean(tail_s)
      else { tail_i <- rets[rets <= thr]; if(length(tail_i) == 0) thr else mean(tail_i) }
    }

    cvar_90 <- compute_cvar(returns_flat, var_90)
    cvar_95 <- compute_cvar(returns_flat, var_95)
    cvar_99 <- compute_cvar(returns_flat, var_99)
    cvar_ok <- (cvar_90 <= var_90) && (cvar_95 <= var_95) && (cvar_99 <= var_99)

    sharpe_monthly <- (ret_mean - rf_monthly) / ret_vol
    sharpe_annual  <- sharpe_monthly * sqrt(periods_per_year)
    annual_return  <- (1 + ret_mean)^periods_per_year - 1
    annual_vol     <- ret_vol * sqrt(periods_per_year)

    list(
      n = n, label = label,
      mean = ret_mean, median = ret_median, vol = ret_vol,
      skewness = ret_skew, excess_kurtosis = ret_excess_kurt,
      var_90 = var_90, var_95 = var_95, var_99 = var_99,
      cvar_90 = cvar_90, cvar_95 = cvar_95, cvar_99 = cvar_99, cvar_ok = cvar_ok,
      sharpe_monthly = sharpe_monthly, sharpe_annual = sharpe_annual,
      annual_return = annual_return, annual_vol = annual_vol,
      pct_positive = mean(returns_flat > 0) * 100,
      pct_negative = mean(returns_flat < 0) * 100,
      max_loss = min(returns_flat), best_return = max(returns_flat)
    )
  }

  # 9b. Compute
  if(length(portfolio_returns_pre_flat) == 0) stop("No pre-stress returns")
  stats_pre <- compute_risk_stats(portfolio_returns_pre_flat, "Pre-Stress")
  if(is_stressed) {
    if(length(portfolio_returns_post_flat) == 0) stop("No post-stress returns")
    stats_post <- compute_risk_stats(portfolio_returns_post_flat, "Post-Stress")
  }

  # Scalar extracts for backward compatibility
  portfolio_mean_pre   <- stats_pre$mean   * 100
  portfolio_median_pre <- stats_pre$median * 100
  portfolio_vol_pre    <- stats_pre$vol    * 100
  portfolio_skew_pre   <- stats_pre$skewness
  portfolio_kurt_pre   <- stats_pre$excess_kurtosis + 3
  var_90_pre  <- stats_pre$var_90  * 100
  var_95_pre  <- stats_pre$var_95  * 100
  var_99_pre  <- stats_pre$var_99  * 100
  cvar_95_pre <- stats_pre$cvar_95 * 100

  if(is_stressed) {
    portfolio_mean_post   <- stats_post$mean   * 100
    portfolio_median_post <- stats_post$median * 100
    portfolio_vol_post    <- stats_post$vol    * 100
    portfolio_skew_post   <- stats_post$skewness
    portfolio_kurt_post   <- stats_post$excess_kurtosis + 3
    var_90_post  <- stats_post$var_90  * 100
    var_95_post  <- stats_post$var_95  * 100
    var_99_post  <- stats_post$var_99  * 100
    cvar_95_post <- stats_post$cvar_95 * 100
  }

  # 9c. Print
  print_stats <- function(s) {
    cat(sprintf("\n  Sample size: %s\n", format(s$n, big.mark = ",")))
    cat(sprintf("\n  RETURN (%s)\n", data_freq_label))
    cat(sprintf("  %-30s : %+.4f%%\n", "Expected Return (mean)",   s$mean   * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "Expected Return (median)", s$median * 100))
    if(abs(s$mean - s$median) * 100 > 0.05)
      cat(sprintf("  %-30s : %+.4f%%  [skew indicator]\n",
                  "Mean − Median", (s$mean - s$median) * 100))
    cat("\n  RETURN (Annualised)\n")
    cat(sprintf("  %-30s : %+.4f%%\n", "Annual Return (geometric)", s$annual_return * 100))
    cat(sprintf("  %-30s : %.4f%%\n",  "Annual Volatility",         s$annual_vol    * 100))
    cat(sprintf("  %-30s : %.4f\n",    "Sharpe Ratio (annual)",     s$sharpe_annual))
    cat(sprintf("\n  RISK (%s)\n", data_freq_label))
    cat(sprintf("  %-30s : %.4f%%\n", "Volatility",       s$vol * 100))
    cat(sprintf("  %-30s : %.4f\n",   "Skewness",         s$skewness))
    cat(sprintf("  %-30s : %.4f\n",   "Excess Kurtosis",  s$excess_kurtosis))
    if(abs(s$skewness) > 1)     cat("  ⚠️ |Skewness| > 1\n")
    if(s$excess_kurtosis > 3)   cat("  ⚠️ Excess kurtosis > 3 — heavy tails\n")
    cat(sprintf("\n  VALUE-AT-RISK (%s, negative = loss)\n", data_freq_label))
    cat(sprintf("  %-30s : %+.4f%%\n", "90% VaR (10th pct)",  s$var_90 * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "95% VaR (5th pct)",   s$var_95 * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "99% VaR (1st pct)",   s$var_99 * 100))
    cat(sprintf("\n  EXPECTED SHORTFALL (%s)\n", data_freq_label))
    cat(sprintf("  %-30s : %+.4f%%\n", "90% CVaR", s$cvar_90 * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "95% CVaR", s$cvar_95 * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "99% CVaR", s$cvar_99 * 100))
    if(!s$cvar_ok) cat("  ⚠️ CVaR > VaR — too few tail observations\n")
    cat("\n  TAIL DIAGNOSTICS\n")
    cat(sprintf("  %-30s : %.2f%%\n",  "% Positive months",  s$pct_positive))
    cat(sprintf("  %-30s : %.2f%%\n",  "% Negative months",  s$pct_negative))
    cat(sprintf("  %-30s : %+.4f%%\n", "Worst single return", s$max_loss    * 100))
    cat(sprintf("  %-30s : %+.4f%%\n", "Best single return",  s$best_return * 100))
  }

  cat("\n========================================================================\n")
  cat(sprintf("PRE-STRESS PORTFOLIO STATISTICS (%s Returns)\n", data_freq_label))
  cat("========================================================================\n")
  print_stats(stats_pre)

  if(is_stressed) {
    cat("\n========================================================================\n")
    cat(sprintf("POST-STRESS PORTFOLIO STATISTICS (%s Returns)\n", data_freq_label))
    cat("========================================================================\n")
    print_stats(stats_post)

    cat("\n========================================================================\n")
    cat("STRESS IMPACT SUMMARY\n")
    cat("========================================================================\n")
    delta_mean   <- (stats_post$mean    - stats_pre$mean)    * 100
    delta_vol    <- (stats_post$vol     - stats_pre$vol)     * 100
    delta_skew   <-  stats_post$skewness - stats_pre$skewness
    delta_kurt   <-  stats_post$excess_kurtosis - stats_pre$excess_kurtosis
    delta_var95  <- (stats_post$var_95  - stats_pre$var_95)  * 100
    delta_var99  <- (stats_post$var_99  - stats_pre$var_99)  * 100
    delta_cvar95 <- (stats_post$cvar_95 - stats_pre$cvar_95) * 100

    cat(sprintf("\n  %-30s %12s %12s %12s\n", "Metric", "Pre-Stress", "Post-Stress", "Change"))
    cat(sprintf("  %s\n", strrep("-", 70)))
    cat(sprintf("  %-30s %11.4f%% %11.4f%% %+11.4f%%\n", "Expected Return",
                stats_pre$mean*100, stats_post$mean*100, delta_mean))
    cat(sprintf("  %-30s %11.4f%% %11.4f%% %+11.4f%%\n", "Volatility",
                stats_pre$vol*100, stats_post$vol*100, delta_vol))
    cat(sprintf("  %-30s %12.4f  %12.4f  %+12.4f\n",  "Skewness",
                stats_pre$skewness, stats_post$skewness, delta_skew))
    cat(sprintf("  %-30s %12.4f  %12.4f  %+12.4f\n",  "Excess Kurtosis",
                stats_pre$excess_kurtosis, stats_post$excess_kurtosis, delta_kurt))
    cat(sprintf("  %-30s %11.4f%% %11.4f%% %+11.4f%%\n", "95% VaR",
                stats_pre$var_95*100, stats_post$var_95*100, delta_var95))
    cat(sprintf("  %-30s %11.4f%% %11.4f%% %+11.4f%%\n", "99% VaR",
                stats_pre$var_99*100, stats_post$var_99*100, delta_var99))
    cat(sprintf("  %-30s %11.4f%% %11.4f%% %+11.4f%%\n", "95% CVaR",
                stats_pre$cvar_95*100, stats_post$cvar_95*100, delta_cvar95))

    if(delta_var95  < -1) cat(sprintf("\n  Warning: VaR95 worsened by %.2f%%\n",  abs(delta_var95)))
    if(delta_cvar95 < -1) cat(sprintf("  Warning: CVaR95 worsened by %.2f%%\n", abs(delta_cvar95)))
    if(delta_vol    >  1) cat(sprintf("  Warning: Vol increased by %.2f%%\n",    delta_vol))

    stress_impact_summary <- data.frame(
      Metric      = c("Expected Return (%)", "Volatility (%)", "Skewness",
                      "Excess Kurtosis", "VaR 95% (%)", "VaR 99% (%)", "CVaR 95% (%)"),
      Pre_Stress  = c(stats_pre$mean*100, stats_pre$vol*100, stats_pre$skewness,
                      stats_pre$excess_kurtosis, stats_pre$var_95*100,
                      stats_pre$var_99*100, stats_pre$cvar_95*100),
      Post_Stress = c(stats_post$mean*100, stats_post$vol*100, stats_post$skewness,
                      stats_post$excess_kurtosis, stats_post$var_95*100,
                      stats_post$var_99*100, stats_post$cvar_95*100),
      Change      = c(delta_mean, delta_vol, delta_skew, delta_kurt,
                      delta_var95, delta_var99, delta_cvar95),
      stringsAsFactors = FALSE
    )
  } else {
    stress_impact_summary <- NULL
  }

  # ===========================================================================
  # 10. PORTFOLIO DIAGNOSTIC PLOTS
  # ===========================================================================

  cat("\n========================================================================\n")
  cat("GENERATING PORTFOLIO DIAGNOSTIC PLOTS\n")
  cat("========================================================================\n")

  portfolio_plots <- list()

  theme_risk <- theme_minimal() +
    theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      axis.text.x   = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  # Plot 1: Asset distributions with fitted overlay
  compute_fitted_density <- function(x_vals, fp, dist) {
    tryCatch({
      if(dist == "rmvnorm") dnorm(x_vals, mean = fp$mean, sd = fp$sd)
      else if(dist == "rmvt") dt((x_vals - fp$mean) / fp$sd, df = fp$df) / fp$sd
      else dsgt(x_vals, mu = fp$mean, sigma = fp$sd,
                lambda = fp$lambda, p = fp$p, q = fp$q)
    }, error = function(e) rep(NA_real_, length(x_vals)))
  }

  asset_df <- do.call(rbind, lapply(seq_len(n_assets), function(i)
    data.frame(Return = returns_numeric[, i] * 100, Asset = asset_names[i])))

  has_fitted_overlay <- !is.null(fitted_params) && length(fitted_params) == n_assets

  overlay_df <- if(has_fitted_overlay) {
    do.call(rbind, Filter(Negate(is.null), lapply(seq_len(n_assets), function(i) {
      fp     <- fitted_params[[i]]
      x_raw  <- returns_numeric[, i]
      x_vals <- seq(min(x_raw) * 1.1, max(x_raw) * 1.1, length.out = 300)
      dens   <- compute_fitted_density(x_vals, fp, sim_returns_dist)
      if(all(is.na(dens))) return(NULL)
      # Density is computed on DECIMAL scale. The histogram x-axis is in PERCENT.
      # By the change-of-variables formula: f_Y(y) = f_X(y/100) / 100
      # So we must divide by 100 to match the histogram's density scale.
      data.frame(x = x_vals * 100, y = dens / 100, Asset = asset_names[i])
    })))
  } else NULL

  p_asset <- ggplot(asset_df, aes(x = Return)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, fill = "lightblue", alpha = 0.6, color = "white") +
    facet_wrap(~Asset, scales = "free") +
    labs(title    = sprintf("Asset %s Return Distributions", data_freq_label),
         subtitle = if(has_fitted_overlay)
           sprintf("Blue: Empirical | Red: Fitted %s", toupper(sim_returns_dist))
           else "Empirical distribution",
         x = sprintf("%s Return (%%)", data_freq_label), y = "Density") +
    theme_risk

  if(!is.null(overlay_df))
    p_asset <- p_asset +
      geom_line(data = overlay_df, aes(x = x, y = y),
                color = "red", linewidth = 0.9, alpha = 0.85, inherit.aes = FALSE)

  portfolio_plots$asset_distributions <- p_asset
  cat("  ✅ Plot 1: Asset distributions\n")

  # Plot 2: Portfolio distribution
  port_mean_pct  <- stats_pre$mean   * 100
  port_vol_pct   <- stats_pre$vol    * 100
  var95_pct      <- stats_pre$var_95 * 100
  var99_pct      <- stats_pre$var_99 * 100
  cvar95_pct     <- stats_pre$cvar_95 * 100
  portfolio_df_pre <- data.frame(Return = portfolio_returns_pre_flat * 100)

  p_port_dist <- ggplot(portfolio_df_pre, aes(x = Return)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60, fill = "lightblue", alpha = 0.6, color = "white") +
    geom_density(color = "darkblue", linewidth = 0.9) +
    geom_vline(xintercept = port_mean_pct,              color = "blue",        linewidth = 1,   linetype = "solid") +
    geom_vline(xintercept = port_mean_pct - port_vol_pct, color = "forestgreen", linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = port_mean_pct + port_vol_pct, color = "forestgreen", linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = var95_pct,                  color = "orange",      linewidth = 0.9, linetype = "dotted") +
    geom_vline(xintercept = var99_pct,                  color = "red",         linewidth = 0.9, linetype = "dotted") +
    annotate("text", x = port_mean_pct, y = Inf,
             label = sprintf("μ=%.2f%%", port_mean_pct),
             vjust = 2, hjust = -0.1, size = 3, color = "blue") +
    annotate("text", x = var95_pct, y = Inf,
             label = sprintf("VaR95\n%.2f%%", var95_pct),
             vjust = 2, hjust = 1.1, size = 2.5, color = "orange") +
    annotate("text", x = var99_pct, y = Inf,
             label = sprintf("VaR99\n%.2f%%", var99_pct),
             vjust = 2, hjust = 1.1, size = 2.5, color = "red") +
    labs(title    = sprintf("Portfolio %s Return Distribution (Pre-Stress)", data_freq_label),
         subtitle = sprintf("Mean:%.3f%%  Vol:%.3f%%  VaR95:%.3f%%  CVaR95:%.3f%%  Skew:%.3f  ExKurt:%.3f",
                            port_mean_pct, port_vol_pct, var95_pct, cvar95_pct,
                            stats_pre$skewness, stats_pre$excess_kurtosis),
         x = sprintf("%s Return (%%)", data_freq_label), y = "Density") +
    theme_risk

  portfolio_plots$portfolio_distribution <- p_port_dist
  cat("  ✅ Plot 2: Portfolio distribution\n")

  # Plot 3: QQ plot
  portfolio_plots$qq_plot <- ggplot(portfolio_df_pre, aes(sample = Return)) +
    stat_qq(alpha = 0.4, size = 0.8, color = "steelblue") +
    stat_qq_line(color = "red", linewidth = 0.9, linetype = "dashed") +
    labs(title    = "Normal QQ Plot — Portfolio Returns (Pre-Stress)",
         subtitle = sprintf("Simulated via %s | Departure = non-normality", toupper(sim_returns_dist)),
         x = "Standard Normal Quantiles", y = "Sample Quantiles (%)") +
    theme_risk
  cat("  ✅ Plot 3: QQ plot\n")

  # Plot 4: Historical correlation heatmap
  build_cor_heatmap <- function(cor_mat, title_str) {
    cor_melt <- reshape2::melt(cor_mat)
    colnames(cor_melt) <- c("Var1", "Var2", "value")
    cor_melt$label <- ifelse(cor_melt$Var1 == cor_melt$Var2, "1.00",
                              sprintf("%.2f", cor_melt$value))
    ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2166ac", high = "#d6604d", mid = "white",
                           midpoint = 0, limits = c(-1, 1), name = "ρ") +
      geom_text(aes(label = label), size = 3, fontface = "bold") +
      labs(title = title_str, x = "", y = "") +
      theme_risk + theme(panel.grid = element_blank())
  }

  portfolio_plots$historical_correlation <- build_cor_heatmap(historical_cor, "Historical Correlation Matrix")
  cat("  ✅ Plot 4: Historical correlation\n")

  # Plot 5: Copula correlation + difference
  if(copula && copula_success && !is.null(copula_cor)) {
    portfolio_plots$copula_correlation <- build_cor_heatmap(
      copula_cor, sprintf("%s Copula Correlation", toupper(copula_type)))

    diff_cor  <- copula_cor - historical_cor
    diff_melt <- reshape2::melt(diff_cor)
    colnames(diff_melt) <- c("Var1", "Var2", "value")
    diff_lim <- max(abs(diff_melt$value), na.rm = TRUE)

    portfolio_plots$copula_vs_historical <- ggplot(diff_melt, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2166ac", high = "#d6604d", mid = "white",
                           midpoint = 0, limits = c(-diff_lim, diff_lim), name = "Δρ") +
      geom_text(aes(label = sprintf("%+.3f", value)), size = 2.8) +
      labs(title    = "Copula − Historical Correlation",
           subtitle = "Positive = copula captures stronger dependence",
           x = "", y = "") +
      theme_risk + theme(panel.grid = element_blank())

    cat("  ✅ Plot 5: Copula correlation + difference\n")
  }

  # Plot 6: Covariance matrix
  sigma_melt <- reshape2::melt(base_sigma)
  colnames(sigma_melt) <- c("Var1", "Var2", "value")
  vol_vec <- sqrt(diag(base_sigma))
  sigma_melt$label <- apply(sigma_melt, 1, function(row) {
    i_name <- as.character(row["Var1"]); j_name <- as.character(row["Var2"])
    val    <- as.numeric(row["value"])
    if(i_name == j_name) sprintf("σ=%.3f%%", sqrt(val) * 100)
    else {
      i_idx <- which(asset_names == i_name); j_idx <- which(asset_names == j_name)
      sprintf("ρ=%.3f", val / (vol_vec[i_idx] * vol_vec[j_idx]))
    }
  })

  portfolio_plots$covariance_matrix <- ggplot(sigma_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#ffffb2", high = "#d73027", name = "Cov") +
    geom_text(aes(label = label), size = 2.8) +
    labs(title    = sprintf("Final Covariance Matrix (%s)", data_freq_label),
         subtitle = "Diagonal: σ | Off-diagonal: ρ",
         x = "", y = "") +
    theme_risk + theme(panel.grid = element_blank())

  cat("  ✅ Plot 6: Covariance matrix\n")

  # ===========================================================================
  # Plot 7: GARCH time-series — conditional vol, returns, unconditional vol,
  #         and 1-step-ahead forecast annotation.
  #         Only produced when cov_estimation = "garch".
  # ===========================================================================

  if(!is.null(garch_result) && cov_estimation == "garch") {

    cat("\n  Building GARCH conditional volatility time-series plots...\n")

    # Build observation dates (fall back to integer index if no date column)
    plot_dates <- if(!is.null(dates_raw) && length(dates_raw) == n_obs)
                    dates_raw
                  else
                    seq_len(n_obs)

    # Assemble long-form data frame for faceting
    ts_rows <- lapply(seq_len(n_assets), function(i) {
      data.frame(
        Date         = plot_dates,
        Return       = returns_numeric[, i] * 100,
        CondVol      = garch_result$cond_vol_mat[, i] * 100,
        UncondVol    = garch_result$uncond_vol[i] * 100,
        ForecastVol  = garch_result$forecast_vol[i] * 100,
        Asset        = asset_names[i],
        stringsAsFactors = FALSE
      )
    })
    ts_df <- do.call(rbind, ts_rows)

    # Summary for subtitle annotation per asset
    summ_df <- data.frame(
      Asset       = asset_names,
      UncondVol   = garch_result$uncond_vol * 100,
      ForecastVol = garch_result$forecast_vol * 100,
      stringsAsFactors = FALSE
    )
    summ_df$Label <- sprintf(
      "Unconditional \u03c3 = %.3f%% | 1-step forecast = %.3f%%",
      summ_df$UncondVol, summ_df$ForecastVol)

    # ── Colour palette ──────────────────────────────────────────────────────
    col_return   <- "#E8998D"    # salmon  — raw returns
    col_cond     <- "#2C5F8A"    # dark blue — conditional vol
    col_uncond   <- "#E63946"    # red — unconditional (long-run) vol hline
    col_forecast <- "#2ECC71"    # green — 1-step forecast point

    # ── Build per-asset panels ──────────────────────────────────────────────
    ts_panels <- lapply(seq_len(n_assets), function(i) {
      asset_nm   <- asset_names[i]
      df_i       <- ts_df[ts_df$Asset == asset_nm, ]
      uncond_val <- garch_result$uncond_vol[i] * 100
      fcast_val  <- garch_result$forecast_vol[i] * 100
      last_date  <- if(inherits(df_i$Date[1], "Date")) max(df_i$Date) else max(df_i$Date)

      p <- ggplot(df_i, aes(x = Date)) +
        # Raw returns (thin, semi-transparent)
        geom_line(aes(y = Return,  colour = "Returns"),
                  linewidth = 0.35, alpha = 0.55) +
        # Conditional volatility (bold)
        geom_line(aes(y = CondVol, colour = "GARCH(1,1)"),
                  linewidth = 0.85) +
        # Unconditional (long-run) volatility — horizontal reference line
        geom_hline(yintercept = uncond_val,
                   colour = col_uncond, linewidth = 0.6, linetype = "dashed") +
        annotate("text",
                 x     = if(inherits(df_i$Date[1], "Date"))
                           min(df_i$Date) + as.numeric(diff(range(df_i$Date))) * 0.02
                         else min(df_i$Date),
                 y     = uncond_val * 1.08,
                 label = sprintf("Unconditional \u03c3 = %.3f%%", uncond_val),
                 colour = col_uncond, size = 2.4, hjust = 0, fontface = "italic") +
        # 1-step-ahead forecast — point at the right edge
        annotate("point",
                 x     = last_date,
                 y     = fcast_val,
                 colour = col_forecast, size = 2.8, shape = 18) +
        annotate("text",
                 x     = last_date,
                 y     = fcast_val * 1.12,
                 label = sprintf("1-step: %.3f%%", fcast_val),
                 colour = col_forecast, size = 2.4, hjust = 1, fontface = "bold") +
        scale_colour_manual(
          name   = NULL,
          values = c("Returns" = col_return, "GARCH(1,1)" = col_cond),
          guide  = guide_legend(override.aes = list(linewidth = c(0.5, 1.2),
                                                    alpha     = c(0.7, 1.0)))
        ) +
        labs(title    = asset_nm,
             subtitle = sprintf("Uncond. \u03c3 = %.3f%% | 1-step fcast = %.3f%%",
                                uncond_val, fcast_val),
             x = NULL, y = sprintf("Return / Vol (%%)")) +
        theme_risk +
        theme(legend.position  = "bottom",
              legend.key.width  = unit(1.2, "cm"),
              plot.title        = element_text(size = 9, face = "bold"),
              plot.subtitle     = element_text(size = 7.5, colour = "gray40"),
              axis.text.x       = element_text(size = 7),
              axis.text.y       = element_text(size = 7))
      p
    })

    # ── Combine all panels onto a single canvas: 4 columns, n_rows rows ────
    n_cols_ts <- min(4L, n_assets)
    n_rows_ts <- ceiling(n_assets / n_cols_ts)

    portfolio_plots$time_series <-
      wrap_plots(ts_panels, ncol = n_cols_ts) +
      plot_annotation(
        title    = sprintf("GJR-GARCH(1,1) Conditional Volatility — %s Returns", data_freq_label),
        subtitle = paste0("Blue = conditional \u03c3 | Salmon = raw return | ",
                          "Red dashed = unconditional \u03c3 | Green \u25c6 = 1-step forecast"),
        theme    = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "gray40")
        )
      )

    cat(sprintf("  \u2705 Plot 7: GARCH time-series (%d asset(s), %d col \u00d7 %d row canvas)\n",
                n_assets, n_cols_ts, n_rows_ts))

  } else {
    if(cov_estimation != "garch")
      cat("  \u2139\ufe0f Plot 7: GARCH time-series skipped (cov_estimation != 'garch')\n")
  }

  # ===========================================================================
  # 11. PAIRWISE COPULA DIAGNOSTIC PLOTS
  # ===========================================================================

  cat("\n========================================================================\n")
  cat("GENERATING PAIRWISE COPULA DIAGNOSTIC PLOTS\n")
  cat(sprintf("  Reference asset: %s\n", reference_asset))
  cat("========================================================================\n")

  if(copula && copula_success && !is.null(latent_Z)) {

    other_assets   <- asset_names[asset_names != reference_asset]
    n_to_plot      <- min(n_pairwise_plots, length(other_assets))

    if(n_to_plot == 0) {
      cat("  ℹ️ No other assets for pairwise plots\n")
    } else {
      assets_to_plot <- other_assets[seq_len(n_to_plot)]
      cat(sprintf("  Pairs: %s vs %s\n", reference_asset, paste(assets_to_plot, collapse = ", ")))

      # Individual plot panels stored as flat list; wrap_plots will arrange in 4 columns:
      # [Copula A | Historical A | Copula B | Historical B | ...]
      all_panels <- list()

      for(k in seq_along(assets_to_plot)) {
        target_asset <- assets_to_plot[k]
        target_idx   <- which(asset_names == target_asset)

        Z_ref    <- latent_Z[, ref_idx];       Z_target <- latent_Z[, target_idx]
        r_ref    <- returns_numeric[, ref_idx]; r_target <- returns_numeric[, target_idx]

        extreme_threshold <- quantile(r_ref, 0.05)
        extreme_idx       <- r_ref <= extreme_threshold
        n_extreme         <- sum(extreme_idx)
        n_extreme_min     <- max(5L, ceiling(n_extreme * 0.5))

        obs_cor     <- cor(r_ref, r_target)
        extreme_cor <- if(n_extreme >= n_extreme_min) cor(r_ref[extreme_idx], r_target[extreme_idx]) else NA
        cop_cor_val <- copula_cor[ref_idx, target_idx]
        td_val      <- if(!is.null(tail_dependence)) tail_dependence[ref_idx, target_idx] else NA

        latent_df <- data.frame(x = Z_ref, y = Z_target, extreme = extreme_idx)
        x_ann     <- quantile(Z_ref, 0.05)
        y_ann_top <- quantile(Z_target, 0.92)
        y_ann_td  <- quantile(Z_target, 0.84)

        p_copula <- ggplot(latent_df, aes(x = x, y = y, color = extreme)) +
          geom_point(alpha = 0.35, size = 0.7) +
          scale_color_manual(values = c("FALSE" = "#4393c3", "TRUE" = "#d6604d")) +
          geom_smooth(method = "lm", se = TRUE, color = "darkgreen",
                      linewidth = 0.7, alpha = 0.15, inherit.aes = FALSE,
                      aes(x = x, y = y), data = latent_df) +
          annotate("text", x = x_ann, y = y_ann_top,
                   label = sprintf("Copula \u03c1 = %.3f", cop_cor_val),
                   hjust = 0, size = 2.5, fontface = "bold", color = "purple") +
          {if(copula_type == "t" && !is.na(td_val))
            annotate("text", x = x_ann, y = y_ann_td,
                     label = sprintf("Tail \u03bb = %.3f", td_val),
                     hjust = 0, size = 2.5, color = "darkorange")} +
          labs(title    = sprintf("COPULA: %s vs %s", reference_asset, target_asset),
               subtitle = "Normal scores | Red: bottom 5%%",
               x = sprintf("%s (Z)", reference_asset),
               y = sprintf("%s (Z)", target_asset)) +
          theme_risk + theme(legend.position = "none",
                             plot.title    = element_text(size = 7.5, face = "bold"),
                             plot.subtitle = element_text(size = 6.5, color = "gray50"))

        observed_df <- data.frame(x = r_ref*100, y = r_target*100, extreme = extreme_idx)
        x_finite    <- observed_df$x[is.finite(observed_df$x)]
        y_finite    <- observed_df$y[is.finite(observed_df$y)]
        x_ann_obs   <- quantile(x_finite, 0.55)
        y_ann_obs1  <- quantile(y_finite, 0.95)
        y_ann_obs2  <- quantile(y_finite, 0.87)
        y_ann_obs3  <- quantile(y_finite, 0.79)

        p_obs <- ggplot(observed_df, aes(x = x, y = y, color = extreme)) +
          geom_point(alpha = 0.35, size = 0.7) +
          scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "#d6604d")) +
          geom_smooth(data = observed_df, aes(x = x, y = y),
                      method = "lm", se = FALSE, color = "steelblue",
                      linewidth = 0.8, inherit.aes = FALSE) +
          {if(n_extreme >= n_extreme_min)
            geom_smooth(data = observed_df[extreme_idx, ], aes(x = x, y = y),
                        method = "lm", se = FALSE, color = "#d6604d",
                        linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)} +
          annotate("text", x = x_ann_obs, y = y_ann_obs1,
                   label = sprintf("Hist \u03c1 = %.3f", obs_cor),
                   hjust = 0, size = 2.5, color = "steelblue") +
          {if(!is.na(extreme_cor))
            annotate("text", x = x_ann_obs, y = y_ann_obs2,
                     label = sprintf("Extreme \u03c1 = %.3f", extreme_cor),
                     hjust = 0, size = 2.5, color = "#d6604d")} +
          annotate("text", x = x_ann_obs, y = y_ann_obs3,
                   label = sprintf("Copula \u03c1 = %.3f", cop_cor_val),
                   hjust = 0, size = 2.2, color = "purple", fontface = "italic") +
          labs(title    = sprintf("HISTORICAL: %s vs %s", reference_asset, target_asset),
               subtitle = "Observed returns",
               x = sprintf("%s %s Return (%%%%)", reference_asset, data_freq_label),
               y = sprintf("%s %s Return (%%%%)", target_asset,    data_freq_label)) +
          theme_risk + theme(legend.position = "none",
                             plot.title    = element_text(size = 7.5, face = "bold"),
                             plot.subtitle = element_text(size = 6.5, color = "gray50"))

        # Append copula panel then historical panel for this pair
        all_panels[[length(all_panels) + 1]] <- p_copula
        all_panels[[length(all_panels) + 1]] <- p_obs
      }

      # 4-column layout: each row = Copula(A) | Historical(A) | Copula(B) | Historical(B)
      portfolio_plots$pairwise_diagnostics <- wrap_plots(all_panels, ncol = 4) +
        plot_annotation(
          title    = sprintf("Pairwise Dependence: %s vs All Assets", reference_asset),
          subtitle = paste0("Each row: Copula Space | Historical Space  |  ",
                            "Copula Space | Historical Space\n",
                            "Red points & dashed line: bottom 5%% of reference asset"),
          theme    = theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 8, color = "gray40")
          )
        )
      cat(sprintf("  \u2705 Plot 11: Pairwise diagnostics (%d pairs, 4-column layout)\n", n_to_plot))
    }

  } else {
    cat("  ℹ️ Pairwise plots skipped")
    if(!copula)            cat(" — copula = FALSE")
    if(copula && !copula_success) cat(" — copula fitting failed")
    if(copula && is.null(latent_Z)) cat(" — latent_Z unavailable")
    cat("\n")
  }

  # ===========================================================================
  # 11b. TAIL DEPENDENCY + COPULA vs HISTORICAL CORRELATION PLOT
  #      Assigned to: result$plots$portfolio_diagnostics$tail_dependency
  # ===========================================================================

  if(copula && copula_success && !is.null(copula_cor)) {

    cat("\n  Building tail dependency & correlation divergence plots...\n")

    # Build all unique asset pairs with their statistics
    #
    # The scatter plot needs to compare LIKE-FOR-LIKE quantities:
    #   - CopRho (copula ρ) comes from Kendall τ inversion: ρ_cop = sin(π/2 × τ)
    #   - HistRho (Pearson ρ) is a different statistic, always ≥ Kendall-based ρ
    #     for non-normal data, causing ALL points to plot below the 45° line.
    #
    # Fix: add KendallRho = sin(π/2 × τ_hist) as the apples-to-apples baseline.
    # The scatter x-axis uses KendallRho; HistRho is shown as a secondary label.
    # Delta/AbsDelta measure |CopRho − KendallRho| (i.e. how much copula fitting
    # diverges from the simple τ-inversion estimate on the raw returns).
    #
    # Note: since psgt() is strictly monotone, Kendall τ is identical whether
    # computed on raw returns or u_fitted.  We use returns_numeric for clarity.
    kendall_tau_mat <- cor(returns_numeric, method = "kendall")
    # Kendall → copula ρ conversion (same formula as warm_start)
    kendall_rho_mat <- if(copula_type == "t") sin(pi / 2 * kendall_tau_mat) else
                         2 * sin(pi / 6 * kendall_tau_mat)   # Spearman → Pearson for Gaussian

    pair_rows <- list()
    for(ii in seq_len(n_assets - 1)) {
      for(jj in seq(ii + 1, n_assets)) {
        td_val      <- if(!is.null(tail_dependence) && !is.na(tail_dependence[ii, jj]))
                         tail_dependence[ii, jj] else NA_real_
        hist_rho    <- historical_cor[ii, jj]          # Pearson ρ (for reference)
        kendall_rho <- kendall_rho_mat[ii, jj]         # sin(π/2 × τ)  (apples-to-apples)
        cop_rho     <- copula_cor[ii, jj]              # fitted copula ρ
        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          Pair       = paste(asset_names[ii], "\u2014", asset_names[jj]),
          Asset_i    = asset_names[ii],
          Asset_j    = asset_names[jj],
          TailDep    = td_val,
          HistRho    = hist_rho,       # kept for reference only
          KendallRho = kendall_rho,    # x-axis: sin(π/2 × τ_hist) — same concept as CopRho
          CopRho     = cop_rho,
          Delta      = cop_rho - kendall_rho,         # deviation from τ-inversion baseline
          AbsDelta   = abs(cop_rho - kendall_rho)
        )
      }
    }
    pair_df <- do.call(rbind, pair_rows)

    # ---- LEFT PLOT: Tail Dependency Bar Chart ----
    # Only meaningful for t-copula; for Gaussian copula λ = 0 by definition
    if(copula_type == "t" && !all(is.na(pair_df$TailDep))) {

      # Risk threshold: λ > 0.2 is commonly used as a "meaningful tail dependence" cutoff
      td_risk_cutoff <- 0.2

      pair_df_td <- pair_df[!is.na(pair_df$TailDep), ]
      pair_df_td <- pair_df_td[order(-pair_df_td$TailDep), ]
      pair_df_td$Pair <- factor(pair_df_td$Pair, levels = rev(pair_df_td$Pair))
      pair_df_td$Risky <- pair_df_td$TailDep >= td_risk_cutoff

      p_tail <- ggplot(pair_df_td, aes(x = Pair, y = TailDep, fill = Risky)) +
        geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
        geom_hline(yintercept = td_risk_cutoff, linetype = "dashed",
                   color = "#d73027", linewidth = 0.9) +
        annotate("text", x = nrow(pair_df_td) * 0.05, y = td_risk_cutoff + 0.01,
                 label = sprintf("\u03bb = %.2f  (risk threshold)", td_risk_cutoff),
                 hjust = 0, size = 2.8, color = "#d73027", fontface = "italic") +
        scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4393c3"),
                          labels = c("TRUE" = sprintf("\u03bb \u2265 %.2f (risky)", td_risk_cutoff),
                                     "FALSE" = sprintf("\u03bb < %.2f", td_risk_cutoff)),
                          name = "Tail Dependence") +
        scale_y_continuous(limits = c(0, max(pair_df_td$TailDep, td_risk_cutoff) * 1.15),
                           labels = scales::number_format(accuracy = 0.01)) +
        coord_flip() +
        labs(title    = "Tail Dependence Coefficients (\u03bb)",
             subtitle = sprintf("t-Copula (df = %.2f) | Dashed = risk threshold (\u03bb \u2265 %.2f)",
                                copula_df_fitted, td_risk_cutoff),
             x = "Asset Pair",
             y = "Tail Dependence \u03bb") +
        theme_risk +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.4, "cm"),
              legend.text     = element_text(size = 7))

    } else {
      # Gaussian copula: asymptotic tail dependence is always zero — show info panel
      p_tail <- ggplot(data.frame(x = 0.5, y = 0.5,
                                   label = paste0(
                                     "Gaussian Copula\n",
                                     "Asymptotic tail dependence = 0\n",
                                     "by construction.\n\n",
                                     "Switch to copula_type = 't'\n",
                                     "to estimate tail dependence.")),
                       aes(x = x, y = y, label = label)) +
        geom_text(size = 4, color = "gray30", hjust = 0.5, vjust = 0.5) +
        labs(title    = "Tail Dependence",
             subtitle = "Not applicable for Gaussian copula") +
        theme_risk + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                           panel.grid = element_blank())
    }

    # ---- RIGHT PLOT: Copula ρ vs Historical ρ Scatter ----
    # Select top diverging pairs to label (highest |Δρ|)
    n_label <- min(5L, nrow(pair_df))
    top_delta_idx <- order(-pair_df$AbsDelta)[seq_len(n_label)]
    pair_df$Label <- ifelse(seq_len(nrow(pair_df)) %in% top_delta_idx, pair_df$Pair, NA_character_)

    rho_range <- range(c(pair_df$KendallRho, pair_df$CopRho), na.rm = TRUE)
    rho_pad   <- diff(rho_range) * 0.1
    rho_lims  <- c(rho_range[1] - rho_pad, rho_range[2] + rho_pad)

    p_rho <- ggplot(pair_df, aes(x = KendallRho, y = CopRho)) +
      # 45-degree reference line: if copula = τ-inversion exactly, points sit here
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "gray50", linewidth = 0.8) +
      annotate("text", x = rho_lims[1] + rho_pad * 0.3,
               y = rho_lims[2] - rho_pad * 0.3,
               label = "Copula stronger", size = 2.5, color = "gray60", fontface = "italic") +
      annotate("text", x = rho_lims[2] - rho_pad * 0.3,
               y = rho_lims[1] + rho_pad * 0.3,
               label = "\u03c4-inversion stronger", size = 2.5, color = "gray60",
               fontface = "italic", hjust = 1) +
      geom_point(aes(color = AbsDelta), size = 2.5, alpha = 0.85) +
      scale_color_gradient(low = "#4393c3", high = "#d73027",
                            name = "|\u0394\u03c1|",
                            labels = scales::number_format(accuracy = 0.001)) +
      {if(requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(
          data = pair_df[!is.na(pair_df$Label), ],
          aes(label = Label),
          size = 2.3, color = "gray20", max.overlaps = 20,
          segment.size = 0.3, segment.color = "gray60",
          box.padding = 0.4, point.padding = 0.3
        )
      } else {
        geom_text(data = pair_df[!is.na(pair_df$Label), ],
                  aes(label = Label),
                  size = 2.3, hjust = -0.1, color = "gray20", check_overlap = TRUE)
      }} +
      scale_x_continuous(limits = rho_lims) +
      scale_y_continuous(limits = rho_lims) +
      labs(title    = "Copula \u03c1 vs \u03c4-Inversion \u03c1 (Kendall baseline)",
           subtitle = paste0(
             "x-axis = sin(\u03c0/2 \u00d7 \u03c4_hist): same concept as copula \u03c1, apples-to-apples\n",
             sprintf("Dashed = 45\u00b0 line (copula = \u03c4-inversion) | Labels: top %d pairs by |\u0394\u03c1|",
                     n_label)),
           x = "Kendall \u03c4-based \u03c1 = sin(\u03c0/2 \u00d7 \u03c4)  [apples-to-apples baseline]",
           y = "Copula Correlation (\u03c1)") +
      theme_risk +
      theme(legend.position = "right",
            legend.key.height = unit(0.6, "cm"))

    # Combine into 1-row, 2-column layout
    portfolio_plots$tail_dependency <- p_tail + p_rho +
      plot_annotation(
        title    = "Copula Diagnostics: Tail Dependence & Correlation Divergence",
        subtitle = sprintf("%s Copula | %d asset pairs",
                           toupper(copula_type), nrow(pair_df)),
        theme    = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40")
        )
      )

    cat("  \u2705 Plot 11b: Tail dependency & correlation divergence\n")

  } else {
    cat("  \u2139\ufe0f  Tail dependency plot skipped — copula = FALSE or fitting failed\n")
  }

  # ===========================================================================
  # 12. STRESS PROPAGATION PLOTS
  # ===========================================================================
  #
  # Four scenarios based on what was shocked:
  #   Scenario A — Endogenous only (portfolio assets shocked, no exogenous)
  #   Scenario B — Exogenous only  (non-portfolio instruments shocked)
  #   Scenario C — Mixed           (both endogenous and exogenous shocked)
  #
  # In all scenarios, a common set of helper data frames is built first,
  # then the appropriate plots are assembled.
  # ===========================================================================

  stressed_plots <- list()

  if(is_stressed) {

    # ── Identify what was shocked ──────────────────────────────────────────────
    has_endo_plot  <- !is.null(asset_shock) && any(asset_shock != 0)
    has_exog_plot  <- has_exogenous && (
      (!is.null(exogenous_shock)            && length(exogenous_shock)            > 0) ||
      (!is.null(exogenous_volatility_shock) && length(exogenous_volatility_shock) > 0)
    )

    shocked_indices      <- if(has_endo_plot) which(asset_shock != 0) else integer(0)
    shocked_count        <- length(shocked_indices)
    shocked_assets_names <- asset_names[shocked_indices]

    # All exog names that have ANY shock (return or vol) — used for display and beta loops
    exog_ret_shocked_names <- if(!is.null(exogenous_shock))            names(exogenous_shock)            else character(0)
    exog_vol_shocked_names <- if(!is.null(exogenous_volatility_shock)) names(exogenous_volatility_shock) else character(0)
    exog_shocked_names     <- unique(c(exog_ret_shocked_names, exog_vol_shocked_names))
    if(!has_exog_plot) exog_shocked_names <- character(0)
    n_exog_s               <- length(exog_shocked_names)

    display_cor <- tryCatch(cov2cor(base_sigma), error = function(e) historical_cor)

    cat("\n========================================================================\n")
    cat("CREATING STRESS PROPAGATION PLOTS\n")
    cat("========================================================================\n")
    scenario_label <- if(has_endo_plot && has_exog_plot) "C: Mixed (endogenous + exogenous)"
                      else if(has_endo_plot)             "A: Endogenous (portfolio assets)"
                      else if(has_exog_plot)             "B: Exogenous (non-portfolio instruments)"
                      else                               "None"
    cat(sprintf("  Scenario             : %s\n", scenario_label))
    if(has_endo_plot)
      cat(sprintf("  Endo shocked (%d)   : %s\n", shocked_count,
                  paste(shocked_assets_names, collapse = ", ")))
    if(has_exog_plot)
      cat(sprintf("  Exog shocked (%d)   : %s\n", n_exog_s,
                  paste(exog_shocked_names, collapse = ", ")))

    # ── COMMON HELPER: portfolio × exogenous beta computation ─────────────────
    # beta_pe[i, e] = Cov(X_port_i, X_exog_e) / Var(X_exog_e)
    # Interpretation: for each 1% move in exog asset e, portfolio asset i moves
    # beta_pe[i,e] % on average (conditional mean propagation coefficient).
    # This is the structural propagation channel regardless of scenario.
    build_exog_beta_df <- function() {
      if(!has_exog_plot || is.null(Sigma_pe_hist) || is.null(Sigma_ee_hist))
        return(NULL)

      # Correlation must be computed consistently: all three terms (covariance,
      # portfolio variance, exogenous variance) must come from the SAME source.
      #
      # Sigma_pe_hist and Sigma_ee_hist are from the historical sample covariance
      # of [returns_numeric, exog_numeric]. The portfolio variance used in the
      # denominator must therefore also be historical: historical_cov[pi, pi].
      #
      # Using base_sigma[pi, pi] (GARCH forward-looking variance) would mix
      # historical covariance with forward-looking variance, producing correlations
      # outside [-1, 1] whenever GARCH vol diverges from historical vol.
      hist_var_port <- diag(historical_cov)   # historical marginal variances

      rows <- list()
      for(en in exog_shocked_names) {
        ei      <- which(exog_names == en)
        var_e   <- Sigma_ee_hist[ei, ei]
        shock_e <- if(!is.null(exogenous_shock) && en %in% names(exogenous_shock))
                     exogenous_shock[en] else 0   # 0 when instrument is vol-shocked only
        vol_m   <- if(!is.null(exogenous_volatility_shock) && en %in% names(exogenous_volatility_shock))
                     exogenous_volatility_shock[en] else 1

        for(pi in seq_len(n_assets)) {
          cov_pe   <- Sigma_pe_hist[pi, ei]
          beta     <- if(var_e > 1e-12) cov_pe / var_e else 0
          prop_ch  <- beta * shock_e * 100
          var_port <- hist_var_port[pi]

          # ρ = Cov(X_port_i, X_exog_e) / sqrt(Var(X_port_i) * Var(X_exog_e))
          # All terms from historical sample covariance → guaranteed ∈ [-1, 1]
          rho <- if(var_e > 1e-12 && var_port > 1e-12)
                   cov_pe / sqrt(var_port * var_e)
                 else NA_real_
          rho <- if(!is.na(rho)) max(-1, min(1, rho)) else NA_real_  # numerical safety

          rows[[length(rows) + 1]] <- data.frame(
            Portfolio_Asset       = asset_names[pi],
            Exog_Asset            = en,
            Beta                  = beta,
            Shock_Size_pct        = shock_e * 100,
            Propagated_Change_pct = prop_ch,
            Vol_Multiplier        = vol_m,
            Correlation           = rho,
            stringsAsFactors      = FALSE
          )
        }
      }
      do.call(rbind, rows)
    }

    exog_beta_df <- build_exog_beta_df()

    # ── COMMON HELPER: mean-change decomposition ───────────────────────────────
    # For each portfolio asset, split Δμ into:
    #   direct_endo     : asset_shock[i] if asset i is directly shocked
    #   prop_from_endo  : propagated effect from endogenous shocked assets
    #   prop_from_exog  : propagated effect from exogenous shocks
    #   total           : mu_stressed[i] - base_mu[i]
    build_decomp_df <- function() {
      total_delta <- (mu_stressed - base_mu) * 100

      # Exogenous contribution per portfolio asset
      exog_contrib <- rep(0, n_assets)
      names(exog_contrib) <- asset_names
      if(has_exog_plot && !is.null(exog_beta_df)) {
        for(en in exog_shocked_names) {
          sub   <- exog_beta_df[exog_beta_df$Exog_Asset == en, ]
          shock_e <- if(!is.null(exogenous_shock) && en %in% names(exogenous_shock))
                       exogenous_shock[en] * 100 else 0
          for(pi in seq_len(n_assets)) {
            row_i <- sub[sub$Portfolio_Asset == asset_names[pi], ]
            if(nrow(row_i) == 1)
              exog_contrib[pi] <- exog_contrib[pi] + row_i$Beta * shock_e
          }
        }
      }

      # Direct endogenous shock for directly shocked assets
      direct_endo <- rep(0, n_assets)
      names(direct_endo) <- asset_names
      if(has_endo_plot)
        direct_endo[shocked_indices] <- asset_shock[shocked_indices] * 100

      # Endogenous propagation (residual for non-shocked assets after exog subtracted)
      prop_endo <- total_delta - direct_endo - exog_contrib

      data.frame(
        Asset          = asset_names,
        Total_Delta    = total_delta,
        Direct_Endo    = direct_endo,
        Prop_From_Endo = prop_endo,
        Prop_From_Exog = exog_contrib,
        Is_Endo_Shocked = asset_names %in% shocked_assets_names,
        stringsAsFactors = FALSE
      )
    }

    decomp_df <- build_decomp_df()

    # ── SHARED PLOT HELPERS ────────────────────────────────────────────────────
    # P_means: Base vs Stressed means dodge bar (used in A, B, C)
    make_means_plot <- function(subtitle_txt = "Conditional mean after shock propagation") {
      means_long <- data.frame(
        Asset    = rep(asset_names, 2),
        Return   = c(base_mu*100, mu_stressed*100),
        Scenario = rep(c("Base", "Stressed"), each = n_assets)
      )
      ggplot(means_long, aes(x = Asset, y = Return, fill = Scenario)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
        scale_fill_manual(values = c("Base" = "steelblue", "Stressed" = "coral")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(aes(label = sprintf("%+.3f%%", Return)),
                  position = position_dodge(width = 0.9),
                  vjust = ifelse(means_long$Return >= 0, -0.4, 1.2), size = 2.5) +
        labs(title = "Expected Returns: Base vs Stressed",
             subtitle = subtitle_txt, x = "", y = sprintf("%s Return (%%)", data_freq_label)) +
        theme_risk
    }

    # P_exog_channel: exogenous propagation betas (used in B, C)
    make_exog_channel_plot <- function() {
      if(is.null(exog_beta_df)) return(NULL)
      # For display, aggregate across exog shocks if multiple
      agg <- aggregate(Propagated_Change_pct ~ Portfolio_Asset, data = exog_beta_df, FUN = sum)
      agg$Positive <- agg$Propagated_Change_pct >= 0

      # If single exog asset, also show individual betas
      if(n_exog_s == 1) {
        en      <- exog_shocked_names[1]
        sub     <- exog_beta_df[exog_beta_df$Exog_Asset == en, ]
        sub$Positive <- sub$Beta >= 0
        p <- ggplot(sub, aes(x = reorder(Portfolio_Asset, -Beta), y = Beta, fill = Positive)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_manual(values = c("TRUE" = "#4393c3", "FALSE" = "#d6604d"), guide = "none") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          geom_text(aes(label = sprintf("%.4f", Beta)),
                    vjust = ifelse(sub$Beta >= 0, -0.4, 1.2), size = 2.8) +
          labs(title    = sprintf("Exogenous Propagation Betas (\u03b2 = Cov / Var(%s))", en),
               subtitle = sprintf("For each 1%% move in %s, portfolio asset moves \u03b2%%", en),
               x = "Portfolio Asset", y = sprintf("\u03b2 w.r.t. %s", en)) +
          theme_risk
      } else {
        # Multiple exog shocks: grouped bar by exog asset
        exog_beta_df$Positive <- exog_beta_df$Beta >= 0
        p <- ggplot(exog_beta_df,
                    aes(x = Portfolio_Asset, y = Beta, fill = Exog_Asset)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          labs(title    = "Exogenous Propagation Betas by Instrument",
               subtitle = "\u03b2 = Cov(portfolio asset, exog asset) / Var(exog asset)",
               x = "Portfolio Asset", y = "\u03b2", fill = "Exog Instrument") +
          theme_risk
      }
      p
    }

    # P_exog_corr: portfolio × exogenous correlation (used in B, C)
    make_exog_corr_plot <- function() {
      if(is.null(exog_beta_df)) return(NULL)
      exog_beta_df$Positive <- exog_beta_df$Correlation >= 0

      if(n_exog_s == 1) {
        en  <- exog_shocked_names[1]
        sub <- exog_beta_df[exog_beta_df$Exog_Asset == en, ]
        sub$Positive <- sub$Correlation >= 0
        ggplot(sub, aes(x = reorder(Portfolio_Asset, -Correlation),
                        y = Correlation, fill = Correlation)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_gradient2(low = "#2166ac", high = "#d73027", mid = "white",
                                midpoint = 0, limits = c(-1, 1), name = "\u03c1") +
          geom_hline(yintercept = c(0, 0.5, -0.5),
                     linetype = c("dashed","dotted","dotted"), color = "gray60") +
          geom_text(aes(label = sprintf("%.3f", Correlation)),
                    vjust = ifelse(sub$Correlation >= 0, -0.4, 1.2), size = 2.8) +
          labs(title    = sprintf("Portfolio \u2014 %s Correlation", en),
               subtitle = "Drives the magnitude of exogenous shock propagation",
               x = "Portfolio Asset", y = sprintf("\u03c1 with %s", en)) +
          theme_risk + theme(legend.position = "none")
      } else {
        ggplot(exog_beta_df,
               aes(x = Portfolio_Asset, y = Correlation, fill = Exog_Asset)) +
          geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          labs(title    = "Portfolio \u00d7 Exogenous Correlations",
               subtitle = "Correlation drives the shock propagation magnitude",
               x = "Portfolio Asset", y = "Correlation \u03c1", fill = "Exog Instrument") +
          theme_risk
      }
    }

    # P_exog_summary: exogenous shock magnitude bar (used in B, C)
    make_exog_summary_plot <- function() {
      if(!has_exog_plot) return(NULL)
      exog_sum_df <- data.frame(
        Asset        = exog_shocked_names,
        Shock_pct    = sapply(exog_shocked_names, function(en)
                         if(!is.null(exogenous_shock) && en %in% names(exogenous_shock))
                           exogenous_shock[en] * 100 else 0),
        Vol_Mult     = sapply(exog_shocked_names, function(en)
          if(!is.null(exogenous_volatility_shock) && en %in% names(exogenous_volatility_shock))
            exogenous_volatility_shock[en] else 1),
        stringsAsFactors = FALSE
      )
      exog_sum_df$Positive <- exog_sum_df$Shock_pct >= 0
      exog_sum_df$Label    <- sprintf("%+.2f%%\n(vol \u00d7%.1f)", exog_sum_df$Shock_pct, exog_sum_df$Vol_Mult)

      ggplot(exog_sum_df, aes(x = reorder(Asset, -abs(Shock_pct)),
                               y = Shock_pct, fill = Positive)) +
        geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
        scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "#d6604d"), guide = "none") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(aes(label = Label),
                  vjust = ifelse(exog_sum_df$Shock_pct >= 0, -0.2, 1.2), size = 3) +
        labs(title    = "Exogenous Instrument Shocks Applied",
             subtitle = "Instruments NOT held in portfolio — used for propagation only",
             x = "Exogenous Instrument", y = "Return Shock (%)") +
        theme_risk
    }

    # ============================================================
    # SCENARIO A — Endogenous shocks only
    # ============================================================
    if(has_endo_plot && !has_exog_plot) {

      if(shocked_count == 1) {

        # A1: Single endogenous shock
        shocked_asset_name <- shocked_assets_names[1]
        shocked_idx        <- shocked_indices[1]
        shock_size_pct     <- asset_shock[shocked_idx] * 100
        var_shocked        <- base_sigma[shocked_idx, shocked_idx]
        covariances        <- base_sigma[shocked_idx, ]
        betas              <- as.vector(covariances / var_shocked)
        names(betas)       <- asset_names
        betas[shocked_idx] <- 1

        prop_data <- data.frame(
          Asset           = factor(asset_names, levels = asset_names[order(-betas)]),
          Beta            = betas,
          Expected_Change = betas * asset_shock[shocked_idx] * 100,
          Base_Mean       = base_mu * 100,
          Stressed_Mean   = mu_stressed * 100,
          Correlation     = display_cor[shocked_idx, ],
          Is_Shocked      = asset_names == shocked_asset_name,
          stringsAsFactors = FALSE
        )

        p1 <- ggplot(prop_data, aes(x = reorder(Asset, -Beta), y = Beta, fill = Is_Shocked)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_manual(values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
                            labels = c("TRUE" = "Shocked", "FALSE" = "Propagated"),
                            name = "") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          geom_text(aes(label = sprintf("%.3f", Beta)),
                    vjust = ifelse(prop_data$Beta >= 0, -0.4, 1.2), size = 2.8) +
          labs(title    = sprintf("Propagation Betas from %s", shocked_asset_name),
               subtitle = sprintf("\u03b2 = Cov(i, %s)/Var(%s) — 1%% shock \u2192 \u03b2%% change per asset",
                                  shocked_asset_name, shocked_asset_name),
               x = "", y = "\u03b2") +
          theme_risk + theme(legend.position = "none")

        p2 <- ggplot(prop_data, aes(x = reorder(Asset, -Expected_Change),
                                     y = Expected_Change, fill = Expected_Change >= 0)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "#d6604d"), guide = "none") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          geom_text(aes(label = sprintf("%+.3f%%", Expected_Change)),
                    vjust = ifelse(prop_data$Expected_Change >= 0, -0.4, 1.2), size = 2.8) +
          labs(title    = sprintf("Propagated Mean Change from %.2f%% Shock", shock_size_pct),
               subtitle = "Conditional mean update via Schur complement",
               x = "", y = "Expected \u0394\u03bc (%)") +
          theme_risk + theme(legend.position = "none")

        p3 <- make_means_plot(
          sprintf("Endogenous shock: %.2f%% to %s", shock_size_pct, shocked_asset_name))

        p4 <- ggplot(prop_data, aes(x = reorder(Asset, -Correlation),
                                     y = Correlation, fill = Correlation)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_gradient2(low = "#2166ac", high = "#d73027", mid = "white",
                                midpoint = 0, limits = c(-1, 1), name = "\u03c1") +
          geom_hline(yintercept = c(-0.5, 0, 0.5),
                     linetype = c("dotted","dashed","dotted"), color = "gray60") +
          geom_text(aes(label = sprintf("%.3f", Correlation)),
                    vjust = ifelse(prop_data$Correlation >= 0, -0.4, 1.2), size = 2.8) +
          labs(title    = sprintf("Portfolio Correlation with %s", shocked_asset_name),
               subtitle = "Higher \u03c1 \u2192 larger conditional mean update",
               x = "", y = "Correlation \u03c1") +
          theme_risk + theme(legend.position = "none")

        stressed_plots$propagation_effects <- (p1 + p2) / (p3 + p4) +
          plot_annotation(
            title    = sprintf("Scenario A: Single Endogenous Shock — %.2f%% to %s",
                               shock_size_pct, shocked_asset_name),
            subtitle = if(propagation) "Propagation via conditional covariance (Schur complement)" else "Propagation: DISABLED",
            theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                             plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
          )
        cat("  \u2705 Plot 12a: Scenario A — single endogenous shock\n")

      } else {

        # A2: Multiple endogenous shocks
        shock_df <- data.frame(
          Asset        = asset_names,
          Direct_Shock = asset_shock * 100,
          Total_Change = (mu_stressed - base_mu) * 100,
          Is_Shocked   = asset_shock != 0,
          stringsAsFactors = FALSE
        )

        p1 <- ggplot(shock_df, aes(x = reorder(Asset, -Total_Change))) +
          geom_bar(aes(y = Total_Change, fill = Is_Shocked), stat = "identity", alpha = 0.8) +
          geom_point(aes(y = Direct_Shock), color = "black", size = 2.5, shape = 18) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "steelblue"),
                            labels = c("TRUE" = "Directly shocked", "FALSE" = "Propagated")) +
          labs(title    = "Total \u0394\u03bc per Portfolio Asset",
               subtitle = "Bars = total | Diamonds = direct shock only",
               x = "", y = "\u0394\u03bc (%)", fill = "") +
          theme_risk

        impact_matrix <- matrix(0, n_assets, n_assets,
                                dimnames = list(asset_names, asset_names))
        for(ii in seq_len(n_assets))
          for(jj in shocked_indices)
            impact_matrix[ii, jj] <- (base_sigma[ii,jj]/base_sigma[jj,jj])*asset_shock[jj]*100

        impact_melt <- reshape2::melt(impact_matrix)
        colnames(impact_melt) <- c("Asset", "Shocked_Asset", "Impact")
        impact_melt <- impact_melt[impact_melt$Shocked_Asset %in% shocked_assets_names, ]

        p2 <- ggplot(impact_melt, aes(x = Shocked_Asset, y = Asset, fill = Impact)) +
          geom_tile(color = "white") +
          scale_fill_gradient2(low = "#d6604d", high = "#4dac26", mid = "white",
                                midpoint = 0, name = "\u0394\u03bc (%)") +
          geom_text(aes(label = sprintf("%+.3f%%", Impact)), size = 2.8) +
          labs(title    = "Endogenous Propagation Impact Matrix",
               subtitle = "Cell[i,j] = \u0394\u03bc in asset i from shock to portfolio asset j",
               x = "Directly Shocked Asset", y = "Receiving Asset") +
          theme_risk

        shock_summary <- data.frame(
          Asset            = shocked_assets_names,
          Direct_Shock_pct = asset_shock[shocked_indices] * 100,
          Vol_Multiplier   = if(!is.null(volatility_shock)) volatility_shock[shocked_indices]
                             else rep(1, shocked_count),
          stringsAsFactors = FALSE
        )
        shock_summary$Label <- sprintf("%+.2f%%\n(vol \u00d7%.1f)",
                                        shock_summary$Direct_Shock_pct,
                                        shock_summary$Vol_Multiplier)

        p3 <- ggplot(shock_summary, aes(x = reorder(Asset, -abs(Direct_Shock_pct)),
                                         y = Direct_Shock_pct,
                                         fill = Direct_Shock_pct >= 0)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "#d6604d"), guide = "none") +
          geom_text(aes(label = Label),
                    vjust = ifelse(shock_summary$Direct_Shock_pct >= 0, -0.2, 1.2), size = 3) +
          labs(title = "Direct Endogenous Shocks Applied",
               subtitle = "Portfolio assets shocked directly",
               x = "", y = "Shock Size (%)") +
          theme_risk + theme(legend.position = "none")

        impact_total <- rowSums(impact_matrix[, shocked_indices, drop = FALSE])
        impact_df    <- data.frame(Asset = asset_names, Total_Impact = impact_total,
                                    Is_Shocked = asset_shock != 0, stringsAsFactors = FALSE)

        p4 <- ggplot(impact_df, aes(x = reorder(Asset, Total_Impact),
                                     y = Total_Impact, fill = Is_Shocked)) +
          geom_bar(stat = "identity", alpha = 0.85) +
          coord_flip() +
          scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "steelblue"),
                            labels = c("TRUE" = "Directly shocked", "FALSE" = "Propagated")) +
          geom_text(aes(label = sprintf("%+.3f%%", Total_Impact)),
                    hjust = ifelse(impact_df$Total_Impact >= 0, -0.1, 1.1), size = 2.8) +
          labs(title = "Total Propagated Impact per Asset",
               x = "", y = "Total \u0394\u03bc (%)", fill = "") +
          theme_risk

        stressed_plots$multi_shock_effects <- (p1 + p2) / (p3 + p4) +
          plot_annotation(
            title    = sprintf("Scenario A: Multiple Endogenous Shocks — %s",
                               paste(shocked_assets_names, collapse = ", ")),
            subtitle = if(propagation) "Propagation via conditional covariance (Schur complement)" else "Propagation: DISABLED",
            theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                             plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
          )
        cat("  \u2705 Plot 12a: Scenario A — multiple endogenous shocks\n")
      }
    }

    # ============================================================
    # SCENARIO B — Exogenous shocks only
    # ============================================================
    if(!has_endo_plot && has_exog_plot) {

      has_ret_shock_b <- length(exog_ret_shocked_names) > 0
      has_vol_shock_b <- length(exog_vol_shocked_names) > 0

      p_exog_summary <- make_exog_summary_plot()
      p_exog_corr    <- make_exog_corr_plot()
      p_exog_beta    <- make_exog_channel_plot()
      p_means        <- if(has_ret_shock_b)
                          make_means_plot("Portfolio means after exogenous shock propagation")
                        else NULL

      b_subtitle <- if(has_ret_shock_b && has_vol_shock_b)
        "Exogenous return shock propagated via Σ_pe Σ_ee⁻¹ δ_e; vol shock widened Σ_pp via LTV addon"
      else if(has_ret_shock_b)
        "Exogenous instruments NOT held in portfolio. Return shock propagated to portfolio means via Schur."
      else
        "Exogenous vol shock widens portfolio covariance via LTV. No mean effect; mean plot omitted."

      b_panels <- Filter(Negate(is.null), list(p_exog_summary, p_exog_corr, p_exog_beta, p_means))
      if(length(b_panels) >= 2) {
        n_cols_b <- min(2L, length(b_panels))
        stressed_plots$exogenous_propagation <-
          wrap_plots(b_panels, ncol = n_cols_b) +
          plot_annotation(
            title    = sprintf("Scenario B: Exogenous Shock%s — %s",
                               if(n_exog_s > 1) "s" else "",
                               paste(exog_shocked_names, collapse = ", ")),
            subtitle = b_subtitle,
            theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                             plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
          )
        cat("  ✅ Plot 12b: Scenario B — exogenous-only propagation
")
      } else {
        cat("  ℹ️ Plot 12b: insufficient panels for Scenario B
")
      }
    }

    # ============================================================
    # SCENARIO C — Mixed: both endogenous and exogenous shocks
    # ============================================================
    if(has_endo_plot && has_exog_plot) {

      # C1: Decomposition waterfall — stacked bar showing three components of Δμ
      decomp_long <- do.call(rbind, list(
        data.frame(Asset = decomp_df$Asset, Value = decomp_df$Direct_Endo,
                   Component = "Direct (endo)", stringsAsFactors = FALSE),
        data.frame(Asset = decomp_df$Asset, Value = decomp_df$Prop_From_Endo,
                   Component = "Propagated from endo", stringsAsFactors = FALSE),
        data.frame(Asset = decomp_df$Asset, Value = decomp_df$Prop_From_Exog,
                   Component = "Propagated from exog", stringsAsFactors = FALSE)
      ))
      decomp_long$Component <- factor(decomp_long$Component,
                                       levels = c("Direct (endo)",
                                                  "Propagated from endo",
                                                  "Propagated from exog"))

      p_decomp <- ggplot(decomp_long, aes(x = Asset, y = Value, fill = Component)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.87) +
        geom_point(data = decomp_df,
                   aes(x = Asset, y = Total_Delta, shape = "Total \u0394\u03bc"),
                   inherit.aes = FALSE, size = 3, color = "black") +
        scale_shape_manual(values = c("Total \u0394\u03bc" = 18), name = "") +
        scale_fill_manual(values = c("Direct (endo)"        = "#d6604d",
                                     "Propagated from endo" = "#f4a582",
                                     "Propagated from exog" = "#4393c3"),
                          name = "\u0394\u03bc Component") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(data = decomp_df,
                  aes(x = Asset, y = Total_Delta,
                      label = sprintf("%+.3f%%", Total_Delta),
                      vjust = ifelse(Total_Delta >= 0, -0.6, 1.4)),
                  inherit.aes = FALSE, size = 2.8, fontface = "bold") +
        labs(title    = "Mean Change Decomposition",
             subtitle = "Stacked: direct shock + endogenous propagation + exogenous propagation",
             x = "Portfolio Asset", y = "Total \u0394\u03bc (%)") +
        theme_risk + theme(legend.position = "bottom")

      # C2: Exogenous propagation channel
      p_exog_beta <- make_exog_channel_plot()

      # C3: Endogenous shock sizes (for reference)
      endo_sum_df <- data.frame(
        Asset        = if(shocked_count > 0) shocked_assets_names else asset_names,
        Shock_pct    = if(shocked_count > 0) asset_shock[shocked_indices]*100 else rep(0, n_assets),
        Vol_Mult     = if(!is.null(volatility_shock) && shocked_count > 0)
                         volatility_shock[shocked_indices] else rep(1, max(shocked_count, 1)),
        stringsAsFactors = FALSE
      )
      endo_sum_df$Label <- sprintf("%+.2f%%\n(vol \u00d7%.1f)",
                                    endo_sum_df$Shock_pct, endo_sum_df$Vol_Mult)
      endo_sum_df$Positive <- endo_sum_df$Shock_pct >= 0

      p_endo_sum <- ggplot(endo_sum_df, aes(x = reorder(Asset, -abs(Shock_pct)),
                                             y = Shock_pct, fill = Positive)) +
        geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
        scale_fill_manual(values = c("TRUE" = "forestgreen", "FALSE" = "#d6604d"), guide = "none") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(aes(label = Label),
                  vjust = ifelse(endo_sum_df$Shock_pct >= 0, -0.2, 1.2), size = 3) +
        labs(title    = "Endogenous (Portfolio) Shocks Applied",
             subtitle = "Assets directly shocked within the portfolio",
             x = "", y = "Return Shock (%)") +
        theme_risk + theme(legend.position = "none")

      # C4: Exogenous shock summary
      p_exog_summary <- make_exog_summary_plot()

      # C5: Resulting base vs stressed means
      p_means <- make_means_plot(
        sprintf("Mixed shock: {%s} endo + {%s} exog",
                paste(shocked_assets_names, collapse=", "),
                paste(exog_shocked_names, collapse=", ")))

      # Layout: row 1 = decomposition (full width) | exog beta
      #         row 2 = endo shocks | exog shocks | means
      stressed_plots$mixed_shock_effects <-
        (p_decomp + p_exog_beta) / (p_endo_sum + p_exog_summary + p_means) +
        plot_annotation(
          title    = sprintf("Scenario C: Mixed Shocks — Endo {%s} + Exog {%s}",
                             paste(shocked_assets_names, collapse=", "),
                             paste(exog_shocked_names,   collapse=", ")),
          subtitle = paste0(
            "Stacked bar shows how each component contributes to the total mean change.\n",
            "Exogenous instruments propagate via \u03a3_pe \u03a3_ee\u207b\u00b9 \u03b4_e; endogenous via \u03a3_ns \u03a3_ss\u207b\u00b9 \u03b4_s"),
          theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                           plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
        )
      cat("  \u2705 Plot 12c: Scenario C — mixed endogenous + exogenous shocks\n")
    }

    if(!has_endo_plot && !has_exog_plot)
      cat("  \u2139\ufe0f No return shocks — propagation plots skipped\n")
  }


  # Pre/post distribution comparison
  if(is_stressed) {

    cat("\n  Building pre/post stress distribution comparison...\n")

    pre_mean  <- stats_pre$mean   * 100;  pre_v95 <- stats_pre$var_95 * 100;  pre_v99 <- stats_pre$var_99 * 100
    post_mean <- stats_post$mean  * 100; post_v95 <- stats_post$var_95 * 100; post_v99 <- stats_post$var_99 * 100

    comparison_df <- data.frame(
      Return   = c(portfolio_returns_pre_flat*100, portfolio_returns_post_flat*100),
      Scenario = rep(c("Pre-Stress","Post-Stress"),
                     times = c(length(portfolio_returns_pre_flat),
                               length(portfolio_returns_post_flat)))
    )

    x_min <- min(pre_v99, post_v99) - 2
    x_max <- max(pre_mean, post_mean) + max(stats_pre$vol, stats_post$vol) * 200
    scenario_colours <- c("Pre-Stress" = "steelblue", "Post-Stress" = "coral")

    p_density <- ggplot(comparison_df, aes(x = Return, fill = Scenario, color = Scenario)) +
      geom_density(alpha = 0.4, linewidth = 0.9) +
      scale_fill_manual(values = scenario_colours) +
      scale_color_manual(values = scenario_colours) +
      geom_vline(xintercept = pre_v95,  color = "steelblue", linetype = "dashed",  linewidth = 0.8, alpha = 0.8) +
      geom_vline(xintercept = pre_v99,  color = "steelblue", linetype = "dotted",  linewidth = 0.8, alpha = 0.8) +
      geom_vline(xintercept = post_v95, color = "coral",     linetype = "dashed",  linewidth = 0.8, alpha = 0.8) +
      geom_vline(xintercept = post_v99, color = "coral",     linetype = "dotted",  linewidth = 0.8, alpha = 0.8) +
      annotate("text", x = pre_v95,  y = Inf, label = sprintf("%.2f%%", pre_v95),
               angle = 90, vjust = 2, hjust = -0.1, size = 2.5, color = "steelblue") +
      annotate("text", x = pre_v99,  y = Inf, label = sprintf("%.2f%%", pre_v99),
               angle = 90, vjust = 2, hjust = -0.1, size = 2.5, color = "steelblue") +
      annotate("text", x = post_v95, y = Inf, label = sprintf("%.2f%%", post_v95),
               angle = 90, vjust = 2, hjust = -0.1, size = 2.5, color = "coral") +
      annotate("text", x = post_v99, y = Inf, label = sprintf("%.2f%%", post_v99),
               angle = 90, vjust = 2, hjust = -0.1, size = 2.5, color = "coral") +
      coord_cartesian(xlim = c(x_min, x_max)) +
      labs(title = "Left Tail Focus", subtitle = "Dashed=VaR95 | Dotted=VaR99",
           x = sprintf("%s Return (%%)", data_freq_label), y = "Density") +
      theme_risk

    p_hist <- ggplot(comparison_df, aes(x = Return, fill = Scenario)) +
      geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.45,
                     position = "identity", color = "white", linewidth = 0.2) +
      geom_density(aes(color = Scenario), linewidth = 0.9, alpha = 0) +
      scale_fill_manual(values = scenario_colours) +
      scale_color_manual(values = scenario_colours) +
      geom_vline(xintercept = pre_mean,  color = "steelblue", linetype = "solid",  linewidth = 1,   alpha = 0.8) +
      geom_vline(xintercept = post_mean, color = "coral",     linetype = "solid",  linewidth = 1,   alpha = 0.8) +
      geom_vline(xintercept = pre_v95,   color = "steelblue", linetype = "dashed", linewidth = 0.7, alpha = 0.7) +
      geom_vline(xintercept = post_v95,  color = "coral",     linetype = "dashed", linewidth = 0.7, alpha = 0.7) +
      labs(title = "Full Distribution", subtitle = "Solid=mean | Dashed=VaR95",
           x = sprintf("%s Return (%%)", data_freq_label), y = "Density") +
      theme_risk

    n_qq <- min(length(portfolio_returns_pre_flat), length(portfolio_returns_post_flat))
    theo <- qnorm(ppoints(n_qq))
    qq_df <- rbind(
      data.frame(Theoretical = theo, Sample = sort(portfolio_returns_pre_flat[seq_len(n_qq)]*100),  Scenario = "Pre-Stress"),
      data.frame(Theoretical = theo, Sample = sort(portfolio_returns_post_flat[seq_len(n_qq)]*100), Scenario = "Post-Stress")
    )

    p_qq <- ggplot(qq_df, aes(x = Theoretical, y = Sample, color = Scenario)) +
      geom_point(alpha = 0.25, size = 0.7) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
      scale_color_manual(values = scenario_colours) +
      labs(title = "QQ Plot Comparison", subtitle = "Departure = non-normality / fat tails",
           x = "Standard Normal Quantiles", y = "Sample Quantiles (%)") +
      theme_risk

    stats_bar_df <- data.frame(
      Metric      = c("Mean Return (%)", "Volatility (%)", "95% VaR (%)", "99% VaR (%)", "95% CVaR (%)"),
      Pre_Stress  = c(stats_pre$mean*100,  stats_pre$vol*100,  stats_pre$var_95*100,  stats_pre$var_99*100,  stats_pre$cvar_95*100),
      Post_Stress = c(stats_post$mean*100, stats_post$vol*100, stats_post$var_95*100, stats_post$var_99*100, stats_post$cvar_95*100)
    )
    stats_long <- reshape2::melt(stats_bar_df, id.vars = "Metric",
                                  variable.name = "Scenario", value.name = "Value")
    stats_long$Scenario <- gsub("_", "-", stats_long$Scenario)

    p_stats <- ggplot(stats_long, aes(x = Metric, y = Value, fill = Scenario)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
      geom_text(aes(label = sprintf("%.2f", Value)),
                position = position_dodge(width = 0.9),
                vjust = ifelse(stats_long$Value >= 0, -0.4, 1.2), size = 2.5) +
      scale_fill_manual(values = scenario_colours) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      labs(title = "Risk Metrics Comparison", x = "", y = "Value (%)") +
      theme_risk

    stressed_plots$portfolio_distribution_comparison <-
      (p_hist + p_density) / (p_qq + p_stats) +
      plot_annotation(
        title = "Portfolio Impact: Pre-Stress vs Post-Stress",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )

    cat("  ✅ Plot 12c: Pre/post distribution comparison\n")
  }

  # ===========================================================================
  # 13. RETURN RESULTS
  # ===========================================================================

  cat("\n========================================================================\n")
  cat("COMPILING RESULTS\n")
  cat("========================================================================\n")

  results <- list(

    # Metadata
    asset_names   = asset_names,
    n_assets      = n_assets,
    n_obs         = n_obs,
    n_weight_sims = n_weight_sims,
    n_sim_returns = n_sim_returns,

    # Configuration
    simulation_params = list(
      sim_returns_dist   = sim_returns_dist,
      cov_estimation     = cov_estimation,
      cov_lambda         = cov_lambda,
      copula             = copula,
      copula_type        = if(copula) copula_type else NULL,
      copula_df_fitted   = if(copula && copula_success && exists("copula_df_fitted")) copula_df_fitted else NULL,
      copula_df_input    = if(copula && copula_type == "t") copula_df else NULL,
      copula_success     = if(copula) copula_success else NULL,
      calibrate_params   = calibrate_params,
      calibration_method = if(calibrate_params) calibration_method else NULL,
      df_joint           = if(sim_returns_dist == "rmvt" && exists("df_joint")) df_joint else NULL,
      asset_lambdas      = if(sim_returns_dist == "rsgt") asset_lambdas else NULL,
      asset_p            = if(sim_returns_dist == "rsgt") asset_p else NULL,
      asset_q            = if(sim_returns_dist == "rsgt") asset_q else NULL,
      propagation        = propagation,
      asset_shock        = asset_shock,
      volatility_shock   = volatility_shock,
      reference_asset    = reference_asset,
      seed               = seed,
      frequency          = data_freq
    ),

    # Calibrated parameters
    calibrated_params = if(calibrate_params) {
      if(sim_returns_dist == "rmvt" && !is.null(calibrated_dfs))
        data.frame(Asset = asset_names, df = round(calibrated_dfs, 4))
      else if(sim_returns_dist == "rsgt" && !is.null(calibrated_lambdas))
        data.frame(Asset  = asset_names,
                   lambda = round(calibrated_lambdas, 4),
                   p      = round(calibrated_ps, 4),
                   q      = round(calibrated_qs, 4))
      else NULL
    } else NULL,

    # Fitted marginals
    fitted_params = fitted_params,

    # Asset empirical summary
    asset_empirical_summary = asset_empirical_summary,

    # Matrices
    historical_correlation = historical_cor,
    historical_covariance  = historical_cov,
    copula_correlation     = if(copula && copula_success) copula_cor else NULL,
    copula_vs_historical   = if(copula && copula_success) (copula_cor - historical_cor) else NULL,
    tail_dependence        = if(copula && copula_type == "t" && copula_success) tail_dependence else NULL,
    final_means            = base_mu,
    final_volatilities     = sqrt(diag(base_sigma)),
    final_covariance       = base_sigma,

    # Exogenous instruments
    exogenous_names                    = if(has_exogenous) exog_names else NULL,
    exogenous_means                    = if(has_exogenous) exog_means else NULL,
    exogenous_portfolio_crosscov       = if(has_exogenous) Sigma_pe_hist else NULL,
    exogenous_covariance               = if(has_exogenous) Sigma_ee_hist else NULL,
    exogenous_asset_empirical_summary  = if(has_exogenous) exogenous_asset_empirical_summary else NULL,
    exogenous_historical_correlation   = if(has_exogenous) exogenous_historical_correlation  else NULL,

    # Copula diagnostics
    ks_test = ks_results_df,

    # Stress
    is_stressed           = is_stressed,
    stressed_means        = if(is_stressed) mu_stressed    else NULL,
    stressed_covariance   = if(is_stressed) sigma_stressed else NULL,
    stress_impact_summary = if(is_stressed) stress_impact_summary else NULL,

    # Simulated returns
    portfolio_returns_pre          = portfolio_returns_pre_flat,
    portfolio_returns_pre_percent  = portfolio_returns_pre_flat * 100,
    portfolio_returns_post         = if(is_stressed) portfolio_returns_post_flat else NULL,
    portfolio_returns_post_percent = if(is_stressed) portfolio_returns_post_flat * 100 else NULL,
    weights_used                   = weights_matrix,

    # Per-asset simulated draws (n_sim_returns x n_assets data frames, first weight set)
    asset_draws        = asset_draws_pre,
    asset_draws_stress = if(is_stressed) asset_draws_post else NULL,

    # Full stats objects (retained for backward compatibility)
    stats_pre  = stats_pre,
    stats_post = if(is_stressed) stats_post else NULL,

    # Consolidated portfolio statistics data frame
    # Rows: Metric names | Columns: Pre_Stress, Post_Stress (NA when not stressed)
    portfolio_stats = {
      metrics <- c("Mean (%)", "Median (%)", "Volatility (%)", "Skewness",
                   "Excess Kurtosis", "Sharpe (Annual)", "Annual Return (%)",
                   "Annual Vol (%)",
                   "VaR 90% (%)", "VaR 95% (%)", "VaR 99% (%)",
                   "CVaR 90% (%)", "CVaR 95% (%)", "CVaR 99% (%)")
      pre_vals <- c(
        stats_pre$mean   * 100, stats_pre$median * 100, stats_pre$vol * 100,
        stats_pre$skewness, stats_pre$excess_kurtosis, stats_pre$sharpe_annual,
        stats_pre$annual_return * 100, stats_pre$annual_vol * 100,
        stats_pre$var_90 * 100, stats_pre$var_95 * 100, stats_pre$var_99 * 100,
        stats_pre$cvar_90 * 100, stats_pre$cvar_95 * 100, stats_pre$cvar_99 * 100
      )
      post_vals <- if(is_stressed) c(
        stats_post$mean   * 100, stats_post$median * 100, stats_post$vol * 100,
        stats_post$skewness, stats_post$excess_kurtosis, stats_post$sharpe_annual,
        stats_post$annual_return * 100, stats_post$annual_vol * 100,
        stats_post$var_90 * 100, stats_post$var_95 * 100, stats_post$var_99 * 100,
        stats_post$cvar_90 * 100, stats_post$cvar_95 * 100, stats_post$cvar_99 * 100
      ) else rep(NA_real_, length(metrics))
      data.frame(Metric = metrics, Pre_Stress = pre_vals,
                 Post_Stress = post_vals, stringsAsFactors = FALSE)
    },

    # Plots
    plots          = list(portfolio_diagnostics = portfolio_plots),
    stressed_plots = if(length(stressed_plots) > 0) stressed_plots else NULL
  )

  cat("\n========================================================================\n")
  cat("SIMULATION COMPLETE\n")
  cat("========================================================================\n")
  cat(sprintf("  Pre-stress samples  : %s\n",
              format(length(portfolio_returns_pre_flat), big.mark = ",")))
  if(is_stressed)
    cat(sprintf("  Post-stress samples : %s\n",
                format(length(portfolio_returns_post_flat), big.mark = ",")))
  cat(sprintf("  Diagnostic plots    : %d\n", length(portfolio_plots)))
  if(length(stressed_plots) > 0)
    cat(sprintf("  Stress plots        : %d\n", length(stressed_plots)))

  return(results)

}
