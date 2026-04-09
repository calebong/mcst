#' Decompose Portfolio Risk into Return and Risk Contributions
#'
#' @description
#' Decomposes a portfolio's return and risk into asset-level contributions
#' using simulation output from \code{\link{portfolio_risk_simulation}}.
#'
#' The function computes:
#' \itemize{
#'   \item \strong{Return contributions}: each asset's weighted central
#'     tendency (mean or median, depending on the fitted distribution).
#'   \item \strong{Standard deviation contributions} via the Euler homogeneity
#'     decomposition: \eqn{\sigma_p = \sum_i w_i \cdot (\boldsymbol{\Sigma}\mathbf{w})_i / \sigma_p}.
#'   \item \strong{Standalone VaR}: per-asset VaR in isolation.
#'   \item \strong{Marginal VaR}: \eqn{\text{MVaR}_i = \beta_i \cdot \text{VaR}_p},
#'     where \eqn{\beta_i} is estimated from simulations (full-distribution or
#'     tail-conditional, depending on \code{var_method}).
#'   \item \strong{Component VaR}: \eqn{\text{CVaR}_i = w_i \cdot \text{MVaR}_i}.
#'     Sums approximately to \eqn{\text{VaR}_p} by Euler's theorem.
#'   \item \strong{Incremental VaR}: the actual change in portfolio VaR when
#'     asset \eqn{i} is removed and remaining weights are renormalised.
#' }
#'
#' When \code{result} contains stressed draws, all measures are computed for
#' both pre-stress and post-stress scenarios and compared side by side.
#'
#' @section Return measure:
#' If \code{result$simulation_params$sim_returns_dist} is \code{"rsgt"}
#' (Skewed Generalised-t), the \strong{median} is used as the first-moment
#' estimator for return contributions, as the distribution may be substantially
#' skewed. For all other distributions (\code{"rmvnorm"}, \code{"rmvt"}),
#' the \strong{mean} is used. This choice propagates to table column names,
#' plot labels, and the \code{central_label} element of the returned list.
#'
#' @section VaR proportion methodology:
#' The \code{Component_VaR_prop_pct} column in \code{var_decomp} and the VaR
#' stacked bar in \code{plots$risk_decomposition} use a beta estimated according
#' to \code{var_method}:
#' \itemize{
#'   \item \code{"simulation"}: full-distribution beta
#'     \eqn{\beta_i = \text{Cov}(R_i, R_p) / \text{Var}(R_p)}. Mathematically
#'     equivalent to the sigma Euler proportion; the std dev and VaR stacked bars
#'     in \code{risk_decomposition} will therefore look identical.
#'   \item \code{"tail"}: tail-conditional beta estimated using only the
#'     simulated observations where \eqn{R_p \le \hat{V}}. Captures which assets
#'     co-move with the portfolio in the crisis tail; the two stacked bars will
#'     show genuinely different compositions.
#' }
#' Proportions are normalised to sum to 100\% in both cases.
#'
#' @param result A named list returned by \code{portfolio_risk_simulation()}.
#'   Must contain \code{result$asset_draws} (pre-stress per-asset simulation
#'   draws, an \eqn{m \times n} data frame) and
#'   \code{result$portfolio_returns_pre} (a numeric vector of length \eqn{m}).
#'   For stressed decomposition, \code{result$asset_draws_stress} and
#'   \code{result$portfolio_returns_post} must also be present.
#' @param weights Named numeric vector of portfolio weights summing to 1.
#'   If \code{NULL} (default), equal weights \eqn{1/n} are assumed for all
#'   \eqn{n} assets. Names must match column names of \code{result$asset_draws}.
#'   Weights that deviate from sum-to-1 by more than \code{1e-6} are silently
#'   renormalised; deviations greater than 5\\% raise an error.
#' @param focus_assets Character vector of asset names to highlight in plots and
#'   printed tables. Only these assets are shown in the filtered output. If
#'   \code{NULL} (default), all assets are included.
#' @param confidence_level Numeric scalar in \eqn{(0,1)}. VaR confidence level.
#'   Default \code{0.95} (i.e. the 5th percentile of simulated returns).
#' @param var_method Character scalar. Method used to estimate beta for marginal
#'   and component VaR. One of:
#'   \describe{
#'     \item{\code{"simulation"} (default)}{Full-distribution beta:
#'       \eqn{\beta_i = \text{Cov}(R_i, R_p) / \text{Var}(R_p)} computed over
#'       all \eqn{m} simulated observations.}
#'     \item{\code{"tail"}}{Tail-conditional beta: same formula but restricted
#'       to observations where \eqn{R_p \le \text{VaR}_p}. Falls back to
#'       full-distribution beta with a warning if fewer than 10 tail observations
#'       are available.}
#'   }
#' @param plot Logical. If \code{TRUE} (default), generate all diagnostic plots
#'   and return them in \code{result$plots}.
#' @param digits Integer. Number of decimal places for printed tables.
#'   Default \code{4}.
#'
#' @return A named list (returned invisibly) with the following elements:
#' \describe{
#'   \item{\code{return_decomp}}{If stressed: a list with elements \code{pre}
#'     and \code{post}, each a data frame with columns \code{Asset},
#'     \code{Weight_pct}, \code{Asset_Mean_pct} (or \code{Asset_Median_pct}),
#'     and \code{Return_Contrib_Mean_pct} (or \code{Return_Contrib_Median_pct}).
#'     If not stressed: a single data frame.}
#'   \item{\code{risk_decomp}}{If stressed: a list with \code{pre} and
#'     \code{post}, each a data frame with columns \code{Asset},
#'     \code{Weight_pct}, \code{Asset_Vol_pct}, \code{MCTR_pct},
#'     \code{Component_Risk_pct}, \code{Component_Risk_prop}.
#'     If not stressed: a single data frame.}
#'   \item{\code{var_decomp}}{If stressed: a list with \code{pre} and
#'     \code{post}, each a data frame with columns \code{Asset},
#'     \code{Weight_pct}, \code{Asset_Vol_pct}, \code{Standalone_VaR_pct},
#'     \code{Beta}, \code{Marginal_VaR_pct}, \code{Component_VaR_pct},
#'     \code{Component_VaR_prop_pct}, \code{Incremental_VaR_pct}.
#'     If not stressed: a single data frame.}
#'   \item{\code{summary}}{Wide-format data frame combining all measures for
#'     all assets and scenarios. Contains all columns from the three tables
#'     above plus \code{Scenario}.}
#'   \item{\code{port_stats}}{If stressed: a list with \code{pre} and
#'     \code{post}, each a named list containing \code{scenario},
#'     \code{n_sim}, \code{central_pct}, \code{central_label}, \code{vol_pct},
#'     \code{var_pct}, \code{cvar_pct}, \code{euler_sum_pct}, and
#'     \code{euler_error_pct}. If not stressed: a single named list.}
#'   \item{\code{plots}}{Named list of ggplot/patchwork objects (only present
#'     when \code{plot = TRUE}):
#'     \describe{
#'       \item{\code{var_comparison}}{Grouped bar chart of all four VaR measures
#'         (standalone, marginal, component, incremental) per asset, faceted by
#'         scenario. Bar values are overlaid as text labels.}
#'       \item{\code{var_heatmap}}{Heatmap of standalone, component, and
#'         incremental VaR per asset across pre- and post-stress scenarios.
#'         Only produced when stressed data is present.}
#'       \item{\code{risk_decomposition}}{2×2 patchwork canvas (only when
#'         stressed) showing: (top-left) stacked bar of std dev composition
#'         by asset; (top-right) delta in std dev share per asset; (bottom-left)
#'         stacked bar of VaR composition by asset; (bottom-right) delta in VaR
#'         share per asset. Colours are generated by \code{make_qual_pal()}.}
#'       \item{\code{summary_canvas}}{Patchwork combining \code{var_comparison}
#'         and \code{var_heatmap} (when stressed) or \code{var_comparison}
#'         alone (when not stressed).}
#'     }
#'   }
#'   \item{\code{weights}}{Named numeric vector of weights used.}
#'   \item{\code{confidence_level}}{Numeric. The confidence level used.}
#'   \item{\code{var_method}}{Character. The beta method used.}
#'   \item{\code{is_stressed}}{Logical. Whether stressed draws were available.}
#'   \item{\code{focus_assets}}{Character vector or \code{NULL}.}
#'   \item{\code{asset_names}}{Character vector of all asset names.}
#'   \item{\code{n_sim}}{Integer. Number of simulations used.}
#'   \item{\code{central_label}}{Character. \code{"Mean"} or \code{"Median"}.}
#'   \item{\code{sim_dist}}{Character. The distribution string read from
#'     \code{result$simulation_params$sim_returns_dist}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Run simulation first
#' result <- portfolio_risk_simulation(returns_df)
#'
#' # Basic decomposition — equal weights, 95% VaR, full-distribution beta
#' decomp <- decompose_portfolio_risk(result)
#' decomp$var_decomp
#' decomp$plots$var_comparison
#' decomp$plots$summary_canvas
#'
#' # Custom weights, tail-conditional beta, focus on two assets
#' w <- c(NVDA=0.15, TSM=0.10, JPM=0.12, AXP=0.10,
#'        LMT=0.10, BA=0.10, UAL=0.08, XOM=0.10, NGG=0.10, PG=0.05)
#' decomp <- decompose_portfolio_risk(
#'   result,
#'   weights          = w,
#'   focus_assets     = c("NVDA", "XOM"),
#'   confidence_level = 0.99,
#'   var_method       = "tail"
#' )
#' decomp$plots$risk_decomposition
#'
#' # Stressed decomposition (result must contain asset_draws_stress)
#' result_stressed <- portfolio_risk_simulation(
#'   returns_df,
#'   asset_shock = list(NVDA = -0.10)
#' )
#' decomp_s <- decompose_portfolio_risk(result_stressed, var_method = "tail")
#' decomp_s$plots$risk_decomposition   # shows pre vs post stress composition
#' decomp_s$plots$var_heatmap
#' }
#'
#' @seealso \code{\link{portfolio_risk_simulation}}
#'
#' @importFrom ggplot2 ggplot aes aes_string geom_col geom_hline geom_text
#'   geom_tile facet_wrap coord_flip position_dodge position_stack
#'   scale_fill_manual scale_fill_gradient2 scale_y_continuous
#'   theme_minimal theme element_text labs
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom grDevices col2rgb rgb2hsv hsv
#' @importFrom stats cov var sd median quantile
#'
#' @export
decompose_portfolio_risk <- function(
    result,
    weights          = NULL,
    focus_assets     = NULL,
    confidence_level = 0.95,
    var_method       = c("simulation", "tail"),
    plot             = TRUE,
    digits           = 4
) {

  var_method <- match.arg(var_method)

  # ---------------------------------------------------------------------------
  # 0. Input validation
  # ---------------------------------------------------------------------------
  if (!is.list(result))
    stop("'result' must be a list returned by portfolio_risk_simulation().")
  if (is.null(result$asset_draws))
    stop("result$asset_draws is NULL. Re-run portfolio_risk_simulation() to generate per-asset draws.")

  asset_draws_pre <- as.data.frame(result$asset_draws)
  asset_names     <- colnames(asset_draws_pre)
  n_assets        <- ncol(asset_draws_pre)
  n_sim           <- nrow(asset_draws_pre)

  is_stressed  <- !is.null(result$asset_draws_stress) &&
                  !is.null(result$portfolio_returns_post)

  asset_draws_post <- if (is_stressed) as.data.frame(result$asset_draws_stress) else NULL
  port_pre         <- as.numeric(result$portfolio_returns_pre)
  port_post        <- if (is_stressed) as.numeric(result$portfolio_returns_post) else NULL

  # ---------------------------------------------------------------------------
  # MOD 3: detect distribution from result$simulation_params
  # Use median for rsgt (skewed), mean otherwise
  # ---------------------------------------------------------------------------
  sim_dist <- tryCatch(result$simulation_params$sim_returns_dist, error = function(e) "rmvnorm")
  if (is.null(sim_dist)) sim_dist <- "rmvnorm"
  use_median   <- sim_dist == "rsgt"
  central_label <- if (use_median) "Median" else "Mean"

  cat(sprintf("  Return central tendency: %s (distribution: %s)\n",
              central_label, sim_dist))

  # Helper: apply correct central tendency to a vector
  central <- function(x) if (use_median) median(x) else mean(x)

  # ---------------------------------------------------------------------------
  # 1. Weights
  # ---------------------------------------------------------------------------
  if (is.null(weights)) {
    w <- setNames(rep(1 / n_assets, n_assets), asset_names)
    cat(sprintf("  decompose_portfolio_risk: equal weights 1/%d applied\n", n_assets))
  } else {
    if (!is.numeric(weights) || length(weights) != n_assets)
      stop(sprintf("'weights' must be a numeric vector of length %d.", n_assets))
    if (!is.null(names(weights))) weights <- weights[asset_names]
    if (abs(sum(weights) - 1) > 0.05)
      stop(sprintf("weights sum to %.4f - must sum to 1.", sum(weights)))
    if (abs(sum(weights) - 1) > 1e-6) weights <- weights / sum(weights)
    w <- setNames(as.numeric(weights), asset_names)
  }

  # ---------------------------------------------------------------------------
  # 2. Focus assets validation
  # ---------------------------------------------------------------------------
  if (!is.null(focus_assets)) {
    bad <- setdiff(focus_assets, asset_names)
    if (length(bad) > 0)
      stop(sprintf("focus_assets not found: %s", paste(bad, collapse = ", ")))
  }

  alpha <- 1 - confidence_level

  # ---------------------------------------------------------------------------
  # 3. Core decomposition function
  # ---------------------------------------------------------------------------
  compute_decomp <- function(asset_mat, port_vec, scenario_label) {

    asset_mat <- as.matrix(asset_mat)
    colnames(asset_mat) <- asset_names

    # ── 3a. Return contributions ────────────────────────────────────────────
    # MOD 3: use median or mean depending on distribution
    ct_assets   <- apply(asset_mat, 2, central)   # per-asset central tendency
    ct_port     <- central(port_vec)              # portfolio central tendency
    ret_contrib <- w * ct_assets                  # weighted contribution

    # MOD 1: return_df keeps only Weight_pct and the central tendency column
    # Column name reflects whether mean or median is used
    ct_col_name <- paste0("Asset_", central_label, "_pct")

    return_df <- data.frame(
      Asset         = asset_names,
      Weight_pct    = w * 100,
      stringsAsFactors = FALSE
    )
    return_df[[ct_col_name]] <- ct_assets * 100
    # Absolute contribution: w_i * central_i (in %)
    return_df[[paste0("Return_Contrib_", central_label, "_pct")]] <- ret_contrib * 100

    # ── 3b. Std dev contributions (Euler decomposition) ─────────────────────
    Sigma   <- cov(asset_mat)
    sigma_p <- sd(port_vec)
    Sigma_w <- as.vector(Sigma %*% w)
    mctr    <- Sigma_w / sigma_p
    ccr     <- w * mctr
    ccr_prop <- ccr / sigma_p

    risk_df <- data.frame(
      Asset               = asset_names,
      Weight_pct          = w * 100,
      Asset_Vol_pct       = sqrt(diag(Sigma)) * 100,
      MCTR_pct            = mctr * 100,
      Component_Risk_pct  = ccr * 100,
      Component_Risk_prop = ccr_prop * 100,
      stringsAsFactors    = FALSE
    )

    # ── 3c. VaR quantities ──────────────────────────────────────────────────
    var_p <- quantile(port_vec, alpha, names = FALSE)

    standalone_var_pct <- apply(asset_mat, 2, function(x)
      quantile(x, alpha, names = FALSE)) * 100

    if (var_method == "simulation") {
      cov_ip <- apply(asset_mat, 2, function(ri) cov(ri, port_vec))
      beta   <- cov_ip / var(port_vec)
    } else {
      tail_idx <- port_vec <= var_p
      if (sum(tail_idx) < 10) {
        warning("Fewer than 10 tail observations; using full-distribution beta.")
        cov_ip <- apply(asset_mat, 2, function(ri) cov(ri, port_vec))
        beta   <- cov_ip / var(port_vec)
      } else {
        port_tail  <- port_vec[tail_idx]
        asset_tail <- asset_mat[tail_idx, , drop = FALSE]
        cov_ip     <- apply(asset_tail, 2, function(ri) cov(ri, port_tail))
        beta       <- cov_ip / var(port_tail)
      }
    }

    marginal_var_pct   <- beta * var_p * 100
    component_var_pct  <- w * marginal_var_pct

    # VaR proportion: which beta to use depends on var_method.
    #   "simulation": full-distribution beta → same formula as sigma Euler proportion,
    #     meaning sigma and VaR bars will look identical in risk_decomposition.
    #     Appropriate when you want a consistent full-distribution view.
    #   "tail": tail-conditional beta (obs where R_p <= VaR_p) → captures which
    #     assets co-move with the portfolio specifically in the crisis tail.
    #     Appropriate when VaR and sigma should tell genuinely different stories.
    # Both are normalised so that asset proportions sum to 100%.
    beta_for_prop <- if (var_method == "tail") {
      tail_idx_var <- port_vec <= var_p
      if (sum(tail_idx_var) >= 10) {
        port_tail_v  <- port_vec[tail_idx_var]
        asset_tail_v <- asset_mat[tail_idx_var, , drop = FALSE]
        cov_tail     <- apply(asset_tail_v, 2, function(ri) cov(ri, port_tail_v))
        cov_tail / var(port_tail_v)
      } else {
        warning("Fewer than 10 tail observations for VaR proportion; using full-distribution beta.")
        beta
      }
    } else {
      beta  # "simulation": full-distribution beta
    }
    raw_var_prop       <- w * beta_for_prop
    prop_sum           <- sum(raw_var_prop)
    component_var_prop <- if (abs(prop_sum) > 1e-10)
      raw_var_prop / prop_sum * 100 else raw_var_prop * 100

    euler_sum_pct <- sum(component_var_pct)

    incremental_var_pct <- vapply(seq_len(n_assets), function(i) {
      w_excl <- w; w_excl[i] <- 0
      w_sum  <- sum(w_excl)
      if (w_sum < 1e-10) return(NA_real_)
      w_excl   <- w_excl / w_sum
      var_excl <- quantile(as.vector(asset_mat %*% w_excl), alpha, names = FALSE)
      (var_p - var_excl) * 100
    }, numeric(1))
    names(incremental_var_pct) <- asset_names

    var_df <- data.frame(
      Asset                  = asset_names,
      Weight_pct             = w * 100,
      Asset_Vol_pct          = sqrt(diag(Sigma)) * 100,
      Standalone_VaR_pct     = standalone_var_pct,
      Beta                   = beta,
      Marginal_VaR_pct       = marginal_var_pct,
      Component_VaR_pct      = component_var_pct,
      Component_VaR_prop_pct = component_var_prop,
      Incremental_VaR_pct    = incremental_var_pct,
      stringsAsFactors       = FALSE
    )

    port_summary <- list(
      scenario        = scenario_label,
      n_sim           = length(port_vec),
      central_pct     = ct_port * 100,         # mean or median of portfolio
      central_label   = central_label,
      vol_pct         = sigma_p * 100,
      var_pct         = var_p * 100,
      cvar_pct        = mean(port_vec[port_vec <= var_p]) * 100,
      euler_sum_pct   = euler_sum_pct,
      euler_error_pct = (euler_sum_pct - var_p * 100) / abs(var_p * 100) * 100
    )

    list(return_df = return_df, risk_df = risk_df,
         var_df = var_df, port_summary = port_summary)
  }

  # ---------------------------------------------------------------------------
  # 4. Run
  # ---------------------------------------------------------------------------
  cat("\n========================================================================\n")
  cat("PORTFOLIO RISK DECOMPOSITION\n")
  cat("========================================================================\n")
  cat(sprintf("  Assets        : %d - %s\n", n_assets, paste(asset_names, collapse = ", ")))
  cat(sprintf("  Simulations   : %s\n", format(n_sim, big.mark = ",")))
  cat(sprintf("  VaR confidence: %.0f%%\n", confidence_level * 100))
  cat(sprintf("  Beta method   : %s\n", var_method))
  cat(sprintf("  Return measure: %s\n", central_label))
  cat(sprintf("  Stressed      : %s\n", ifelse(is_stressed, "YES", "NO")))
  if (!is.null(focus_assets))
    cat(sprintf("  Focus assets  : %s\n", paste(focus_assets, collapse = ", ")))

  res_pre  <- compute_decomp(asset_draws_pre, port_pre, "Pre-Stress")
  res_post <- if (is_stressed)
    compute_decomp(asset_draws_post, port_post, "Post-Stress") else NULL

  # ---------------------------------------------------------------------------
  # 5. Print tables
  # ---------------------------------------------------------------------------
  print_section <- function(res) {
    ps <- res$port_summary
    cat(sprintf("\n--- %s ---\n", ps$scenario))
    cat(sprintf("  %s: %+.4f%%  Vol: %.4f%%  VaR%.0f: %.4f%%  CVaR: %.4f%%\n",
                ps$central_label, ps$central_pct, ps$vol_pct,
                confidence_level * 100, ps$var_pct, ps$cvar_pct))
    cat(sprintf("  Euler: sum(CVaR)=%.4f%% vs VaR=%.4f%% (error=%.2f%%)\n\n",
                ps$euler_sum_pct, ps$var_pct, ps$euler_error_pct))

    show <- function(df, title) {
      cat(title, "\n")
      d <- if (!is.null(focus_assets)) df[df$Asset %in% focus_assets, ] else df
      print(d, digits = digits, row.names = FALSE)
      cat("\n")
    }
    show(res$return_df, sprintf("Return Contributions (%s-based):", central_label))
    show(res$risk_df,   "Std Dev Contributions (Euler):")
    show(res$var_df,    sprintf("VaR Decomposition (%.0f%%):", confidence_level * 100))
  }

  print_section(res_pre)
  if (is_stressed) print_section(res_post)

  # ---------------------------------------------------------------------------
  # 6. Build combined summary
  # Dynamic column name for the central tendency (mean vs median)
  # ---------------------------------------------------------------------------
  ct_col        <- paste0("Asset_", central_label, "_pct")
  rc_col        <- paste0("Return_Contrib_", central_label, "_pct")

  build_summary <- function(rp, rk, vr, scenario_label) {
    df <- data.frame(
      Scenario               = scenario_label,
      Asset                  = rp$Asset,
      Weight_pct             = rp$Weight_pct,
      Asset_Vol_pct          = rk$Asset_Vol_pct,
      MCTR_pct               = rk$MCTR_pct,
      Component_Risk_pct     = rk$Component_Risk_pct,
      Component_Risk_prop    = rk$Component_Risk_prop,
      Standalone_VaR_pct     = vr$Standalone_VaR_pct,
      Beta                   = vr$Beta,
      Marginal_VaR_pct       = vr$Marginal_VaR_pct,
      Component_VaR_pct      = vr$Component_VaR_pct,
      Component_VaR_prop_pct = vr$Component_VaR_prop_pct,
      Incremental_VaR_pct    = vr$Incremental_VaR_pct,
      stringsAsFactors       = FALSE
    )
    # Attach central tendency and contribution columns with dynamic names
    df[[ct_col]] <- rp[[ct_col]]
    df[[rc_col]] <- rp[[rc_col]]
    df
  }

  summary_pre  <- build_summary(res_pre$return_df, res_pre$risk_df,
                                res_pre$var_df, "Pre-Stress")
  summary_post <- if (is_stressed)
    build_summary(res_post$return_df, res_post$risk_df,
                  res_post$var_df, "Post-Stress") else NULL
  summary_all  <- if (is_stressed) rbind(summary_pre, summary_post) else summary_pre

  # ---------------------------------------------------------------------------
  # 7. Plots
  # ---------------------------------------------------------------------------
  plots <- list()

  if (plot) {

    theme_d <- theme_minimal(base_size = 10) +
      theme(plot.title      = element_text(face = "bold", size = 11, hjust = 0),
            plot.subtitle   = element_text(size = 8.5, color = "gray40"),
            axis.text.x     = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")

    sc_cols <- c("Pre-Stress" = "steelblue", "Post-Stress" = "coral")

    pdf <- if (!is.null(focus_assets))
      summary_all[summary_all$Asset %in% focus_assets, ] else summary_all

    # ── Plot: All four VaR types faceted by scenario, with bar labels ─────────
    vtdf <- do.call(rbind, lapply(unique(summary_all$Scenario), function(sc) {
      sub <- summary_all[summary_all$Scenario == sc, ]
      if (!is.null(focus_assets)) sub <- sub[sub$Asset %in% focus_assets, ]
      rbind(
        data.frame(Asset=sub$Asset, Scenario=sc, VaR_Type="Standalone",
                   Value=sub$Standalone_VaR_pct, stringsAsFactors=FALSE),
        data.frame(Asset=sub$Asset, Scenario=sc, VaR_Type="Marginal",
                   Value=sub$Marginal_VaR_pct, stringsAsFactors=FALSE),
        data.frame(Asset=sub$Asset, Scenario=sc, VaR_Type="Component",
                   Value=sub$Component_VaR_pct, stringsAsFactors=FALSE),
        data.frame(Asset=sub$Asset, Scenario=sc, VaR_Type="Incremental",
                   Value=sub$Incremental_VaR_pct, stringsAsFactors=FALSE)
      )
    }))
    vtdf$VaR_Type <- factor(vtdf$VaR_Type,
                            levels = c("Standalone","Marginal","Component","Incremental"))
    plots$var_comparison <- ggplot(
      vtdf, aes(x = Asset, y = Value, fill = VaR_Type)) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(
        aes(label = sprintf("%.2f", Value),
            vjust = ifelse(Value >= 0, -0.35, 1.25)),
        position = position_dodge(width = 0.9),
        size = 2.2, color = "gray20") +
      scale_fill_manual(values = c("Standalone"="#d6604d","Marginal"="#f4a582",
                                   "Component"="#4393c3","Incremental"="#2166ac")) +
      facet_wrap(~Scenario, ncol = 1) +
      labs(title = sprintf("VaR Measures by Asset (%.0f%%)", confidence_level * 100),
           subtitle = "Standalone: solo | Marginal: beta x VaR_p | Component: w x MVaR | Incremental: VaR_p - VaR_excl",
           x = "Asset", y = "VaR (%)") + theme_d

    # ── var_heatmap: works for both stressed and unstressed ─────────────────
    # When unstressed: shows only Pre-Stress columns (single scenario)
    # When stressed: shows both Pre-Stress and Post-Stress columns
    hdf <- do.call(rbind, lapply(unique(summary_all$Scenario), function(sc) {
      sub <- summary_all[summary_all$Scenario == sc, ]
      if (!is.null(focus_assets)) sub <- sub[sub$Asset %in% focus_assets, ]
      data.frame(
        Asset    = rep(sub$Asset, 3),
        Measure  = rep(c("Standalone","Component","Incremental"), each = nrow(sub)),
        Scenario = sc,
        Value    = c(sub$Standalone_VaR_pct, sub$Component_VaR_pct, sub$Incremental_VaR_pct),
        stringsAsFactors = FALSE
      )
    }))
    hdf$Panel <- paste(hdf$Scenario, hdf$Measure, sep = "\n")
    heatmap_title <- if (is_stressed) "VaR Heatmap: Pre vs Post Stress"
                     else sprintf("VaR Heatmap (%.0f%% Confidence)", confidence_level * 100)
    plots$var_heatmap <- ggplot(hdf, aes(x = Panel, y = Asset, fill = Value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#d6604d", mid = "white", high = "#4393c3",
                           midpoint = 0, name = "VaR (%)") +
      geom_text(aes(label = sprintf("%.3f", Value)), size = 2.8) +
      labs(title = heatmap_title, x = "", y = "Asset") +
      theme_d + theme(axis.text.x = element_text(angle = 30, hjust = 1))

    # ── Internal helper plots for summary_canvas (not added to plots list) ───
    # Return contributions
    p_ret_internal <- ggplot(
      pdf, aes_string(x = sprintf("reorder(Asset, -abs(%s))", rc_col),
                      y = rc_col, fill = "Scenario")) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(aes_string(label = sprintf('sprintf("%%+.3f%%%%", %s)', rc_col)),
                position = position_dodge(width = 0.9),
                vjust = ifelse(pdf[[rc_col]] >= 0, -0.4, 1.2), size = 2.5) +
      scale_fill_manual(values = sc_cols) +
      labs(title    = sprintf("Return Contributions (%s-based)", central_label),
           subtitle = sprintf("w_i x %s_i", tolower(central_label)),
           x = "Asset", y = sprintf("Return Contribution (%s) (%%)", central_label)) +
      theme_d

    # Component std dev
    p_risk_internal <- ggplot(
      pdf, aes(x = reorder(Asset, -Component_Risk_pct),
               y = Component_Risk_pct, fill = Scenario)) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(aes(label = sprintf("%.3f%%", Component_Risk_pct)),
                position = position_dodge(width = 0.9),
                vjust = -0.4, size = 2.5) +
      scale_fill_manual(values = sc_cols) +
      labs(title    = "Component Std Dev (Euler)",
           subtitle = "w_i x (Sigma w)_i / sigma_p",
           x = "Asset", y = "Component Std Dev (%)") +
      theme_d

    # Component VaR horizontal bar
    pdf$VaR_Dir <- ifelse(pdf$Component_VaR_pct <= 0, "Risk", "Hedge")
    p_cvar_internal <- ggplot(
      pdf, aes(x = reorder(Asset, Component_VaR_pct),
               y = Component_VaR_pct,
               fill = interaction(VaR_Dir, Scenario))) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(aes(label = sprintf("%.3f%%", Component_VaR_pct)),
                position = position_dodge(width = 0.9),
                hjust = ifelse(pdf$Component_VaR_pct <= 0, 1.1, -0.1), size = 2.5) +
      scale_fill_manual(
        values = c("Risk.Pre-Stress"   = "#d6604d", "Risk.Post-Stress"  = "#b22222",
                   "Hedge.Pre-Stress"  = "#4393c3", "Hedge.Post-Stress" = "#2166ac"),
        name = "") +
      coord_flip() +
      labs(title    = sprintf("Component VaR (%.0f%%)", confidence_level * 100),
           subtitle = "Red = risk | Blue = diversifier",
           x = "Asset", y = "Component VaR (%)") +
      theme_d

    # % of portfolio VaR
    p_varprop_internal <- ggplot(
      pdf, aes(x = Asset, y = Component_VaR_prop_pct, fill = Scenario)) +
      geom_col(position = "dodge", alpha = 0.85) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_text(aes(label = sprintf("%.1f%%", Component_VaR_prop_pct)),
                position = position_dodge(width = 0.9),
                vjust = ifelse(pdf$Component_VaR_prop_pct >= 0, -0.4, 1.2), size = 2.5) +
      scale_fill_manual(values = sc_cols) +
      labs(title    = "% of Portfolio VaR by Asset",
           subtitle = "w_i x beta_i x 100 | Positive = risk contributor | Negative = diversifier",
           x = "Asset", y = "% of Portfolio VaR") +
      theme_d

    # ── summary_canvas: restored 2x2 grid matching original layout ──────────
    plots$summary_canvas <- wrap_plots(
      list(p_ret_internal, p_risk_internal,
           p_cvar_internal, p_varprop_internal), ncol = 2) +
      plot_annotation(
        title    = sprintf("Portfolio Risk Decomposition (VaR %.0f%%)", confidence_level * 100),
        subtitle = sprintf("%d assets | %s simulations | beta: %s | return: %s",
                           n_assets, format(n_sim, big.mark = ","), var_method, central_label),
        theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
                         plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"))
      )

    # ── Qualitative palette generator (hoisted: used in both stressed + unstressed paths) ──
    make_qual_pal <- function(asset_names_sorted) {
      base_hex <- c("#E63946","#2196F3","#FF9800","#4CAF50","#9C27B0",
                    "#00BCD4","#FF5722","#3F51B5","#CDDC39","#009688")
      n        <- length(asset_names_sorted)
      n_base   <- length(base_hex)
      colours  <- character(n)
      for (k in seq_len(n)) {
        base_idx <- ((k - 1L) %% n_base) + 1L
        cycle    <- (k - 1L) %/% n_base
        if (cycle == 0L) {
          colours[k] <- base_hex[base_idx]
        } else {
          rgb_vals   <- col2rgb(base_hex[base_idx]) / 255
          hsv_vals   <- rgb2hsv(rgb_vals[1], rgb_vals[2], rgb_vals[3])
          new_h      <- (hsv_vals[1] + cycle * 0.042) %% 1
          new_s      <- pmin(hsv_vals[2] * (1 - cycle * 0.08), 1)
          new_v      <- pmax(hsv_vals[3] * (1 - cycle * 0.12), 0.3)
          colours[k] <- hsv(new_h, new_s, new_v)
        }
      }
      setNames(colours, asset_names_sorted)
    }

    # ── risk_decomposition 2x2 canvas ───────────────────────────────────────
    # When stressed: shows pre vs post stacked bars + delta bars
    # When not stressed: shows pre-only stacked bar (single bar) + note panels
    if (is_stressed) {

      # Build long-form stacked data for std dev and VaR
      risk_stack <- rbind(
        data.frame(Scenario="Pre-Stress",  Asset=res_pre$risk_df$Asset,
                   Component_Risk_prop=res_pre$risk_df$Component_Risk_prop,
                   stringsAsFactors=FALSE),
        data.frame(Scenario="Post-Stress", Asset=res_post$risk_df$Asset,
                   Component_Risk_prop=res_post$risk_df$Component_Risk_prop,
                   stringsAsFactors=FALSE)
      )
      risk_stack$Scenario <- factor(risk_stack$Scenario,
                                    levels = c("Pre-Stress","Post-Stress"))

      var_stack <- rbind(
        data.frame(Scenario="Pre-Stress",  Asset=res_pre$var_df$Asset,
                   Component_VaR_prop_pct=res_pre$var_df$Component_VaR_prop_pct,
                   stringsAsFactors=FALSE),
        data.frame(Scenario="Post-Stress", Asset=res_post$var_df$Asset,
                   Component_VaR_prop_pct=res_post$var_df$Component_VaR_prop_pct,
                   stringsAsFactors=FALSE)
      )
      var_stack$Scenario <- factor(var_stack$Scenario,
                                   levels = c("Pre-Stress","Post-Stress"))

      # Delta data for std dev and VaR
      risk_delta <- merge(
        res_pre$risk_df[,  c("Asset","Component_Risk_prop")],
        res_post$risk_df[, c("Asset","Component_Risk_prop")],
        by = "Asset", suffixes = c("_pre","_post")
      )
      risk_delta$Delta <- risk_delta$Component_Risk_prop_post -
                          risk_delta$Component_Risk_prop_pre

      var_delta <- merge(
        res_pre$var_df[,  c("Asset","Component_VaR_prop_pct")],
        res_post$var_df[, c("Asset","Component_VaR_prop_pct")],
        by = "Asset", suffixes = c("_pre","_post")
      )
      var_delta$Delta <- var_delta$Component_VaR_prop_pct_post -
                         var_delta$Component_VaR_prop_pct_pre

      # Apply focus_assets filter if set
      if (!is.null(focus_assets)) {
        risk_stack <- risk_stack[risk_stack$Asset %in% focus_assets, ]
        var_stack  <- var_stack[var_stack$Asset %in% focus_assets, ]
        risk_delta <- risk_delta[risk_delta$Asset %in% focus_assets, ]
        var_delta  <- var_delta[var_delta$Asset %in% focus_assets, ]
      }

      # Top-left: stacked bar — std dev composition
      n_pal     <- length(unique(risk_stack$Asset))
      qual_pal  <- make_qual_pal(sort(unique(risk_stack$Asset)))

      p_rd_stack <- ggplot(
        risk_stack,
        aes(x = Scenario, y = Component_Risk_prop, fill = Asset)) +
        geom_col(position = "stack", alpha = 0.88, width = 0.55) +
        geom_text(aes(label = sprintf("%.1f%%", Component_Risk_prop)),
                  position = position_stack(vjust = 0.5), size = 2.5, color = "white",
                  fontface = "bold") +
        scale_fill_manual(values = qual_pal) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(title    = "Std Dev Composition by Asset",
             subtitle = "Component_Risk_prop | Bars sum to 100%",
             x = "", y = "% of Portfolio Sigma") +
        theme_d +
        theme(legend.position = "right")

      # Top-right: delta std dev composition
      p_rd_delta <- ggplot(
        risk_delta,
        aes(x = reorder(Asset, Delta), y = Delta, fill = Delta > 0)) +
        geom_col(alpha = 0.85) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(aes(label = sprintf("%+.2f%%", Delta)),
                  hjust = ifelse(risk_delta$Delta > 0, -0.1, 1.1), size = 2.8) +
        coord_flip() +
        scale_fill_manual(
          values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
          labels = c("TRUE" = "Risk share increased", "FALSE" = "Risk share reduced"),
          name   = "") +
        labs(title    = "Change in Std Dev Share (Stressed - Unstressed)",
             subtitle = "Post minus Pre | Positive = larger risk share under stress",
             x = "Asset", y = "Delta Component_Risk_prop (pp)") +
        theme_d

      # Bottom-left: stacked bar — VaR composition
      n_pal_v   <- length(unique(var_stack$Asset))
      qual_pal_v <- make_qual_pal(sort(unique(var_stack$Asset)))

      p_vd_stack <- ggplot(
        var_stack,
        aes(x = Scenario, y = Component_VaR_prop_pct, fill = Asset)) +
        geom_col(position = "stack", alpha = 0.88, width = 0.55) +
        geom_text(aes(label = sprintf("%.1f%%", Component_VaR_prop_pct)),
                  position = position_stack(vjust = 0.5), size = 2.5, color = "white",
                  fontface = "bold") +
        scale_fill_manual(values = qual_pal_v) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(title    = sprintf("VaR %.0f%% Composition by Asset", confidence_level * 100),
             subtitle = sprintf("Beta method: %s | Proportions normalised to 100%%", var_method),
             x = "", y = sprintf("%% of Portfolio VaR %.0f%%", confidence_level * 100)) +
        theme_d +
        theme(legend.position = "right")

      # Bottom-right: delta VaR composition
      p_vd_delta <- ggplot(
        var_delta,
        aes(x = reorder(Asset, Delta), y = Delta, fill = Delta > 0)) +
        geom_col(alpha = 0.85) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        geom_text(aes(label = sprintf("%+.2f%%", Delta)),
                  hjust = ifelse(var_delta$Delta > 0, -0.1, 1.1), size = 2.8) +
        coord_flip() +
        scale_fill_manual(
          values = c("TRUE" = "#d6604d", "FALSE" = "#4393c3"),
          labels = c("TRUE" = "VaR share increased", "FALSE" = "VaR share reduced"),
          name   = "") +
        labs(title    = sprintf("Change in VaR %.0f%% Share (Stressed - Unstressed)",
                                confidence_level * 100),
             subtitle = "Post minus Pre | Positive = larger VaR share under stress",
             x = "Asset",
             y = sprintf("Delta Component_VaR_prop (pp, VaR %.0f%%)", confidence_level * 100)) +
        theme_d

      plots$risk_decomposition <- wrap_plots(
        list(p_rd_stack, p_rd_delta, p_vd_stack, p_vd_delta), ncol = 2) +
        plot_annotation(
          title    = sprintf(
            "Risk Decomposition: Std Dev and VaR %.0f%% Composition (Pre vs Post Stress)",
            confidence_level * 100),
          subtitle = sprintf(
            "%d assets | %s simulations | Beta method: %s | Return measure: %s",
            n_assets, format(n_sim, big.mark = ","), var_method, central_label),
          theme    = theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40")
          )
        )

    } else {
      # Non-stressed: single-bar stacked composition (one bar per chart)
      risk_stack_ns <- data.frame(
        Scenario            = "Pre-Stress",
        Asset               = res_pre$risk_df$Asset,
        Component_Risk_prop = res_pre$risk_df$Component_Risk_prop,
        stringsAsFactors    = FALSE
      )
      var_stack_ns <- data.frame(
        Scenario               = "Pre-Stress",
        Asset                  = res_pre$var_df$Asset,
        Component_VaR_prop_pct = res_pre$var_df$Component_VaR_prop_pct,
        stringsAsFactors       = FALSE
      )
      if (!is.null(focus_assets)) {
        risk_stack_ns <- risk_stack_ns[risk_stack_ns$Asset %in% focus_assets, ]
        var_stack_ns  <- var_stack_ns[var_stack_ns$Asset  %in% focus_assets, ]
      }
      qual_pal_ns  <- make_qual_pal(sort(unique(risk_stack_ns$Asset)))
      qual_pal_v_ns <- make_qual_pal(sort(unique(var_stack_ns$Asset)))

      p_rd_stack_ns <- ggplot(
        risk_stack_ns, aes(x = Scenario, y = Component_Risk_prop, fill = Asset)) +
        geom_col(position = "stack", alpha = 0.88, width = 0.35) +
        geom_text(aes(label = sprintf("%.1f%%", Component_Risk_prop)),
                  position = position_stack(vjust = 0.5), size = 2.5,
                  color = "white", fontface = "bold") +
        scale_fill_manual(values = qual_pal_ns) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(title = "Std Dev Composition by Asset",
             subtitle = "Component_Risk_prop | Bar sums to 100%",
             x = "", y = "% of Portfolio Sigma") +
        theme_d + theme(legend.position = "right")

      p_vd_stack_ns <- ggplot(
        var_stack_ns, aes(x = Scenario, y = Component_VaR_prop_pct, fill = Asset)) +
        geom_col(position = "stack", alpha = 0.88, width = 0.35) +
        geom_text(aes(label = sprintf("%.1f%%", Component_VaR_prop_pct)),
                  position = position_stack(vjust = 0.5), size = 2.5,
                  color = "white", fontface = "bold") +
        scale_fill_manual(values = qual_pal_v_ns) +
        scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
        labs(title = sprintf("VaR %.0f%% Composition by Asset", confidence_level * 100),
             subtitle = sprintf("Beta method: %s | Bar sums to 100%%", var_method),
             x = "", y = sprintf("%% of Portfolio VaR %.0f%%", confidence_level * 100)) +
        theme_d + theme(legend.position = "right")

      # For non-stressed, use an info panel instead of delta bars
      make_info_panel <- function(msg) {
        ggplot(data.frame(x = 0.5, y = 0.5, label = msg),
               aes(x = x, y = y, label = label)) +
          geom_text(size = 3.5, color = "gray50", hjust = 0.5, vjust = 0.5) +
          labs(title = "", x = "", y = "") +
          theme_d + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                          panel.grid = element_blank())
      }

      p_note <- make_info_panel(
        "Delta bars not available:\nno stressed scenario provided.\nRe-run with asset_shock or\nexogenous_shock to enable.")

      plots$risk_decomposition <- wrap_plots(
        list(p_rd_stack_ns, p_note, p_vd_stack_ns, p_note), ncol = 2) +
        plot_annotation(
          title    = sprintf(
            "Risk Decomposition: Std Dev and VaR %.0f%% Composition (Pre-Stress Only)",
            confidence_level * 100),
          subtitle = sprintf(
            "%d assets | %s simulations | Beta method: %s | Return measure: %s",
            n_assets, format(n_sim, big.mark = ","), var_method, central_label),
          theme    = theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40")
          )
        )
    }

    cat(sprintf("\n  Generated %d plots.\n", length(plots)))
  }

  cat("\n========================================================================\n")
  cat("DECOMPOSITION COMPLETE\n")
  cat("========================================================================\n")

  invisible(list(
    return_decomp = if(is_stressed) list(pre=res_pre$return_df, post=res_post$return_df)
                    else res_pre$return_df,
    risk_decomp   = if(is_stressed) list(pre=res_pre$risk_df, post=res_post$risk_df)
                    else res_pre$risk_df,
    var_decomp    = if(is_stressed) list(pre=res_pre$var_df, post=res_post$var_df)
                    else res_pre$var_df,
    summary       = summary_all,
    port_stats    = if(is_stressed) list(pre=res_pre$port_summary, post=res_post$port_summary)
                    else res_pre$port_summary,
    plots         = plots,
    weights       = w,
    confidence_level = confidence_level,
    var_method    = var_method,
    is_stressed   = is_stressed,
    focus_assets  = focus_assets,
    asset_names   = asset_names,
    n_sim         = n_sim,
    central_label = central_label,
    sim_dist      = sim_dist
  ))
}
