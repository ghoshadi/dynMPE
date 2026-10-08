# Baseline 2: the naive long-run LLR. One cross-sectional local linear regression at t = 0 of
# each discounted sum Gamma^R_{i,0}, R in {Y, A}, on Z_{i,0}; the estimate is the ratio of the
# two jumps. Units with Z_{i,0} = -Inf drop out. The bandwidth follows the fuzzy-RD version of
# the Imbens-Kalyanaraman rule: a sharp pilot on Gamma^Y_{i,0} gives a pilot ratio, and the sharp
# rule is then applied to Gamma^Y_{i,0} - ratio * Gamma^A_{i,0}.

naive_first_period <- function(panel) {
  d0 <- panel[panel$t == 0 & is.finite(panel$Z_it), , drop = FALSE]
  if (!nrow(d0)) stop("No unit has a finite Z_it at t = 0.")
  d0
}

naive_fit <- function(d0, c0, h, K, vcov_type) {
  zc <- d0$Z_it - c0
  w <- kernel_function(K)(zc / h)
  keep <- w > 0
  if (!any(zc[keep] >= 0) || !any(zc[keep] < 0))
    stop("The bandwidth h = ", format(h, digits = 4), " has observations on only one side of the threshold (naive).")
  fit <- cluster_fit(rd_design(zc[keep]), cbind(d0$Gamma_Y, d0$Gamma_A)[keep, ], w[keep], d0$i[keep], vcov_type)
  c(jump_ratio(fit), n_obs = sum(keep), n_clusters = length(unique(d0$i[keep])))
}

naive_bandwidth <- function(d0, c0, K, ik_kernel, vcov_type) {
  h_pilot <- bandwidth_or_silverman(d0$Gamma_Y, d0$Z_it, c0, ik_kernel, "naive pilot")
  ratio <- naive_fit(d0, c0, h_pilot, K, vcov_type)$est
  # A pilot ratio that dominates the residual makes the bandwidth track its noise; keep the pilot then.
  if (abs(ratio) * stats::sd(d0$Gamma_A) > stats::sd(d0$Gamma_Y)) {
    message(sprintf("The pilot ratio (%.3g) dominates the residual; the pilot bandwidth %.4g is used.", ratio, h_pilot))
    return(list(bandwidth = h_pilot, rule = "fuzzy_ik_pilot", h_pilot = h_pilot, ratio_pilot = ratio))
  }
  h <- bandwidth_or_silverman(d0$Gamma_Y - ratio * d0$Gamma_A, d0$Z_it, c0, ik_kernel, "naive")
  list(bandwidth = h, rule = "fuzzy_ik", h_pilot = h_pilot, ratio_pilot = ratio)
}

naive_ik_bandwidth <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, gamma, c0, K = "uniform",
                               vcov_type = c("CR3", "CR2", "CR1", "CR0"), ik_kernel = NULL) {
  d0 <- naive_first_period(discounted_gamma(as_rd_panel(df, Z, Y, A, c0 = c0), gamma))
  naive_bandwidth(d0, c0, K, if (is.null(ik_kernel)) K else ik_kernel, match.arg(vcov_type))
}

naive_rd <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, gamma, c0, h = NULL, K = "uniform",
                     vcov_type = c("CR3", "CR2", "CR1", "CR0"), alpha = 0.05, bias_sd_ratio = 0.5,
                     ik_kernel = NULL) {
  vcov_type <- match.arg(vcov_type)
  d0 <- naive_first_period(discounted_gamma(as_rd_panel(df, Z, Y, A, c0 = c0), gamma))
  bw_rule <- "supplied"
  if (is.null(h)) {
    bw <- naive_bandwidth(d0, c0, K, if (is.null(ik_kernel)) K else ik_kernel, vcov_type)
    h <- bw$bandwidth
    bw_rule <- bw$rule
  }
  stopifnot(length(h) == 1, is.finite(h), h > 0)
  r <- naive_fit(d0, c0, h, K, vcov_type)
  cbind(data.frame(method = "naive", gamma = gamma, c0 = c0, h = h, bw_rule = bw_rule),
        confidence_intervals(r$est, r$se, alpha, bias_sd_ratio),
        data.frame(jump_Y = r$jump_Y, jump_A = r$jump_A, se_jump_Y = r$se_jump_Y, se_jump_A = r$se_jump_A,
                   t_denom = r$t_denom, n_obs = r$n_obs, n_clusters = r$n_clusters))
}
