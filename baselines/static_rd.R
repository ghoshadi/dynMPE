# Baseline 1: the standard local linear regression of Y_it on Z_it over all unit-periods, with
# period fixed effects and kernel weights K(|Z_it - c| / h), but neither the discounted sums nor
# the gamma^t weights of the proposed estimator. It does not depend on gamma and targets the
# one-period jump at the threshold. With fuzzy = TRUE it reports the ratio of the jumps in Y and A.
# The default bandwidth is the sharp Imbens-Kalyanaraman rule on the pairs (Y_it, Z_it).

static_rd <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, c0, h = NULL, K = "uniform",
                      time_fe = TRUE, fuzzy = FALSE, vcov_type = c("CR3", "CR2", "CR1", "CR0"),
                      alpha = 0.05, bias_sd_ratio = 0.5, ik_kernel = NULL) {
  vcov_type <- match.arg(vcov_type)
  kernel <- kernel_function(K)
  panel <- as_rd_panel(df, Z, Y, A, c0 = c0)
  finite <- is.finite(panel$Z_it)
  if (!any(finite)) stop("No finite Z_it values.")
  bw_rule <- "supplied"
  if (is.null(h)) {
    h <- bandwidth_or_silverman(panel$Y_it[finite], panel$Z_it[finite], c0,
                                if (is.null(ik_kernel)) K else ik_kernel, "static")
    bw_rule <- "sharp_ik"
  }
  stopifnot(length(h) == 1, is.finite(h), h > 0)
  zc <- panel$Z_it - c0
  w <- ifelse(finite, kernel(zc / h), 0)
  keep <- w > 0
  if (!any(zc[keep] >= 0) || !any(zc[keep] < 0))
    stop("The bandwidth h = ", format(h, digits = 4), " has observations on only one side of the threshold (static).")
  Ymat <- if (fuzzy) cbind(panel$Y_it, panel$A_it)[keep, ] else panel$Y_it[keep]
  fit <- cluster_fit(rd_design(zc[keep], panel$t[keep], time_fe), Ymat, w[keep], panel$i[keep], vcov_type)
  r <- if (fuzzy) jump_ratio(fit) else {
    se <- sqrt(sum(fit$influence[[1]][, 2]^2))
    list(est = fit$coef[2, 1], se = se, jump_Y = fit$coef[2, 1], jump_A = NA_real_, se_jump_Y = se,
         se_jump_A = NA_real_, t_denom = NA_real_)
  }
  cbind(data.frame(method = "static", gamma = NA_real_, c0 = c0, h = h, bw_rule = bw_rule),
        confidence_intervals(r$est, r$se, alpha, bias_sd_ratio),
        data.frame(jump_Y = r$jump_Y, jump_A = r$jump_A, se_jump_Y = r$se_jump_Y, se_jump_A = r$se_jump_A,
                   t_denom = r$t_denom, n_obs = sum(keep), n_clusters = length(unique(panel$i[keep]))))
}
