# Baseline 3: the recursive estimator of Cellini, Ferreira and Rothstein (2010, QJE 125(1)),
# with their equation numbers. For each lag m, delta_y_m and delta_a_m are the jumps at the
# threshold in Y_{i,t0+m} and A_{i,t0+m} given Z_{i,t0} (their eq. 6; delta_a_0 = 1). Their
# eqs. (8)-(11) invert delta_y_m = sum_{h=0}^{m} delta_a_h theta_{m-h} to the effect theta_m of
# A_{i,t0} alone; we add the discounted sum sum_m gamma^m theta_m.
#
# spec = "global_poly" is their specification: ordinary least squares on all focal units with
# a polynomial of order poly_order (3 in their main tables) common to both sides. spec =
# "local_linear" is the local analogue, with its own IK bandwidth for each lag. Each unit has one
# focal period t0 = focal_t, so their pooled eq. (7) reduces to eq. (6) fitted lag by lag.
# Standard errors use the joint cluster-robust covariance of all 2 Tn - 1 jumps.

cellini_lag_fits <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, c0, Tn,
                             spec = c("global_poly", "local_linear"), poly_order = 3, h = NULL,
                             K = "uniform", vcov_type = c("CR3", "CR2", "CR1", "CR0"), ik_kernel = NULL,
                             balanced = TRUE, focal_t = 0) {
  spec <- match.arg(spec)
  vcov_type <- match.arg(vcov_type)
  kernel <- kernel_function(K)
  if (is.null(ik_kernel)) ik_kernel <- K
  stopifnot(length(Tn) == 1, Tn >= 1, Tn == round(Tn), poly_order >= 1, poly_order == round(poly_order))
  panel <- as_rd_panel(df, Z, Y, A, c0 = c0)
  d0 <- panel[panel$t == focal_t & is.finite(panel$Z_it), , drop = FALSE]
  if (!nrow(d0)) stop("No unit has a finite Z_it at the focal period t = ", focal_t, ".")
  lookup <- function(column, m) panel[[column]][match(paste(d0$i, focal_t + m), paste(panel$i, panel$t))]
  y_lags <- lapply(seq_len(Tn) - 1, function(m) lookup("Y_it", m))
  a_lags <- lapply(seq_len(Tn) - 1, function(m) lookup("A_it", m))
  if (balanced) {
    complete <- Reduce(`&`, lapply(c(y_lags, a_lags), function(v) !is.na(v)))
    if (!any(complete)) stop("No unit is observed through lag ", Tn - 1, "; lower Tn or set balanced = FALSE.")
    if (!all(complete))
      message("cellini: ", sum(!complete), " of ", nrow(d0), " focal units dropped for follow-up shorter than ", Tn, " periods.")
    d0 <- d0[complete, , drop = FALSE]
    y_lags <- lapply(y_lags, `[`, complete)
    a_lags <- lapply(a_lags, `[`, complete)
  }
  v <- d0$Z_it - c0
  if (!any(v >= 0) || !any(v < 0)) stop("The focal units lie on only one side of the threshold (cellini).")
  clusters <- sort(unique(d0$i))

  fit_one <- function(y, hm) {
    use <- !is.na(y)
    if (spec == "global_poly") {
      X <- cbind(1, v >= 0, outer(v, seq_len(poly_order), `^`))
      w <- rep(1, length(y))
      hm <- NA_real_
    } else {
      if (is.null(hm)) hm <- bandwidth_or_silverman(y[use], d0$Z_it[use], c0, ik_kernel, "cellini")
      X <- rd_design(v)
      w <- kernel(v / hm)
      use <- use & w > 0
    }
    if (!any(v[use] >= 0) || !any(v[use] < 0) || sum(use) <= ncol(X)) stop("A lag regression has too few usable observations.")
    fit <- cluster_fit(X[use, , drop = FALSE], y[use], w[use], d0$i[use], vcov_type)
    influence <- numeric(length(clusters))
    influence[match(rownames(fit$influence[[1]]), as.character(clusters))] <- fit$influence[[1]][, 2]
    list(beta = fit$coef[2, 1], influence = influence, h = hm, n = sum(use))
  }

  h_in <- if (is.null(h)) vector("list", 2 * Tn - 1) else as.list(rep(h, length.out = 2 * Tn - 1))
  y_fits <- lapply(seq_len(Tn), function(k) fit_one(y_lags[[k]], h_in[[k]]))
  a_fits <- lapply(seq_len(Tn - 1), function(m) fit_one(a_lags[[m + 1]], h_in[[Tn + m]]))
  U <- sapply(c(y_fits, a_fits), `[[`, "influence")
  list(delta_y = vapply(y_fits, `[[`, 1, "beta"), delta_a = vapply(a_fits, `[[`, 1, "beta"),
       Sigma = crossprod(matrix(U, nrow = length(clusters))), Tn = Tn, spec = spec,
       poly_order = if (spec == "global_poly") poly_order else NA, c0 = c0, focal_t = focal_t,
       se_delta_y = sqrt(colSums(matrix(U, nrow = length(clusters))[, seq_len(Tn), drop = FALSE]^2)),
       h_lag = vapply(y_fits, `[[`, 1, "h"), n_lag = vapply(y_fits, `[[`, 1, "n"),
       n_clusters = length(clusters), n_focal = nrow(d0))
}

# Their eqs. (8)-(11), with the gradient of each theta_m with respect to
# (delta_y_0, ..., delta_y_{Tn-1}, delta_a_1, ..., delta_a_{Tn-1}).
cellini_recursion <- function(delta_y, delta_a) {
  Tn <- length(delta_y)
  unit <- function(k) replace(numeric(2 * Tn - 1), k, 1)
  theta <- numeric(Tn)
  grad <- vector("list", Tn)
  for (m in seq_len(Tn) - 1) {
    theta[m + 1] <- delta_y[m + 1]
    grad[[m + 1]] <- unit(m + 1)
    for (s in seq_len(m)) {
      theta[m + 1] <- theta[m + 1] - delta_a[s] * theta[m - s + 1]
      grad[[m + 1]] <- grad[[m + 1]] - theta[m - s + 1] * unit(Tn + s) - delta_a[s] * grad[[m - s + 1]]
    }
  }
  list(theta_Cel = theta, grad = grad)
}

cellini_lag_profile <- function(fits, alpha = 0.05) {
  rec <- cellini_recursion(fits$delta_y, fits$delta_a)
  se <- vapply(rec$grad, function(g) sqrt(max(drop(crossprod(g, fits$Sigma %*% g)), 0)), 1)
  z <- cv_normal(alpha)
  data.frame(m = seq_len(fits$Tn) - 1, spec = fits$spec, theta_Cel = rec$theta_Cel, se_Cel = se,
             Cel_ci_l = rec$theta_Cel - z * se, Cel_ci_u = rec$theta_Cel + z * se,
             delta_a = c(1, fits$delta_a), h_lag = fits$h_lag, n_lag = fits$n_lag)
}

# One row per gamma; the per-lag estimates are attached as attr(<result>, "cellini").
cellini_rd <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, gamma, c0, Tn,
                       spec = c("global_poly", "local_linear"), poly_order = 3, h = NULL, K = "uniform",
                       vcov_type = c("CR3", "CR2", "CR1", "CR0"), alpha = 0.05, bias_sd_ratio = 0.5,
                       ik_kernel = NULL, balanced = TRUE, focal_t = 0, fits = NULL) {
  if (is.null(fits))
    fits <- cellini_lag_fits(df, Z, Y, A, c0 = c0, Tn = Tn, spec = match.arg(spec), poly_order = poly_order,
                             h = h, K = K, vcov_type = match.arg(vcov_type), ik_kernel = ik_kernel,
                             balanced = balanced, focal_t = focal_t)
  rec <- cellini_recursion(fits$delta_y, fits$delta_a)
  global <- fits$spec == "global_poly"
  out <- do.call(rbind, lapply(gamma, function(g) {
    w <- g^(seq_len(fits$Tn) - 1)
    grad <- Reduce(`+`, Map(`*`, w, rec$grad))
    cbind(data.frame(method = "Cellini2010", gamma = g, c0 = fits$c0,
                     h = if (global) NA_real_ else stats::median(fits$h_lag),
                     bw_rule = if (global) paste0("global_poly_g", fits$poly_order) else "per_lag_ik"),
          confidence_intervals(sum(w * rec$theta_Cel), sqrt(max(drop(crossprod(grad, fits$Sigma %*% grad)), 0)),
                               alpha, bias_sd_ratio),
          data.frame(jump_Y = fits$delta_y[1], jump_A = NA_real_, se_jump_Y = fits$se_delta_y[1],
                     se_jump_A = NA_real_, t_denom = NA_real_, n_obs = fits$n_focal, n_clusters = fits$n_clusters))
  }))
  attr(out, "cellini") <- cellini_lag_profile(fits, alpha)
  out
}
