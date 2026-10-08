if (!requireNamespace("dynMPE", quietly = TRUE)) remotes::install_github("ghoshadi/dynMPE")

# Shared internals of the baseline estimators in static_rd.R, naive_rd.R and cellini_rd.R:
# kernels, the panel format, discounted sums, weighted least squares with unit-clustered
# standard errors (the same CR0-CR3 as the dynMPE package), the sharp Imbens-Kalyanaraman
# bandwidth (dynMPE::IK_bandwidth) and the critical values. Source this file before the estimator
# files. Kernels are given by name: "triangular", "uniform" or "epanechnikov".
#
# Data: a long data.frame with columns i (unit), t (period, 0, 1, 2, ... without gaps),
# Z_it (running variable, -Inf where the rule does not apply), Y_it and optionally A_it
# (built as 1{Z_it >= c0} if absent); or n x T matrices Z, Y and optionally A.

kernel_function <- function(kernel = c("triangular", "uniform", "epanechnikov"))
  switch(match.arg(kernel), triangular = function(u) pmax(1 - abs(u), 0),
         uniform = function(u) 0.5 * (abs(u) <= 1), epanechnikov = function(u) 0.75 * pmax(1 - u^2, 0))

# The bias-aware value solves P(|N(b, 1)| <= cv) = 1 - alpha; b = 1/2 at an MSE-optimal h.
cv_normal <- function(alpha = 0.05) stats::qnorm(1 - alpha / 2)
cv_bias_aware <- function(alpha = 0.05, bias_sd_ratio = 0.5)
  sqrt(stats::qchisq(1 - alpha, df = 1, ncp = bias_sd_ratio^2))

confidence_intervals <- function(est, se, alpha = 0.05, bias_sd_ratio = 0.5) {
  z <- cv_normal(alpha)
  zb <- cv_bias_aware(alpha, bias_sd_ratio)
  data.frame(est = est, se = se, ci_l = est - z * se, ci_u = est + z * se,
             ci_l_ba = est - zb * se, ci_u_ba = est + zb * se, cv = z, cv_ba = zb)
}

RD_SCHEMA <- c("method", "gamma", "c0", "h", "bw_rule", "est", "se", "ci_l", "ci_u", "ci_l_ba",
               "ci_u_ba", "cv", "cv_ba", "jump_Y", "jump_A", "se_jump_Y", "se_jump_A", "t_denom",
               "n_obs", "n_clusters")

as_rd_panel <- function(df = NULL, Z = NULL, Y = NULL, A = NULL, c0 = NULL) {
  if (!is.null(df)) {
    if (!is.null(Z) || !is.null(Y) || !is.null(A)) stop("Supply either `df` or the matrices, not both.")
    need <- c("i", "t", "Z_it", "Y_it")
    if (!all(need %in% names(df))) stop("`df` must contain columns: ", paste(need, collapse = ", "), ".")
    out <- data.frame(i = df$i, t = as.integer(df$t), Z_it = as.numeric(df$Z_it), Y_it = as.numeric(df$Y_it))
    if ("A_it" %in% names(df)) out$A_it <- as.numeric(df$A_it)
  } else {
    if (is.null(Z) || is.null(Y)) stop("Supply both `Z` and `Y` matrices.")
    Z <- as.matrix(Z); Y <- as.matrix(Y)
    if (!identical(dim(Z), dim(Y))) stop("`Z` and `Y` must have the same dimensions.")
    out <- data.frame(i = rep(seq_len(nrow(Z)), times = ncol(Z)), t = rep(seq_len(ncol(Z)) - 1, each = nrow(Z)),
                      Z_it = as.numeric(Z), Y_it = as.numeric(Y))
    if (!is.null(A)) {
      A <- as.matrix(A)
      if (!identical(dim(A), dim(Z))) stop("`A` must have the same dimensions as `Z`.")
      out$A_it <- as.numeric(A)
    }
  }
  out$Z_it[is.na(out$Z_it)] <- -Inf
  if (any(!is.finite(out$Y_it))) stop("Y_it must be finite.")
  if (is.null(out$A_it)) {
    if (is.null(c0)) stop("No A_it supplied and no `c0` to build 1{Z_it >= c0}.")
    out$A_it <- as.numeric(is.finite(out$Z_it) & out$Z_it >= c0)
  }
  if (any(!is.finite(out$A_it))) stop("A_it must be finite.")
  out <- out[order(out$i, out$t), , drop = FALSE]
  rownames(out) <- NULL
  bad <- vapply(split(out$t, out$i), function(t) t[1] != 0 || any(diff(t) != 1), TRUE)
  if (any(bad)) stop("t must run 0, 1, 2, ... without gaps within each unit; offending units: ",
                     paste(utils::head(names(bad)[bad], 5), collapse = ", "), ".")
  out
}

# Adds Gamma_Y and Gamma_A: Gamma_t = R_t + gamma Gamma_{t+1} within each unit.
discounted_gamma <- function(panel, gamma) {
  panel <- panel[order(panel$i, panel$t), , drop = FALSE]
  first <- c(TRUE, panel$i[-1] != panel$i[-nrow(panel)])
  lengths <- diff(c(which(first), nrow(panel) + 1))
  if (min(lengths) < max(lengths))
    warning("Units are observed for different numbers of periods (", min(lengths), " to ", max(lengths),
            "); each discounted sum runs to the unit's own last period.", call. = FALSE)
  has_next <- c(!first[-1], FALSE)
  G <- cbind(panel$Y_it, panel$A_it)
  for (k in sort(unique(panel$t[has_next]), decreasing = TRUE)) {
    r <- which(panel$t == k & has_next)
    G[r, ] <- G[r, ] + gamma * G[r + 1, , drop = FALSE]
  }
  panel$Gamma_Y <- G[, 1]
  panel$Gamma_A <- G[, 2]
  panel
}

# Local linear design [1, D, Z - c, D (Z - c), period effects]; the jump is the coefficient on D.
rd_design <- function(zc, period = NULL, time_fe = FALSE) {
  d <- as.numeric(zc >= 0)
  X <- cbind(1, d, zc, d * zc)
  levels <- sort(unique(period))
  if (time_fe && length(levels) > 1) X <- cbind(X, outer(period, levels[-1], "==") * 1)
  X
}

# Weighted least squares on standardized columns; a tiny ridge only if the design is singular.
weighted_fit <- function(X, Y, w) {
  sw <- sqrt(w)
  Xw <- X * sw
  scale <- sqrt(colSums(Xw^2))
  if (any(scale == 0)) stop("A regressor is identically zero within the bandwidth.")
  Xs <- t(t(Xw) / scale)
  A <- crossprod(Xs)
  if (rcond(A) < 1e-12) {
    warning("The design is singular; a ridge of 1e-8 (on standardized columns) was added.", call. = FALSE)
    diag(A) <- diag(A) + 1e-8
  }
  bread <- chol2inv(chol(A))
  coef <- bread %*% crossprod(Xs, as.matrix(Y) * sw)
  list(Xs = Xs, bread = bread, scale = scale, coef = coef / scale, resid = as.matrix(Y) * sw - Xs %*% coef)
}

# influence[[q]][g, j] is cluster g's contribution to coefficient j for response q, so variances
# and covariances, also across responses and regressions, are sums of products.
cluster_fit <- function(X, Y, w, cluster, vcov = c("CR3", "CR2", "CR1", "CR0")) {
  vcov <- match.arg(vcov)
  f <- weighted_fit(X, Y, w)
  E <- f$resid
  if (vcov %in% c("CR2", "CR3")) {
    if (!anyDuplicated(cluster)) {
      leverage <- rowSums((f$Xs %*% f$bread) * f$Xs)
      E <- E / if (vcov == "CR3") 1 - leverage else sqrt(1 - leverage)
    } else for (rows in split(seq_len(nrow(E)), cluster)) {
      Xg <- f$Xs[rows, , drop = FALSE]
      M <- diag(length(rows)) - Xg %*% f$bread %*% t(Xg)
      if (rcond(M) < 1e-10) {
        warning("A cluster has leverage 1; its ", vcov, " adjustment uses a ridge of 1e-8.", call. = FALSE)
        diag(M) <- diag(M) + 1e-8
      }
      E[rows, ] <- if (vcov == "CR3") solve(M, E[rows, , drop = FALSE]) else {
        e <- eigen(M, symmetric = TRUE)
        e$vectors %*% (e$values^(-1 / 2) * crossprod(e$vectors, E[rows, , drop = FALSE]))
      }
    }
  }
  G <- length(unique(cluster))
  correction <- if (vcov == "CR1") sqrt(G / (G - 1) * (nrow(X) - 1) / (nrow(X) - ncol(X))) else 1
  influence <- lapply(seq_len(ncol(E)), function(q)
    correction * t(t(rowsum(f$Xs * E[, q], cluster) %*% f$bread) / f$scale))
  list(coef = f$coef, influence = influence)
}

# Ratio of the jumps (coefficient 2) in responses 1 and 2 of a cluster_fit, with its
# delta-method standard error and the t-statistic of the denominator.
jump_ratio <- function(fit) {
  jump <- fit$coef[2, 1:2]
  num <- fit$influence[[1]][, 2]
  den <- fit$influence[[2]][, 2]
  if (abs(jump[2]) < .Machine$double.eps) stop("The jump in the denominator is numerically zero.")
  est <- jump[1] / jump[2]
  se_den <- sqrt(sum(den^2))
  list(est = est, se = sqrt(sum(((num - est * den) / jump[2])^2)), jump_Y = jump[1], jump_A = jump[2],
       se_jump_Y = sqrt(sum(num^2)), se_jump_A = se_den, t_denom = jump[2] / se_den)
}

silverman_bandwidth <- function(x, c0) {
  x <- (x - c0)[is.finite(x)]
  1.84 * stats::sd(x) * length(x)^(-1 / 5)
}

bandwidth_or_silverman <- function(y, x, c0, kernel, what = "") {
  tryCatch(dynMPE::IK_bandwidth(y, x, c0, kernel)$bandwidth, error = function(e) {
    message("The sharp IK bandwidth", if (nzchar(what)) paste0(" (", what, ")"), " failed (",
            conditionMessage(e), "); Silverman's rule was used instead.")
    silverman_bandwidth(x, c0)
  })
}
