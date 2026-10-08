kernels <- list(
  triangular = function(u) pmax(1 - abs(u), 0),
  uniform = function(u) 0.5 * (abs(u) <= 1),
  epanechnikov = function(u) 0.75 * pmax(1 - u^2, 0))

# xi1: bias constant of the one-sided local linear kernel; CK: constant of the IK bandwidth.
kernel_constants <- cbind(
  triangular = c(xi1 = -1 / 10, CK = 480^(1 / 5)),
  uniform = c(xi1 = -1 / 6, CK = 144^(1 / 5)),
  epanechnikov = c(xi1 = -11 / 95, CK = (284160 / 847)^(1 / 5)))

# Sorts by unit and period and computes Gamma_t = R_t + gamma Gamma_{t+1} within each unit.
make_panel <- function(data, outcome, running_var, time_index, unit_index, threshold,
                       treatment, gamma) {
  absent <- setdiff(c(outcome, running_var, time_index, unit_index, treatment), names(data))
  if (length(absent)) stop("Columns not found in data: ", paste(absent, collapse = ", "))
  o <- order(data[[unit_index]], data[[time_index]])
  unit <- data[[unit_index]][o]
  t <- data[[time_index]][o]
  z <- as.numeric(data[[running_var]][o])
  y <- as.numeric(data[[outcome]][o])
  a <- if (is.null(treatment)) as.numeric(z >= threshold) else as.numeric(data[[treatment]][o])
  if (anyNA(z) || any(z == Inf)) stop("`", running_var, "` must be finite, or -Inf where the rule does not apply.")
  if (!all(is.finite(y))) stop("`", outcome, "` must be finite.")
  if (!all(is.finite(a))) stop("`", treatment, "` must be finite.")
  first <- c(TRUE, unit[-1] != unit[-length(unit)])
  if (any(t[first] != 0) || any(diff(t)[!first[-1]] != 1))
    stop("`", time_index, "` must run 0, 1, 2, ... without gaps within each unit.")
  lengths <- diff(c(which(first), length(t) + 1))
  if (min(lengths) < max(lengths))
    warning("Units are observed for different numbers of periods (", min(lengths), " to ",
            max(lengths), "); each discounted sum runs to the unit's own last period.", call. = FALSE)
  has_next <- c(!first[-1], FALSE)
  G <- cbind(y, a)
  for (k in sort(unique(t[has_next]), decreasing = TRUE)) {
    r <- which(t == k & has_next)
    G[r, ] <- G[r, ] + gamma * G[r + 1, , drop = FALSE]
  }
  list(z = z, t = t, cluster = cumsum(first), gy = G[, 1], ga = G[, 2], n_units = sum(first))
}

design <- function(x, t, time_fe) {
  d <- as.numeric(x >= 0)
  X <- cbind(1, d, x, d * x)
  levels <- sort(unique(t))
  if (time_fe && length(levels) > 1) X <- cbind(X, outer(t, levels[-1], "==") * 1)
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
    warning("The local design is singular; a ridge of 1e-8 (on standardized columns) was added.",
            call. = FALSE)
    diag(A) <- diag(A) + 1e-8
  }
  bread <- chol2inv(chol(A))
  coef <- bread %*% crossprod(Xs, as.matrix(Y) * sw)
  list(Xs = Xs, bread = bread, scale = scale, coef = coef / scale,
       resid = as.matrix(Y) * sw - Xs %*% coef)
}

# influence[[q]][g, j] is cluster g's contribution to coefficient j for response q.
cluster_fit <- function(X, Y, w, cluster, vcov) {
  f <- weighted_fit(X, Y, w)
  E <- f$resid
  if (vcov %in% c("CR2", "CR3")) {
    for (rows in split(seq_len(nrow(E)), cluster)) {
      Xg <- f$Xs[rows, , drop = FALSE]
      M <- diag(length(rows)) - Xg %*% f$bread %*% t(Xg)
      if (rcond(M) < 1e-10) {
        warning("A cluster has leverage 1; its ", vcov, " adjustment uses a ridge of 1e-8.",
                call. = FALSE)
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

ratio_fit <- function(panel, threshold, gamma, h, kernel, vcov, time_fe) {
  x <- panel$z - threshold
  w <- ifelse(is.finite(x), kernels[[kernel]](x / h), 0) * gamma^panel$t
  keep <- w > 0
  if (!any(x[keep] >= 0) || !any(x[keep] < 0))
    stop("Bandwidth h = ", format(h, digits = 4), " has observations on only one side of the threshold.")
  fit <- cluster_fit(design(x[keep], panel$t[keep], time_fe), cbind(panel$gy, panel$ga)[keep, ],
                     w[keep], panel$cluster[keep], vcov)
  jump <- fit$coef[2, ]
  num <- fit$influence[[1]][, 2]
  den <- fit$influence[[2]][, 2]
  if (jump[2] == 0) stop("The jump in the discounted treatment is exactly zero.")
  estimate <- jump[1] / jump[2]
  list(estimate = estimate, se = sqrt(sum(((num - estimate * den) / jump[2])^2)),
       numerator = c(estimate = jump[1], se = sqrt(sum(num^2))),
       denominator = c(estimate = jump[2], se = sqrt(sum(den^2))), n_obs = sum(keep))
}

# Sharp Imbens-Kalyanaraman bandwidth for y on x - threshold, as in plrd::IK_bandwidth.
ik_bandwidth <- function(y, x, kernel) {
  n <- length(y)
  right <- x >= 0
  left <- !right
  if (n < 10 || !any(left) || !any(right)) stop("too few observations on one side of the threshold")
  h1 <- 1.84 * stats::sd(x) * n^(-1 / 5)
  near <- list(left & x >= -h1, right & x <= h1)
  n1 <- vapply(near, sum, 1)
  if (any(n1 <= 1)) stop("too few observations near the threshold")
  sigma2 <- vapply(near, function(i) stats::var(y[i]), 1)
  f <- sum(n1) / (2 * n * h1)
  m3 <- 6 * stats::lm.fit(cbind(1, right, x, x^2, x^3), y)$coefficients[[5]]
  if (is.na(m3)) stop("the cubic pilot regression is rank-deficient")
  h2 <- 7200^(1 / 7) * (sigma2 / (f * m3^2))^(1 / 7) * c(sum(left), sum(right))^(-1 / 7)
  near <- list(left & x >= -h2[1], right & x <= h2[2])
  n2 <- vapply(near, sum, 1)
  if (any(n2 <= 2)) stop("too few observations for the curvature pilot")
  m2 <- vapply(near, function(i) 2 * stats::lm.fit(cbind(1, x[i], x[i]^2), y[i])$coefficients[[3]], 1)
  if (anyNA(m2)) stop("the quadratic pilot regression is rank-deficient")
  r <- 2160 * sigma2 / (n2 * h2^4)
  h <- kernel_constants[["CK", kernel]] * (sum(sigma2) / (f * ((m2[2] - m2[1])^2 + sum(r))) / n)^(1 / 5)
  if (!is.finite(h) || h <= 0) stop("the bandwidth is not positive and finite")
  h
}

# Procedure (dynamic-ik-bandwidth) of Ghosh and Wager (2025), with n the number of units.
dynamic_bandwidth <- function(panel, threshold, gamma, kernel, vcov, time_fe) {
  x <- panel$z - threshold
  finite <- is.finite(x)
  n <- panel$n_units
  h_pilot <- tryCatch(ik_bandwidth(panel$gy[finite], x[finite], kernel), error = function(e) {
    warning("The sharp IK pilot bandwidth failed (", conditionMessage(e),
            "); Silverman's rule was used instead.", call. = FALSE)
    1.84 * stats::sd(x[finite]) * sum(finite)^(-1 / 5)
  })
  pilot <- ratio_fit(panel, threshold, gamma, h_pilot, kernel, vcov, time_fe)
  V <- n * h_pilot * (pilot$se * pilot$denominator[["estimate"]])^2
  resid <- panel$gy - pilot$estimate * panel$ga

  wt <- gamma^panel$t
  pool <- finite & wt > 0
  if (!any(pool & x >= 0) || !any(pool & x < 0))
    stop("The dynamic bandwidth needs observations on both sides of the threshold.")
  h1 <- 1.84 * stats::sd(x[finite]) * n^(-1 / 5)
  f <- sum(wt[finite & abs(x) <= h1]) / (2 * h1 * sum(wt[finite]))
  if (!is.finite(f) || f <= 0) stop("The discounted density at the threshold is not positive.")
  mid <- pool & x >= stats::median(x[pool & x < 0]) & x <= stats::median(x[pool & x >= 0])
  if (sum(mid) <= 5) stop("Too few observations for the cubic curvature pilot.")
  cubic <- weighted_fit(cbind(1, x[mid] >= 0, x[mid], x[mid]^2 / 2, x[mid]^3 / 6), resid[mid], wt[mid])
  sigma2 <- sum(cubic$resid^2) / sum(wt[mid])
  share <- c(sum(wt[pool & x >= 0]), sum(wt[pool & x < 0])) / sum(wt[pool])
  h2 <- 7200^(1 / 7) * (sigma2 / (f * cubic$coef[5]^2))^(1 / 7) * (n * share)^(-1 / 7)
  if (!all(is.finite(h2) & h2 > 0)) stop("The curvature pilot bandwidths are not positive and finite.")

  k <- kernels[[kernel]]
  right <- finite & x >= 0 & x <= h2[1]
  left <- finite & x < 0 & x >= -h2[2]
  kw <- numeric(length(x))
  kw[right] <- k(x[right] / h2[1])
  kw[left] <- k(x[left] / h2[2])
  keep <- kw * wt > 0
  if (!any(keep & x >= 0) || !any(keep & x < 0))
    stop("The curvature pilot windows miss one side of the threshold.")
  X <- design(x[keep], panel$t[keep], time_fe)
  X <- cbind(X, x[keep]^2 / 2, X[, 2] * x[keep]^2 / 2)
  curvature <- cluster_fit(X, resid[keep], (kw * wt)[keep], panel$cluster[keep], vcov)
  B <- curvature$coef[ncol(X), 1]
  R <- sum(curvature$influence[[1]][, ncol(X)]^2)
  h <- (V / (kernel_constants[["xi1", kernel]]^2 * (B^2 + 3 * R)))^(1 / 5) * n^(-1 / 5)
  if (!is.finite(h) || h <= 0) stop("The dynamic bandwidth is not positive and finite.")
  h
}

#' Print a dynMPE object
#' @param x dynMPE object
#' @param digits number of digits to print
#' @param ... Additional arguments passed to print methods.
#' @export
print.dynMPE <- function(x, digits = 4, ...) {
  cat("Dynamic marginal policy effect\n")
  cat(sprintf("Threshold: %g   gamma: %g   period fixed effects: %s\n",
              x$threshold, x$gamma, if (x$time_fe) "yes" else "no"))
  cat(sprintf("Bandwidth: %s (%s)   Kernel: %s   Variance: %s, clustered by unit\n",
              format(x$bandwidth, digits = digits), x$bandwidth_rule, x$kernel, x$vcov))
  cat(sprintf("Units: %d   Unit-periods within the bandwidth: %d\n\n", x$n_units, x$n_obs))
  print(summary(x), digits = digits, ...)
  cat(sprintf("\n%g%% CI: [%s, %s]   (bias-aware critical value %.6f)\n", 100 * (1 - x$alpha),
              format(x$ci[[1]], digits = digits), format(x$ci[[2]], digits = digits), x$critical_value))
  invisible(x)
}

#' dynMPE summary
#' @param object dynMPE object
#' @param ... Additional arguments (currently ignored).
#' @export
summary.dynMPE <- function(object, ...) {
  rows <- rbind(effect = c(object$estimate, object$se),
                numerator = object$numerator, denominator = object$denominator)
  data.frame(Estimate = rows[, 1], `Std. Error` = rows[, 2], `t value` = rows[, 1] / rows[, 2],
             check.names = FALSE)
}
