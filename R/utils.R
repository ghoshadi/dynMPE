kernels <- list(
  triangular = function(u) pmax(1 - abs(u), 0),
  uniform = function(u) 0.5 * (abs(u) <= 1),
  epanechnikov = function(u) 0.75 * pmax(1 - u^2, 0))

# Bias constant of the one-sided local linear equivalent kernel.
kernel_xi1 <- c(triangular = -1 / 10, uniform = -1 / 6, epanechnikov = -11 / 95)

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

#' Compute MSE-optimal Imbens-Kalyanaraman bandwidth for a sharp RD.
#'
#' This convenience function computes weights using the Imbens-Kalyanaraman bandwidth procedure.
#' The code does exactly what the MATLAB code available on the author's website does.
#'
#' @param Y The outcomes.
#' @param X The running variable.
#' @param threshold The threshold.
#' @param kernel The kernel type used to construct weights within the bandwidth.
#'
#' @return A list containing the sample weights along with optimal bandwidth.
#'
#' @references Imbens, G., and Kalyanaraman, K. (2012).
#'  Optimal Bandwidth Choice for the Regression Discontinuity Estimator.
#'  The Review of Economic Studies, 79(3).
#'
#' @examples
#' set.seed(42)
#' n = 1000; threshold = 0
#' X = runif(n, -1, 1)
#' W = as.numeric(X >= threshold)
#' Y = (1 + 2*W)*(1 + X^2) + 1 / (1 + exp(X)) + rnorm(n, sd = .5)
#' out = IK_bandwidth(Y, X, threshold)
#'
#' @export
IK_bandwidth <- function(Y, X, threshold,
                         kernel = c("triangular", "uniform", "epanechnikov")) {
  if (length(Y) != length(X)) stop("'Y' and 'X' must have the same length.")
  if (threshold >= max(X) || threshold <= min(X))
    stop("RD threshold is outside the running variable range.")

  kernel <- match.arg(kernel)
  x <- X - threshold; n <- length(x)
  left <- x < 0; right <- !left

  # Density and conditional variances at the threshold
  h1 <- 1.84 * stats::sd(x) * n^(-1/5)
  i.min <- left & x >= -h1
  i.plus <- right & x <= h1
  n1 <- c(sum(i.min), sum(i.plus))
  if (any(n1 <= 1)) stop("Insufficient observations near discontinuity.")

  sigma2 <- c(stats::var(Y[i.min]), stats::var(Y[i.plus]))
  fc <- sum(n1) / (2 * n * h1)

  # Pilot third derivative and bandwidths for second derivatives
  m3 <- 6 * unname(stats::lm.fit(cbind(1, right, x, x^2, x^3), Y)$coefficients[5])
  if (is.na(m3))
    stop("The IK cubic pilot regression is rank-deficient.")

  h2 <- 7200^(1/7) *
    (sigma2 / (fc * m3^2))^(1/7) *
    c(sum(left), sum(right))^(-1/7)

  i.min <- left & x >= -h2[1]
  i.plus <- right & x <= h2[2]
  n2 <- c(sum(i.min), sum(i.plus))
  if (any(n2 <= 2)) stop("Insufficient observations near discontinuity.")

  m2 <- c(
    2 * unname(stats::lm.fit(cbind(1, x[i.min], x[i.min]^2), Y[i.min])$coefficients[3]),
    2 * unname(stats::lm.fit(cbind(1, x[i.plus], x[i.plus]^2), Y[i.plus])$coefficients[3])
  )

  if (is.na(m2[1]))
    stop("The IK quadratic pilot regression is rank-deficient below the threshold.")
  if (is.na(m2[2]))
    stop("The IK quadratic pilot regression is rank-deficient above the threshold.")

  # Regularization and optimal bandwidth
  r <- 2160 * sigma2 / (n2 * h2^4)
  CK <- switch(kernel,
               triangular   = 480^(1/5),
               uniform      = 144^(1/5),
               epanechnikov = (284160 / 847)^(1/5)
  )

  h.opt <- CK *
    (sum(sigma2) / (fc * ((m2[2] - m2[1])^2 + sum(r))))^(1/5) *
    n^(-1/5)

  if (h.opt <= 0)
    stop("The calculated IK bandwidth is not positive.")

  if (!is.finite(h.opt))
    stop("The calculated IK bandwidth is not finite.")

  # Kernel weights, normalized to sum to one
  u <- abs(x / h.opt)
  weights <- switch(kernel,
                    triangular   = pmax(1 - u, 0),
                    uniform      = as.numeric(u <= 1),
                    epanechnikov = pmax(1 - u^2, 0)
  )

  list(
    bandwidth = unname(h.opt),
    weights = weights / sum(weights)
  )
}

# Procedure (dynamic-ik-bandwidth) of Ghosh and Wager (2025), with n the number of units.
dynamic_bandwidth <- function(panel, threshold, gamma, kernel, vcov, time_fe) {
  x <- panel$z - threshold
  finite <- is.finite(x)
  n <- panel$n_units
  h_pilot <- tryCatch(IK_bandwidth(panel$gy[finite], panel$z[finite], threshold, kernel)$bandwidth, error = function(e) {
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
  h <- (V / (kernel_xi1[[kernel]]^2 * (B^2 + 3 * R)))^(1 / 5) * n^(-1 / 5)
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
