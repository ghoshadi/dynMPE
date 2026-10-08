#' Dynamic marginal policy effect of a threshold rule
#'
#' Estimates the dynamic marginal policy effect of moving the cutoff of a
#' threshold rule that is applied to the same units period after period, by
#' the twice-discounted local linear regression of Ghosh and Wager (2025). The
#' discounted outcome and treatment sums \eqn{\Gamma^Y_{i,t}} and
#' \eqn{\Gamma^A_{i,t}} are each regressed on the running variable near the
#' cutoff over all unit-periods, with weights \eqn{\gamma^t K(|Z_{i,t} - c|/h)}
#' and period fixed effects; the estimate is the ratio of the two jumps.
#'
#' @param data A data frame in long format, one row per unit and period.
#' @param outcome Name of the outcome column \eqn{Y_{i,t}}.
#' @param running_var Name of the running-variable column \eqn{Z_{i,t}}. Use
#'   \code{-Inf} in periods where the threshold rule does not apply.
#' @param time_index Name of the period column, running 0, 1, 2, ... without
#'   gaps within each unit.
#' @param unit_index Name of the unit identifier column; standard errors are
#'   clustered by unit.
#' @param threshold The cutoff \eqn{c}.
#' @param treatment Name of the treatment column \eqn{A_{i,t}}, for a fuzzy
#'   design. If \code{NULL}, the design is sharp: \eqn{A_{i,t} = 1\{Z_{i,t} \ge c\}}.
#' @param gamma Discount factor in \eqn{[0, 1]}.
#' @param alpha The confidence interval has level \eqn{1 - \alpha}.
#' @param h Bandwidth. If \code{NULL}, it is chosen by the dynamic
#'   Imbens-Kalyanaraman procedure of Ghosh and Wager (2025).
#' @param kernel \code{"triangular"}, \code{"uniform"} or \code{"epanechnikov"}.
#' @param vcov Cluster-robust variance estimator: \code{"CR3"}, \code{"CR2"},
#'   \code{"CR1"} or \code{"CR0"}.
#' @param time_fe Whether to include period fixed effects.
#'
#' @return An object of class \code{"dynMPE"}: a list with the \code{estimate},
#'   its standard error \code{se}, the confidence interval \code{ci}, the jumps
#'   in the discounted outcome (\code{numerator}) and treatment
#'   (\code{denominator}) with their standard errors, the \code{bandwidth} and
#'   how it was chosen (\code{bandwidth_rule}), the \code{critical_value}, and
#'   the settings used. The interval is \code{estimate +/- critical_value * se},
#'   where the critical value is the \eqn{1 - \alpha} quantile of the folded
#'   normal \eqn{|N(1/2, 1)|} (Armstrong and Kolesar, 2020).
#'
#' @references Ghosh, A. and Wager, S. (2025). Non-parametric causal inference
#'   in dynamic thresholding designs. arXiv:2512.15244.
#'
#'   Armstrong, T. B. and Kolesar, M. (2020). Simple and honest confidence
#'   intervals in nonparametric regression. Quantitative Economics, 11, 1-39.
#'
#' @examples
#' set.seed(1)
#' n <- 500; periods <- 8
#' sim <- do.call(rbind, lapply(seq_len(n), function(i) {
#'   z <- numeric(periods); y <- numeric(periods); x <- 100 + 9 * rnorm(1)
#'   for (t in seq_len(periods)) {
#'     z[t] <- x
#'     x <- 100 + 0.9 * (x - 100) - (x >= 110) * 0.1 * (x - 100) + 4 * rnorm(1)
#'     y[t] <- -x
#'   }
#'   data.frame(id = i, period = 0:(periods - 1), z = z, y = y)
#' }))
#' fit <- dynMPE(sim, "y", "z", "period", "id", threshold = 110, gamma = 0.8)
#' fit
#' summary(fit)
#' @export
dynMPE <- function(data, outcome, running_var, time_index, unit_index, threshold,
                   treatment = NULL, gamma = 1, alpha = 0.05, h = NULL,
                   kernel = c("triangular", "uniform", "epanechnikov"),
                   vcov = c("CR3", "CR2", "CR1", "CR0"), time_fe = TRUE) {
  kernel <- match.arg(kernel)
  vcov <- match.arg(vcov)
  stopifnot(is.data.frame(data), length(threshold) == 1, is.finite(threshold),
            length(alpha) == 1, alpha > 0, alpha < 1,
            length(gamma) == 1, gamma >= 0, gamma <= 1)
  panel <- make_panel(data, outcome, running_var, time_index, unit_index, threshold,
                      treatment, gamma)
  rule <- if (is.null(h)) "dynamic IK" else "supplied"
  if (is.null(h)) h <- dynamic_bandwidth(panel, threshold, gamma, kernel, vcov, time_fe)
  stopifnot(length(h) == 1, is.finite(h), h > 0)
  fit <- ratio_fit(panel, threshold, gamma, h, kernel, vcov, time_fe)
  cv <- sqrt(stats::qchisq(1 - alpha, df = 1, ncp = 1 / 4))
  structure(list(
    estimate = fit$estimate, se = fit$se,
    ci = c(lower = fit$estimate - cv * fit$se, upper = fit$estimate + cv * fit$se),
    numerator = fit$numerator, denominator = fit$denominator,
    bandwidth = h, bandwidth_rule = rule, critical_value = cv, alpha = alpha,
    gamma = gamma, threshold = threshold, kernel = kernel, vcov = vcov,
    time_fe = time_fe, n_units = panel$n_units, n_obs = fit$n_obs, call = match.call()),
    class = "dynMPE")
}
