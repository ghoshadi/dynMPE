# Checks that every baseline runs on one small replication of the autoregressive design and
# agrees with itself. Run from this folder: Rscript smoke_test.R

for (f in c("shared_helpers.R", "static_rd.R", "naive_rd.R", "cellini_rd.R")) source(f)

pass <- function(msg) cat("  [ok]  ", msg, "\n", sep = "")
chk <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); pass(msg) }
hdr <- function(msg) cat("\n== ", msg, " ==\n", sep = "")

set.seed(20260914)
c0 <- 110; Tn <- 12; n <- 600
X <- matrix(NA_real_, n, Tn + 1)
X[, 1] <- 100 + 4 / sqrt(1 - 0.9^2) * rnorm(n)
for (t in seq_len(Tn)) {
  xc <- X[, t] - 100
  X[, t + 1] <- 100 + t + 0.9 * xc - (X[, t] >= c0) * 0.1 * pmax(xc, 0) + 4 * rnorm(n)
}
mats <- list(Z = X[, 1:Tn], Y = -X[, -1], A = (X[, 1:Tn] >= c0) * 1)
panel <- as_rd_panel(Z = mats$Z, Y = mats$Y, A = mats$A)
gammas <- c(0.5, 0.8, 1)

hdr("1. as_rd_panel")
chk(isTRUE(all.equal(panel, as_rd_panel(df = panel))), "the data.frame round trip is a fixed point")
chk(nrow(panel) == n * Tn && all(range(panel$t) == c(0, Tn - 1)), "n T rows, t runs from 0")
chk(inherits(try(as_rd_panel(df = panel[panel$t != 3, ]), silent = TRUE), "try-error"), "a gap in t is rejected")
chk(isTRUE(all.equal(as_rd_panel(Z = mats$Z, Y = mats$Y, c0 = c0)$A_it, panel$A_it)),
    "A defaults to 1{Z >= c0}")

hdr("2. Discounted sums")
pg <- discounted_gamma(panel, gamma = 0.8)
u1 <- pg[pg$i == 1, ]
chk(abs(u1$Gamma_Y[1] - sum(0.8^(seq_len(Tn) - 1) * u1$Y_it)) < 1e-10, "Gamma^Y_{1,0} is the discounted sum")
chk(abs(tail(u1$Gamma_Y, 1) - tail(u1$Y_it, 1)) < 1e-12, "Gamma at the last period equals Y")

hdr("3. Critical values")
root <- uniroot(function(v) pnorm(v - 0.5) - pnorm(-v - 0.5) - 0.95, c(0, 10), tol = 1e-12)$root
chk(abs(cv_bias_aware(0.05) - root) < 1e-6, sprintf("cv_bias_aware = %.6f", cv_bias_aware()))
chk(abs(cv_bias_aware(0.05, 0) - cv_normal(0.05)) < 1e-8, "bias_sd_ratio = 0 gives the normal quantile")

hdr("4. The estimators run, under every variance estimator")
res <- c(list(static_rd(df = panel, c0 = c0)), lapply(gammas, function(g) naive_rd(df = panel, gamma = g, c0 = c0)))
cel <- cellini_rd(df = panel, gamma = gammas, c0 = c0, Tn = 7)
tbl <- do.call(rbind, lapply(c(res, list(cel)), `[`, RD_SCHEMA))
print(tbl[, c("method", "gamma", "h", "bw_rule", "est", "se", "ci_l", "ci_u", "t_denom")], row.names = FALSE, digits = 4)
chk(all(is.finite(tbl$est)) && all(tbl$se > 0), "estimates and standard errors are finite and positive")
for (v in c("CR0", "CR1", "CR2", "CR3")) {
  se <- c(static_rd(df = panel, c0 = c0, h = 3, vcov_type = v)$se, naive_rd(df = panel, gamma = 0.8, c0 = c0, h = 8, vcov_type = v)$se,
          cellini_rd(df = panel, gamma = 0.8, c0 = c0, Tn = 7, vcov_type = v)$se)
  chk(all(is.finite(se) & se > 0), paste(v, "runs for every estimator"))
}

hdr("5. Cluster-robust variances agree with direct formulas")
w <- kernel_function("triangular")((panel$Z_it - c0) / 3)
keep <- w > 0
X5 <- rd_design(panel$Z_it[keep] - c0, panel$t[keep], TRUE)
fit <- cluster_fit(X5, panel$Y_it[keep], w[keep], panel$i[keep], "CR0")
bread <- solve(crossprod(X5 * sqrt(w[keep])))
scores <- rowsum(X5 * w[keep] * drop(panel$Y_it[keep] - X5 %*% fit$coef), panel$i[keep])
chk(abs(sum(fit$influence[[1]][, 2]^2) - (bread %*% crossprod(scores) %*% bread)[2, 2]) < 1e-10,
    "CR0 matches the sandwich formula")

hdr("6. Matrix and data.frame inputs agree")
chk(isTRUE(all.equal(static_rd(Z = mats$Z, Y = mats$Y, A = mats$A, c0 = c0), res[[1]])), "static_rd")
chk(isTRUE(all.equal(naive_rd(Z = mats$Z, Y = mats$Y, A = mats$A, gamma = 0.5, c0 = c0), res[[2]])), "naive_rd")
chk(isTRUE(all.equal(cellini_rd(Z = mats$Z, Y = mats$Y, A = mats$A, gamma = gammas, c0 = c0, Tn = 7), cel)), "cellini_rd")

hdr("7. Confidence intervals")
r <- res[[1]]
chk(abs((r$ci_u - r$est) / r$se - cv_normal()) < 1e-10, "ci_u uses the normal quantile")
chk(abs((r$ci_u_ba - r$est) / r$se - cv_bias_aware()) < 1e-10, "ci_u_ba uses the bias-aware value")

hdr("8. Naive bandwidth")
nbw <- naive_ik_bandwidth(df = panel, gamma = 0.8, c0 = c0)
chk(nbw$bandwidth > 0, sprintf("bandwidth %.4f (pilot %.4f, rule %s)", nbw$bandwidth, nbw$h_pilot, nbw$rule))
chk(abs(naive_rd(df = panel, gamma = 0.8, c0 = c0, h = nbw$bandwidth)$est - res[[3]]$est) < 1e-12,
    "supplying the selected bandwidth reproduces the default fit")

hdr("9. Cellini: the recursion inverts their eq. (5) and the delta method is right")
info <- attr(cel, "cellini")
fits <- cellini_lag_fits(df = panel, c0 = c0, Tn = 7)
recon <- vapply(seq_len(7), function(k) sum(c(1, fits$delta_a)[seq_len(k)] * info$theta_Cel[k:1]), 1)
chk(max(abs(recon - fits$delta_y)) < 1e-10, "eq. (5) reconstructs delta_y from theta and delta_a")
target <- function(p) sum(0.8^(0:6) * cellini_recursion(p[1:7], p[-(1:7)])$theta_Cel)
p0 <- c(fits$delta_y, fits$delta_a)
num <- vapply(seq_along(p0), function(k) {
  d <- max(1e-6 * abs(p0[k]), 1e-7)
  (target(replace(p0, k, p0[k] + d)) - target(replace(p0, k, p0[k] - d))) / (2 * d)
}, 1)
chk(abs(sqrt(drop(crossprod(num, fits$Sigma %*% num))) - cel$se[cel$gamma == 0.8]) / cel$se[cel$gamma == 0.8] < 1e-4,
    "the reported standard error matches a numerical Jacobian")
chk(isTRUE(all.equal(cellini_rd(gamma = gammas, c0 = c0, Tn = 7, fits = fits), cel)), "precomputed `fits` reproduce the call")
loc <- cellini_rd(df = panel, gamma = 0.8, c0 = c0, Tn = 7, spec = "local_linear", K = "triangular")
chk(is.finite(loc$est) && loc$se > 0, "the local linear specification runs")

hdr("10. Failure modes are errors")
chk(inherits(try(static_rd(df = panel, c0 = 1e6, h = 1), silent = TRUE), "try-error"), "no data near the threshold")
chk(inherits(try(cellini_rd(df = panel, gamma = 0.8, c0 = c0, Tn = Tn + 50), silent = TRUE), "try-error"),
    "Tn beyond the panel")
chk(inherits(try(as_rd_panel(Z = mats$Z, Y = mats$Y), silent = TRUE), "try-error"), "no A and no c0")

cat("\nALL SMOKE TESTS PASSED\n")
