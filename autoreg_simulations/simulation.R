# Autoregressive simulation for Ghosh and Wager (arXiv:2512.15244).
if (!requireNamespace("dynMPE", quietly = TRUE)) remotes::install_github("ghoshadi/dynMPE")
library(dynMPE)
for (f in c("shared_helpers", "static_rd", "naive_rd", "cellini_rd"))
  source(file.path("..", "baselines", paste0(f, ".R")))

# Z_{t+1} = mu + delta (t+1) + rho (Z_t - mu) - A_t e_t + sigma eps_t,  A_t = 1{Z_t >= c0},
# Y_t = -Z_{t+1} + nu_t,  e_t = tau (Z_t - mu)_+ in Settings 1-2 and tau (c0 - mu) in Setting 3.
dgp <- list(c0 = 110, horizon = 12, mu = 100, sigma = 4, rho = 0.9, tau = 0.1, noise = 1)
settings <- data.frame(setting = 1:3, delta = c(0, 1, 0), const = c(FALSE, FALSE, TRUE))

effect <- function(s, zc, p) if (s$const) p$tau * (p$c0 - p$mu) else p$tau * pmax(zc, 0)

simulate_z <- function(n, s, thresh, p = dgp) {
  z <- matrix(NA_real_, n, p$horizon + 1)
  z[, 1] <- p$mu + p$sigma / sqrt(1 - p$rho^2) * rnorm(n)
  for (t in seq_len(p$horizon)) {
    zc <- z[, t] - p$mu
    z[, t + 1] <- p$mu + s$delta * t + p$rho * zc - (z[, t] >= thresh) * effect(s, zc, p) + p$sigma * rnorm(n)
  }
  z
}

# tau_RD by a central difference in the threshold; the per-lag jumps at Z_0 = c0 by setting
# A_0 to 1 and to 0 with the same innovations.
oracle <- function(s, gammas, n_oracle = 1e6, cores = 1, h = 0.25, batch = 25000, seed = 1e6, p = dgp) {
  sizes <- diff(unique(c(seq(0, n_oracle, by = batch), n_oracle)))
  w <- outer(0:(p$horizon - 1), gammas, function(t, g) g^t)
  sums <- function(f) Reduce(`+`, parallel::mclapply(seq_along(sizes), f, mc.cores = cores)) / n_oracle
  value <- sums(function(k) {
    paths <- lapply(p$c0 + c(-h, h), function(thresh) {
      set.seed(seed + k)
      z <- simulate_z(sizes[k], s, thresh, p)
      list(y = -z[, -1] %*% w, a = (z[, -(p$horizon + 1)] >= thresh) %*% w)
    })
    dy <- paths[[2]]$y - paths[[1]]$y
    da <- paths[[2]]$a - paths[[1]]$a
    rbind(colSums(dy), colSums(da), colSums(dy^2), colSums(da^2), colSums(dy * da))
  })
  jump <- sums(function(k) {
    set.seed(seed + 2e6 + k)
    eps <- matrix(rnorm(sizes[k] * p$horizon), sizes[k], p$horizon)
    paths <- lapply(1:0, function(a0) {
      z <- rep(p$c0, sizes[k]); a <- rep(a0, sizes[k]); y <- a_path <- matrix(0, sizes[k], p$horizon)
      for (t in seq_len(p$horizon)) {
        a_path[, t] <- a
        zc <- z - p$mu
        z <- p$mu + s$delta * t + p$rho * zc - a * effect(s, zc, p) + p$sigma * eps[, t]
        y[, t] <- -z
        a <- z >= p$c0
      }
      list(y = y, a = a_path)
    })
    rbind(colSums(paths[[1]]$y - paths[[2]]$y), colSums(paths[[1]]$a - paths[[2]]$a))
  })
  stopifnot(jump[2, 1] == 1)
  lags <- data.frame(setting = s$setting, m = 0:(p$horizon - 1), delta_y = jump[1, ], delta_a = jump[2, ],
                     theta_Cel = cellini_recursion(jump[1, ], jump[2, -1])$theta_Cel)
  ratio <- value[1, ] / value[2, ]
  var_ratio <- value[3, ] - value[1, ]^2 + ratio^2 * (value[4, ] - value[2, ]^2) -
    2 * ratio * (value[5, ] - value[1, ] * value[2, ])
  list(lags = lags, targets = data.frame(
    setting = s$setting, gamma = gammas, tau_RD = ratio,
    se_tau_RD = sqrt(pmax(var_ratio, 0) / n_oracle) / abs(value[2, ]),
    tau_PE = p$tau * (p$c0 - p$mu),
    tau_naive = drop(t(w) %*% lags$delta_y) / drop(t(w) %*% lags$delta_a),
    tau_Cel = drop(t(w) %*% lags$theta_Cel)))
}

fit_columns <- c("est", "se", "ci_l", "ci_u", "ci_l_ba", "ci_u_ba", "h", "bw_rule")
lag_columns <- c("m", "theta_Cel", "Cel_ci_l", "Cel_ci_u")

proposed <- function(d, gamma, c0, kernel, vcov) {
  f <- dynMPE(d, "Y_it", "Z_it", "t", "i", c0, gamma = gamma, kernel = kernel, vcov = vcov)
  q <- qnorm(0.975)
  data.frame(est = f$estimate, se = f$se, ci_l = f$estimate - q * f$se, ci_u = f$estimate + q * f$se,
             ci_l_ba = f$ci[["lower"]], ci_u_ba = f$ci[["upper"]], h = f$bandwidth, bw_rule = f$bandwidth_rule)
}

fit_rows <- function(method, gamma, fit) {
  if (is.character(fit))
    return(data.frame(method, gamma, as.list(setNames(rep(NA_real_, length(fit_columns)), fit_columns)), err = fit))
  data.frame(method, gamma, fit[fit_columns], err = NA_character_)
}

# The smaller samples are the first n units of the largest.
fit_seed <- function(seed, s, n_grid, gammas, vcov = "CR0", kernel = "triangular", seed_base = 7e6, p = dgp) {
  set.seed(seed_base + seed)
  z <- simulate_z(max(n_grid), s, p$c0, p)
  zt <- z[, -(p$horizon + 1)]
  full <- as_rd_panel(Z = zt, Y = -z[, -1] + p$noise * matrix(rnorm(length(zt)), nrow(zt)), A = (zt >= p$c0) * 1)
  attempt <- function(expr) tryCatch(suppressMessages(expr), error = conditionMessage)
  by_n <- lapply(n_grid, function(n) {
    d <- full[full$i <= n, ]
    cel <- attempt(cellini_rd(df = d, gamma = gammas, c0 = p$c0, Tn = p$horizon, vcov_type = vcov))
    rows <- rbind(
      fit_rows("static", gammas, attempt(static_rd(df = d, c0 = p$c0, K = kernel, vcov_type = vcov))),
      fit_rows("Cellini2010", gammas, if (is.character(cel)) cel else cel[cel$method == "Cellini2010", ]),
      do.call(rbind, lapply(gammas, function(g) rbind(
        fit_rows("dynamic", g, attempt(proposed(d, g, p$c0, kernel, vcov))),
        fit_rows("naive", g, attempt(naive_rd(df = d, gamma = g, c0 = p$c0, K = kernel, vcov_type = vcov)))))))
    list(rows = cbind(n = n, rows), lags = if (is.data.frame(cel)) cbind(n = n, attr(cel, "cellini")[lag_columns]))
  })
  stack <- function(part) cbind(setting = s$setting, seed = seed, do.call(rbind, lapply(by_n, `[[`, part)))
  list(rows = stack("rows"), lags = stack("lags"))
}

method_labels <- c(static = "Standard LLR", dynamic = "Proposed LLR (1.96)", dynamic_ba = "Proposed LLR (2.18)",
                   naive = "Naive long-run LLR", Cellini2010 = "Cellini et al. (2010)")
own_target <- c(static = "tau_PE", dynamic = "tau_RD", dynamic_ba = "tau_RD", naive = "tau_naive",
                Cellini2010 = "tau_Cel")
target_names <- c("tau_RD", "tau_PE", "tau_naive", "tau_Cel")

# Scoring "rd" evaluates every method against tau_RD, "own" each against its own target.
score <- function(raw, targets, reps) {
  ba <- transform(raw[raw$method == "dynamic", ], method = "dynamic_ba", ci_l = ci_l_ba, ci_u = ci_u_ba)
  d <- merge(rbind(raw, ba), targets, by = c("setting", "gamma"))
  d$ok <- is.finite(d$est + d$se + d$ci_l + d$ci_u)
  values <- as.matrix(d[target_names])
  do.call(rbind, lapply(c("rd", "own"), function(scoring) {
    d$target_name <- if (scoring == "rd") "tau_RD" else unname(own_target[d$method])
    d$target <- values[cbind(seq_len(nrow(d)), match(d$target_name, target_names))]
    do.call(rbind, lapply(split(d, list(d$setting, d$gamma, d$method, d$n), drop = TRUE), function(g) {
      o <- g[g$ok, ]
      data.frame(scoring, setting = g$setting[1], gamma = g$gamma[1], method = g$method[1],
                 n = g$n[1], target_name = g$target_name[1], target = g$target[1],
                 coverage = mean(o$ci_l <= o$target & o$target <= o$ci_u),
                 width_mean = mean(o$ci_u - o$ci_l), width_med = median(o$ci_u - o$ci_l),
                 bias = mean(o$est) - g$target[1], rmse = sqrt(mean((o$est - g$target[1])^2)),
                 reps_ok = nrow(o), reps_failed = reps - nrow(o))
    }))
  }))
}

score_lags <- function(lags, truth) {
  d <- merge(lags, truth, by = c("setting", "m"), suffixes = c("", "_true"))
  do.call(rbind, lapply(split(d, list(d$setting, d$n, d$m), drop = TRUE), function(g) data.frame(
    setting = g$setting[1], n = g$n[1], m = g$m[1],
    cov_Cel = mean(g$Cel_ci_l <= g$theta_Cel_true & g$theta_Cel_true <= g$Cel_ci_u, na.rm = TRUE),
    medw_Cel = median(g$Cel_ci_u - g$Cel_ci_l, na.rm = TRUE),
    bias_Cel = mean(g$theta_Cel, na.rm = TRUE) - g$theta_Cel_true[1])))
}

print_tables <- function(scores, lag_scores, reps) {
  cell <- function(cv, w) if (is.finite(cv)) sprintf("%5.1f %8.2f", 100 * cv, w) else sprintf("%14s", "--")
  table <- function(row_labels, first, cols, cov, width) {
    w <- max(nchar(c(row_labels, first)))
    cat(sprintf("%-*s", w, ""), sprintf("   %14s", cols), "\n", sprintf("%-*s", w, first),
        rep(sprintf("   %5s %8s", "cvrg", "width"), length(cols)), "\n", sep = "")
    for (k in seq_along(row_labels))
      cat(sprintf("%-*s", w, row_labels[k]), paste0("   ", mapply(cell, cov[k, ], width[k, ])), "\n", sep = "")
  }
  grid <- function(d, by, keys, across, value) tapply(d[[value]], list(d[[by]], d[[across]]), `[`, 1)[keys, , drop = FALSE]
  for (scoring in c("rd", "own")) for (s in unique(scores$setting)) for (g in sort(unique(scores$gamma))) {
    d <- scores[scores$scoring == scoring & scores$setting == s & scores$gamma == g, ]
    keys <- intersect(names(method_labels), d$method)
    lab <- if (scoring == "rd") method_labels[keys] else paste0(method_labels[keys], " [for ", own_target[keys], "]")
    cat(sprintf("\nSetting %d, gamma = %g, %s, %d seeds\n", s, g,
                if (scoring == "rd") "every method vs tau_RD" else "each method vs its own target", reps))
    table(lab, "method", paste("n =", sort(unique(d$n))),
          grid(d, "method", keys, "n", "coverage"), grid(d, "method", keys, "n", "width_med"))
    u <- d[!duplicated(d$target_name), ]
    cat("targets: ", paste(sprintf("%s = %.4f", u$target_name, u$target), collapse = ", "), "\n", sep = "")
  }
  l <- lag_scores[lag_scores$n == max(lag_scores$n), ]
  m <- as.character(sort(unique(l$m)))
  cat(sprintf("\nCellini et al. (2010) per-lag estimates, n = %d, %d seeds\n", max(l$n), reps))
  table(m, "lag", paste("Setting", sort(unique(l$setting))),
        grid(l, "m", m, "setting", "cov_Cel"), grid(l, "m", m, "setting", "medw_Cel"))
}

center <- function(x, w) formatC(paste0(strrep(" ", (w - nchar(x)) %/% 2), x), width = -w)

print_targets <- function(targets) {
  for (s in sort(unique(targets$setting))) {
    d <- targets[targets$setting == s, ]
    d <- d[order(d$gamma), ]
    cat(sprintf("\nTargets, Setting %d\n%-9s", s, ""), sprintf("   %12s", paste("gamma =", d$gamma)), "\n", sep = "")
    for (k in target_names) cat(sprintf("%-9s", k), sprintf("   %12.4f", d[[k]]), "\n", sep = "")
  }
}

at_a_glance <- function(scores) {
  n_max <- max(scores$n)
  d <- scores[scores$n == n_max, ]
  d <- merge(d[d$scoring == "rd", c("setting", "gamma", "method", "coverage", "width_med")],
             d[d$scoring == "own", c("setting", "gamma", "method", "coverage")], by = c("setting", "gamma", "method"))
  ss <- sort(unique(d$setting))
  w <- max(nchar(method_labels))
  for (g in sort(unique(d$gamma))) {
    cat(sprintf("\ngamma = %g, n = %d, cvrg against tau_RD, cvrg* against own target\n", g, n_max))
    cat(sprintf("%-*s", w, ""), paste0("   ", center(paste("Setting", ss), 19)), "\n",
        sprintf("%-*s", w, "method"), rep(sprintf("   %5s %5s %7s", "cvrg", "cvrg*", "width"), length(ss)), "\n", sep = "")
    for (m in names(method_labels)) {
      y <- d[d$gamma == g & d$method == m, ][match(ss, d$setting[d$gamma == g & d$method == m]), ]
      cat(sprintf("%-*s", w, method_labels[[m]]),
          sprintf("   %5.1f %5.1f %7.2f", 100 * y$coverage.x, 100 * y$coverage.y, y$width_med), "\n", sep = "")
    }
  }
}
