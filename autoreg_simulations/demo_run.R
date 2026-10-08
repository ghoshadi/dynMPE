# A small version of the simulation, saved to demo_results/. The full results (1000 seeds,
# n = 500, 1000, 2000, 5000, same seeds) are in results/; print_tables.R prints them.
# Run from this folder.
source("simulation.R")

reps <- 100
n_grid <- c(500, 1000)
gammas <- c(0.5, 0.8, 1)
out <- "demo_results"
cores <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores() - 1)

truth <- lapply(split(settings, settings$setting), oracle, gammas = gammas, cores = cores)
fits <- parallel::mclapply(seq_len(nrow(settings) * reps), function(k)
  fit_seed((k - 1) %% reps + 1, settings[(k - 1) %/% reps + 1, ], n_grid, gammas), mc.cores = cores)

bind <- function(x, part) do.call(rbind, lapply(x, `[[`, part))
res <- list(targets = bind(truth, "targets"), targets_lags = bind(truth, "lags"),
            raw = bind(fits, "rows"), raw_lags = bind(fits, "lags"))
res$scores <- score(res$raw, res$targets, reps)
res$lag_scores <- score_lags(res$raw_lags, res$targets_lags)
dir.create(out, showWarnings = FALSE)
for (f in names(res)) write.csv(res[[f]], file.path(out, paste0(f, ".csv")), row.names = FALSE)
print_tables(res$scores, res$lag_scores, reps)
