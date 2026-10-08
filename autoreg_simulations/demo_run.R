# A small version of the simulation (gamma = 0.8 only), saved to demo_results/. The full results
# (1000 seeds, n = 500, 1000, 2000, 5000, gamma = 0.5, 0.8, 1) are in results/; print_tables.R
# prints them. Run from this folder.
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply", repos = "https://cloud.r-project.org")
source("simulation.R")

reps <- 100
n_grid <- c(500, 1000)
gammas <- 0.8
out <- "demo_results"
cores <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores() - 1)

truth <- lapply(split(settings, settings$setting), oracle, gammas = gammas, cores = cores)
fits <- unlist(lapply(seq_len(nrow(settings)), function(s) {
  cat(sprintf("Setting %d\n", s))
  pbmcapply::pbmclapply(seq_len(reps), fit_seed, s = settings[s, ], n_grid = n_grid, gammas = gammas,
                        mc.cores = cores, ignore.interactive = TRUE)
}), recursive = FALSE)

bind <- function(x, part) do.call(rbind, lapply(x, `[[`, part))
res <- list(targets = bind(truth, "targets"), targets_lags = bind(truth, "lags"),
            raw = bind(fits, "rows"), raw_lags = bind(fits, "lags"))
res$scores <- score(res$raw, res$targets, reps)
res$lag_scores <- score_lags(res$raw_lags, res$targets_lags)
dir.create(out, showWarnings = FALSE)
for (f in names(res)) write.csv(res[[f]], file.path(out, paste0(f, ".csv")), row.names = FALSE)
print_tables(res$scores, res$lag_scores, reps)
