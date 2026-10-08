# A small version of the autoreg simulations, outputs saved to demo_results/. 
# The full results (1000 seeds, n = 500, 1000, 2000, 5000, gamma = 0.5, 0.8, 1) 
# are in results/ and print_tables.R prints them. Run from this folder.
rm(list=ls())
if (!requireNamespace("pbmcapply", quietly = TRUE)) install.packages("pbmcapply")
source("simulation.R")

reps <- 100
n_grid <- c(500, 1000)
gammas <- 0.8
out.dir <- "demo_results"
cores <- max(1, parallel::detectCores() - 1)

truth <- lapply(split(settings, settings$setting), 
                oracle, 
                gammas = gammas, 
                cores = cores)
fits <- unlist(
  lapply(seq_len(nrow(settings)), 
         function(s) {
           cat(sprintf("Setting %d\n", s))
           pbmcapply::pbmclapply(seq_len(reps), 
                                 fit_seed, 
                                 s = settings[s, ], 
                                 n_grid = n_grid, 
                                 gammas = gammas,
                                 mc.cores = cores, 
                                 ignore.interactive = TRUE)
           }), 
  recursive = FALSE)

bind <- function(x, part) do.call(rbind, lapply(x, `[[`, part))
res <- list(targets = bind(truth, "targets"), 
            targets_lags = bind(truth, "lags"),
            raw = bind(fits, "rows"), 
            raw_lags = bind(fits, "lags"))
res$scores <- score(res$raw, res$targets, reps)
res$lag_scores <- score_lags(res$raw_lags, res$targets_lags)
dir.create(out.dir, showWarnings = FALSE)
for (f in names(res)) write.csv(res[[f]], file.path(out.dir, paste0(f, ".csv")), row.names = FALSE)
print_tables(res$scores, res$lag_scores, reps)
