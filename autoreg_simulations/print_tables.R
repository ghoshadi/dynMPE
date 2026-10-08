# Tables of the full results: Settings 1-3, 
# 1000 seeds, n = 500, 1000, 2000, 5000, read from results/.
# Run from this folder.
source("simulation.R")

full <- lapply(c(scores = "scores", 
                 lag_scores = "lag_scores", 
                 targets = "targets"),
               function(f){
                 read.csv(file.path("results", paste0(f, ".csv")))
                 }
               )
print_targets(full$targets)
print_tables(full$scores, full$lag_scores, 1000)
at_a_glance(full$scores)
