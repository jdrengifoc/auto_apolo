# Based on the output files of any of `apolo_mh_v\\d.R` files, the code `MIGUEL`
# create `incomplete_individuals.RData`. This code, split the missing
# individuals in `n_files`, a parameter subject to calibration.
library(magrittr)

pp <- readRDS('After/ABC4SFAv2/Data/Inputs/incomplete_individuals.RData')

names(pp) <- gsub('scenario', '', names(pp))

firms_to_repeat <- rep(0, 16)
names(firms_to_repeat) <- names(pp)

for (scenario in names(pp)) {
  names(pp[[scenario]]) <- paste0('sim', 1L:100L)
  n_firms <- 0L
  for (sim in names(pp[[scenario]])) {
    current_n_firms <- length(pp[[scenario]][[sim]])
    n_firms <- n_firms + current_n_firms
    if (current_n_firms == 0L) {
      pp[[scenario]][[sim]] <- NULL
    }
  }
  firms_to_repeat[[scenario]] <- n_firms
}

# Estimate duration
n_files <- 16 # Calibrate
(horas_por_firma <- 13.5 * 24 / (50 * 100))
iters_required4success <- 1 / (1 - firms_to_repeat / (50 * 100))
sum(firms_to_repeat * iters_required4success * horas_por_firma) / (24 * n_files)

(iters_per_file <- floor(sum(firms_to_repeat * iters_required4success) / n_files))


ultimate_pp <- list()
idx_file <- 1L
n_iters_in_current_file <- 0

for (scenario in names(pp)) {
  iters_required <- iters_required4success[[scenario]]
  for (sim in names(pp[[scenario]])) {
    firms <- pp[[scenario]][[sim]]
    n_firms <- length(firms)
    n_required_iters <- n_firms * iters_required
    
    if (n_required_iters +  n_iters_in_current_file <= iters_per_file) {
      ultimate_pp[[sprintf('file%d', idx_file)]][[scenario]][[sim]] <- firms
      n_iters_in_current_file <- n_iters_in_current_file + n_required_iters
    } else {
      # Pasar al siguiente archivo y guardar.
      n_iters_in_current_file <- 0
      idx_file <- idx_file + 1L
      ultimate_pp[[sprintf('file%d', idx_file)]][[scenario]][[sim]] <- firms
    }
  }
}

saveRDS(ultimate_pp, 'After/ABC4SFAv2/Data/Inputs/ultimate_MH_exp.RData')

