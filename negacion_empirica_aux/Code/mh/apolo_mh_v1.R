library(dplyr)
source('Code/Fixed/functions.R')
# This code make metropolis-hasting with different tuns one per estimate.

# Parameters
S <- 1e4
S_step_accept <- 100
scenarios <- c('s15')
tun <- 1.8


# The game is on,  not wait for it needed ---------------------------------
inputs_filepath <- 'Data/Inputs/BK_exp_simData.RData'
results_folder <- 'Data/Outputs/Outputs_exp/'

df_scenarios <- readxl::read_excel('Data/Inputs/Badunenko&Kumbhakar.xlsx')
inputs <- readRDS(inputs_filepath)

X <- inputs$X
resultsMH <- list()
for (scenario in scenarios) {
  output_file <-sprintf("Data/Outputs/MH%s_tun%.1f.RData", scenario, tun) 
  resultsMH[[scenario]] <- list()
  # Get real parameters
  n_periods <- inputs[[scenario]][['params']][['n']]
  t_periods <- dim(X)[1] / n_periods
  variances <- df_scenarios %>% dplyr::filter(ID_inverted == scenario) %>% 
    dplyr::select(starts_with('sigma'))
  betas <- inputs[[scenario]][['params']][['beta']]
  
  # Save in list.
  resultsMH[[scenario]][['variances']] <- variances
  resultsMH[[scenario]][['betas']] <- betas
  
  # Get abc results.
  results_filepath <- list.files(results_folder, full.names = T,
                                 pattern = sprintf('%s\\.', scenario))
  abc_results <- readRDS(results_filepath)
  
  sims <- names(inputs[[scenario]])
  sims <- sims[grep('^sim\\d+$', sims)]
  sims <- paste0('sim', 1:100)
  idx_firms <- rep(1:n_periods, each = t_periods)
  for (sim in sims) {
    # Compute inputs for initial.
    y <- inputs[[scenario]][[sim]][['y']]
    ## Real.
    perturbations <- y - X%*%betas
    
    # ABC.
    abc_estimates <- abc_results[[scenario]][[sim]][['ABCpostChain']] %>%
      apply(2, mean)
    abc_betas <- lm(y ~ X[, 2L:3L])
    abc_betas <- c(abc_estimates[['B0']], 
                   abc_betas$coefficients[-1L] %>% unname)
    residuals <- y - X%*%abc_betas
    abc_variances <- abc_estimates[2L:5L]
    
    # Save in list.
    resultsMH[[scenario]][[sim]] <- list()
    resultsMH[[scenario]][[sim]][['perturbations']] <- perturbations
    resultsMH[[scenario]][[sim]][['residuals']] <- residuals
    resultsMH[[scenario]][[sim]][['abc_betas']] <- abc_betas
    resultsMH[[scenario]][[sim]][['abc_variances']] <- abc_variances
    
    EF_abc <- NULL
    EF_reales <- NULL
    
    for (firm in 1:n_periods) {
      sprintf('Began scenario %s - simulation %s - firm %d',
              scenario, sim, firm) %>% print
      tic <- Sys.time()
      current_perturbations <- perturbations[firm == idx_firms]
      current_residuals <- residuals[firm == idx_firms]
      inputs_MH_reales <- inputsMH(current_perturbations, t_periods, variances)
      inputs_MH_abc <- inputsMH(current_residuals, t_periods, abc_variances)
      
      draws_reales <- NULL
      draws_abc <- NULL
      for (s in 1:S) {
        inputs_MH_reales <- uDrawMH(inputs_MH_reales, current_perturbations, t_periods, variances, tun)
        inputs_MH_abc <- uDrawMH(inputs_MH_abc, current_residuals, t_periods, abc_variances, tun)
        if (s %% S_step_accept == 0) {
          draws_reales <- rbind(draws_reales, inputs_MH_reales$u)
          draws_abc <- rbind(draws_abc, inputs_MH_abc$u)
        }
      }
      ef_reales <- colMeans(apply(-draws_reales, 2, exp))
      ef_abc <- colMeans(apply(-draws_abc, 2, exp))
      EF_reales <- rbind(EF_reales, ef_reales)
      EF_abc <- rbind(EF_abc, ef_abc)
      sprintf('End scenario %s - simulation %s - firm %d.', 
              scenario, sim, firm) %>% print
      (Sys.time() - tic) %>% print
      
    }
    
    resultsMH[[scenario]][[sim]][['EF_reales']] <- EF_reales
    resultsMH[[scenario]][[sim]][['EF_abc']] <- EF_abc
    saveRDS(resultsMH, output_file)
  }
}
saveRDS(resultsMH, output_file)
