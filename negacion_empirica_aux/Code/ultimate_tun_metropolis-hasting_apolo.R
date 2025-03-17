file <- 'file1'
tun <- 1.8
results_filename <- sprintf('Data/Outputs/ultimate_MH_%s_tun%1f.RData', file, tun)


# Dependencies ------------------------------------------------------------
library(dplyr)
inputsMH <- function(residuals, t_periods, variances) {
  sv <- variances[["sigma_v"]]
  sa <- variances[["sigma_v0"]]
  su <- variances[["sigma_u"]]
  seta <- variances[["sigma_u0"]]
  
  A <- cbind(1, diag(t_periods))
  it <- rep(1, t_periods)
  Sigma <- sv^2 * diag(t_periods) + sa^2*it%*%t(it) 
  V <- diag(c(seta^2, rep(su^2, t_periods)))
  iSigma <- solve(Sigma)
  Lambda <- solve(solve(V) + t(A)%*%iSigma%*%A)
  R <- Lambda%*%t(A)%*%iSigma
  # Proposal. Must have the same support (0, inf).
  u <- tmvtnorm::rtmvnorm(1, mean = c(R%*%residuals), sigma = Lambda,
                          lower = rep(0, length = t_periods + 1), 
                          upper = rep(Inf, length = t_periods + 1), 
                          algorithm = 'gibbs')
  return(list(u = u, A = A, Sigma = Sigma, V = V, R = R, 
              Lambda = Lambda, iSigma = iSigma))
}
uDrawMH <- function(inputs_MH, residuals, t_periods, variances, tun = 1) {
  # Unpack parameters.
  u <- inputs_MH$u
  A <- inputs_MH$A
  Sigma <- inputs_MH$Sigma
  iSigma <- inputs_MH$iSigma
  V <- inputs_MH$V
  R <- inputs_MH$R
  Lambda <- tun *  inputs_MH$Lambda
  seta <- variances[["sigma_u0"]]
  su <- variances[["sigma_u"]]
  
  uc <- tmvtnorm::rtmvnorm(1, mean = c(R%*%residuals), sigma = Lambda,
                           lower=rep(0, length = t_periods + 1), 
                           upper=rep(Inf, length = t_periods + 1), 
                           algorithm = 'gibbs')
  quc <- tmvtnorm::dtmvnorm(uc, mean = c(R%*%residuals), sigma = Lambda,
                            lower = rep(0, length = t_periods + 1), 
                            upper = rep(Inf, length = t_periods + 1))
  qu <- tmvtnorm::dtmvnorm(u, mean = c(R%*%residuals), sigma = Lambda,
                           lower = rep(0, length = t_periods + 1), 
                           upper = rep(Inf, length = t_periods + 1))
  fuc <- -0.5*t(residuals - A%*%c(uc))%*%iSigma%*%(residuals - A%*%c(uc)) - 
    t(c(seta, rep(su, t_periods)))%*%t(uc)
  fu <- -0.5*t(residuals - A%*%c(u))%*%iSigma%*%(residuals - A%*%c(u)) - 
    t(c(seta, rep(su, t_periods)))%*%t(u)
  # Criterio de transición. Cómo actualizamos la cadena.
  alpha <- min(exp(fuc-fu) * qu/quc, 1, na.rm = T)
  unif <- runif(1, 0, 1)
  if(unif < alpha){
    unew <- uc
    accept <- 1
  } else{
    unew <- u
    accept <- 0
  }
  inputs_MH$u <- unew
  inputs_MH$accept <- c(inputs_MH$accept, accept)
  return(inputs_MH)
}

# The game is on,  not wait for it needed ---------------------------------
inputs_filepath <- 'Data/Inputs/BK_exp_simData.RData'
results_folder <- 'Data/Outputs/Outputs_exp/'

df_scenarios <- readxl::read_excel('Data/Inputs/Badunenko&Kumbhakar.xlsx')
inputs <- readRDS(inputs_filepath)

ultimate_pp_filepath <- 'Data/Inputs/ultimate_MH_exp.RData'
ultimate_pp <- readRDS(ultimate_pp_filepath)
ultimate_pp <- ultimate_pp[[file]]

S <- 1e4
S_step_accept <- 100
scenarios <- names(ultimate_pp)

X <- inputs$X
resultsMH <- list()
for (scenario in scenarios) {
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
  
  sims <- names(ultimate_pp[[scenario]])
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
    
    firms <- ultimate_pp[[scenario]][[sim]]
    for (firm in firms) {
      ef_reales <- NA
      ef_abc <- NA
      while ( any(c(0, NA) %in% c(ef_abc, ef_reales)) ) {
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
      }
      
      EF_reales <- rbind(EF_reales, ef_reales)
      EF_abc <- rbind(EF_abc, ef_abc)
      sprintf('End scenario %s - simulation %s - firm %d.', 
              scenario, sim, firm) %>% print
      (Sys.time() - tic) %>% print
      
    }
    
    resultsMH[[scenario]][[sim]][['EF_reales']] <- EF_reales
    resultsMH[[scenario]][[sim]][['EF_abc']] <- EF_abc
    saveRDS(resultsMH, results_filename)
  }
}
saveRDS(resultsMH, results_filename)
