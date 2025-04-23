
# Requirements ------------------------------------------------------------
setwd('auto_apolo')
rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))
library(mvtnorm)
library(dplyr)
library(stringr)

# Globals -----------------------------------------------------------------
FOLDER_INPUT <- "abc4sfa/1_simulation_resources/Data/Inputs/"
INPUT_DATA <- readRDS(file.path(FOLDER_INPUT, "BK_simData.RData"))

FOLDER_GIBBS <- '../../ABC4SFA/After/ABC4SFAv2/Data/Outputs/Gibs_experiment'
FOLDER_ABC <- 'abc4sfa/1_simulation_resources/Data/Outputs/Simulation hf'
FOLDER_OUTPUT <- "abc4sfa/4_error_measurements/data/output/simulations/hf"
dir.create(FOLDER_OUTPUT, recursive = T)

# ABC: Compute efficiency estimates ---------------------------------------




#' It is also a common practice in ABC to perform a regression adjustment after 
#' retained draws (Beaumont et al., 2002; Leuenberger and Wegmann, 2010;
#' Sisson et al., 2018)
# Preallocate.
# Without adjustment.
scenarios <- paste0('s', 1:16)

efficiency_abc <- list()

scenario <- scenarios[1]
sim <- 'sim1'
# For each scenario.
for (scenario in scenarios) {
  # Read results.
  outputData <- readRDS(sprintf("%s/BK_%s.RData", FOLDER_ABC, scenario))
  
  # Extract real betas and covariates.
  X <- INPUT_DATA[['X']]
  betas <- INPUT_DATA[[scenario]]$params$beta
  # Extract simulations realized.
  sims <- names(outputData[[scenario]])
  
  # Fix inverted sigmas in experiments.
  sigmas <- get_real_sigmas(INPUT_DATA, scenario)
  Ts <- INPUT_DATA[[scenario]]$params$t
  n <- INPUT_DATA[[scenario]]$params$n
  theta <- list(n = n, t = Ts)
  
  efficiency_abc[[scenario]] <- list()
  
  for (sim in sims) {
    tic <- Sys.time()
    message(scenario, sim)
    y <- INPUT_DATA[[scenario]][[sim]]$y
    est_ABCparams <- apply(outputData[[scenario]][[sim]]$ABCpostChain, 2, mean)
    RegY <- lm(y ~ .,  data = data.frame(y, X[,2:3]))
    SumStatY <- getSampleSumStats(theta = theta, y = y, X = X)
    est_betas <- c(est_ABCparams[1], RegY$coefficients[-1])
    est_sigmas <- get_abc_sigmas(est_ABCparams)
    PostDraws <- outputData[[scenario]][[sim]]$ABCpostChain
    K <- dim(PostDraws)[1]
    SumStatPar <- matrix(0, K, 5)
    for (k in 1:K){
      ResSim <- simRes(theta = theta, prior = PostDraws[k, 1:5])#, dist = "halfnormal")
      SumStatPar[k, ] <- getSumStats(theta = theta, resS = ResSim)
    }
    
    TI_real <- ColombiExpectation(n, Ts, y, X, betas, sigmas, p = 1)
    TI_est <- ColombiExpectation(n, Ts, y, X, est_betas, est_sigmas, p = 1)
    
    ls0 <- list(
      TE_real = t(sapply(TI_real, inv))[, 1:(Ts+1)],
      TE_est = t(sapply(TI_est, inv))[, 1:(Ts+1)]
    )
    efficiency_abc[[scenario]][[sim]] <- ls0
    print(Sys.time() - tic)
  }
}


efficiency_abc %>% 
  saveRDS(
    sprintf(
      "%s/efficiency_abc_hf_%s.rds",
      FOLDER_OUTPUT, Sys.Date()
    )
  )

# Gibbs: Efficiency --------------------------------------------------------


media <- function(x) {
  mean(x, na.rm = TRUE)
}

# Parameters
p <- 1L

# Read data.
files <- list.files(FOLDER_GIBBS, full.names = T)

efficiency_abc <- list()
for (file in files) {
  data <- readRDS(file)
  sims <- names(data)
  scenario <- stringr::str_extract(file, "s\\d+")
  # Extract params.
  params <- INPUT_DATA[[scenario]]$params
  n <- params$n
  t <- params$t
  X <- INPUT_DATA$X
  betas <- INPUT_DATA[[scenario]]$params$beta
  sigmas <- get_real_sigmas(INPUT_DATA, scenario)
  
  efficiency_abc[[scenario]] <- list()
  for (sim in sims[1:5]) {
    tic <- Sys.time()
    message("Scenario ", scenario, ". Simulation: ", sim)
    # Read sim related input/outputs.
    y <- INPUT_DATA[[scenario]][[sim]]$y
    gibbs_results <- data[[sim]][['postChain']]
    sigmas_gibbs <- get_gibbs_sigmas(gibbs_results)
    betas_gibbs <- apply(gibbs_results$Thetachain, 1L, mean)
    
    TI_real <- ColombiExpectation(n, t, y, X, betas, sigmas, p)
    TI_gibbs <- ColombiExpectation(n, t, y, X, betas_gibbs, sigmas_gibbs, p)
    
    ls0 <- list(
      TE_real = t(sapply(TI_real, inv))[, 1:(t+1)],
      TE_est = t(sapply(TI_gibbs, inv))[, 1:(t+1)]
    )
    efficiency_abc[[scenario]][[sim]] <- ls0
    print(Sys.time() - tic)
  }
}
# Save data.
scenarios <- stringr::str_extract(files, "s[0-9]+")
(df <- data.frame(
  scenario = scenarios,
  RelativeBiasPerm = RelativeBias[,1],
  RelativeBiasTrans = RelativeBias[,2],
  RelativeBiasAll= RelativeBias[,3],
  UpperBiasPerm = UpperBias[,1],
  UpperBiasTrans = UpperBias[,2],
  UpperBiasAll = UpperBias[,3],
  PearsonCorrCoefPerm = PearsonCorrCoef[,1],
  PearsonCorrCoefTrans = PearsonCorrCoef[,2],
  PearsonCorrCoefAll = PearsonCorrCoef[,3],
  RelativeMSEPerm = RelativeMSE[,1],
  RelativeMSETrans = RelativeMSE[,2],
  RelativeMSEAll = RelativeMSE[,3]
))

tibble(
  scenario = rep(scenarios, each = 3),
  efficiency_type = rep(c('Persistent', 'Transient', 'Overall'), times = length(scenarios)),
  relative_bias = c(t(RelativeBias)),
  pearson_correlation = c(t(PearsonCorrCoef)),
  relative_MSE = c(t(RelativeMSE))
) %>% 
  pivot_longer(
    cols = c('relative_bias', 'pearson_correlation', 'relative_MSE'),
    names_to = 'metric_type', values_to = 'metric_value'
  ) %>%
  fix_inverted_ids
  writexl::write_xlsx(
    sprintf("%s/efficiency_metrics_gibbs_hf_%s.xlsx", FOLDER_OUTPUT, Sys.Date())
  )



# ABC + Gibbs: Sigmas -------------------------------------------------------------

scenarios <- paste0('s', c(1, 14))

files_gibbs <- list.files(FOLDER_GIBBS, full.names = T)
files_abc <- list.files(FOLDER_ABC, full.names = T)

df <- NULL
scenario <- scenarios[1]
for (scenario in scenarios) {
  # Ground
  sigmas <- get_real_sigmas(INPUT_DATA, scenario)
  
  # ABC.
  file_abc <- files_abc[str_extract(files_abc, 's\\d+') == scenario]
  abc_data <- readRDS(file_abc)
  sigmas_est_abc <- sapply(
    abc_data[[scenario]], function(sim) {get_abc_sigmas(sim)}
  )
  relative_bias_abc <- apply(sigmas_est_abc / sigmas - 1, 1, mean)
  root_mse_abc <- apply((sigmas_est_abc - sigmas)^2, 1, mean) %>% sqrt
  
  # Gibbs
  file_gibbs <- files_gibbs[str_extract(files_gibbs, 's\\d+') == scenario]
  if (!identical(file_gibbs, character(0))) {
    gibbs_data <- readRDS(file_gibbs)
    
    sigmas_est_gibbs <- sapply(
      gibbs_data, function(sim) {get_gibbs_sigmas(sim$postChain)}
    )
    relative_bias_gibbs <- apply(sigmas_est_gibbs / sigmas - 1, 1, mean)
    root_mse_gibbs <- apply((sigmas_est_gibbs - sigmas)^2, 1, mean) %>% sqrt
  } else {
    relative_bias_gibbs <- NULL
    root_mse_gibbs <- NULL
  }
  
  # Save results.
  df0 <- left_join(
    tibble(
      scenario = rep(scenario, length(sigmas)),
      sigmas = paste0('sigma_', c('v', 'v0',	'u', 'u0')),
      relative_bias_abc,
      root_mse_abc
    ),
    tibble(
      scenario = rep(scenario, length(sigmas)),
      sigmas = paste0('sigma_', c('v', 'v0',	'u', 'u0')),
      relative_bias_gibbs,
      root_mse_gibbs
    ),
    by = c('scenario', 'sigmas')
  )
  
  df <- bind_rows(df, df0)
}

df %>% 
  fix_inverted_ids %>% 
  writexl::write_xlsx(
    sprintf(
      "%s/sigmas_metrics_hf_%s.xlsx",
      FOLDER_OUTPUT, Sys.Date()
    )
  )

