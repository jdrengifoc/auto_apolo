
# Requirements ------------------------------------------------------------
rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))
library(mvtnorm)
library(dplyr)

# Globals -----------------------------------------------------------------
FOLDER_INPUT <- "abc4sfa/1_simulation_resources/Data/Inputs/"
INPUT_DATA <- readRDS(file.path(FOLDER_INPUT, "BK_simData.RData"))

FOLDER_GIBBS <- '../../ABC4SFA/After/ABC4SFAv2/Data/Outputs/Gibs_experiment'
FOLDER_ABC <- 'abc4sfa/1_simulation_resources/Data/Outputs/Simulation hf'
FOLDER_OUTPUT <- "abc4sfa/4_error_measurements/data/output/simulations/hf"
dir.create(FOLDER_OUTPUT, recursive = T)

# ABC: Efficiency --------------------------------------------------------




#' It is also a common practice in ABC to perform a regression adjustment after 
#' retained draws (Beaumont et al., 2002; Leuenberger and Wegmann, 2010;
#' Sisson et al., 2018)
# Preallocate.
# Without adjustment.
ID <- paste0('s', 1:16)
n_pp <- length(ID)
RelativeBias <- matrix(0, n_pp, 3)
UpperBias <- matrix(0, n_pp, 3)
PearsonCorrCoef <- matrix(0, n_pp, 3)
RelativeMSE <- matrix(0, n_pp, 3)


s <- 1
# For each scenario.
for (id in ID){
  # Read results.
  outputData <- readRDS(sprintf("%s/BK_%s.RData", FOLDER_ABC, id))
  
  # Extract real betas and covariates.
  X <- INPUT_DATA[['X']]
  betas <- INPUT_DATA[[id]]$params$beta
  # Extract simulations realized.
  sims <- names(outputData[[1]])
  nsim <- length(sims)
  # Fix inverted sigmas in experiments.
  sigmas <- sigmas <- get_real_sigmas(INPUT_DATA, id)
  Ts <- INPUT_DATA[[id]]$params$t
  n <- INPUT_DATA[[id]]$params$n
  theta <- list(n = n, t = Ts)
  aux_idx <- c(1, 2:(Ts+1), 2:(Ts+1))
  RBs <- matrix(0, nsim, 3)
  UBs <- matrix(0, nsim, 3)
  RMSEs <- matrix(0, nsim, 3)
  PCCs <- matrix(0, nsim, 3)
  l <- 1
  for (sim in sims){
    tic <- Sys.time()
    
    y <- INPUT_DATA[[id]][[sim]]$y
    # est_ABCparams <- outputData[[id]][[sim]]$ABCpostChain[1,]
    est_ABCparams <- apply(outputData[[id]][[sim]]$ABCpostChain, 2, mean)
    RegY <- lm(y ~ .,  data = data.frame(y, X[,2:3]))
    SumStatY <- getSampleSumStats(theta = theta, y = y, X = X)
    est_betas <- c(est_ABCparams[1],
                   RegY$coefficients[-1])
    est_sigmas <- get_abc_sigmas(est_ABCparams)
    PostDraws <- outputData[[id]][[sim]]$ABCpostChain
    K <- dim(PostDraws)[1]
    SumStatPar <- matrix(0, K, 5)
    for (k in 1:K){
      ResSim <- simRes(theta = theta, prior = PostDraws[k, 1:5])#, dist = "halfnormal")
      SumStatPar[k, ] <- getSumStats(theta = theta, resS = ResSim)
    }
    
    TI_real <- ColombiExpectation(n, Ts, y, X, betas, sigmas, p = 1)
    TE_real <- lapply(TI_real, inv)
    TI_est <- ColombiExpectation(n, Ts, y, X, est_betas, est_sigmas, p = 1)
    TE_est <- lapply(TI_est, inv)
    
    mean_TI_est <- apply(as.matrix(do.call(rbind, TI_est)), 2, mean, na.rm = T)
    mean_TI_real <- apply(as.matrix(do.call(rbind, TI_real)), 2, mean, na.rm = T)
    
    mean_TE_est <- apply(as.matrix(do.call(rbind, TE_est)), 2, mean, na.rm = T)
    mean_TE_real <- apply(as.matrix(do.call(rbind, TE_real)), 2, mean, na.rm = T)
    
    RB <- 0
    UB <- 0
    RMSE <- 0
    PCC_cov <- 0
    PCC_sd1 <- 0
    PCC_sd2 <- 0
    inan_real <- which(is.nan(c(unlist(lapply(TE_real, sum)))))
    inan_est <- which(is.nan(c(unlist(lapply(TE_est, sum)))))
    inan <- union(inan_real,inan_est)
    ids <- 1:n
    idUs <- setdiff(ids, inan) 
    nr <- length(idUs)
    for (i in idUs){
      TEi_est <- TE_est[[i]]
      TEi_real <- TE_real[[i]]
      RB <- RB + (TEi_est / TEi_real - 1)
      (UB <- UB + as.integer(TEi_est > TEi_real))
      PCC_cov <- PCC_cov + (TEi_est - mean_TE_est) * (TEi_real - mean_TE_real)
      PCC_sd1 <- PCC_sd1 + (TEi_est - mean_TE_est)^2
      PCC_sd2 <- PCC_sd2 + (TEi_real - mean_TE_real)^2
      RMSE <- RMSE + (TEi_est / TEi_real - 1)^2
    }
    RBs[l, ] <- c(RB[1] / nr, sum(RB[2:7]) / (nr*Ts), sum(RB[8:13]) / (nr*Ts))
    UBs[l, ] <- c(UB[1] / nr, sum(UB[2:7]) / (nr*Ts), sum(UB[8:13]) / (nr*Ts))
    PCC_cov <- c(PCC_cov[1], sum(PCC_cov[2:7]), sum(PCC_cov[8:13]))
    PCC_sd1 <- c(PCC_sd1[1], sum(PCC_sd1[2:7]), sum(PCC_sd1[8:13]))
    PCC_sd2 <- c(PCC_sd2[1], sum(PCC_sd2[2:7]), sum(PCC_sd2[8:13]))
    PCCs[l, ] <- c((PCC_cov[1] / (sqrt(PCC_sd1[1]) * sqrt(PCC_sd2[1]))),
                   (PCC_cov[2] / (sqrt(PCC_sd1[2]) * sqrt(PCC_sd2[2]))),
                   (PCC_cov[3] / (sqrt(PCC_sd1[3]) * sqrt(PCC_sd2[3]))))
    RMSEs[l, ] <- c(RMSE[1] / nr, sum(RMSE[2:7]) / (nr*Ts), sum(RMSE[8:13]) / (nr*Ts))
    
    l <- l + 1
    print(paste(id, ': ', sim, 'Elapsed time:', Sys.time() - tic, 'seconds'))
  }
  RelativeBias[s,] <- colMeans(RBs, na.rm = TRUE)
  UpperBias[s,] <- colMeans(UBs, na.rm = TRUE)
  PearsonCorrCoef[s,] <- colMeans(PCCs, na.rm = TRUE)
  RelativeMSE[s,] <- colMeans(RMSEs, na.rm = TRUE)
  
  
  s <- s + 1
}

(df <- data.frame(scenario = ID,
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

df %>% 
  select(-starts_with('UpperBias')) %>% 
  left_join(
    readxl::read_excel('abc4sfa/1_simulation_resources/Data/Inputs/Badunenko&Kumbhakar.xlsx') %>%
      select(scenario = ID, ID_inverted),
    by = 'scenario'
  ) %>% relocate(ID_inverted) %>% 
  writexl::write_xlsx(
    sprintf(
      "%s/efficiency_metrics_abc_hf_%s.xlsx",
      FOLDER_OUTPUT, Sys.Date()
    )
  )


# Gibbs: Efficiency --------------------------------------------------------


media <- function(x){
  mean(x, na.rm = TRUE)
}

# Parameters
p <- 1L

# Read data.
files <- list.files(FOLDER_GIBBS, full.names = T)


# Initialize 
n_files <- length(files)
RelativeBias <- matrix(0, n_files, 3)
UpperBias <- matrix(0, n_files, 3)
PearsonCorrCoef <- matrix(0, n_files, 3)
RelativeMSE <- matrix(0, n_files, 3)

s <- 1
for (file in files) {
  data <- readRDS(file)
  sims <- names(data)
  nsim <- length(sims)
  scenario <- stringr::str_extract(file, "s[0-9]+")
  # Extract params.
  params <- INPUT_DATA[[scenario]]$params
  n <- params$n
  t <- params$t
  X <- INPUT_DATA$X
  betas <- INPUT_DATA[[scenario]]$params$beta
  sigmas <- get_real_sigmas(INPUT_DATA, scenario)
  
  # Preallocate
  RBs <- matrix(0, nsim, 3)
  UBs <- matrix(0, nsim, 3)
  RMSEs <- matrix(0, nsim, 3)
  PCCs <- matrix(0, nsim, 3)
  l <- 1L
  for (sim in sims) {
    message("Scenario ", scenario, ". Simulation: ", sim)
    # Read sim related input/outputs.
    y <- INPUT_DATA[[scenario]][[sim]]$y
    gibbs_results <- data[[sim]][['postChain']]
    sigmas_gibbs <- get_gibbs_sigmas(gibbs_results)
    betas_gibbs <- apply(gibbs_results$Thetachain, 1L, mean)
    
    # Compute real technical efficiency.
    TI_real <- ColombiExpectation(n, t, y, X, betas, sigmas, p)
    TE_real <- lapply(TI_real, inv)
    
    # Compute Gibbs' estimates.
    TI_gibbs <- ColombiExpectation(n, t, y, X, betas_gibbs, sigmas_gibbs, p)
    TE_est <- lapply(TI_gibbs, inv)
    apply(as.matrix(do.call(rbind, TE_est)), 2, function(x) any(is.na(x))) %>% any
    # Compute metrics.
    mean_TE_est <- apply(as.matrix(do.call(rbind, TE_est)), 2, media)
    mean_TE_real <- apply(as.matrix(do.call(rbind, TE_real)), 2, media)
    
    RB <- 0
    UB <- 0
    RMSE <- 0
    PCC_cov <- 0
    PCC_sd1 <- 0
    PCC_sd2 <- 0
    inan_real <- which(is.nan(c(unlist(lapply(TE_real, sum)))))
    inan_est <- which(is.nan(c(unlist(lapply(TE_est, sum)))))
    inan <- union(inan_real,inan_est)
    iInf_real <- which(is.infinite(c(unlist(lapply(TE_real, sum)))))
    iInf_est <- which(is.infinite(c(unlist(lapply(TE_est, sum)))))
    iInf <- union(iInf_real, iInf_est)
    ids <- 1:n
    idUs <- setdiff(ids, inan) %>% setdiff(iInf)
    nr <- length(idUs)
    for (i in idUs){
      TEi_est <- TE_est[[i]]
      TEi_real <- TE_real[[i]]
      RB <- RB + (TEi_est / TEi_real - 1)
      (UB <- UB + as.integer(TEi_est > TEi_real))
      PCC_cov <- PCC_cov + (TEi_est - mean_TE_est) * (TEi_real - mean_TE_real)
      PCC_sd1 <- PCC_sd1 + (TEi_est - mean_TE_est)^2
      PCC_sd2 <- PCC_sd2 + (TEi_real - mean_TE_real)^2
      RMSE <- RMSE + (TEi_est / TEi_real - 1)^2
    }
    RBs[l, ] <- c(RB[1] / nr, sum(RB[2:7]) / (nr*t), sum(RB[8:13]) / (nr*t))
    UBs[l, ] <- c(UB[1] / nr, sum(UB[2:7]) / (nr*t), sum(UB[8:13]) / (nr*t))
    PCC_cov <- c(PCC_cov[1], sum(PCC_cov[2:7]), sum(PCC_cov[8:13]))
    PCC_sd1 <- c(PCC_sd1[1], sum(PCC_sd1[2:7]), sum(PCC_sd1[8:13]))
    PCC_sd2 <- c(PCC_sd2[1], sum(PCC_sd2[2:7]), sum(PCC_sd2[8:13]))
    PCCs[l, ] <- c((PCC_cov[1] / (sqrt(PCC_sd1[1]) * sqrt(PCC_sd2[1]))),
                   (PCC_cov[2] / (sqrt(PCC_sd1[2]) * sqrt(PCC_sd2[2]))),
                   (PCC_cov[3] / (sqrt(PCC_sd1[3]) * sqrt(PCC_sd2[3]))))
    RMSEs[l, ] <- c(RMSE[1] / nr, sum(RMSE[2:7]) / (nr*t), sum(RMSE[8:13]) / (nr*t))
    l <- l + 1
  }
  
  RelativeBias[s,] <- colMeans(RBs, na.rm = TRUE)
  UpperBias[s,] <- colMeans(UBs, na.rm = TRUE)
  PearsonCorrCoef[s,] <- colMeans(PCCs, na.rm = TRUE)
  RelativeMSE[s,] <- colMeans(RMSEs, na.rm = TRUE)
  s <- s + 1
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
  left_join(
    readxl::read_excel(file.path(FOLDER_INPUT, 'Badunenko&Kumbhakar.xlsx')) %>%
      select(scenario = ID, ID_inverted),
    by = 'scenario'
  ) %>% select(-scenario) %>% relocate(scenario = ID_inverted) %>% 
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
  sigmas <- INPUT_DATA[[scenario]]$params$sigma
  
  # ABC.
  file_abc <- files_abc[str_extract(files_abc, 's\\d+') == scenario]
  abc_data <- readRDS(file_abc)
  sigmas_est_abc <- sapply(
    abc_data[[scenario]], 
    function(sim) {get_abc_sigmas(sim$sumary$mean_ABC)}
    )
  relative_bias_abc <- apply(sigmas_est_abc / sigmas - 1, 1, mean)
  root_mse_abc <- apply((sigmas_est_abc - sigmas)^2, 1, mean)
  
  # Gibbs
  file_gibbs <- files_gibbs[str_extract(files_gibbs, 's\\d+') == scenario]
  if (!identical(file_gibbs, character(0))) {
    gibbs_data <- readRDS(file_gibbs)
    
    sigmas_est_gibbs <- sapply(
      gibbs_data, function(sim) {get_gibbs_sigmas(sim$postChain)}
    )
    relative_bias_gibbs <- apply(sigmas_est_gibbs / sigmas - 1, 1, mean)
    root_mse_gibbs <- apply((sigmas_est_gibbs - sigmas)^2, 1, mean)
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
  left_join(
    readxl::read_excel('abc4sfa/1_simulation_resources/Data/Inputs/Badunenko&Kumbhakar.xlsx') %>%
      select(scenario = ID, ID_inverted),
    by = 'scenario'
  ) %>% relocate(ID_inverted) %>% 
  writexl::write_xlsx(
    sprintf(
      "%s/sigmas_metrics_hf_%s.xlsx",
      FOLDER_OUTPUT, Sys.Date()
    )
  )

