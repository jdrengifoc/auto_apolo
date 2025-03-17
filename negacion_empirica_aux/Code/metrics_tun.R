# Count NAs -----------------------------------------------------
setwd("/Users/juanrengifo101/OneDrive - Universidad EAFIT/_Publicaciones/ABC4SFA")

library(dplyr)
root_folder <- 'After/ABC4SFAv2/Data/Outputs'

exps <- paste0('exp', 0:4)
#sims <- paste0('sim', 1:100)

df <- NULL
for (exp in exps) {
  if (exp == 'exp0') {
    exp_files <- file.path(root_folder, 'MH') %>% 
      list.files(full.names = T)
    tuns <- rep(1, length(exp_files))
  } else {
    exp_files <- file.path(root_folder, sprintf('MH_tun_%s', exp)) %>% 
      list.files(full.names = T)
    tuns <- sub(".*tun([0-9]+\\.[0-9]+).*", "\\1", exp_files) %>% as.numeric
  }
  for (exp_file in exp_files) {
    data <- readRDS(exp_file)
    scenario <- names(data)
    sims <- names(data[[scenario]])
    sims <- sims[grepl('sim', sims)]
    tun <- tuns[exp_files == exp_file]
    for (sim in sims) {
      aux <- (data[[scenario]][[sim]][['EF_reales']] + 
                data[[scenario]][[sim]][['EF_abc']]) %>% 
        apply(1, function(x) all(!is.na(x))) %>% unname
      df0 <- tibble(trai = exp, tun, scenario, sim, id = 1L:50L, not_nas = aux)
      df <- rbind(df, df0)
    }
  }
}

df %>% 
  #filter(trai != 'exp4') %>% 
  #filter(trai %in% paste0('exp', 2)) %>% 
  group_by(scenario, sim, id) %>% summarise(not_nas = max(not_nas)) %>% 
  group_by(scenario, sim) %>% summarise(not_nas = mean(not_nas)) %>% 
  group_by(scenario) %>% summarise(not_nas = mean(not_nas)) %>% 
  arrange(not_nas)

# Join data ---------------------------------------------------------------

root_folder <- 'After/ABC4SFAv2/Data/Outputs'


base_folder <- file.path(root_folder, 'MH_final')
new_folder <- file.path(root_folder, 'MH_final_este_si')
scenarios <- paste0('s', 1:16)
exps <- paste0('exp', 1:4)

# Para cada `scenario` y `sim`, determinar un escenario base y a partir de este
# buscar cuales filas tienen NAs por la estimaciones de la eficiencia real y 
# por ABC.
for (scenario in scenarios) {
  # Leer datos base.
  file <- file.path(base_folder, sprintf('MH_%s.RData', scenario))
  base_data <- readRDS(file)[[scenario]]
  sims <- names(base_data)
  sims <- sims[grepl('sim', sims)]
  for (sim in sims) {
    ids_real_to_search <- (base_data[[sim]][['EF_reales']]) %>% 
      apply(1, function(x) any(is.na(x))) %>% unname %>% which
    ids_abc_to_search <- (base_data[[sim]][['EF_abc']]) %>% 
      apply(1, function(x) any(is.na(x))) %>% unname %>% which
    
    ids_to_search <- union(ids_real_to_search, ids_abc_to_search)
    # Si no hay nada para buscar siga.
    if ( length(ids_to_search) == 0) {
      next
    }
    # Si hay algún NA, ver que archivos pueden tener la solución.
    new_files <- list.files(root_folder, recursive = T, full.names = T,
                            pattern = sprintf('MH%s_', scenario))
    # Por cada archivo explorar de manera independiente las estimaciones reales 
    for (new_file in new_files) {
      # Si no hay nada para buscar, entonces siga.
      if ( length(ids_real_to_search) == 0) {
        break
      }
      
      new_data <- readRDS(new_file)[[scenario]][[sim]][['EF_reales']]
      ids_new_without_na <- (new_data %>% apply(1, function(x) all(!is.na(x))) %>%
          which %>% unname)
      # Determinar que filas se van a reemplazar y con qué se va a reemplazar.
      if (length(ids_new_without_na) > length(ids_real_to_search)) {
        ids_selected_to_replace <- sample(ids_new_without_na, length(ids_real_to_search),
                                          replace = F)
        ids_tobe_replaced <- ids_real_to_search
      } else {
        ids_selected_to_replace <- ids_new_without_na
        ids_tobe_replaced <- ids_real_to_search[1:length(ids_selected_to_replace)]
      }
      # Actualizar la base y los ids a buscar.
      base_data[[sim]][['EF_reales']][ids_tobe_replaced,] <- new_data[ids_selected_to_replace,]
      ids_real_to_search <- setdiff(ids_real_to_search, ids_tobe_replaced)
    }
    # Por cada archivo explorar de manera independiente las estimaciones abc 
    for (new_file in new_files) {
      # Si no hay nada para buscar, entonces siga.
      if ( length(ids_abc_to_search) == 0) {
        break
      }
      
      new_data <- readRDS(new_file)[[scenario]][[sim]][['EF_abc']]
      ids_new_without_na <- (new_data %>% apply(1, function(x) all(!is.na(x))) %>%
                               which %>% unname)
      # Determinar que filas se van a reemplazar y con qué se va a reemplazar.
      if (length(ids_new_without_na) > length(ids_abc_to_search)) {
        ids_selected_to_replace <- sample(ids_new_without_na, length(ids_abc_to_search),
                                          replace = F)
        ids_tobe_replaced <- ids_abc_to_search
      } else {
        ids_selected_to_replace <- ids_new_without_na
        ids_tobe_replaced <- ids_abc_to_search[1:length(ids_selected_to_replace)]
      }
      # Actualizar la base y los ids a buscar.
      base_data[[sim]][['EF_abc']][ids_tobe_replaced,] <- new_data[ids_selected_to_replace,]
      ids_real_to_search <- setdiff(ids_abc_to_search, ids_tobe_replaced)
    }
  }
    new_file <- file.path(new_folder, sprintf('MH_%s.RData', scenario))
    new_data <- list()
    new_data[[scenario]] <- base_data
    saveRDS(new_data, new_file)
}



# Validar join ------------------------------------------------------------

setwd("/Users/juanrengifo101/OneDrive - Universidad EAFIT/_Publicaciones/ABC4SFA")


df <- NULL
exp_files <- file.path(new_folder) %>% list.files(full.names = T)

for (exp_file in exp_files) {
  data <- readRDS(exp_file)
  scenario <- names(data)
  sims <- names(data[[scenario]])
  sims <- sims[grepl('sim', sims)]
  for (sim in sims) {
    aux <- (data[[scenario]][[sim]][['EF_reales']] + 
              data[[scenario]][[sim]][['EF_abc']]) %>% 
      apply(1, function(x) all(!is.na(x))) %>% unname
    df0 <- tibble(scenario, sim, id = 1L:50L, not_nas = aux)
    df <- rbind(df, df0)
  }
}

df %>% 
  group_by(scenario, sim, id) %>% summarise(not_nas = max(not_nas)) %>% 
  group_by(scenario, sim) %>% summarise(not_nas = mean(not_nas)) %>% 
  group_by(scenario) %>% summarise(not_nas = mean(not_nas)) %>% 
  arrange(not_nas)
# Compute metrics ---------------------------------------------------------


setwd('After/ABC4SFAv2')
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')
library(mvtnorm)

media <- function(x){
  mean(x, na.rm = TRUE)
}
# Parameters
p <- 1L

# Read data.
folder_data <- 'Data/Outputs/Gibs_experiment'
files <- list.files(folder_data, full.names = T)
inputData <- readRDS("Data/Inputs/BK_simData.RData")

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
  params <- inputData[[scenario]]$params
  n <- params$n
  t <- params$t
  X <- inputData$X
  betas <- inputData[[scenario]]$params$beta
  sigmas <- rev(inputData[[scenario]]$params$sigma)
  
  # Preallocate
  RBs <- matrix(0, nsim, 3)
  UBs <- matrix(0, nsim, 3)
  RMSEs <- matrix(0, nsim, 3)
  PCCs <- matrix(0, nsim, 3)
  l <- 1L
  for (sim in sims) {
    # Read sim related input/outputs.
    y <- inputData[[scenario]][[sim]]$y
    gibbs_results <- data[[sim]][['postChain']]
    sigmas_gibbs <- c(
      mean(gibbs_results$sigmaeta2chain), 
      mean(gibbs_results$sigmau2chain),
      mean(gibbs_results$sigmaalpha2chain),
      mean(gibbs_results$sigmaepsilon2chain)
    )
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
(df <- data.frame(scenario = scenarios,
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
writexl::write_xlsx(df, 'Data/Outputs/gibbs_metrics.xlsx')
