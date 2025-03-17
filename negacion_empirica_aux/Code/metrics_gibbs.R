#setwd('After/ABC4SFAv2')
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
