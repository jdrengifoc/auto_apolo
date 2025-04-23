rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))
library(mvtnorm)
library(stringr)

# Results halfnormal-normal GTRE ---------------------------------------------
FOLDER <- 'abc4sfa/1_simulation_resources'
ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_Ystats.RData'))
inputData <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_simData.RData'))
files <- list.files(
  file.path(FOLDER, 'Data/Outputs'), full.names = T, pattern = '2025-04-06'
  )
output_file <- sprintf(
  "abc4sfa/4_error_measurements/data/output/simulation_hn_misspecification_exp_%s.xlsx",
  Sys.Date()
)
#' It is also a common practice in ABC to perform a regression adjustment after 
#' retained draws (Beaumont et al., 2002; Leuenberger and Wegmann, 2010;
#' Sisson et al., 2018)
# Preallocate.
# Without adjustment.
ID <- str_extract(files, 's\\d+')
n_ids <- length(ID)
RelativeBias <- matrix(0, n_ids, 3)
UpperBias <- matrix(0, n_ids, 3)
PearsonCorrCoef <- matrix(0, n_ids, 3)
RelativeMSE <- matrix(0, n_ids, 3)
# With adjustment.
RelativeBiasAdj <- matrix(0, n_ids, 3)
UpperBiasAdj <- matrix(0, n_ids, 3)
PearsonCorrCoefAdj <- matrix(0, n_ids, 3)
RelativeMSEAdj <- matrix(0, n_ids, 3)

tic <- Sys.time()
s <- 1
# For each scenario.
for (id in ID){
  # Read results.
  outputData <- readRDS(files[grep(sprintf("%s_", id), files)])
  
  # Extract real betas and covariates.
  X <- inputData[['X']]
  betas <- inputData[[id]]$params$beta
  # Extract simulations realized.
  sims <- names(outputData[[id]])
  nsim <- length(sims)
  # Fix inverted sigmas in experiments.
  sigmas <- rev(inputData[[id]]$params$sigma)
  Ts <- inputData[[id]]$params$t
  n <- inputData[[id]]$params$n
  theta <- list(n = n, t = Ts)
  aux_idx <- c(1, 2:(Ts+1), 2:(Ts+1))
  RBs <- matrix(0, nsim, 3)
  UBs <- matrix(0, nsim, 3)
  RMSEs <- matrix(0, nsim, 3)
  PCCs <- matrix(0, nsim, 3)
  RBsAdj <- matrix(0, nsim, 3)
  UBsAdj <- matrix(0, nsim, 3)
  RMSEsAdj <- matrix(0, nsim, 3)
  PCCsAdj <- matrix(0, nsim, 3)
  l <- 1
  for (sim in sims){
    y <- inputData[[id]][[sim]]$y
    # est_ABCparams <- outputData[[id]][[sim]]$ABCpostChain[1,]
    est_ABCparams <- apply(outputData[[id]][[sim]]$ABCpostChain, 2, mean)
    RegY <- lm(y ~ .,  data = data.frame(y, X[,2:3]))
    SumStatY <- getSampleSumStats(theta = theta, y = y, X = X)
    est_betas <- c(est_ABCparams[1],
                   RegY$coefficients[-1])
    est_sigmas <- est_ABCparams[2:5]
    PostDraws <- outputData[[id]][[sim]]$ABCpostChain
    K <- dim(PostDraws)[1]
    SumStatPar <- matrix(0, K, 5)
    for (k in 1:K){
      ResSim <- simRes(theta = theta, prior = PostDraws[k, 1:5])#, dist = "halfnormal")
      SumStatPar[k, ] <- getSumStats(theta = theta, resS = ResSim)
    }
    PostDrawsAdj <- matrix(0, K, 5)
    for (p in 1:5){
      RegABCp <- lm(log(PostDraws[, p])~SumStatPar)
      ABCp <- exp(c(c(1,SumStatY)%*%RegABCp$coefficients)+RegABCp$residuals)
      PostDrawsAdj[, p] <- ABCp
    }
    ABCAdj <- colMeans(PostDrawsAdj)
    est_betasAdj <- c(ABCAdj[1],
                      RegY$coefficients[-1]) 
    est_sigmasAdj <- ABCAdj[2:5] 
    
    TI_real <- ColombiExpectation(n, Ts, y, X, betas, sigmas, p = 1)
    TE_real <- lapply(TI_real, inv)
    TI_est <- ColombiExpectation(n, Ts, y, X, est_betas, est_sigmas, p = 1)
    TE_est <- lapply(TI_est, inv)
    
    mean_TI_est <- apply(as.matrix(do.call(rbind, TI_est)), 2, mean, na.rm = T)
    mean_TI_real <- apply(as.matrix(do.call(rbind, TI_real)), 2, mean, na.rm = T)
    
    mean_TE_est <- apply(as.matrix(do.call(rbind, TE_est)), 2, mean, na.rm = T)
    mean_TE_real <- apply(as.matrix(do.call(rbind, TE_real)), 2, mean, na.rm = T)
    
    TI_estAdj <- ColombiExpectation(n, Ts, y, X, est_betasAdj, est_sigmasAdj, p = 1)
    TE_estAdj <- lapply(TI_estAdj, inv)
    
    mean_TI_estAdj <- apply(as.matrix(do.call(rbind, TI_estAdj)), 2, mean, na.rm = T)
    
    mean_TE_estAdj <- apply(as.matrix(do.call(rbind, TE_estAdj)), 2, mean, na.rm = T)
    
    
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
    
    RBAdj <- 0
    UBAdj <- 0
    RMSEAdj <- 0
    PCC_covAdj <- 0
    PCC_sd1Adj <- 0
    PCC_sd2Adj <- 0
    inan_estAdj <- which(is.nan(c(unlist(lapply(TE_estAdj, sum)))))
    inan <- union(inan_real,inan_estAdj)
    ids <- 1:n
    idUs <- setdiff(ids, inan) 
    nr <- length(idUs)
    for (i in idUs){
      TEi_est <- TE_estAdj[[i]]
      TEi_real <- TE_real[[i]]
      RBAdj <- RBAdj + (TEi_est / TEi_real - 1)
      (UBAdj <- UBAdj + as.integer(TEi_est > TEi_real))
      PCC_covAdj <- PCC_covAdj + (TEi_est - mean_TE_estAdj) * (TEi_real - mean_TE_real)
      PCC_sd1Adj <- PCC_sd1Adj + (TEi_est - mean_TE_estAdj)^2
      PCC_sd2Adj <- PCC_sd2Adj + (TEi_real - mean_TE_real)^2
      RMSEAdj <- RMSEAdj + (TEi_est / TEi_real - 1)^2
    }
    RBsAdj[l, ] <- c(RBAdj[1] / nr, sum(RBAdj[2:7]) / (nr*Ts), sum(RBAdj[8:13]) / (nr*Ts))
    UBsAdj[l, ] <- c(UBAdj[1] / nr, sum(UBAdj[2:7]) / (nr*Ts), sum(UBAdj[8:13]) / (nr*Ts))
    PCC_covAdj <- c(PCC_covAdj[1], sum(PCC_covAdj[2:7]), sum(PCC_covAdj[8:13]))
    PCC_sd1Adj <- c(PCC_sd1Adj[1], sum(PCC_sd1Adj[2:7]), sum(PCC_sd1Adj[8:13]))
    PCC_sd2Adj <- c(PCC_sd2Adj[1], sum(PCC_sd2Adj[2:7]), sum(PCC_sd2Adj[8:13]))
    PCCsAdj[l, ] <- c((PCC_covAdj[1] / (sqrt(PCC_sd1Adj[1]) * sqrt(PCC_sd2Adj[1]))),
                      (PCC_covAdj[2] / (sqrt(PCC_sd1Adj[2]) * sqrt(PCC_sd2Adj[2]))),
                      (PCC_covAdj[3] / (sqrt(PCC_sd1Adj[3]) * sqrt(PCC_sd2Adj[3]))))
    RMSEsAdj[l, ] <- c(RMSEAdj[1] / nr, sum(RMSEAdj[2:7]) / (nr*Ts), sum(RMSEAdj[8:13]) / (nr*Ts))
    
    l <- l + 1
    print(paste(id, ': ', sim, 'Elapsed time:', Sys.time() - tic, 'seconds'))
  }
  RelativeBias[s,] <- colMeans(RBs, na.rm = TRUE)
  UpperBias[s,] <- colMeans(UBs, na.rm = TRUE)
  PearsonCorrCoef[s,] <- colMeans(PCCs, na.rm = TRUE)
  RelativeMSE[s,] <- colMeans(RMSEs, na.rm = TRUE)
  
  RelativeBiasAdj[s,] <- colMeans(RBsAdj, na.rm = TRUE)
  UpperBiasAdj[s,] <- colMeans(UBsAdj, na.rm = TRUE)
  PearsonCorrCoefAdj[s,] <- colMeans(PCCsAdj, na.rm = TRUE)
  RelativeMSEAdj[s,] <- colMeans(RMSEsAdj, na.rm = TRUE)
  
  
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

df %>% mutate(across(-scenario, ~ round(100 * ., 2)))
writexl::write_xlsx(df, output_file)

# No sirve.
# (dfAdj <- data.frame(scenario = ID,
#                      RelativeBiasPermAdj = RelativeBiasAdj[,1],
#                      RelativeBiasTransAdj = RelativeBiasAdj[,2],
#                      RelativeBiasAllAdj = RelativeBiasAdj[,3],
#                      UpperBiasPermAdj = UpperBiasAdj[,1],
#                      UpperBiasTransAdj = UpperBiasAdj[,2],
#                      UpperBiasAllAdj = UpperBiasAdj[,3],
#                      PearsonCorrCoefPermAdj = PearsonCorrCoefAdj[,1],
#                      PearsonCorrCoefTransAdj = PearsonCorrCoefAdj[,2],
#                      PearsonCorrCoefAllAdj = PearsonCorrCoefAdj[,3],
#                      RelativeMSEPermAdj = RelativeMSEAdj[,1],
#                      RelativeMSETransAdj = RelativeMSEAdj[,2],
#                      RelativeMSEAllAdj = RelativeMSEAdj[,3]
# ))
# writexl::write_xlsx(dfAdj, 'Results_sims_final_real_meanNewV1Adj.xlsx')

