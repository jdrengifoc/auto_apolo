# source('requirements.R')
rm(list = ls())
set.seed(010101)
library('mvtnorm')

source("C:/ANDRES/Portatil/ABC/ABC_SFA/JuanDavid/ABC4SFA_Hassan/Code/Fixed/functions.R")

# Function ----------------------------------------------------------------
#'  Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
#'  A = −p × [1Ti ITi] is a matrix of dimension Ti × (Ti + 1),  where 1Ti is the
#' column vector of length Ti, and ITi is the identity matrix of dimension T.
get_A <- function(T_i, p){ 
  A <- -p * cbind(rep(1,T_i), diag(T_i))
  return(A)
}

#'  Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
#'  Composed error term of the joint density function.
get_epsiloni <- function(y, X, beta){
  epsiloni <- y - X%*%beta
  return(epsiloni)
}

#' Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
#' Diagonal? matrix of dimension T_i + 1, 
#' where its diagonal is [sigma_u0, sigma_u_1, ..., sigma_u_T_i].
get_V <- function(sigma_u0, sigma_u, T_i){
  V <- diag(c(sigma_u0^2, rep(sigma_u^2, T_i)))
  return(V)
}

#' Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
get_Sigma <- function(sigma_v, sigma_v0, T_i){
  Sigma <- sigma_v^2 * diag(T_i) + sigma_v0^2 * rep(1, T_i)%*%t(rep(1,T_i))
  return(Sigma)
}

#' Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
get_Lambda <- function(V, A, Sigma){
  Lambda <- solve(solve(V) + t(A)%*%solve(Sigma)%*%A)
  return(Lambda)
}

#' Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
get_R <- function(Lambda, A, Sigma){
  R <- Lambda%*%t(A)%*%solve(Sigma)
  return(R)
}

#' Probability that a q-variate normal variable of expected value \pmb\mu
#' and variance \Sigma belongs to the positive orthant.
probPosOrthant <- function(Mean, Sigma) {
  Mean <- as.vector(Mean)
  Phi_bar <- pmvnorm(lower = 0, upper = Inf, mean = Mean, sigma = Sigma)
  return(Phi_bar)
}

#' Estimation of technical efficiency of Colombi et al. (2014).
ColombiExpectation <- function(n, Ts, y, X, betas, sigmas, p = 1){
  # Change balanced notation to unbalanced panel data notation.
  if (length(Ts) == 1 & length(Ts) < n) {
    Ts <- rep(Ts, n)
  }
  # Unwraps variances.
  sigma_v <- sigmas[1]
  sigma_v0 <- sigmas[2]
  sigma_u <- sigmas[3]
  sigma_u0 <- sigmas[4]
  # Preallocate.
  ExpectedTIs <- vector(mode = "list", length = length(ts))
  epsilonis <- get_epsiloni(y, X, betas)
  
  # For each individual.
  idx <- 0
  for (i in 1:n){
    # Get indexes for individual i.
    T_i <- Ts[i]
    I <- (1:T_i) + idx
    
    A <- get_A(T_i, p)
    epsiloni <- epsilonis[I]
    V <- get_V(sigma_u0, sigma_u, T_i)
    Sigma <- get_Sigma(sigma_v, sigma_v0, T_i)
    Lambda <- get_Lambda(V, A, Sigma)
    R <- get_R(Lambda, A, Sigma)
    
    tBold <- diag(T_i+1)
    ExpectedTI_i <- rep(0, T_i+1)
    for (t_i in 1:(T_i+1)){
      ExpectedTI_i[t_i] <- probPosOrthant(R%*%epsiloni + Lambda%*%tBold[t_i,], Lambda) /
        probPosOrthant(R%*%epsiloni, Lambda) *
        exp(t(tBold[t_i,])%*%R%*%epsiloni + 0.5*t(tBold[t_i,])%*%Lambda%*%tBold[t_i,])
    }
    ExpectedTIs[[i]] <- ExpectedTI_i
    # Update individual index.
    idx <- T_i + idx
  }
  return(ExpectedTIs)
}

inv <- function(x){
  xinv <- 1/x
  overallEf <- xinv[1]*xinv[2:length(x)]
  return(c(xinv,overallEf))
}

media <- function(x){
  xmedia <- mean(x, na.rm = TRUE)
  return(xmedia)
}

#' A = −p × [1Ti ITi] is a matrix of dimension Ti × (Ti + 1), 
#' where 1Ti is the column vector of length Ti, and ITi is the
#' identity matrix of dimension Ti.
#' p is 1 if c production frontier (TI) p is -1 if cost function (TE).
#' @param ts vector de periodos por individuo,


# Results -----------------------------------------------------------------

ID <- paste0('s', 1:16)
RelativeBias <- matrix(0, 16, 3)
UpperBias <- matrix(0, 16, 3)
PearsonCorrCoef <- matrix(0, 16, 3)
RelativeMSE <- matrix(0, 16, 3)

RelativeBiasAdj <- matrix(0, 16, 3)
UpperBiasAdj <- matrix(0, 16, 3)
PearsonCorrCoefAdj <- matrix(0, 16, 3)
RelativeMSEAdj <- matrix(0, 16, 3)

tic <- Sys.time()
s <- 1
for (id in ID){
  outputData <- readRDS(sprintf("BK_%s.RData", id))
  inputData <- readRDS("BK_simData.RData")
  
  X <- inputData[['X']]
  
  sims <- names(outputData[[1]])
  nsim <- length(sims)
  betas <- inputData[[id]]$params$beta
  # Inverted sigmas in experiments.
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
      ResSim <- simRes(theta = theta, prior = PostDraws[k, 1:5], dist = "halfnormal")
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
    
    mean_TI_est <- apply(as.matrix(do.call(rbind, TI_est)), 2, media)
    mean_TI_real <- apply(as.matrix(do.call(rbind, TI_real)), 2, media)
    
    mean_TE_est <- apply(as.matrix(do.call(rbind, TE_est)), 2, media)
    mean_TE_real <- apply(as.matrix(do.call(rbind, TE_real)), 2, media)
    
    TI_estAdj <- ColombiExpectation(n, Ts, y, X, est_betasAdj, est_sigmasAdj, p = 1)
    TE_estAdj <- lapply(TI_estAdj, inv)
    
    mean_TI_estAdj <- apply(as.matrix(do.call(rbind, TI_estAdj)), 2, media)
    
    mean_TE_estAdj <- apply(as.matrix(do.call(rbind, TE_estAdj)), 2, media)

    
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
# writexl::write_xlsx(df, 'Results_sims_final_real_meanNewV1.xlsx')

(dfAdj <- data.frame(scenario = ID,
                  RelativeBiasPermAdj = RelativeBiasAdj[,1],
                  RelativeBiasTransAdj = RelativeBiasAdj[,2],
                  RelativeBiasAllAdj = RelativeBiasAdj[,3],
                  UpperBiasPermAdj = UpperBiasAdj[,1],
                  UpperBiasTransAdj = UpperBiasAdj[,2],
                  UpperBiasAllAdj = UpperBiasAdj[,3],
                  PearsonCorrCoefPermAdj = PearsonCorrCoefAdj[,1],
                  PearsonCorrCoefTransAdj = PearsonCorrCoefAdj[,2],
                  PearsonCorrCoefAllAdj = PearsonCorrCoefAdj[,3],
                  RelativeMSEPermAdj = RelativeMSEAdj[,1],
                  RelativeMSETransAdj = RelativeMSEAdj[,2],
                  RelativeMSEAllAdj = RelativeMSEAdj[,3]
))
writexl::write_xlsx(dfAdj, 'Results_sims_final_real_meanNewV1Adj.xlsx')

# pp1 <- as.matrix(do.call(rbind, TI_est))[,1]
# pp2 <- as.matrix(do.call(rbind, TI_real))[,1]
# plot(pp1, pp2, xlim = c(1, 1.5), ylim = c(1, 1.5))
# abline(coef = c(0,1))
