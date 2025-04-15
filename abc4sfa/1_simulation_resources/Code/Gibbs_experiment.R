# Mejor s1 y peor s13-14 invertido, para halfnormal con código Miguel 
distribution <- c('hfn', 'exp')[1L]
scenario <- c('s1', 's13')[1]

output_file <- sprintf(
  'Data/Outputs/Gibs_experiment/Gibs_experiment_%s_%s.RData', 
  distribution, scenario)

# Load libraries ----------------------------------------------------------
# libs <- c('invgamma', 'MASS', 'spdep', 'numDeriv',
#           'tmvtnorm', 'progress', 'coda', 'dplyr',
#           'future', 'furrr')
# for (lib in libs) {
#   install.packages(lib)
# }

library(LearnBayes)
library(truncnorm)
library(invgamma)
library(MASS)
library(spdep)
library(numDeriv)
library(tmvtnorm)
#library(tidyverse)
library(progress)
library(coda)
library(dplyr)
library(future)
library(furrr)


# Secondary functions -----------------------------------------------------

create_skewed <- function(Ri, sigmaepsilon2, sigmadelta2, lambdadelta) {
  function(x) {
    A <- exp(- sum((Ri-x)**2) / (2 * sigmaepsilon2)
             - x^2 / sigmadelta2) * pnorm(lambdadelta * x / sqrt(sigmadelta2))
    return(drop(A))
  }
}

create_log_skew <- function(Ri, sigmaepsilon2, sigmadelta2, lambdadelta){
  function(x){
    p <-  -sum((Ri-x)**2)/(2*sigmaepsilon2)-x**2/sigmadelta2 +
      pnorm(lambdadelta*x/sqrt(sigmadelta2), log.p = T)
    
    return(p)
  }
}

etasampling <- function(Di, Sigma, t) {
  Di <- matrix(Di, length(Di), 1)  
  phi2 <- sigmaeta2 * solve(1 + sigmaeta2 * matrix(1, 1, t) %*%
                              solve(Sigma) %*% matrix(1, t, 1))
  m <- phi2 * matrix(1, 1, t) %*% solve(Sigma) %*% Di
  eta <- rtruncnorm(1, 0, Inf, m, sd = sqrt(phi2))
  return(eta)
}

usampling <- function(u_bar, sigmau2, sigmaepsilon2) {
  V <- (sigmau2*sigmaepsilon2)/(sigmaepsilon2 + sigmau2)
  u <- rtruncnorm(1, 0, Inf, u_bar, sqrt(V))
  return(u)
}

# deltasampling1 <- function(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta) {
#   Ri <- matrix(Ri, length(Ri), 1)
#   skewed <- create_skewed(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
#   log_skewed <- function(x) log(skewed(x))
#   
#   Mode <- grid[which.max(sapply(grid, skewed))]
#   s2 <- -1 / hessian(skewed, x = Mode)
#   delta <- rnorm(1, mean = Mode, sd = sqrt(s2))
#   
#   prob <- exp(log_skewed(delta) - log_skewed(Mode)) * exp(-(delta - Mode)^2 / (2 * s2))
#   return(list(delta, prob))
# }

deltasampling1 <- function(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta) {
  Ri <- matrix(Ri, length(Ri), 1)
  skewed <- create_skewed(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
  log_skewed <- create_log_skew(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
  log_skewed_min <- function(x) -log_skewed(x)
  
  Mode <- optim(-0.5, log_skewed_min, method = 'Brent', 
                upper = 15, lower = -15)$par
  
  s2 <- - 1/hessian(skewed, x = Mode)
  delta <- rnorm(1, mean = Mode, sd = sqrt(s2))
  
  prob <- exp(log_skewed(delta) - log_skewed(Mode)) * exp(-(delta - Mode)^2 / (2 * s2))
  return(list(delta, prob))
}

get_Xi <- function(i, X, N, t){
  X[seq(i, N*t, by = N), ]
}



# Main function -----------------------------------------------------------

# Inputs:
## data: lista con matriz de diseño "X", indices i y t "loc", y vector dependiente "y"
#
## model: cost function o production function
get_posterior <- function(data, model = "cost", burnin=10000, 
                          n_samples=25000, thinning=10) {
  start_time <- Sys.time()
  # Prepare data matrices
  X <- data$X %>% as.matrix()
  y <- data$y
  N <- data$loc %>% dplyr::select(id) %>% unique() %>% nrow() %>% as.numeric()
  t <- data$loc %>% dplyr::select(time) %>% unique() %>% nrow() %>% as.numeric()
  Zt <- X
  # beta0: vector de coeficientes estimados por MCO
  beta0 <- solve(t(X)%*%X)%*%t(X)%*%y %>% as.vector
  # beta0 <- matrix(beta0[-1], length(beta0)-1, 1)
  # Ztm_intercept <- postmatrix(X[,-1]%*%beta0, N, t, byrow = FALSE)
  k1 <- length(beta0)
  k2 <- N
  k3 <- N * t
  
  # Initialize parameters
  sigmaepsilon2 <- 0.15^2
  sigmaeta2 <- 0.28^2
  sigmaalpha2 <- 0.12^2
  sigmau2 <- 0.43^2
  
  p <- ifelse(model == 'cost', 1, -1)
  Sigma <- tcrossprod(matrix(1, t, 1), matrix(1, t, 1)) * as.numeric(sigmaalpha2) +
    as.numeric(sigmaepsilon2) * diag(t)
  # B01 <- 1000
  # b01 <- 0
  B01 <- solve(diag(k1) * 10^10)
  b01 <- matrix(0, k1, 1)
  etai <- abs(rnorm(N, 0, sqrt(sigmaeta2)))
  eta <- matrix(etai, N, t, byrow = FALSE)
  alphai <- matrix(rnorm(N, 0, sqrt(sigmaalpha2)), N, 1)
  alpha <- matrix(alphai, N, t, byrow = FALSE)
  deltai <- etai + p*alphai
  delta <- eta + p*alpha
  um <- matrix(abs(rnorm(N*t, 0, sqrt(sigmau2))), N, t)
  u <- as.vector(t(um))
  
  # Tener mucho cuidado con como esta organizada la data en cuanto a t e id
  ym <- matrix(y, N, t, byrow = FALSE)
  nbar <- 1
  nsigma <- N * t + nbar
  nsigma1 <- N + nbar
  qbar <- 10^(-4)
  
  
  # Storage for chains
  Thetachain <- matrix(0, k1, n_samples)
  deltachain <- matrix(0, k2, n_samples)
  etachain <- matrix(0, k2, n_samples)
  alphachain <- matrix(0, k2, n_samples)
  uchain <- matrix(0, k3, n_samples)
  epsilonchain <- matrix(0, k3, n_samples)
  sigmau2chain <- matrix(0, 1, n_samples)
  sigmaalpha2chain <- matrix(0, 1, n_samples)
  sigmaeta2chain <- matrix(0, 1, n_samples)
  sigmaepsilon2chain <- matrix(0, 1, n_samples)
  sigmadelta2chain <- matrix(0, 1, n_samples)
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Sampling [:bar] :percent in :elapsed",
    total = n_samples,
    width = 60
  )
  
  Xi <- lapply(1:N, function(i) get_Xi(i, X, N, t))
  
  # Main Gibbs sampling loop
  plan(multisession)
  
  results <- future_map(1:n_samples, .progress = TRUE, .options = furrr_options(seed = TRUE), function(i) {
    tryCatch({
      
      # Y <- ym - Ztm_intercept - delta - um
      # B2 <- 1/(matrix(1, 1, t)%*%solve(Sigma)%*%matrix(1,t,1) + B01)
      # b3 <- B2*(sum((sapply(1:N, function(i) matrix(1, 1, t)%*%solve(Sigma)%*%Y[i,]))) +
      #             b01/B01)
      # intercept <- rnorm(1, b3, sqrt(B2))
      # theta <- rbind(intercept, beta0)
      
      Ym <- ym - p*delta - p*um
      XSX <- lapply(Xi, function(X){
        t(X)%*%solve(Sigma)%*%X
        })
      
      sumXSX <- Reduce('+', XSX)
      B2 <- solve(sumXSX + B01)
      
      XSY <- lapply(1:N, function(i){
        t(Xi[[i]])%*%solve(Sigma)%*%Ym[i,]
      })
      
      sumXSY <- Reduce('+', XSY)
      b3 <- B2%*%(sumXSY + B01%*%b01)
      theta <- t(rmnorm(1, b3, B2))
      
      Ztm <- matrix(Zt %*% theta, N, t, byrow = FALSE)
      sigmau2 <- (qbar + crossprod(u, u)) / rchisq(1, nsigma)
      R <- ym - (Ztm + p*um)
      
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      lambdadelta <- sqrt(sigmaeta2 / sigmaalpha2)
      
      candp <- t(do.call(cbind, lapply(1:N, function(i) {
        Ri <- as.matrix(R[i,])
        deltasampling1(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta)
      })))
      
      deltai <- matrix(sapply(1:N, function(i) {
        ifelse(rbinom(1, 1, min(1, unlist(candp[i, 2]))) == 1, 
               as.numeric(candp[i, 1]), as.numeric(deltai[i]))}), N, 1)
      
      epsilonhatm <- ym - (Ztm + delta + p*um)
      epsilonhat <- as.vector(epsilonhatm)
      
      sigmaepsilon2 <- 1 / rchisq(1, nsigma) * (qbar + crossprod(epsilonhat, epsilonhat))
      utilde <- as.vector(ym - (Ztm + delta))
      u <- sapply(1:(N*t), function(i) usampling(utilde[i], sigmau2, sigmaepsilon2))
      um <- matrix(u, N, t, byrow = FALSE)
      
      D <- ym - (Ztm + p*um)
      etai <- matrix(sapply(1:N, function(i) etasampling(D[i, ], Sigma, t)), N, 1)
      sigmaeta2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(etai, etai)))
      eta <- matrix(etai, N, t, byrow = FALSE)
      
      delta <- matrix(deltai, N, t, byrow = FALSE)
      
      alphai <- as.numeric(deltai) - p*as.numeric(etai)
      alpha <- matrix(alphai, N, t, byrow = FALSE)
      sigmaalpha2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(alphai, alphai)))
      
      Sigma <- tcrossprod(matrix(1, t, 1), matrix(1, t, 1)) * as.numeric(sigmaalpha2) + as.numeric(sigmaepsilon2) * diag(t)
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      
      # Check for NaNs
      if (any(is.nan(c(theta, deltai, etai, alphai, u, sigmaepsilon2, sigmaeta2, sigmaalpha2)))) {
        return(NULL)  # Skip to the next iteration
      } else {
        # Store samples
        list(
          theta = theta,
          deltai = deltai,
          etai = etai,
          alphai = alphai,
          u = u,
          epsilon = epsilonhat,
          sigmau2 = sigmau2,
          sigmaepsilon2 = sigmaepsilon2,
          sigmaeta2 = sigmaeta2,
          sigmaalpha2 = sigmaalpha2,
          sigmadelta2 = sigmadelta2
        )
      }
      
    }, error = function(e) {
      cat("Error at iteration:", i, "\n")
      cat("Error message:", e$message, "\n")
      NULL  # Skip to the next iteration
    }, warning = function(w) {
      cat("Warning at iteration:", i, "\n")
      cat("Warning message:", w$message, "\n")
      NULL  # Skip to the next iteration
    })
  })
  
  # Process results
  valid_results <- results[!sapply(results, is.null)]
  for (i in 1:length(valid_results)) {
    res <- valid_results[[i]]
    Thetachain[, i] <- res$theta
    deltachain[, i] <- res$deltai
    etachain[, i] <- res$etai
    alphachain[, i] <- res$alphai
    uchain[, i] <- res$u
    epsilonchain[, i] <- res$epsilon
    sigmau2chain[, i] <- res$sigmau2
    sigmaepsilon2chain[, i] <- res$sigmaepsilon2
    sigmaeta2chain[, i] <- res$sigmaeta2
    sigmaalpha2chain[, i] <- res$sigmaalpha2
    sigmadelta2chain[, i] <- res$sigmadelta2
  }
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  thinned_indices <- seq(1, length(valid_results), by = thinning)
  
  list(
    Thetachain = Thetachain[, thinned_indices],
    deltachain = deltachain[, thinned_indices],
    etachain = etachain[, thinned_indices],
    alphachain = alphachain[, thinned_indices],
    uchain = uchain[, thinned_indices],
    epsilonchain = epsilonchain[, thinned_indices],
    sigmau2chain = sigmau2chain[, thinned_indices],
    sigmaepsilon2chain = sigmaepsilon2chain[, thinned_indices],
    sigmaeta2chain = sigmaeta2chain[, thinned_indices],
    sigmaalpha2chain = sigmaalpha2chain[, thinned_indices],
    sigmadelta2chain = sigmadelta2chain[, thinned_indices],
    elapsed_time = elapsed_time
  )
}


# run ---------------------------------------------------------------------
# Read data.
if (distribution == 'hfn') {
  ABC_sample <- readRDS('Data/Inputs/BK_simData.RData')
  ABC_inputs <- readRDS('Data/Inputs/BK_Ystats.RData')[[scenario]][['params']]
} else {
  ABC_sample <- readRDS('Data/Inputs/BK_exp_simData.RData')
  ABC_inputs <- readRDS('Data/Inputs/BK_exp_Ystats.RData')[[scenario]][['params']]
}

# Process data for `get_posterior()`
data_base <- list(
  X = ABC_sample$X,
  loc = tibble(
    id = rep(1:ABC_inputs$n, each = ABC_inputs$t),
    time = rep(1:ABC_inputs$t, times = ABC_inputs$n)
    ))

# Initialize.
results <- list()
sims <- names(ABC_sample[[scenario]])
sims <- sims[grepl('sim', sims)]

# Wait for it,
for (sim in sims) {
  data <- c(data_base, list(y = ABC_sample[[scenario]][[sim]]$y))
  postChain <- get_posterior(data, model = "production", 
                             n_samples = 25e3, burnin = 10e3, thinning = 10)
  results[[sim]][['postChain']] <- postChain
}

# Save results.
saveRDS(results, output_file)
