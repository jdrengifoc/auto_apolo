# Load libraries ----------------------------------------------------------

library(LearnBayes)
library(truncnorm)
library(invgamma)
library(MASS)
library(spdep)
library(numDeriv)
library(tmvtnorm)
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

etasampling <- function(Di, Sigma, t, sigmaeta2) {
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

deltasampling1 <- function(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta) {
  Ri <- matrix(Ri, length(Ri), 1)
  skewed <- create_skewed(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
  log_skewed <- create_log_skew(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
  log_skewed_min <- function(x) -log_skewed(x)
  
  Mode <- optim(-0.5, log_skewed_min, method = 'Brent', 
                upper = 50, lower = -50)$par
  
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
get_posterior <- function(data, model = "cost", burnin=10, 
                          n_samples=250, thinning=2) {
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
  sigmaepsilon2 <- 0.058^2
  sigmaeta2 <- 0.094^2
  sigmaalpha2 <- 0.163^2
  sigmau2 <- 0.217^2
  
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
  
  n_keep <- floor((n_samples - burnin) / thinning)
  # Storage for chains
  sigmau2chain    <- numeric(n_keep)
  sigmaepsilon2chain <- numeric(n_keep)
  sigmaeta2chain  <- numeric(n_keep)
  sigmaalpha2chain<- numeric(n_keep)
  sigmadelta2chain<- numeric(n_keep)
  
  Xi <- lapply(1:N, function(i) get_Xi(i, X, N, t))
  
  start_time <- Sys.time()
  
  results <- vector("list", n_samples)
  keep <- 0L
  
  for (i in seq_len(n_samples)) {
    
    res <- tryCatch({
      
      ## ——— cuerpo original de la función anónima en future_map ———
      
      Ym <- ym - p*delta - p*um
      XSX <- lapply(Xi, function(X){
        t(X) %*% solve(Sigma) %*% X
      })
      sumXSX <- Reduce('+', XSX)
      B2 <- solve(sumXSX + B01)
      
      XSY <- lapply(1:N, function(i){
        t(Xi[[i]]) %*% solve(Sigma) %*% Ym[i,]
      })
      sumXSY <- Reduce('+', XSY)
      b3    <- B2 %*% (sumXSY + B01 %*% b01)
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
               as.numeric(candp[i, 1]), as.numeric(deltai[i]))
      }), N, 1)
      
      epsilonhatm <- ym - (Ztm + delta + p*um)
      epsilonhat  <- as.vector(epsilonhatm)
      
      sigmaepsilon2 <- 1 / rchisq(1, nsigma) * (qbar + crossprod(epsilonhat, epsilonhat))
      utilde <- as.vector(ym - (Ztm + delta))
      u      <- sapply(1:(N*t), function(i)
        usampling(utilde[i], sigmau2, sigmaepsilon2))
      um <- matrix(u, N, t, byrow = FALSE)
      
      D    <- ym - (Ztm + p*um)
      etai <- matrix(sapply(1:N, function(i)
        etasampling(D[i, ], Sigma, t, sigmaeta2)), N, 1)
      sigmaeta2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(etai, etai)))
      eta <- matrix(etai, N, t, byrow = FALSE)
      
      delta <- matrix(deltai, N, t, byrow = FALSE)
      
      alphai <- as.numeric(deltai) - p*as.numeric(etai)
      alpha <- matrix(alphai, N, t, byrow = FALSE)
      sigmaalpha2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(alphai, alphai)))
      
      Sigma      <- tcrossprod(matrix(1, t, 1), matrix(1, t, 1)) *
        as.numeric(sigmaalpha2) +
        as.numeric(sigmaepsilon2) * diag(t)
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      
      # Si hay NaNs, devuelve NULL para saltar esta iteración
      if (any(is.nan(c(theta, deltai, etai, alphai, u,
                       sigmaepsilon2, sigmaeta2, sigmaalpha2)))) {
        NULL
      } else {
        list(
          theta         = theta,
          deltai        = deltai,
          etai          = etai,
          alphai        = alphai,
          u             = u,
          epsilon       = epsilonhat,
          sigmau2       = sigmau2,
          sigmaepsilon2 = sigmaepsilon2,
          sigmaeta2     = sigmaeta2,
          sigmaalpha2   = sigmaalpha2,
          sigmadelta2   = sigmadelta2
        )
      }
      
    }, error   = function(e) NULL,
    warning = function(w) NULL)
    
    if (is.null(res)) next
    
    # guarda sólo cada thinning-ésima iteración
    if (i > burnin && ((i - burnin) %% thinning) == 0L) {
      keep <- keep + 1L
      sigmau2chain[keep]         <- res$sigmau2
      sigmaepsilon2chain[keep]   <- res$sigmaepsilon2
      sigmaeta2chain[keep]       <- res$sigmaeta2
      sigmaalpha2chain[keep]     <- res$sigmaalpha2
    }
  }
  
  elapsed_time <- Sys.time() - start_time
  
  # retorno final con cadenas ya "thinned"
  list(
    sigmau2chain       = sigmau2chain,
    sigmaepsilon2chain = sigmaepsilon2chain,
    sigmaeta2chain     = sigmaeta2chain,
    sigmaalpha2chain   = sigmaalpha2chain,
    elapsed_time       = elapsed_time
  )
}