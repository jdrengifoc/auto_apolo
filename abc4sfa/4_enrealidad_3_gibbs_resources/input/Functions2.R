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
library(progress)

create_skewed <- function(Ri, sigmaepsilon2, sigmadelta2, lambdadelta) {
  function(x) {
    A <- exp(-crossprod(Ri - x * matrix(1, length(Ri[, 1]), 1),
                        Ri - x * matrix(1, dim(Ri)[1], 1)) / (2 * sigmaepsilon2)
             - x^2 / sigmadelta2) * pnorm(-lambdadelta * x / sqrt(sigmadelta2))
    return(drop(A))
  }
}

# Gibbs sampler for eta+
etasampling <- function(Di, Sigma, phi2, t) {
  Di <- matrix(Di, length(Di), 1)  # Ensure Di is a column vector
  m <- phi2 * sum(solve(Sigma, Di))
  eta <- rtruncnorm(1, 0, Inf, m, sd = sqrt(phi2))  # Posterior conditional eta+
  return(eta)
}

# Gibbs sampler for u+
usampling <- function(ui, sigmaepsilon2, sigmau2) {
  u <- rtruncnorm(1, 0, Inf, ui, sd = sqrt((sigmaepsilon2 * sigmau2) / (sigmaepsilon2 + sigmau2)))
  return(u)
}

# Gibbs sampler for delta
deltasampling1 <- function(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta) {
  Ri <- matrix(Ri, length(Ri), 1)
  skewed <- create_skewed(Ri, sigmaepsilon2, sigmadelta2, lambdadelta)
  log_skewed <- function(x) log(skewed(x))
  
  Mode <- grid[which.max(sapply(grid, skewed))]
  s2 <- -1 / hessian(skewed, x = Mode)
  delta <- rnorm(1, mean = Mode, sd = sqrt(s2))
  
  prob <- exp(log_skewed(delta) - log_skewed(Mode)) * exp(-(delta - Mode)^2 / (2 * s2))
  return(list(delta, prob))
}

get_posterior_time_limit <- function(
  data, model = 'cost', burnin_rate=0.3, time_limit_seconds = 60,
  max_na_iterations = 5, fixed_beta = TRUE, thinning = 10
  ) {
  # library(readr)
  library(mvtnorm)  # for rmnorm
  p <- ifelse(model == 'cost', 1, -1)
  # Prepare data matrices
  X <- as.matrix(data$X)
  # X <- as.matrix(cbind(1, X))
  y <- data$y
  N <- data$loc %>% dplyr::select(id) %>% unique() %>% nrow() %>% as.numeric()
  t <- data$loc %>% dplyr::select(time) %>% unique() %>% nrow() %>% as.numeric()
  Zt <- X
  
  # Initialize parameters
  sigmaepsilon2 <- 0.1^2
  sigmaeta2     <- 0.1^2
  sigmaalpha2   <- 0.1^2
  sigmau2       <- 0.1^2
  # beta0: vector de coeficientes estimados por MCO.
  b01 <- solve(t(X)%*%X)%*%t(X)%*%y %>% as.vector
  k <- length(b01)
  # Estimacion de cov por MCO.
  B01 <- as.numeric((1/(N*t-k))*crossprod(y - X%*%b01)) * solve(t(X) %*% X)
  
  etai <- abs(rnorm(N, 0, sqrt(sigmaeta2)))
  eta <- matrix(etai, N, t, byrow = TRUE)
  alphai <- matrix(rnorm(N, 0, sqrt(sigmaalpha2)), N, 1)
  alpha <- matrix(alphai, N, t, byrow = TRUE)
  deltai <- alphai + p*etai
  delta <- alpha + p*eta
  
  um <- matrix(abs(rnorm(N * t, 0, sqrt(sigmau2))), N, t)
  u <- as.vector(t(um))
  
  ym <- matrix(y, N, t, byrow = TRUE)
  
  nbar   <- 1
  nsigma <- N * t + nbar
  nsigma1 <- N + nbar
  qbar   <- 10^(-4)
  grid   <- seq(-5, 5, 0.05)
  
  # State tracking for reverting when NAs occur
  last_valid_state <- list(
    theta = b01,
    deltai = deltai,
    etai = etai,
    alphai = alphai,
    u = u,
    sigmaepsilon2 = sigmaepsilon2,
    sigmaeta2 = sigmaeta2,
    sigmaalpha2 = sigmaalpha2,
    sigmau2 = sigmau2,
    sigmadelta2 = sigmaalpha2 + sigmaeta2
  )
  consecutive_na_count <- 0
  
  results_list <- list()
  

  valid_iters <- 0
  start_time_SMC <- Sys.time()
  while (as.numeric(difftime(Sys.time(), start_time_SMC, units = "secs")) < time_limit_seconds) {
    print(valid_iters)
    tryCatch({
      Y <- as.vector(t(ym - delta - p*um))
      B2 <- chol2inv(chol((as.numeric(1/sigmaepsilon2) * crossprod(X, X) + B01)))
      b3 <- B2 %*% (as.numeric(1/sigmaepsilon2) * crossprod(X, Y) + B01 %*% b01)
      if (fixed_beta) {
        pp <- b01
        pp[1] <- b3[1]
        b3 <- pp
        pp <- B01
        pp[1, 1] <- B2[1, 1]
        B2 <- pp
      }
      theta <- t(rmnorm(1, b3, B2))
      if (fixed_beta) {
        theta[-1] <- last_valid_state$theta[-1]
      }

      Ztm <- matrix(Zt %*% theta, N, t, byrow = TRUE)
      sigmau2 <- (qbar + crossprod(u, u)) / rchisq(1, nsigma)
      R <- t(ym - (Ztm + p*um))
      
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      lambdadelta <- sqrt(sigmau2 / sigmaepsilon2)
      
      # For each observation, propose a new deltai using your custom function
      candp <- t(do.call(cbind, lapply(1:N, function(i) {
        Ri <- as.matrix(R[, i])
        deltasampling1(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta)
      })))
      
      # Update deltai with a Metropolis step
      deltai <- matrix(sapply(1:N, function(i) {
        ifelse(rbinom(1, 1, min(1, unlist(candp[i, 2]))) == 1, 
               as.numeric(candp[i, 1]), as.numeric(deltai[i]))
      }), N, 1)
      
      utilde <- as.vector(t(ym - (Ztm + delta)))
      u <- sapply(1:(N * t), function(i) usampling(utilde[i], sigmaepsilon2, sigmau2))
      
      um <- matrix(u, N, t, byrow = TRUE)
      epsilonhatm <- ym - (Ztm + delta + p*um)
      epsilonhat <- as.vector(t(epsilonhatm))
      
      sigmaepsilon2 <- 1 / rchisq(1, nsigma) * (qbar + crossprod(epsilonhat, epsilonhat))
      
      Sigma <- tcrossprod(matrix(1, t, 1), matrix(1, t, 1)) * as.numeric(sigmaalpha2) +
        as.numeric(sigmaepsilon2) * diag(t)
      D <- t(ym - (Ztm + p*um))
      phi2 <- sigmaeta2 * chol2inv(chol(1 + sigmaeta2 * matrix(1, 1, t) %*%
                                          chol2inv(chol(Sigma)) %*% matrix(1, t, 1)))
      
      etai <- matrix(sapply(1:N, function(i) etasampling(D[, i], Sigma, phi2, t)), N, 1)
      sigmaeta2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(etai, etai)))
      eta <- matrix(etai, N, t, byrow = FALSE)
      
      delta <- matrix(deltai, N, t, byrow = FALSE)
      
      alphai <- as.numeric(deltai) - p*as.numeric(etai)
      alpha <- matrix(alphai, N, t, byrow = FALSE)
      sigmaalpha2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(alphai, alphai)))
      
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      
      # Check for NaNs in key variables
      if (any(is.nan(c(theta, deltai, etai, alphai, u, sigmaepsilon2, sigmaeta2, sigmaalpha2)))) {
        consecutive_na_count <- consecutive_na_count + 1
	print("NA")
	print(consecutive_na_count)
        if (consecutive_na_count > max_na_iterations) {
          cat("Reverting to last valid state at iteration:", valid_iters, "\n")
          theta    <- last_valid_state$theta + rnorm(length(last_valid_state$theta), 0, 1e-3)
          deltai   <- last_valid_state$deltai + rnorm(length(last_valid_state$deltai), 0, 1e-3)
          etai     <- last_valid_state$etai + rnorm(length(last_valid_state$etai), 0, 1e-3)
          alphai   <- last_valid_state$alphai + rnorm(length(last_valid_state$alphai), 0, 1e-3)
          u        <- last_valid_state$u + rnorm(length(last_valid_state$u), 0, 1e-3)
          sigmaepsilon2 <- last_valid_state$sigmaepsilon2
          sigmaeta2   <- last_valid_state$sigmaeta2
          sigmaalpha2 <- last_valid_state$sigmaalpha2
          sigmau2     <- last_valid_state$sigmau2
          sigmadelta2 <- last_valid_state$sigmadelta2
          consecutive_na_count <- 0
        }
      } else {
        consecutive_na_count <- 0
        last_valid_state <- list(
          theta = theta,
          deltai = deltai,
          etai = etai,
          alphai = alphai,
          u = u,
          sigmaepsilon2 = sigmaepsilon2,
          sigmaeta2 = sigmaeta2,
          sigmaalpha2 = sigmaalpha2,
          sigmau2 = sigmau2,
          sigmadelta2 = sigmadelta2
        )
        valid_iters <- valid_iters + 1
      }
      
      if (valid_iters %% thinning == 0) {
        results_list[[valid_iters]] <- list(
          theta = theta,
          deltai = deltai,
          etai = etai,
          alphai = alphai,
          u = u,
          sigmau2 = sigmau2,
          sigmaepsilon2 = sigmaepsilon2,
          sigmaeta2 = sigmaeta2,
          sigmaalpha2 = sigmaalpha2,
          sigmadelta2 = sigmadelta2
        )
      }
    }, error = function(e) {
      cat("Error at iteration:", valid_iters, "\n")
      cat("Error message:", e$message, "\n")
    }, warning = function(w) {
      cat("Warning at iteration:", valid_iters, "\n")
      cat("Warning message:", w$message, "\n")
    })
  }

  # Remove any iterations that failed
  results <- results_list[!sapply(results_list, is.null)]
  
  extract_results <- function(results, param) {
    do.call(cbind, lapply(results, function(res) res[[param]]))
  }
  burnin <- floor(length(results_list) * burnin_rate)
  list(
    beta = extract_results(results, "theta")[, (burnin + 1):length(results)],
    delta = extract_results(results, "deltai")[, (burnin + 1):length(results)],
    eta = extract_results(results, "etai")[, (burnin + 1):length(results)],
    alpha = extract_results(results, "alphai")[, (burnin + 1):length(results)],
    u = extract_results(results, "u")[, (burnin + 1):length(results)],
    sigmau2 = extract_results(results, "sigmau2")[(burnin + 1):length(results)],
    sigmaepsilon2 = extract_results(results, "sigmaepsilon2")[(burnin + 1):length(results)],
    sigmaeta2 = extract_results(results, "sigmaeta2")[(burnin + 1):length(results)],
    sigmaalpha2 = extract_results(results, "sigmaalpha2")[(burnin + 1):length(results)],
    sigmadelta2 = extract_results(results, "sigmadelta2")[(burnin + 1):length(results)]
  )
}

get_posterior <- function(data, model = 'cost', burnin=5e3, n_samples = 15000, 
                                     max_na_iterations = 5, fixed_beta = TRUE) {
  #library(readr)
  library(mvtnorm)  # for rmnorm
  p <- ifelse(model == 'cost', 1, -1)
  # Prepare data matrices
  X <- as.matrix(data$X)
  X <- as.matrix(cbind(1, X))
  y <- data$y
  N <- data$loc %>% dplyr::select(id) %>% unique() %>% nrow() %>% as.numeric()
  t <- data$loc %>% dplyr::select(time) %>% unique() %>% nrow() %>% as.numeric()
  Zt <- X
  
  # Initialize parameters
  sigmaepsilon2 <- 0.1^2
  sigmaeta2     <- 0.1^2
  sigmaalpha2   <- 0.1^2
  sigmau2       <- 0.1^2
  # beta0: vector de coeficientes estimados por MCO.
  b01 <- solve(t(X)%*%X)%*%t(X)%*%y %>% as.vector
  k <- length(b01)
  # Estimacion de cov por MCO.
  B01 <- (1/(N*t-k))*crossprod(y - X%*%b01) * solve(t(X) %*% X)
  
  etai <- abs(rnorm(N, 0, sqrt(sigmaeta2)))
  eta <- matrix(etai, N, t, byrow = TRUE)
  alphai <- matrix(rnorm(N, 0, sqrt(sigmaalpha2)), N, 1)
  alpha <- matrix(alphai, N, t, byrow = TRUE)
  deltai <- alphai + p*etai
  delta <- alpha + p*eta
  
  um <- matrix(abs(rnorm(N * t, 0, sqrt(sigmau2))), N, t)
  u <- as.vector(t(um))
  
  ym <- matrix(y, N, t, byrow = TRUE)
  
  nbar   <- 1
  nsigma <- N * t + nbar
  nsigma1 <- N + nbar
  qbar   <- 10^(-4)
  grid   <- seq(-5, 5, 0.05)
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Sampling [:bar] :percent in :elapsed",
    total = n_samples,
    width = 60
  )
  
  # State tracking for reverting when NAs occur
  last_valid_state <- list(
    theta = b01,
    deltai = deltai,
    etai = etai,
    alphai = alphai,
    u = u,
    sigmaepsilon2 = sigmaepsilon2,
    sigmaeta2 = sigmaeta2,
    sigmaalpha2 = sigmaalpha2,
    sigmau2 = sigmau2,
    sigmadelta2 = sigmaalpha2 + sigmaeta2
  )
  consecutive_na_count <- 0
  
  results_list <- vector("list", n_samples)
  
  for (i in 1:n_samples) {
    tryCatch({
      Y <- as.vector(t(ym - delta - p*um))
      B2 <- chol2inv(chol((as.numeric(1/sigmaepsilon2) * crossprod(X, X) + B01)))
      b3 <- B2 %*% (as.numeric(1/sigmaepsilon2) * crossprod(X, Y) + B01 %*% b01)
      if (fixed_beta) {
        pp <- b01
        pp[1] <- b3[1]
        b3 <- pp
        pp <- B01
        pp[1, 1] <- B2[1, 1]
        B2 <- pp
      }
      theta <- t(rmnorm(1, b3, B2))
      if (fixed_beta) {
        theta[-1] <- last_valid_state$theta[-1]
      }

      Ztm <- matrix(Zt %*% theta, N, t, byrow = TRUE)
      sigmau2 <- (qbar + crossprod(u, u)) / rchisq(1, nsigma)
      R <- t(ym - (Ztm + p*um))
      
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      lambdadelta <- sqrt(sigmau2 / sigmaepsilon2)
      
      # For each observation, propose a new deltai using your custom function
      candp <- t(do.call(cbind, lapply(1:N, function(i) {
        Ri <- as.matrix(R[, i])
        deltasampling1(Ri, grid, sigmaepsilon2, sigmadelta2, lambdadelta)
      })))
      
      # Update deltai with a Metropolis step
      deltai <- matrix(sapply(1:N, function(i) {
        ifelse(rbinom(1, 1, min(1, unlist(candp[i, 2]))) == 1, 
               as.numeric(candp[i, 1]), as.numeric(deltai[i]))
      }), N, 1)
      
      utilde <- as.vector(t(ym - (Ztm + delta)))
      u <- sapply(1:(N * t), function(i) usampling(utilde[i], sigmaepsilon2, sigmau2))
      
      um <- matrix(u, N, t, byrow = TRUE)
      epsilonhatm <- ym - (Ztm + delta + p*um)
      epsilonhat <- as.vector(t(epsilonhatm))
      
      sigmaepsilon2 <- 1 / rchisq(1, nsigma) * (qbar + crossprod(epsilonhat, epsilonhat))
      
      Sigma <- tcrossprod(matrix(1, t, 1), matrix(1, t, 1)) * as.numeric(sigmaalpha2) +
        as.numeric(sigmaepsilon2) * diag(t)
      D <- t(ym - (Ztm + p*um))
      phi2 <- sigmaeta2 * chol2inv(chol(1 + sigmaeta2 * matrix(1, 1, t) %*%
                                          chol2inv(chol(Sigma)) %*% matrix(1, t, 1)))
      
      etai <- matrix(sapply(1:N, function(i) etasampling(D[, i], Sigma, phi2, t)), N, 1)
      sigmaeta2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(etai, etai)))
      eta <- matrix(etai, N, t, byrow = FALSE)
      
      delta <- matrix(deltai, N, t, byrow = FALSE)
      
      alphai <- as.numeric(deltai) - p*as.numeric(etai)
      alpha <- matrix(alphai, N, t, byrow = FALSE)
      sigmaalpha2 <- 1 / (rchisq(1, nsigma1) / (qbar + crossprod(alphai, alphai)))
      
      sigmadelta2 <- sigmaalpha2 + sigmaeta2
      
      # Check for NaNs in key variables
      if (any(is.nan(c(theta, deltai, etai, alphai, u, sigmaepsilon2, sigmaeta2, sigmaalpha2)))) {
        # Record the iteration where NA occurred
        last_valid_state$iteration <- i
        consecutive_na_count <- consecutive_na_count + 1
        if (consecutive_na_count > max_na_iterations) {
          cat("Reverting to last valid state at iteration:", i, "\n")
          theta    <- last_valid_state$theta + rnorm(length(last_valid_state$theta), 0, 1e-3)
          deltai   <- last_valid_state$deltai + rnorm(length(last_valid_state$deltai), 0, 1e-3)
          etai     <- last_valid_state$etai + rnorm(length(last_valid_state$etai), 0, 1e-3)
          alphai   <- last_valid_state$alphai + rnorm(length(last_valid_state$alphai), 0, 1e-3)
          u        <- last_valid_state$u + rnorm(length(last_valid_state$u), 0, 1e-3)
          sigmaepsilon2 <- last_valid_state$sigmaepsilon2
          sigmaeta2   <- last_valid_state$sigmaeta2
          sigmaalpha2 <- last_valid_state$sigmaalpha2
          sigmau2     <- last_valid_state$sigmau2
          sigmadelta2 <- last_valid_state$sigmadelta2
          consecutive_na_count <- 0
        }
      } else {
        consecutive_na_count <- 0
        last_valid_state <- list(
          theta = theta,
          deltai = deltai,
          etai = etai,
          alphai = alphai,
          u = u,
          sigmaepsilon2 = sigmaepsilon2,
          sigmaeta2 = sigmaeta2,
          sigmaalpha2 = sigmaalpha2,
          sigmau2 = sigmau2,
          sigmadelta2 = sigmadelta2
        )
      }
      
      results_list[[i]] <- list(
        theta = theta,
        deltai = deltai,
        etai = etai,
        alphai = alphai,
        u = u,
        sigmau2 = sigmau2,
        sigmaepsilon2 = sigmaepsilon2,
        sigmaeta2 = sigmaeta2,
        sigmaalpha2 = sigmaalpha2,
        sigmadelta2 = sigmadelta2
      )
    }, error = function(e) {
      cat("Error at iteration:", i, "\n")
      cat("Error message:", e$message, "\n")
      results_list[[i]] <- NULL
    }, warning = function(w) {
      cat("Warning at iteration:", i, "\n")
      cat("Warning message:", w$message, "\n")
      results_list[[i]] <- NULL
    })
    pb$tick()
  }
  
  # Remove any iterations that failed
  results <- results_list[!sapply(results_list, is.null)]
  
  extract_results <- function(results, param) {
    do.call(cbind, lapply(results, function(res) res[[param]]))
  }
  
  list(
    beta = extract_results(results, "theta")[, (burnin + 1):length(results)],
    delta = extract_results(results, "deltai")[, (burnin + 1):length(results)],
    eta = extract_results(results, "etai")[, (burnin + 1):length(results)],
    alpha = extract_results(results, "alphai")[, (burnin + 1):length(results)],
    u = extract_results(results, "u")[, (burnin + 1):length(results)],
    sigmau2 = extract_results(results, "sigmau2")[(burnin + 1):length(results)],
    sigmaepsilon2 = extract_results(results, "sigmaepsilon2")[(burnin + 1):length(results)],
    sigmaeta2 = extract_results(results, "sigmaeta2")[(burnin + 1):length(results)],
    sigmaalpha2 = extract_results(results, "sigmaalpha2")[(burnin + 1):length(results)],
    sigmadelta2 = extract_results(results, "sigmadelta2")[(burnin + 1):length(results)]
  )
}
