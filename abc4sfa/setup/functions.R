
# Simulations -------------------------------------------------------------
# Simulate data.
#' Simulate Badunenko & Kumbhakar regressors.
#' @param theta Parameter list, must include the sample size (n) and (t) as well
#' as the B&K parameters (b).
#' @param fig If TRUE, then plot the simulated data set. The default is FALSE.
#' @return X Regresors (includes intercept) of dimension n x p that distribute
#' f_b(x) = 1 / (b-1) exp(log(b)-x) for 0 < x < log(b).
genData_BK <- function(theta, fig = F){
  # Parameters.
  n <- theta$n
  t <- theta$t
  b <- theta$b
  # Generate data.
  x1 <- -log(1 - (b[1]-1)/b[1] * runif(n*t))
  x2 <- -log(1 - (b[2]-1)/b[2] * runif(n*t))
  X <- cbind(rep(1, n*t), x1, x2)
  # Plot data.
  if (fig) {
    hist(x1, probability = T, col = "#9D388F")
    curve(1/(b[1]-1)*exp(log(b[1]) - x),  0, log(b[1]),
          add = T, lwd = 4, xlab = "x1",
          main = paste("Histogram X1 (b =", b[1]))
    hist(x2, probability = T, col = "#9D388F")
    curve(1/(b[2]-1)*exp(log(b[2]) - x),  0, log(b[2]),
          add = T, lwd = 4, xlab = "x2",
          main = paste("Histogram X2 (b =", b[2]))
  }
  return(X)
}


#' GTRE normal-halfnormal/exponential model. The theta$sigmas are the standar
#' deviation following the notation of Badunenki&Kumbakar.
#' @param X Regressors matrix of dimension n x p.
#' @param theta Parameter list, must include beta, the sample size and the
#' variances of the error terms.
#' @param dist Inefficiency's distribution name; support capital letters.
#' Halfnoralm: "halfnormal", "half", "hn"; default.
#' Exponential: "exponential", "exp", "e".
#' @param name Models name.
#' Linear model: "ln", "linear; default.
#' Cobb-Douglas: "cobb-douglas", "cb".
#' @return the response variable of the GTRE model.
SFM <- function(X, theta, name = "ln"){
  # Parameter treatment.
  name = tolower(name)
  n = theta$n
  t = theta$t
  dist <- ifelse(is.null(theta$model), 'hn', tolower(theta$model))
  
  # Simulate inefficiency.
  if (dist %in% c("halfnormal", "half", "hn")){
    #u <- fdrtool::rhalfnorm(n*t, theta = sqrt(pi/2)/theta$sigma[3])
    u <- fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*theta$sigma[3]^2)))
    # Permanent inefficiency.
    u0 <- fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*theta$sigma[4]^2)))
    #u0 <- fdrtool::rhalfnorm(n, theta = sqrt(pi/2)/theta$sigma[4])
  }
  else if (dist %in% c("exponential", "exp", "e")){
    # Transitory inefficiency.
    u <- rexp(n*t, rate = 1/theta$sigma[3])
    # Permanent inefficiency.
    u0 <- rexp(n, rate = 1/theta$sigma[4])
  }
  u0 <- rep(u0, each = T)
  
  # Measurement error.
  v <- rnorm(n*t, 0, theta$sigma[1])
  # Unobserved heterogeneity.
  v0 <- rnorm(n, 0, theta$sigma[2])
  v0 <- rep(v0, each = T)
  
  # Model.
  model <- function(X, beta, v, v0, u, u0, name = "ln"){
    if (name %in% c("ln", "linear")) {
      y <- X%*%theta$beta + v + v0 - u - u0
    } 
    else if (name %in% c("cobb-douglas", "cb")) {
      if (all(X[,1]==1)) {
        # With intercept, viz. technology.
        y <- exp(beta[1]) * apply(t(t(X[,-1])^theta$beta[-1]),1, sum) * exp(v + v0 - u - u0)
      } else {
        # Without intercept, viz. technology.
        y <- apply(t(t(X)^theta$beta), 1, sum) * exp(v + v0 - u - u0)
      }
    }
    else {
      warning("Non supported model name.")
    }
    return(as.vector(y))
  }
  y <- model(X, theta$beta, v, v0, u, u0, name)
  return(y)
}
# Unsed?
SFM_all <- function(X, theta, name = "ln"){
  # Parameter treatment.
  name = tolower(name)
  n = theta$n
  t = theta$t
  dist <- ifelse(is.null(theta$model), 'hn', tolower(theta$model))
  
  # Simulate inefficiency.
  if (dist %in% c("halfnormal", "half", "hn")){
    # Transitory inefficiency.
    #u <- fdrtool::rhalfnorm(n*t, theta = sqrt(pi/2)/theta$sigma[3])
    u <- fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*theta$sigma[3]^2)))
    # Permanent inefficiency.
    u0 <- fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*theta$sigma[4]^2)))
    #u0 <- fdrtool::rhalfnorm(n, theta = sqrt(pi/2)/theta$sigma[4])
  }
  else if (dist %in% c("exponential", "exp", "e")){
    # Transitory inefficiency.
    u <- rexp(n*t, rate = 1/theta$sigma[3])
    # Permanent inefficiency.
    u0 <- rexp(n, rate = 1/theta$sigma[4])
  }
  u0 <- rep(u0, each = T)
  
  # Measurement error.
  v <- rnorm(n*t, 0, theta$sigma[1])
  # Unobserved heterogeneity. 
  v0 <- rnorm(n, 0, theta$sigma[2])
  v0 <- rep(v0, each = T)
  
  # Model.
  model <- function(X, beta, v, v0, u, u0, name = "ln"){
    if (name %in% c("ln", "linear")) {
      y <- X%*%theta$beta + v + v0 - u - u0
    } 
    else if (name %in% c("cobb-douglas", "cb")) {
      if (all(X[,1]==1)) {
        # With intercept, viz. technology.
        y <- exp(beta[1]) * apply(t(t(X[,-1])^theta$beta[-1]),1, sum) * exp(v + v0 - u - u0)
      } else {
        # Without intercept, viz. technology.
        y <- apply(t(t(X)^theta$beta), 1, sum) * exp(v + v0 - u - u0)
      }
    }
    else {
      warning("Non supported model name.")
    }
    return(as.vector(y))
  }
  y <- model(X, theta$beta, v, v0, u, u0, name)
  return(list(y = y, u0 = u0, u = u, v0 = v0, v = v))
}

#' Simulate the SFA residuals.
#' @param theta Parameter list, the sample size.
#' @param prior Atomic vector with the intercept B0, sigma_v, sigma_v0, sigma_u,
#' sigma_u0.
#' @param dist Inefficiency's distribution name; support capital letters.
#' Halfnoralm: "halfnormal", "half", "hn"; default.
#' Exponential: "exponential", "exp", "e".
#' @return The residuals r = beta0 + v + v0 - (u - mean(u)) - (u0 - mean(u0)).
simRes <- function(theta, prior) {
  # Parameter treatment.
  n = theta$n
  t = theta$t
  dist <- ifelse(is.null(theta$model), 'hn', tolower(theta$model))
  
  # Simulate inefficiency.
  if (dist %in% c("halfnormal", "half", "hn")){
    # Transitory inefficiency.
    # u <- fdrtool::rhalfnorm(n*t, theta = sqrt(pi/2)/prior[4])
    u <- fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*prior[4]^2)))
    # Permanent inefficiency.
    # u0 <- fdrtool::rhalfnorm(n, theta = sqrt(pi/2)/prior[5])
    u0 <- fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*prior[5]^2)))
  }
  else if (dist %in% c("exponential", "exp", "e")){
    # Transitory inefficiency.
    u <- rexp(n*t, rate = 1/prior[4])
    # Permanent inefficiency.
    u0 <- rexp(n, rate = 1/prior[5])
  }
  u0 <- rep(u0, each = t)
  
  # Measurement error.
  v <- rnorm(n*t, 0, prior[2])
  # Unobserved heterogeneity.
  v0 <- rnorm(n, 0, prior[3])
  v0 <- rep(v0, each = t)
  res <- prior[1] + v + v0 - (u - mean(u)) - (u0 - mean(u0))
  return(res)
}

#' Compute the summary statistics of ABC for the simulated residuals.
#' @param theta Parameter list, must include beta, the sample size and the
#' variances of the error terms.
#' @param resS `simRes()` residuals.
#' @return Summary statistics of the residuals: m1, individual variance,
#' idiosyncratic variance, m3, m4. Notice that the variance is unbiased.
getSumStats <- function(theta, resS){
  # Parameters treatment.
  t = theta$t
  n = theta$n
  # Mean of the residuals.
  m1 <- mean(resS)
  # Total unbiased variance of the residuals.
  s2TotS <- var(resS)*(t*n - 1)/(t*n - length(theta$beta))
  # Central residuals with its "ID".
  ResAPs <- cbind(rep(1:n, each = t), rep(1:t, times = n), resS - m1)
  # Initialize the sum to compute the individual and idiosyncratic variance.
  Sum_it <- 0
  for(i in 1:n){
    for(j in 1:(t-1)){
      for(l in (j+1):t){
        Res_it <- ResAPs[which(ResAPs[,1]==i), 3]
        Sum_it <- Sum_it + Res_it[j] * Res_it[l]
      }
    }
  }
  # Individual unbiased variance of the residuals.
  s2IndS <- Sum_it / (n*t*(t-1)/2 - length(theta$beta))
  # Idiosincratic variance (V_tot = V_ind + V_idio).
  s2IdioS <- s2TotS - s2IndS
  # Third and fourth central moments.
  m3 <- mean((resS-m1)^3)
  m4 <- mean((resS-m1)^4)
  # Return the statistis.
  return(c(m1, s2IndS, s2IdioS, m3, m4))
}

#' Compute the summary statistics of ABC for the OLS residuals, viz. the sample
#' statistics..
#' @param theta Parameter list, must include beta, the sample size and the
#' variances of the error terms.
#' @param y output variable.
#' @param X input matrix of dimensions n x p.
#' @return Summary statistics of the residuals: beta0/m1, individual variance,
#' idiosyncratic variance, m3, m4. Notice that the variance is unbiased.
getSampleSumStats <- function(theta, y, X){
  # Parameters treatment.
  t = theta$t
  n = theta$n
  # OLS residuals.
  reg <- lm(y ~ X[,2] + X[,3])
  resOLS <- reg$res
  # OLS residuals' unbiased variance.
  OLSvar <- sum(resOLS^2)/(n*t - length(theta$beta))
  # OLS residuals with its "ID".
  ResAP <- cbind(1:n, 1:t, resOLS)
  # Initialize the sum to compute the individual and idiosyncratic variance.
  Sum_it <- 0
  for(i in 1:n){
    for(j in 1:(t-1)){
      for(l in (j+1):t){
        Res_it <- ResAP[which(ResAP[,1]==i),3]
        Sum_it <- Sum_it + Res_it[j] * Res_it[l]
      }
    }
  }
  # Individual variance with small sample correction of the residuals.
  s2IndOLS <- Sum_it / (n*t*(t-1)/2 - length(theta$beta))
  # Idiosincratic variance (V_tot = V_ind + V_idio).
  s2IdioOLS <- OLSvar - s2IndOLS
  # Return the statistis.
  Ystats <- c(reg$coef[[1]], s2IndOLS, s2IdioOLS,
              mean(resOLS^3), mean(resOLS^4))
  return(as.vector(Ystats))
}

#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' @param theta Parameter list, must include FALTA.
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @return Return the best `a` * `S` priors and its square euclidean distance.
ABCpanelSFA <- function(theta, SumStatsY) {
  # Parameter treatment.
  S <- theta$S
  
  # Selected samples.
  SelS <- round(S*theta$delta)
  
  # Simulate the priors.
  priors <- sapply(1:length(theta$lb), function(i) {runif(S, theta$lb[i], theta$ub[i])})
  # Simulate the residuals for each row of the prior combination.
  z <- apply(priors, 1, function(prior) {
    simRes(theta = theta, prior = prior)})
  # Get the summary statistics of the residuals.
  SumStatsZs <- sapply(1:S, function(x) {
    getSumStats(theta, z[,x])
  })
  # Compute the square of the euclidean distance of the simulated summary statistics against
  # the sample statistics.
  EucDist <- sapply(1:S, function(x) {
    sum((SumStatsZs[, x]-SumStatsY)^2)
  })
  # Sort the distances.
  I <- order(EucDist)
  selPrior <- priors[I, ][1:SelS, ]
  selDist <- EucDist[I][1:SelS]
  # Return the best priors and its euclidean distance.
  ABCpostChain <- cbind(SelPrior, selDist)
  colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                              'sigma_u0', 'EucDist')
  return(ABCpostChain)
}

#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' @param theta Parameter list, must include FALTA and the cluster..
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @return Return the best `a` * `S` priors and its distance.
parallel_ABCpanelSFA <- function(theta, SumStatsY) {
  # Parameter treatment.
  S <- theta$S
  cluster <- theta$cluster
  # Simulate the priors B0, sigma_v, sigma_v0, sigma_u, sigma_u0.
  priors <- parSapply(cluster, 1:length(theta$lb), function(i) {runif(S, theta$lb[i], theta$ub[i])})
  # Simulate the residuals for each row of the prior combination.
  z <- parApply(cluster, priors, 1, function(prior) {
    simRes(theta = theta, prior = prior)})
  # Get the summary statistics of the residuals.
  SumStatsZs <- parSapply(cluster, 1:S, function(x) {
    getSumStats(theta, z[,x])
  })
  # Compute the square euclidean distance of the simulated summary statistics against
  # the sample statistics.
  EucDist <- parSapply(cluster, 1:S, function(x) {
    sum((SumStatsZs[, x]-SumStatsY)^2)
  })
  # Sort the distances.
  I <- order(EucDist)
  EucDistOrd <- EucDist[I]
  OrdPrior <- priors[I, ]
  selDist <- EucDistOrd[1:round(S * theta$delta)]
  SelPrior <- OrdPrior[1:round(S * theta$delta), ]
  # Return the best priors and its euclidean distance.
  ABCpostChain <- cbind(SelPrior, selDist)
  colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                              'sigma_u0', 'EucDist')
  return(ABCpostChain)
}

#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' In the urge to optimize memory use, this rutine divide the simulation in smaller
#' simulations -"batches"-.
#' @param theta Parameter list, must include FALTA and the cluster.
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @param maxSelS Number of maximum combinations of parameters that will be returned.
#' @param SBatch Number of combination for each batch.
#' @return Return the best `a` * `S` priors and its square euclidean distance.
batchABCpanelSFA <- function(theta, SumStatsY, maxSelS = 100, SBatch = 5e5) {
  # Parameter treatment.
  S <- theta$S
  
  # The selected samples are maximum maxSelS.
  SelS <- min(round(S*theta$delta), maxSelS)
  if (S <= SBatch){
    if (theta$parallel){
      ABCpostChain <- parallel_ABCpanelSFA(theta, SumStatsY)
    } else {
      ABCpostChain <- ABCpanelSFA(theta, SumStatsY)
    }
    
  } else {
    # S for each "batch" (iteration).
    SBatches <- c(rep(SBatch, floor(S / SBatch)), S %% SBatch)
    SBatches <- SBatches[SBatches != 0]
    # Preallocation.
    EucDist <- rep(Inf, SelS)
    priors <- matrix(rep(0, SelS*5), ncol = 5)
    thetaBatch <- theta
    
    # Make ABC for each batch and update the best distances.
    cont <- 0
    for (s in SBatches){
      cat(paste('Batch', cont, '/', length(SBatches) , ' time: ', Sys.time(), '\n'))
      # Modify the delta to return up to maxSelS combinations.
      thetaBatch$delta <- min(1, maxSelS / s)
      # Modify the S.
      thetaBatch$S <- s
      # Run ABC for the batch.
      if (theta$parallel){
        batchABCpostChain <- parallel_ABCpanelSFA(thetaBatch, SumStatsY)
      } else {
        batchABCpostChain <- ABCpanelSFA(thetaBatch, SumStatsY)
      }
      BatchDist <- batchABCpostChain[,6]
      BatchPriors <- batchABCpostChain[,1:5]
      
      # There is a smaller distance?
      existSmaller <- BatchDist[1] < EucDist[SelS]
      # Substitute the best combinations until there is no new better combination.
      while (existSmaller && length(BatchDist) > 0) {
        # Sustitute the smaller new distance with the largest old distance.
        EucDist[SelS] <- BatchDist[1]
        priors[SelS,] <- BatchPriors[1,]
        # Sort by distance.
        I <- order(EucDist)
        EucDist <- EucDist[I][1:SelS]
        priors <- priors[I, ][1:SelS, ]
        # Delete distance.
        BatchDist <- BatchDist[-1]
        # Update control.
        existSmaller <- BatchDist[1] < EucDist[SelS]
      }
      cont <- cont + 1
    }
    
    ABCpostChain <- cbind(priors, EucDist)
    colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                                'sigma_u0', 'EucDist')
  }
  return(ABCpostChain)
}

#' Function that returns the number of samples that must be simulated to
#' guarantee O-big convergency.
#' @param n: Number of samples (dot<t_i, n_i>)
#' @param d: Number of parameters to be estimated.
#' @param S_best: number of best samples to be selected.
getNumSamplesToGuaranteeConvergence <- function(n, d = 5, S_best = 100){
  delta <- n^(-d/2) * log(n) ^ (-1)
  # Number of simulation
  S <- S_best / delta
  return(S)
}



# Empirical data ----------------------------------------------------------
#' Simulate the SFA residuals.
#' @param theta Parameter list, the sample size.
#' @param prior Atomic vector with the intercept B0, sigma_v, sigma_v0, sigma_u,
#' sigma_u0.
#' @param dist Inefficiency's distribution name; support capital letters.
#' Halfnoralm: "halfnormal", "half", "hn"; default.
#' Exponential: "exponential", "exp", "e".
#' @return The residuals r = beta0 + v + v0 - (u - mean(u)) - (u0 - mean(u0)).
simResEmpirical <- function(theta, prior) {
  repEachVec <- function(x, eachs){
    # Preallocate.
    y <- rep(0, sum(eachs))
    n <- length(x)
    # Check dimensions.
    if (n == length(eachs)) {
      idx1 <- 1
      idx2 <- eachs[1]
      for (i in 1:n) {
        each <- eachs[i]
        y[idx1:idx2] <- x[i]
        if (i != n){
          idx1 <- idx1 + each
          idx2 <- idx2 + eachs[i+1]
        }
      }
    } else {
      stop("`x` and `eachs` must have the same length")
    }
    return(y)
  }
  # Parameter treatment.
  n <- theta$n
  ts <- theta$ts
  dist <- ifelse(is.null(theta$model), 'hn', tolower(theta$model))
  
  # Simulate inefficiency.
  if (dist %in% c("halfnormal", "half", "hn")){
    # Transitory inefficiency.
    # u <- fdrtool::rhalfnorm(sum(ts), theta = sqrt(pi/2)/prior[4])
    u <- fdrtool::rhalfnorm(n*ts, theta = sqrt((pi-2)/(2*prior[4]^2)))
    # Permanent inefficiency.
    #u0 <- fdrtool::rhalfnorm(n, theta = sqrt(pi/2)/prior[5])
    u0 <- fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*prior[5]^2)))
  }
  else if (dist %in% c("exponential", "exp", "e")){
    # Transitory inefficiency.
    u <- rexp(sum(ts), rate = 1/prior[4])
    # Permanent inefficiency.
    u0 <- rexp(n, rate = 1/prior[5])
  }
  u0 <- repEachVec(u0, eachs = ts)
  
  # Measurement error.
  v <- rnorm(sum(ts), 0, prior[2])
  # Unobserved heterogeneity.
  v0 <- rnorm(n, 0, prior[3])
  v0 <- repEachVec(v0, eachs = ts)
  
  res <- prior[1] + v + v0 - (u - mean(u)) - (u0 - mean(u0))
  return(res)
}

#' Compute the summary statistics of ABC for the simulated residuals.
#' @param theta Parameter list, must include beta, the sample size and the
#' variances of the error terms.
#' @param resS `simRes()` residuals.
#' @return Summary statistics of the residuals: m1, individual variance,
#' idiosyncratic variance, m3, m4. Notice that the variance is unbiased.
getSumStatsEmpirical <- function(theta, resS){
  # Parameters treatment.
  ts = theta$ts
  n = theta$n
  # Mean of the residuals.
  m1 <- mean(resS)
  # Total unbiased variance of the residuals.
  s2TotS <- var(resS)*(sum(ts) - 1) / (sum(ts) - length(names(theta$X)))
  # Central residuals with its "ID".
  ResAPs <- cbind(theta$loc, resS - m1)
  # Initialize the sum to compute the individual and idiosyncratic variance.
  Sum_it <- 0
  for(i in 1:n){
    t <- length(ResAPs[ResAPs$id==i,3])
    for(j in 1:(t-1)){
      for(l in (j+1):t){
        Res_it <- ResAPs[ResAPs[,1]==i,3]
        Sum_it <- Sum_it + Res_it[j] * Res_it[l]
      }
    }
  }
  
  # Individual unbiased variance of the residuals.
  #### WARNING@
  s2IndS <- Sum_it / (sum(ts))#*(t-1)/2 - length(names(data$X)))
  # Idiosincratic variance (V_tot = V_ind + V_idio).
  s2IdioS <- s2TotS - s2IndS
  # Third and fourth central moments.
  m3 <- mean((resS-m1)^3)
  m4 <- mean((resS-m1)^4)
  # Return the statistis.
  return(c(m1, s2IndS, s2IdioS, m3, m4))
}

#' Compute the summary statistics of ABC for the OLS residuals, viz. the sample
#' statistics..
#' @param theta Parameter list, must include beta, the sample size and the
#' variances of the error terms.
#' @param y output variable.
#' @param X input matrix of dimensions n x p.
#' @return Summary statistics of the residuals: beta0/m1, individual variance,
#' idiosyncratic variance, m3, m4. Notice that the variance is unbiased.
getSampleSumStatsEmpirical <- function(theta, y, X){
  # Parameters treatment.
  ts <- theta$ts
  n <- theta$n
  loc <- theta$loc
  
  # OLS residuals.
  if (is.null(loc)){
    reg <- lm(y ~ ., data.frame(y, X))
  } else {
    formulae <- paste('y ~', paste(names(data.frame(X)), collapse = ' + '))
    reg <- plm(formulae, data.frame(y, loc, X), model="random",random.method = 'amemiya',
               index = c("id", "time"))  
  }
  
  resOLS <- reg$res
  # OLS residuals' unbiased variance.
  OLSvar <- sum(resOLS^2)/(sum(ts) - length(names(X)))
  # OLS residuals with its "ID".
  ResAP <- cbind(theta$loc, resOLS)
  # Initialize the sum to compute the individual and idiosyncratic variance.
  Sum_it <- 0
  for(i in 1:n){
    t <- dim(ResAP[ResAP$id==i,])[1]
    for(j in 1:(t-1)){
      for(l in (j+1):t){
        Res_it <- ResAP[ResAP[,1]==i,3]
        Sum_it <- Sum_it + Res_it[j] * Res_it[l]
      }
    }
  }
  # Individual variance with small sample correction of the residuals.
  ##### WARNING #####
  s2IndOLS <- Sum_it / (sum(ts)*(t-1)/2 - length(names(X)))
  # Idiosincratic variance (V_tot = V_ind + V_idio).
  s2IdioOLS <- OLSvar - s2IndOLS
  # Return the statistis.
  Ystats <- c(reg$coef[[1]], s2IndOLS, s2IdioOLS,
              mean(resOLS^3), mean(resOLS^4))
  return(as.vector(Ystats))
  
}
#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' @param theta Parameter list, must include FALTA.
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @return Return the best `a` * `S` priors and its square euclidean distance.
ABCpanelSFAEmpirical <- function(theta, SumStatsY) {
  # Parameter treatment.
  S <- theta$S
  
  # Selected samples.
  SelS <- theta$S_best
  
  # Simulate the priors.
  priors <- sapply(1:length(theta$lb), function(i) {runif(S, theta$lb[i], theta$ub[i])})
  # Simulate the residuals for each row of the prior combination.
  z <- apply(priors, 1, function(prior) {
    simResEmpirical(theta = theta, prior = prior)})
  # Get the summary statistics of the residuals.
  SumStatsZs <- sapply(1:S, function(x) {
    getSumStatsEmpirical(theta, z[,x])
  })
  # Compute the square of the euclidean distance of the simulated summary statistics against
  # the sample statistics.
  EucDist <- sapply(1:S, function(x) {
    sum((SumStatsZs[, x]-SumStatsY)^2)
  })
  # Sort the distances.
  I <- order(EucDist)
  selPrior <- priors[I, ][1:SelS, ]
  selDist <- EucDist[I][1:SelS]
  # Return the best priors and its euclidean distance.
  ABCpostChain <- cbind(selPrior, selDist)
  colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                              'sigma_u0', 'EucDist')
  return(ABCpostChain)
}

#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' @param theta Parameter list, must include FALTA and the cluster..
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @return Return the best `a` * `S` priors and its distance.
parallel_ABCpanelSFAEmpirical <- function(theta, SumStatsY) {
  # Parameter treatment.
  S <- theta$S
  cluster <- theta$cluster
  # Simulate the priors B0, sigma_v, sigma_v0, sigma_u, sigma_u0.
  priors <- parSapply(cluster, 1:length(theta$lb), function(i) {runif(theta$S, theta$lb[i], theta$ub[i])})
  # Simulate the residuals for each row of the prior combination.
  z <- parApply(cluster, priors, 1, function(prior) {
    simResEmpirical(theta = theta, prior = prior)})
  # Get the summary statistics of the residuals.
  # SumStatsZs <- sapply(1:S, function(x) {
  #   getSumStatsEmpirical(theta, z[,x])
  # })
  # 
  SumStatsZs <- parSapply(cluster, 1:S, function(x) {
    getSumStatsEmpirical(theta, z[,x])
  })
  
  SumStatsZs <- parSapply(cluster, 1:S, function(x) {
    getSumStatsEmpirical(theta, z[,x])
  })
  
  # Compute the square euclidean distance of the simulated summary statistics against
  # the sample statistics.
  EucDist <- parSapply(cluster, 1:S, function(x) {
    sum((SumStatsZs[, x]-SumStatsY)^2)
  })
  # Sort the distances.
  I <- order(EucDist)
  EucDistOrd <- EucDist[I]
  OrdPrior <- priors[I, ]
  selDist <- EucDistOrd[1:round(S * theta$delta)]
  SelPrior <- OrdPrior[1:round(S * theta$delta), ]
  # Return the best priors and its euclidean distance.
  ABCpostChain <- cbind(SelPrior, selDist)
  colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                              'sigma_u0', 'EucDist')
  return(ABCpostChain)
}

#' Approximate Bayesian Computation for stochastic frontier panel data model (GTRE).
#' In the urge to optimize memory use, this rutine divide the simulation in smaller
#' simulations -"batches"-.
#' @param theta Parameter list, must include FALTA and the cluster.
#' @param SumStatsY Summary statistics of the sample calculated with `getSampleSumStats()`
#' @param maxSelS Number of maximum combinations of parameters that will be returned.
#' @param SBatch Number of combination for each batch.
#' @return Return the best `a` * `S` priors and its square euclidean distance.
batchABCpanelSFAEmpirical <- function(theta, SumStatsY, SBatch = 1e5) {
  # Parameter treatment.
  S <- theta$S
  maxSelS <- theta$S_best  
  # The selected samples are maximum maxSelS.
  SelS <- min(S, maxSelS)
  
  if (S <= SBatch){
    if (theta$parallel){
      ABCpostChain <- parallel_ABCpanelSFAEmpirical(theta, SumStatsY)
    } else {
      ABCpostChain <- ABCpanelSFAEmpirical(theta, SumStatsY)
    }
    
  } else {
    # S for each "batch" (iteration).
    SBatches <- c(rep(SBatch, floor(S / SBatch)), S %% SBatch)
    SBatches <- SBatches[SBatches != 0]
    # Preallocation.
    EucDist <- rep(Inf, SelS)
    priors <- matrix(rep(0, SelS*5), ncol = 5)
    thetaBatch <- theta
    
    # Make ABC for each batch and update the best distances.
    cont <- 0
    for (s in SBatches){
      cat(paste('Batch', cont, '/', length(SBatches) , ' time: ', Sys.time(), '\n'))
      # Modify the delta to return up to maxSelS combinations.
      thetaBatch$delta <- min(1, maxSelS / s)
      # Modify the S.
      thetaBatch$S <- s
      # Run ABC for the batch.
      if (theta$parallel){
        batchABCpostChain <- parallel_ABCpanelSFAEmpirical(thetaBatch, SumStatsY)
      } else {
        batchABCpostChain <- ABCpanelSFAEmpirical(thetaBatch, SumStatsY)
      }
      BatchDist <- batchABCpostChain[,6]
      BatchPriors <- batchABCpostChain[,1:5]
      
      # There is a smaller distance?
      existSmaller <- BatchDist[1] < EucDist[SelS]
      # Substitute the best combinations until there is no new better combination.
      while (existSmaller && length(BatchDist) > 0) {
        # Sustitute the smaller new distance with the largest old distance.
        EucDist[SelS] <- BatchDist[1]
        priors[SelS,] <- BatchPriors[1,]
        # Sort by distance.
        I <- order(EucDist)
        EucDist <- EucDist[I][1:SelS]
        priors <- priors[I, ][1:SelS, ]
        # Delete distance.
        BatchDist <- BatchDist[-1]
        # Update control.
        existSmaller <- BatchDist[1] < EucDist[SelS]
      }
      cont <- cont + 1
    }
    
    ABCpostChain <- cbind(priors, EucDist)
    colnames(ABCpostChain) <- c('B0', 'sigma_v', 'sigma_v0', 'sigma_u',
                                'sigma_u0', 'EucDist')
  }
  return(ABCpostChain)
}

# Error measurement -------------------------------------------------------
#'  Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
#'  A = −p × [1Ti ITi] is a matrix of dimension Ti × (Ti + 1),  where 1Ti is the
#' column vector of length Ti, and ITi is the identity matrix of dimension T.
get_A <- function(T_i, p) { 
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
#' @param p Production 1; Cost -1
ColombiExpectation <- function(n, Ts, y, X, betas, sigmas, p = 1){
  # Change balanced notation to unbalanced panel data notation.
  if (length(Ts) == 1 & length(Ts) < n) {
    Ts <- rep(Ts, n)
  }
  # Unwrap variances.
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
    # Preallocate.
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
  return(c(xinv,  overallEf))
}

# Parser
get_real_sigmas <- function(input_data, scenario) {
  sigmas <- input_data[[scenario]]$params$sigma
  
  return(sigmas)
}

get_gibbs_sigmas <- function(gibbs_results) {
  sigmas_gibbs <- c(
    mean(sqrt(gibbs_results$sigmaepsilon2chain)),
    mean(sqrt(gibbs_results$sigmaalpha2chain)),
    mean(sqrt(gibbs_results$sigmau2chain)),
    mean(sqrt(gibbs_results$sigmaeta2chain))
  ) %>% rev
  
  return(sigmas_gibbs)
}

get_abc_sigmas <- function(sim_abc) {
  sigmas_abc <- sim_abc$sumary$mean_ABC[2:5]
  
  return(sigmas_abc)
}

cor_na_rm <- function(x1, x2, na.rm = T) {
  n <- sum(!is.na(x1 + x2))
  mean(
    (x1 - mean(x1, na.rm = na.rm)) * (x2 - mean(x2, na.rm = na.rm)),
    na.rm = na.rm
  ) / 
    (sd(x1, na.rm = na.rm) * sd(x2, na.rm = na.rm)) * n / (n - 1)
}


get_metrics_efficiency <- function(efficiencies_ls, name_real, name_est) {
  
  sapply(
    efficiencies_ls,
    function(sim) {
      TE_real <- sim[[name_real]] %>% 
        apply(2, function(x) {
          x[x == 0] <- mean(x[x != 0])
          x
        })
      TE_est <- sim[[name_est]]
      t <- dim(TE_real)[2] - 1
      
      TE_real <- cbind(TE_real, TE_real[, 1] * TE_real[, 2:(t+1)])
      TE_est <- cbind(TE_est, TE_est[, 1] * TE_est[, 2:(t+1)])
      
      ts_trans <- 2:(t+1)
      ts_overall <- (t+2):(2*t+1)
      # Compute overall efficiency
      EF_real <- cbind(
        TE_real, TE_real[, 1] * TE_real[, 2:(t+1)]
      )
      EF_est <- cbind(
        TE_est, TE_est[, 1] * TE_est[, 2:dim(TE_est)[2]]
      )
      
      EF_perm_est <- EF_est[, 1]
      EF_perm_real <- EF_real[, 1]
      idx_clean <- !is.na(EF_perm_est + EF_perm_real) & 
        !is.infinite(EF_perm_est + EF_perm_real)
      EF_perm_est <- EF_perm_est[idx_clean]
      EF_perm_real <- EF_perm_real[idx_clean]
      
      
      EF_trans_est <- c(EF_est[, ts_trans])
      EF_trans_real <- c(EF_real[, ts_trans])
      idx_clean <- !is.na(EF_trans_est + EF_trans_real) &
        !is.infinite(EF_trans_est + EF_trans_real)
      EF_trans_est <- EF_trans_est[idx_clean]
      EF_trans_real <- EF_trans_real[idx_clean]
      
      EF_overall_est <- c(EF_est[, ts_overall])
      EF_overall_real <- c(EF_real[, ts_overall])
      idx_clean <- !is.na(EF_overall_est + EF_overall_real) &
        !is.infinite(EF_overall_est + EF_overall_real)
      EF_overall_est <- EF_overall_est[idx_clean]
      EF_overall_real <- EF_overall_real[idx_clean]
      
      c(
        relative_bias_perm = mean(EF_perm_est / EF_perm_real - 1),
        relative_bias_trans = mean(EF_trans_est / EF_trans_real - 1),
        relative_bias_overall = mean(EF_overall_est / EF_overall_real - 1),
        
        cor_perm = cor(EF_perm_est, EF_perm_real),
        cor_trans = cor(EF_trans_est, EF_trans_real),
        cor_overall = cor(EF_overall_est, EF_overall_real),
        
        rootmse_perm = mean((EF_perm_est - EF_perm_real)^2) %>% sqrt,
        rootmse_trans = mean((EF_trans_est - EF_trans_real)^2) %>% sqrt,
        rootmse_overall = mean((EF_overall_est - EF_overall_real)^2) %>% sqrt
      )
    }
  ) %>% apply(1, function(x) {mean(x, na.rm = T)})
}
# Metropolis-Hasting ------------------------------------------------------
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
  Lambda <- tun * inputs_MH$Lambda
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


# Bendito seas ------------------------------------------------------------

fix_inverted_ids <- function(df) {
  df %>% 
    left_join(
      readxl::read_excel('abc4sfa/1_simulation_resources/Data/Inputs/Badunenko&Kumbhakar.xlsx') %>%
        select(scenario = ID, ID_inverted),
      by = 'scenario'
    ) %>% select(-scenario) %>% relocate(scenario = ID_inverted)
}


