
id <- 's4'
sim <- 'sim2'
outputData <- readRDS(sprintf("Data/Outputs/BK_%s.RData", id))
inputData <- readRDS("Data/Inputs/BK_simData.RData")
X <- readRDS("Data/Inputs/BK_simData.RData")[['X']]
y <- inputData[[id]][[sim]]$y
est_ABCparams <- outputData[[id]][[sim]]$ABCpostChain[1,]
est_betas <- c(est_ABCparams[1],
               lm(y ~ .,  data = data.frame(y, X[,2:3]))$coefficients[-1])
est_sigmas <- est_ABCparams[2:5]

betas <- inputData[[id]]$params$beta
sigmas <- inputData[[id]]$params$sigma
Ts <- inputData[[id]]$params$t
n <- inputData[[id]]$params$n

TI_real <- ColombiExpectation(n, Ts, y, X, betas, sigmas, p = 1)
TI_est <- ColombiExpectation(n, Ts, y, X, est_betas, est_sigmas, p = 1)
w
epsilon <- TI_est[[2]] - TI_real[[2]]
mean(epsilon^2)
(Ts-1)*var(epsilon)/Ts + mean(epsilon)^2



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
  V <- diag(c(sigma_u0, rep(sigma_u, T_i)))
  return(V)
}

#' Auxiliar for estimation of technical efficiency of Colombi et al. (2014).
get_Sigma <- function(sigma_v, sigma_v0, T_i){
  Sigma <- sigma_v * diag(T_i) + sigma_v0 * rep(1, T_i)%*%t(rep(1,T_i))
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
    ExpectedTI_i <- rep(0, T_i)
    for (t_i in 1:T_i){
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

#' A = −p × [1Ti ITi] is a matrix of dimension Ti × (Ti + 1), 
#' where 1Ti is the column vector of length Ti, and ITi is the
#' identity matrix of dimension Ti.
#' p is 1 if c production frontier (TI) p is -1 if cost function (TE).
#' @param ts vector de periodos por individuo,
