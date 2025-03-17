#' Tareas
#' 1. Corregir B0.
#' 2. Pensarlo para datos de panel (between).
#' 3. Buscar datos de Kumba.
#' 4. Añadir random effects (pooling). Listo pero no cambia nada :p
#' 5. Calcular las ineficiencias individuales
#' Preguntas/Resultados
#' 1. ¿Corregir adentro o afuera?
#' 2. No sirvió.
#' 3. Semilisto
#' 4. No se cambió porque da igual y es más lento.
#' 5a. ¿cómo calcular sigma? ¿suma (pp. 4) o asíntótico (pp.11)}
#'    - ¿es igual en panel data?
#' 5b. ¿cómo calcular sigmau y sigmav?
#' 5c. ¿cómo calcular los residuales?

rm(list = ls())
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')

SFM_all <- function(X, theta, dist = "halfnormal", name = "ln"){
  # Parameter treatment.
  dist = tolower(dist)
  name = tolower(name)
  n = theta$n
  t = theta$t
  
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

# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 50, t = 6,           # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 1e3,                 # Number of priors.
              delta = 100 / 1e3,             # Selected priors proportion.
              nsim = 1,                # Number of simulations.
              ID = '1.2',              # Experiment ID.
              parallel = T,            # Parallel computation?
              batch = F
)

n <- theta$n; k <- theta$k;

# Begin parallel cluster.
theta$cluster <- makeCluster(detectCores(), type = "SOCK")
registerDoParallel(theta$cluster)
clusterExport(theta$cluster, list("parallel_ABCpanelSFA", "getSampleSumStats",
                                  "getSumStats", "simRes", "SFM", "genData_BK",
                                  "theta"))

X <- genData_BK(theta)
sfm <- SFM_all(X, theta)
##
sd(sfm$v)
##
y <- sfm$y
Ystats <- getSampleSumStats(theta, y, X)

clusterExport(theta$cluster, list("X", "y", "Ystats"))
ABCpostChain <- parallel_ABCpanelSFA(theta, Ystats)
# Fix intercept.
fix_intercept <- function(theta, beta0, sigma_u0, sigma_u){
  if (theta$t > 1){
    # Panel data.
    B0_fix <- beta0 + sqrt(2 / pi) * (sigma_u0 + sigma_u)
  } else {
    # Cross-sectional data.
    B0_fix <- beta0 + sqrt(2 / pi) * (sigma_u0)
  }
  return(B0_fix)
}
fix_intercept(theta, ABCpostChain[1,1], ABCpostChain[1,5], ABCpostChain[1,4])


# Mean efficiency
meanTE <- function(sigma){
  return(2*(1 - pnorm(sigma)) * exp(sigma^2 / 2))
}
meanTE(ABCpostChain[1,4:5])

# Individual efficiency.
indTE <- function(sigma_u, sigma_v, epsilon){
  sigma <- sqrt(sigma_u^2 + sigma_v^2)
  sigma_aux <- sigma_u * sigma_v / sigma
  # Varies across observations.
  mu_aux <- -epsilon * (sigma_u / sigma)^2
  return((1 - pnorm(1 - mu_aux/sigma_aux)) / (1 - pnorm(- mu_aux/sigma_aux)) *
         exp(-mu_aux + 0.5*sigma_aux^2))
  
}

# Residuals.
reg <- lm(y ~ X[,2] + X[,3])
B0_est <- fix_intercept(theta, ABCpostChain[1,1], ABCpostChain[1,5],
                        ABCpostChain[1,4])
beta_OLS <- as.numeric(reg$coefficients[-1])
beta_est <- c(B0_est, beta_OLS)
data.frame(real = theta$beta, estimated = beta_est,
           error = theta$beta - beta_est)
est_res <- y - X%*%beta_est

## Different ways to compute residuals.
real_res <- y - X%*%theta$beta
real_res2 <- sfm$v + sfm$v0 - sfm$u - sfm$u0
real_res2 <-  theta$beta[1] + sfm$v + sfm$v0 - (sfm$u - mean(sfm$u)) - (sfm$u0 - mean(sfm$u0))
##

# Transient and permanent residuals.
real_resP <- sfm$v0 - sfm$u0
real_resT <- sfm$v - sfm$u

# Technical efficiency.
est_OTE <- indTE(est_res, sum(ABCpostChain[1,1:2]),  sum(ABCpostChain[1,3:4]))


real_PTE <- indTE(real_resP, theta$sigma[2],  theta$sigma[4])
real_TTE <- indTE(real_resT, theta$sigma[1],  theta$sigma[3])
real_OTE <- real_PTE * real_TTE

real_OTE1 <- indTE(real_res2, sum(theta$sigma[1:2]),  sum(theta$sigma[3:4]))
real_OTE2 <- indTE(real_res2, prod(theta$sigma[1:2]),  prod(theta$sigma[3:4]))

plot(real_OTE, real_OTE2, xlim = c(0, 1), ylim = c(0, 1), main = "Overall")
abline(a=0, b=1) 

summary(real_OTE)

meanTE(sum(theta$sigma[3:4]))
meanTE(sum(theta$sigma[3])) * meanTE(sum(theta$sigma[4]))

# Exploring panel data ----------------------------------------------------


getSampleSumStats_rand <- function(theta, y, X){
  # Parameters treatment.
  t = theta$t
  n = theta$n
  if (t > 1) {
    # Random effects residuals.
    DataPanel <- data.frame(cbind(rep(1:n, each = t), rep(1:t, times = n), y, X[,2:3]))
    names(DataPanel) <- c("id", "time", "y", "x1", "x2")
    reg <- plm(y ~ x1 + x2, data = DataPanel, index = c("id", "time"),
               model = "random")
  } else {
    reg <- lm(y ~ X[, 2] + X[, 3])
  }
  resOLS <- reg$residuals
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
SFM_all <- function(X, theta, dist = "halfnormal", name = "ln"){
  # Parameter treatment.
  dist = tolower(dist)
  name = tolower(name)
  n = theta$n
  t = theta$t
  
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


# Simulate data
X <- genData_BK(theta, F)
sfm <- SFM_all(X, theta)
y <- SFM(X, theta)
y <- sfm$y
epsilon <- sfm$v0 + sfm$v - sfm$u + sfm$u0


# Get data summary statistics.
getSampleSumStats(theta, y, X)
getSampleSumStats_rand(theta, y, X)

### OLS vs random effect - Estimation
DataPanel <- data.frame(cbind(rep(1:theta$n, each = theta$t), rep(1:theta$t, times = theta$n), y, X[,2:3]))
colnames(DataPanel) <- c("id", "time", "y", "x1", "x2")

reg <- lm(y ~ X[, 2] + X[, 3])
reg_random <- plm(y ~ x1 + x2, data = DataPanel, index = c("id", "time"), model = "random")
reg_between <- plm(y ~ x1 + x2, data = DataPanel, index = c("id", "time"), model = "between")
reg_within <- plm(y ~ x1 + x2, data = DataPanel, index = c("id", "time"), model = "within")



res_ols <- reg$residuals
res_rand <- as.numeric(reg_random$residuals)
res_bet <- reg_between$residuals
res_bet <- as.numeric(rep(res_bet, each = theta$t))
res_wit <- as.numeric(reg_within$residuals)

epsilon_it <- res_wit - res_bet
epsilon_0 <- res_rand - epsilon_it
epsilon_0 <- matrix(epsilon_0, nrow = theta$t)
epsilon_0 <- apply(epsilon_0, 2, mean)

plot(sort(epsilon), sort(res_rand))
plot(sort(sfm$v0 - sfm$u0), sort(epsilon_0))
plot(sort(sfm$v - sfm$u), sort(epsilon_it))


### OLS vs random effect - Computation time.
tic <- Sys.time()
for (i in 1:100) {
  X <- genData_BK(theta, F)
  y <- SFM(X, theta)
  getSampleSumStats_rand(theta, y, X)
}
Sys.time() - tic

tic <- Sys.time()
for (i in 1:100) {
  X <- genData_BK(theta, F)
  y <- SFM(X, theta)
  getSampleSumStats(theta, y, X)
}
Sys.time() - tic
###
