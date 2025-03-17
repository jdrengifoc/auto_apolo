###### Metropolis-Hastings: Normal-exponential model ######
rm(list = ls())
SFpanel <- function(theta, dist, X, beta) {
  # theta <- EnsOLS$par
  v <- rnorm(n*T, 0, theta[2])
  a <-  rnorm(n, 0, theta[3])
  if(dist == "halfnormal"){
    eta <- fdrtool::rhalfnorm(n, theta = sqrt(pi / 2) / theta[4]) # Fix theta[5]=1
    u <- fdrtool::rhalfnorm(n*T, theta = sqrt(pi / 2) / theta[5]) # Fix theta[5]=1
    Resz <- theta[1] + v + rep(a, each = T) - (u - mean(u)) 
    - (rep(eta, each = T) - mean(eta))
  }else{
    eta <- rexp(n, rate = 1/theta[4]) 
    u <- rexp(n*T, rate = 1/theta[5])
    Resz <- theta[1] + v + rep(a, each = T) - (u - mean(u)) 
    - (rep(eta, each = T) - mean(eta))
  }
  y <- X%*%beta + Resz
  return(y)
}

# Panel dimension and example data.
T <- 5; n <- 1
X <- cbind(1, matrix(runif(n*T), T, 2))
sa <- 0.04; sv <- 0.04; seta <- 0.04; su <- 0.04
betas <- c(0.3, 0.4, 0.6)
y <- SFpanel(theta = c(0.3, c(sv, sa, seta, su)), X = X, beta = betas, dist = "halfnormal")

get_epsiloni <- function(y, X, beta){
  epsiloni <- y - X%*%beta
  return(epsiloni)
}

r <- get_epsiloni(y = y, X = X[1:T,], beta = betas)

# Inputs
T <- 5
A <- cbind(1, diag(T))
r <- get_epsiloni(y = y, X = X[1:T,], beta = betas) 
sa <- 0.04; sv <- 0.04; seta <- 0.04; su <- 0.04
it <- rep(1, T)
Sigma <- sv^2 * diag(T) + sa^2*it%*%t(it) 
V <- diag(c(seta^2, rep(su^2, T)))
iSigma <- solve(Sigma)
tA <- t(A)
Lambda <- solve(solve(V) + tA%*%iSigma%*%A)
R <- Lambda%*%tA%*%iSigma
# Proposal. Must have the same support (0, inf).
# permanete, transitorias en cada T
set.seed(1)
u <- tmvtnorm::rtmvnorm(1, mean = c(R%*%r), sigma = Lambda,
                        lower = rep(0, length = T+1), 
                        upper = rep(Inf, length = T+1))

uDrawMH <- function(u, A, r, Sigma, V) {
  # iSigma <- solve(Sigma)
  # tA <- t(A)
  # Lambda <- solve(solve(V)+tA%*%iSigma%*%A)
  # R <- Lambda%*%tA%*%iSigma
  uc <- tmvtnorm::rtmvnorm(1, mean = c(R%*%r), sigma = Lambda,
                lower=rep(0, length = T+1), 
                upper=rep(Inf, length = T+1))
  quc <- tmvtnorm::dtmvnorm(uc, mean = c(R%*%r), sigma = Lambda,
                           lower = rep(0, length = T+1), 
                           upper = rep(Inf, length = T+1))
  qu <- tmvtnorm::dtmvnorm(u, mean = c(R%*%r), sigma = Lambda,
                            lower = rep(0, length = T+1), 
                            upper = rep(Inf, length = T+1))
  fuc <- -0.5*t(r - A%*%c(uc))%*%iSigma%*%(r - A%*%c(uc)) - t(c(seta, rep(su, T)))%*%t(uc)
  fu <- -0.5*t(r - A%*%c(u))%*%iSigma%*%(r - A%*%c(u)) - t(c(seta, rep(su, T)))%*%t(u)
  # Criterio de transición. Cómo actualizamos la cadena.
  alpha <- min(exp(fuc-fu) * qu/quc, 1)
  unif <- runif(1, 0, 1)
  if(unif < alpha){
    unew <- uc
    accept <- 1
  }else{
    unew <- u
    accept <- 0
  }
  return(list(draws = unew, accept = accept))
}
##
for (scenario in scenarios) {
  X
  varianzas
  inputs_MH_reales
  # samplear por empresa
  for (empresa in empresas) {
    inputs_MH
    for (s in Ss) {
      update_MH # cada 100
    }
    # S / 100 draws
    ef <- colMeans(apply(-UniDraws, 2, exp))
  }
  EF_reales # n x t+1
  for (sim in sims) {
    y
    perturbaciones
    abc_varianzas
    residuales
    # samplear por empresa
    for (empresa in empresas) {
      inputs_MH
      for (s in Ss) {
        update_MH # cada 100
      }
      # S / 100 draws
      ef <- colMeans(apply(-UniDraws, 2, exp))
    }
    EF # n x t+1
  }
}
##
uDrawMH(u = u, A = A, r = r, Sigma = Sigma, V = V)

S <- 100
Draws <- matrix(0, S, T + 1)
Accepted <- rep(0, S)
for(s in 1:S){
  #' u: inicial
  #' A: constante
  #' r: residuales o perturbaciones.
  #' Sigmas: estimadas o poblacionales.
  #' V: estimadas o poblacionales.
  Res <- uDrawMH(u = u, A = A, r = r, Sigma = Sigma, V = V)
  u <- Res$draws
  Draws[s,] <- u 
  Accepted[s] <- Res$accept
}
Draws
# Tasa de aceptacion.
mean(Accepted)

# Mixing properties.
UniDraws <- apply(Draws, 2, unique)
EF <- colMeans(apply(-UniDraws, 2, exp))
EF
colMeans(UniDraws)
var(UniDraws)
chol(var(UniDraws))
t(R%*%r)
chol(Sigma)
chol(V)


abc_outputs <- readRDS('After/ABC4SFAv2/Data/Outputs/BK_s13.RData')
