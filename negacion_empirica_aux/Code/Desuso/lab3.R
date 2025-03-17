source('Code/Fixed/requirements.R')

get_metrics <- function(data, scenario, rand_mean = T, best = F){
  set.seed(0001111)
  # Parmeters
  params <- data$params
  sigma_u <- params$sigma[3]
  sigma_u0 <- params$sigma[4] 
  nsim <- params$nsim
  nsim <- length(names(scenario))
  n <- params$n
  t <- params$t
  
  # Initialize.
  RB <- rep(0, 3); UB <- RB; PCC <- RB; RMSE <- RB;
  SIGMA_U <- rep(0, nsim)
  SIGMA_U0 <- rep(0, nsim)
  DIST <- rep(0, nsim)
  MEAN_SIGMA_U <- rep(0, nsim)
  MEAN_SIGMA_U0 <- rep(0, nsim)
  
  for (i in 1:nsim){
    name <- names(scenario)[i]
    nsamples <- dim(scenario[[name]]$ABCpostChain)[1]
    # Initiallize
    est_u <- NULL; est_u0 <- est_u
    if (rand_mean){
      for (j in 1:nsamples) {
        if (best) {
          est_sigma_u <- sigma_u
          est_sigma_u0 <- sigma_u0
        } else {
          est_sigma_u <- scenario[[name]]$ABCpostChain[j, 4]
          est_sigma_u0 <- scenario[[name]]$ABCpostChain[j, 5]
        }
        
        # Estimated transitory inefficiency.
        est_u <- cbind(
          est_u, fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*est_sigma_u^2)))
        )
        # Estimated permanent inefficiency.
        est_u0 <- cbind(
          est_u0, rep(
            fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*est_sigma_u0^2))), each = t)
        )
      }
      # Mean by rows.
      est_u <- apply(est_u, 1, mean)
      est_u0 <- apply(est_u0, 1, mean)
    } else {
      if (best) {
        est_sigma_u <- sigma_u
        est_sigma_u0 <- sigma_u0
      } else {
        est_sigma_u <- mean(scenario[[name]]$ABCpostChain[, 4])
        est_sigma_u0 <- mean(scenario[[name]]$ABCpostChain[, 5])
      }
      # Estimated transitory inefficiency.
      est_u <- fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*est_sigma_u^2)))
      # Estimated permanent inefficiency.
      est_u0 <- rep(
        fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*est_sigma_u0^2))), each = t)
    }
    
    
    u <- data[[name]]$u
    u0 <- rep(data[[name]]$u0, each = t)
    
    est_TE <- exp(-est_u)
    est_PE <- exp(-est_u0)
    est_OE <- est_PE * est_TE
    
    TE <- exp(-u)
    PE <- exp(-u0)
    OE <- PE * TE
    
    RB <- RB + c(mean(est_PE / PE - 1), mean(est_TE / TE - 1), mean(est_OE / OE - 1))
    UB <- UB + c(mean(est_PE > PE), mean(est_TE > TE), mean(est_OE > OE))
    # PCC <- PCC + c(cor(est_PE, PE), cor(est_TE, TE), cor(est_OE, OE))
    RMSE <- RMSE + c(mean((est_PE / PE - 1)^2), mean((est_TE / TE - 1)^2),
                                                     mean((est_OE / OE - 1)^2))
  }
  
  # df <- data.frame(RelativeBias = RB, UpperBias = UB, PearsonCor = PCC) / nsim
  df <- data.frame(Inefficiency = c('Persistent', 'Transient', 'Overall'),
                   RelativeBias = RB / nsim, UpperBias = UB / nsim,
                   RelativeMSE = sqrt(RMSE / nsim))
  return(df)
  
}

est_results <- NULL; best_results <- NULL; results <- NULL
ID <- paste0('s', 1:16)
# s1: 76, s3:  69
for (id in ID){
  scenario <- readRDS(sprintf("Data/Outputs/BK_%s.RData", id))
  data <- readRDS('Data/Inputs/BK_simData.RData')[[names(scenario)]]
  scenario <- scenario[[1]]
  
  
  est_results <- rbind(
    est_results, cbind(id, get_metrics(data, scenario, rand_mean = T, best = F))
  )
  
  best_results <- rbind(
    best_results, cbind(id, get_metrics(data, scenario, rand_mean = F, best = T))
  )
  results <- rbind(
    results, 
    cbind(type = 'best', id, get_metrics(data, scenario, rand_mean = F, best = T)),
    cbind(type = 'est', id, get_metrics(data, scenario, rand_mean = T, best = F))
  )
}

writexl::write_xlsx(results, 'Results_sims.xlsx')
data$params$sigma[3:4]

# Lab ---------------------------------------------------------------------


# Numero de parámetros
params <- data$params
sigma_u <- params$sigma[3]
sigma_u0 <- params$sigma[4] 
nsim <- params$nsim
n <- params$n
t <- params$t

# Preallocation.

RB <- 0; UB <- 0; PCC <- 0; PCC1 <- 0
SIGMA_U <- rep(0, nsim)
SIGMA_U0 <- rep(0, nsim)
DIST <- rep(0, nsim)
MEAN_SIGMA_U <- rep(0, nsim)
MEAN_SIGMA_U0 <- rep(0, nsim)

# Mean sigmas.
for (i in 1:nsim){
  name <- names(scenario)[i]
  est_sigma_u <- scenario[[name]]$ABCpostChain[1, 4]
  est_sigma_u0 <- scenario[[name]]$ABCpostChain[1, 5]
  
  SIGMA_U[i] <- est_sigma_u
  SIGMA_U0[i] <- est_sigma_u0
  DIST[i] <- scenario[[name]]$ABCpostChain[1, 6]
  MEAN_SIGMA_U[i] <- mean(scenario[[name]]$ABCpostChain[, 4])
  MEAN_SIGMA_U0[i] <- mean(scenario[[name]]$ABCpostChain[, 5])
  est_sigma_u <- MEAN_SIGMA_U[i]
  est_sigma_u0 <- MEAN_SIGMA_U0[i]
  # Estimated transitory inefficiency.
  est_u <- fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*sigma_u^2)))
  # Estimated permanent inefficiency.
  est_u0 <- rep(
    fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*sigma_u0^2))), each = t)
  
  
  u <- data[[name]]$u
  u0 <- rep(data[[name]]$u0, each = t)
  
  est_TE <- exp(-est_u)
  #est_TE <- est_u
  est_PE <- exp(-est_u0)
  #est_PE <- est_u
  est_OE <- est_PE * est_TE
  
  TE <- exp(-u)
  # TE <- u
  PE <- exp(-u0)
  # PE <- u0
  OE <- PE * TE
  
  RB <- RB + mean(est_OE / OE - 1)
  UB <- UB + mean(est_OE > OE)
  # PCC <- PCC + sum((est_OE - mean(est_OE)) * (OE - mean(OE))) / sqrt(sum((est_OE - mean(est_OE))^2) * sum((OE - mean(OE))^2))
  PCC <- PCC + cor(est_OE, OE)
  PCC1 <- PCC1 + cor(sort(est_OE), sort(OE))
  # PCC2[i] <- cor(sort(est_OE), sort(OE))
}
# Mean nsim randoms.
for (i in 1:nsim){
  name <- names(scenario)[i]
  nsamples <- dim(scenario[[name]]$ABCpostChain)[1]
  
  est_u <- NULL; est_u0 <- est_u
  for (j in 1:nsamples) {
    est_sigma_u <- scenario[[name]]$ABCpostChain[j, 4]
    est_sigma_u0 <- scenario[[name]]$ABCpostChain[j, 5]
    # Estimated transitory inefficiency.
    est_u <- cbind(
      est_u, fdrtool::rhalfnorm(n*t, theta = sqrt((pi-2)/(2*est_sigma_u^2)))
      )
    
    # Estimated permanent inefficiency.
    est_u0 <- cbind(
      est_u0, rep(
        fdrtool::rhalfnorm(n, theta = sqrt((pi-2)/(2*est_sigma_u0^2))), each = t)
      )
  }
  est_u <- apply(est_u, 1, mean)
  est_u0 <- apply(est_u0, 1, mean)
  
  u <- data[[name]]$u
  u0 <- rep(data[[name]]$u0, each = t)
  
  est_TE <- exp(-est_u)
  #est_TE <- est_u
  est_PE <- exp(-est_u0)
  #est_PE <- est_u
  est_OE <- est_PE * est_TE
    
  TE <- exp(-u)
  # TE <- u
  PE <- exp(-u0)
  # PE <- u0
  OE <- PE * TE
  
  RB <- RB + mean(est_OE / OE - 1)
  UB <- UB + mean(est_OE > OE)
  # PCC <- PCC + sum((est_OE - mean(est_OE)) * (OE - mean(OE))) / sqrt(sum((est_OE - mean(est_OE))^2) * sum((OE - mean(OE))^2))
  PCC <- PCC + cor(est_OE, OE)
  PCC1 <- PCC1 + cor(sort(est_OE), sort(OE))
  # PCC2[i] <- cor(sort(est_OE), sort(OE))
}

c(RelativeBias = RB, UpperBias = UB, PearsonCor = PCC, SortPearsonCor = PCC1) / nsim

plot(sort(OE), sort(est_OE))
plot(PCC2)
