rm(list = ls())
set.seed(010101)
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')

# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = NULL,            # v, v0[a], u, u0[eta].
              n = 50, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              nsim = 100,              # Number of simulations.
              ID = NULL                # Experiment ID.
)

filename <- 'BK'
# Initialize data to be save.
ABC_data <- list()
ABC_Ystats <- list()

sigmas <- readxl::read_excel('Data/Inputs/Badunenko&Kumbhakar.xlsx')


X <- genData_BK(theta, F)
ABC_data[['X']] <- X

for (j in 1:nrow(sigmas)){
  theta$ID <- sigmas[[j,1]]
  theta$sigma <- as.numeric(sigmas[j,2:5])
  ABC_data[[theta$ID]][['params']] <- theta
  ABC_Ystats[[theta$ID]][['params']] <- theta
  
  for (i in 1:theta$nsim){
    sfm <- SFM_all(X, theta)
    Ystats <- getSampleSumStats(theta, sfm$y, X)
    
    ABC_data[[theta$ID]][[paste0('sim', i)]]<- sfm
    ABC_Ystats[[theta$ID]][[paste0('sim', i)]][['Ystats']] <- Ystats
  }
}

saveRDS(ABC_data, sprintf('Data/Inputs/%s_simData.RData', filename))
saveRDS(ABC_Ystats, sprintf('Data/Inputs/%s_Ystats.RData', filename))

data <- readRDS('Data/ABC_sampledata.RData')
