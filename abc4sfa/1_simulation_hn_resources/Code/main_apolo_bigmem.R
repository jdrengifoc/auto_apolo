rm(list = ls())
set.seed(010101)
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')

# Read inputs.
ABC_inputs <- readRDS('Data/Inputs/BK_Ystats.RData')

# Parameters.
ID <- 's16'
theta <- list(S = 1e7,                # Number of priors.
              delta = 100 / 1e7,      # Selected priors proportion.
              parallel = T,           # Parallel computation?
              batch = T               # Batch computation?
)
theta <- c(ABC_inputs[[ID]]$params, theta)
theta$nsim <- 100 # MAX 100.

filename <- sprintf("Data/Outputs/BK_%s.RData", theta$ID)

source('Code/apolo_sim_low_memory.R')
