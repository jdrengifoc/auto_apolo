rm(list = ls())
set.seed(010101)

FOLDER <- 'auto_apolo/negacion_empirica'
source(file.path(FOLDER, 'Code/Fixed/requirements.R'))
source(file.path(FOLDER, 'Code/Fixed/functions.R'))

# Read inputs.
ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_Ystats.RData'))
# ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_exp_Ystats.RData'))

# Parameters.
ID <- s1
theta <- list(S = 1e7,                # Number of priors.
              delta = 1000 / 1e7,      # Selected priors proportion.
              parallel = T,           # Parallel computation?
              batch = T               # Batch computation?
)
theta <- c(ABC_inputs[[ID]]$params, theta)
theta$nsim <- 100
theta$model <- 'exp'

filename <- file.path(FOLDER, sprintf("Data/Outputs/BK_hn_%s_misspecification_exp.RData", theta$ID))
# filename <- file.path(FOLDER, sprintf("Data/Outputs/BK_exp_%s.RData", theta$ID))

source(file.path(FOLDER, 'Code/apolo_sim_low_memory.R'))
