rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))

FOLDER <- 'abc4sfa/1_simulation_resources'
# Read inputs.
# ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_Ystats.RData'))
ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_exp_Ystats.RData'))

# Parameters.
ID <- "__id__"
theta <- list(
    S = 1e7,                # Number of priors.
    delta = 1000 / 1e7,      # Selected priors proportion.
    parallel = T,           # Parallel computation?
    batch = T               # Batch computation?
    )
theta <- c(ABC_inputs[[ID]]$params, theta)
theta$nsim <- 100
theta$model <- 'hn'

filename <- file.path(
    FOLDER, 
    sprintf("Data/Outputs/BK_exp_%s_misspecification_hn_%s.RData", theta$ID, Sys.Date())
    )
# filename <- file.path(FOLDER, sprintf("Data/Outputs/BK_exp_%s.RData", theta$ID))

source(file.path(FOLDER, 'Code/apolo_sim_low_memory.R'))
