rm(list = ls())
set.seed(010101)
source('Code/requirements.R')
source('Code/functions.R')

ABC_data <- readRDS('Data/Inputs/ABC_samplestats.RData')
# Parameters.
ID <- 's1'                      # Experiment ID
theta <- list(S = 10000,         # Number of priors.
              delta = 0.01,      # Selected priors proportion.
              parallel = T,     # Parallel computation?
              batch = F
)

theta <- c(ABC_data[[ID]]$params, theta)
filename <- "Data/Outputs/ABC_test1.RData"
source('Code/apolo_sim.R')

