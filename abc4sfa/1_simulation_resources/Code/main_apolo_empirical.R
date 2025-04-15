rm(list = ls())
set.seed(010101)
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')

# Read inputs.
# ABC_inputs <- readRDS('Data/Inputs/BK_Ystats.RData')
empiricalData <- readRDS('Data/Inputs/BK_empiricalData_new.RData')

#' Parameters.
#' "swissRailWays", "spainDairy", "usElectricity", "indonesianRice"
#' "unbalanced", "balanced", "unbalanced", "balanced"
exp_ID <- 'swissRailWays'

data <- empiricalData[[exp_ID]]

S_best <- 1000
# aux_S <- getNumSamplesToGuaranteeConvergence(dim(data$loc)[1], 5, 1000)
aux_S <- 1e8

ts <- table(data$loc$id)
names(ts) <- NULL

theta <- list(S = aux_S,                # Number of priors.
              delta = S_best / aux_S,      # Selected priors proportion.
              S_best = S_best,
              parallel = T,           # Parallel computation?
              batch = T,               # Batch computation?
              ID = exp_ID,
              n = length(unique(data$loc$id)),
              ts = ts,
              b = c(2, 10),
              lb = c(0.00, 0.01, 0.01, 0.01, 0.01),
              ub = c(2.0, 0.5, 0.5, 0.5, 0.5)
)

filename <- sprintf("Data/Outputs/BK_%s.RData", theta$ID)

source('Code/apolo_sim_empirical.R')




