rm(list = ls())
set.seed(010101)
source('Code/Fixed/requirements.R')
source('Code/Fixed/functions.R')
# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 50, t = 6,          # Sample size.
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

filename <- "Data/Outputs/ABC_par1M.RData"
filename <- "Data/Outputs/JAJA.RData"
source('Code/apolo_sim.R')

# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 50, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 2e3,                 # Number of priors.
              delta = 100/1e4,             # Selected priors proportion.
              nsim = 1,                # Number of simulations.
              ID = '1.2',              # Experiment ID.
              parallel = T,            # Parallel computation?
              batch = T
)

filename <- "Data/Outputs/ABC_test_local.RData"
source('Code/apolo_sim.R')


