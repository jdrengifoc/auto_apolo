rm(list = ls())
set.seed(010101)
source('Code/requirements.R')
source('Code/functions.R')
# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 100, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 100,                 # Number of priors.
              delta = 0.1,             # Selected priors proportion.
              nsim = 100,                # Number of simulations.
              ID = '1.1',              # Experiment ID.
              parallel = F,            # Parallel computation?
              batch = F
)

filename <- "Data/ABC_100.RData"
source('Code/apolo_sim.R')

# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 100, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 100000,                 # Number of priors.
              delta = 0.1,             # Selected priors proportion.
              nsim = 1000,                # Number of simulations.
              ID = '1.1',              # Experiment ID.
              parallel = F,            # Parallel computation?
              batch = F
)

filename <- "Data/ABC_1k.RData"
source('Code/apolo_sim.R')


# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 100, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 100000,                 # Number of priors.
              delta = 0.1,             # Selected priors proportion.
              nsim = 10000,                # Number of simulations.
              ID = '1.1',              # Experiment ID.
              parallel = F,            # Parallel computation?
              batch = F
)

filename <- "Data/ABC_10k.RData"
source('Code/apolo_sim.R')

# Parameters.
theta <- list(beta = c(0.3, 0.4, 0.6), # B0[a], B1[alpha], B2[1-alpha].
              sigma = rep(0.04, 4),    # v, v0[a], u, u0[eta].
              n = 100, t = 6,          # Sample size.
              b = c(2, 10),            # DGP Badunenko & Kumbhakar parameter.
              lb = c(0, rep(0.01, 4)), # Lower bounds of the priors.
              ub = c(2, rep(0.5, 4)),  # Upper bounds of the priors.
              S = 100000,                 # Number of priors.
              delta = 0.1,             # Selected priors proportion.
              nsim = 10000,                # Number of simulations.
              ID = '1.1',              # Experiment ID.
              parallel = F,            # Parallel computation?
              batch = F
)

filename <- "Data/ABC_100k.RData"
source('Code/apolo_sim.R')


pp <- readRDS('Data_apolo/ABC_par500k.RData')
