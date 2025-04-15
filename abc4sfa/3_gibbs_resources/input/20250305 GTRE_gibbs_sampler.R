# Read the data and run the sampler using the correct data variable
# data <- read_csv('C:/Users/migue/Universidad EAFIT/ABC4SFA/USbanksPanel.csv')
# Ensure beta0 is defined appropriately (e.g., a vector of zeros of the right length)
postChainBanks <- get_posterior(data, n_samples = 15000, burnin = 5000)
saveRDS(postChainBanks, 'postChainBanks2.RData')
