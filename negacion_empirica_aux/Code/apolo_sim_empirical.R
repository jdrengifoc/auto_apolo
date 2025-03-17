# For parallel computation.
if (theta$parallel) {
  n <- theta$n; #k <- theta$k;
  
  # Begin parallel cluster.
  theta$cluster <- makeCluster(detectCores(), type = "SOCK")
  registerDoParallel(theta$cluster)
  clusterExport(theta$cluster, list("parallel_ABCpanelSFAEmpirical",
                                    "getSampleSumStatsEmpirical", "repEachVec",
                              "getSumStatsEmpirical", "simResEmpirical",
                              "theta", "n", "data")) # , "k"))
}

# Initialize data to be save.
if (file.exists(filename)){
  ABC_results <- readRDS(filename)
} else {
  ABC_results <- list()
}


# Initialize experiment path and save parameters.
if (theta$parallel){
  path <- paste0(theta$ID, '.p')
} else {
  path <- theta$ID
}
ABC_results[[path]] <- list()
ABC_results[[path]][['params']] <- theta




# Message to the user.
cat(paste('Empirical aplication ', theta$ID, ' has began at:', Sys.time(), '\n'))

# Get data
X <- data$X
theta$X <- X
ABC_results[[path]][['X']] <- X
y <- data$y
theta$y <- y
loc <- data$loc
theta$loc <- loc

# Get data summary statistics.
Ystats <- getSampleSumStatsEmpirical(theta, y, X)
# For parallel computation.
if (theta$parallel) {
  # Register current simulation data.
  clusterExport(theta$cluster, list("X", "y", "Ystats", "theta"))
}

# Run ABC4SFA algorithm.
if (theta$batch) {
  ABCpostChain <- batchABCpanelSFAEmpirical(theta, Ystats)
} else if (theta$parallel) {
  ABCpostChain <- parallel_ABCpanelSFA(theta, Ystats)
} else{
  ABCpostChain <- ABCpanelSFA(theta, Ystats)
}

# Compute
mean_ABC <- apply(ABCpostChain, 2, mean)
median_ABC <- apply(ABCpostChain, 2, quantile, probs = 0.5)
lb_ABC <- apply(ABCpostChain, 2, quantile, probs = 0.025)
ub_ABC <- apply(ABCpostChain, 2, quantile, probs = 0.975)

# Save simulation results.
ABC_results[[path]][['ABCpostChain']] <- ABCpostChain
ABC_results[[path]][['sumary']] <- data.frame(mean_ABC, median_ABC, lb_ABC, ub_ABC)
ABC_results[[path]][['y']] <- y

# Message to the user.
cat(paste('Empirical aplication ', theta$ID, ' has concluded at:', Sys.time(), '\n'))
  
# Save data.
saveRDS(ABC_results, filename)

# For parallel computation.
if (theta$parallel) {
  # End parallel cluster.
  stopCluster(theta$cluster)
}
