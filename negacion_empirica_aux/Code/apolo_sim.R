# For parallel computation.
if (theta$parallel) {
  n <- theta$n; k <- theta$k;
  
  # Begin parallel cluster.
  theta$cluster <- makeCluster(detectCores(), type = "SOCK")
  registerDoParallel(theta$cluster)
  clusterExport(theta$cluster, list("parallel_ABCpanelSFA", "getSampleSumStats",
                              "getSumStats", "simRes", "SFM", "genData_BK",
                              "theta", "n", "k"))
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


X <- genData_BK(theta, F)
ABC_results[[path]][['X']] <- X
# Run nsim simulations.
for (i in 1:theta$nsim){
  # Message to the user.
  cat(paste('Simulation', i, 'has began at:', Sys.time(), '\n'))
  
  # Simulate data.
  y <- SFM(X, theta)
  # Get data summary statistics.
  Ystats <- getSampleSumStats(theta, y, X)
  # For parallel computation.
  if (theta$parallel) {
    # Register current simulation data.
    clusterExport(theta$cluster, list("X", "y", "Ystats"))
  }
  
  # Run ABC4SFA algorithm.
  if (theta$batch) {
    ABCpostChain <- batchABCpanelSFA(theta, Ystats)
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

  # Simulation name.
  sim_name <- paste0('sim', i)
  
  # Save simulation results.
  ABC_results[[path]][[sim_name]][['ABCpostChain']] <- ABCpostChain
  ABC_results[[path]][[sim_name]][['sumary']] <- data.frame(mean_ABC, median_ABC, lb_ABC, ub_ABC)
  ABC_results[[path]][[sim_name]][['y']] <- y
  
  # Message to the user.
  cat(paste('Simulation', i, 'has concluded at:', Sys.time(), '\n'))
  # For large simulations save data in each iteration.
  if (theta$S > 1e5){
    saveRDS(ABC_results, filename)
  }
}
# Save data.
saveRDS(ABC_results, filename)

# For parallel computation.
if (theta$parallel) {
  # End parallel cluster.
  stopCluster(theta$cluster)
}
