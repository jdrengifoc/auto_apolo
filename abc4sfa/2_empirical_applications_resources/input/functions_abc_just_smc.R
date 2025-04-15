library(furrr)
library(future)
library(progressr)
library(plm)
library(fdrtool)
library(EasyABC)
library(coda)
library(magrittr)
library(foreach)

# Get time obs for each individual
get_Ti <- function(data) {
  ids <- unique(data$id)
  N <- length(ids)
  Ti <- vector("list", N)
  
  i <- 1
  for(id in ids){
    Ti[[i]] <- sort(data$time[data$id == id])
    i <- i + 1
  }
  
  return(Ti)
}

# Define the function to read and process the data
read_and_prepare_data <- function(data_path, formula_spec) {
  data <- read.csv(data_path, header = TRUE)
  attach(data)
  formula <- as.formula(formula_spec)
  list(data = data, formula = formula)
}

# Function to calculate OLS regression
calculate_ols <- function(data, formula) {
  OLSreg <- lm(formula, data = data)
  summary(OLSreg)
  OLSreg
}

# Function to calculate between regression
calculate_between_regression <- function(data, formula) {
  regBetween <- plm::plm(formula, data = data, index = c("id", "time"), model = "between")
  summary(regBetween)
  regBetween
}

# Function to generate synthetic data for SFA
SFpanel <- function(theta, dist = 'halfnormal', model = 'cost', Ti) {
  N <- Ti %>% length()
  p <- ifelse(model == 'cost', 1, -1)
  N.T <- Ti %>% unlist() %>% length()
  v <- rnorm(N.T, 0, theta[2])
  a <- rnorm(N, 0, theta[3])
  if (dist == "halfnormal") {
    eta <- fdrtool::rhalfnorm(N, theta = sqrt(pi / 2) / theta[4])
    u <- fdrtool::rhalfnorm(N.T, theta = sqrt(pi / 2) / theta[5])
    rep_a <- sapply(1:N, function(i) rep(a[i], each = length(Ti[[i]]))) %>% unlist %>% as.vector()
    rep_eta <- sapply(1:N, function(i) rep(eta[i], each = length(Ti[[i]]))) %>% unlist() %>% as.vector()
    Resz <- theta[1] + v + rep_a + p * (u - mean(u)) + p * (rep_eta - mean(eta))
  } else {
    eta <- rexp(N, rate = 1 / theta[4])
    u <- rexp(N.T, rate = 1 / theta[5])
    rep_a <- sapply(1:N, function(i) rep(a[i], each = length(Ti[[i]]))) %>% unlist() %>% as.vector()
    rep_eta <- sapply(1:N, function(i) rep(eta[i], each = length(Ti[[i]]))) %>% unlist() %>% as.vector()
    Resz <- theta[1] + v + rep_a + p * (u - mean(u)) + p * (rep_eta - mean(eta))
  }
  return(Resz)
}

# Function to compute summary statistics for OLS residuals
SumStatOLS <- function(resS, Ti, K) {
  N.T <- Ti %>% unlist() %>% length()
  N <- length(Ti)
  m1 <- mean(resS)
  s2TotS <- var(resS) * (N.T - 1) / (N.T - K)
  rep_n <- sapply(1:N, function(i) rep(i, each = length(Ti[[i]]))) %>% unlist() %>% as.vector()
  rep_t <- Ti %>% unlist()
  ResAPs <- cbind(rep_n, rep_t, resS - m1)
  SUMA <- NULL
  for(i in 1:N){
    observed_years <- Ti[[i]]
    Resit <- ResAPs[ResAPs[, 1] == i, 3]
    names(Resit) <- observed_years
    
    for(t in observed_years){
      for(s in observed_years){
        if(s > t){
          if(!is.na(Resit[as.character(t)]) && !is.na(Resit[as.character(s)])){
            Sumit <- Resit[as.character(t)] * Resit[as.character(s)]
            # Append the result to SUMA
            SUMA <- c(SUMA, Sumit)
          }
        }
      }
    }
  }
  valid_pairs_count <- 0
  for(i in 1:N){
    observed_years <- Ti[[i]]
    num_years <- length(observed_years)
    valid_pairs_count <- valid_pairs_count + num_years * (num_years - 1) / 2
  }
  s2IndS <- sum(SUMA) / (valid_pairs_count - K)
  s2IdioS <- s2TotS - s2IndS
  m3 <- mean((resS - m1)^3)
  m4 <- mean((resS - m1)^4)
  return(c(m1, s2IndS, s2IdioS, m3, m4))
}

# Moment Matching con los plims de las estimaciones
get_MMols <- function(data, OLSreg, Ti, regBetween, model = 'cost') {
  function(theta) {
    p <- ifelse(model == 'cost', 1, -1)
    N <- length(Ti)
    N.T <- length(unlist(Ti))
    K <- length(OLSreg$coefficients)
    resOLS <- OLSreg$residuals
    OLSvar <- sum(resOLS^2) / (N.T - K)
    rep_n <- sapply(1:N, function(i) rep(i, each = length(Ti[[i]]))) %>% unlist() %>% as.vector()
    rep_t <- Ti %>% unlist() %>% as.vector()
    ResAPs <- cbind(rep_n, rep_t, resOLS)
    
    SUMA <- NULL
    for(i in 1:N){
      observed_years <- Ti[[i]]
      Resit <- ResAPs[ResAPs[, 1] == i, 3]
      names(Resit) <- observed_years
      
      for(t in observed_years){
        for(s in observed_years){
          if(s > t){
            if(!is.na(Resit[as.character(t)]) && !is.na(Resit[as.character(s)])){
              Sumit <- Resit[as.character(t)] * Resit[as.character(s)]
              # Append the result to SUMA
              SUMA <- c(SUMA, Sumit)
            }
          }
        }
      }
    }
    
    valid_pairs_count <- 0
    for(i in 1:N){
      observed_years <- Ti[[i]]
      num_years <- length(observed_years)
      valid_pairs_count <- valid_pairs_count + num_years * (num_years - 1) / 2
    }
    
    s2IndOLS <- sum(SUMA) / (valid_pairs_count - K)
    s2IdioOLS <- OLSvar - s2IndOLS
    
    B0 <- theta[1]
    sige <- theta[2]
    siga <- theta[3]
    sigeta <- theta[4]
    sigu <- theta[5]
    
    if (sige < 0 || siga < 0 || sigeta < 0 || sigu < 0) {
      return(Inf)
    }
    
    Eq1 <- OLSreg$coef[1] - p*(B0 + 2/pi * (sigu^2 + sigeta^2))
    Eq2 <- s2IdioOLS - (sige^2 + sigu^2 * (1 - 2/pi))
    Eq3 <- s2IndOLS - (siga^2 + sigeta^2 * (1 - 2/pi))
    
    Eq5 <- sum((regBetween$res)^2) / (length(regBetween$res) - length(regBetween$coef)) - 
      (sige^2 / mean(sapply(Ti, length)) + siga^2 + (1 - 2/pi) * (sigeta^2 + sigu^2 / mean(sapply(Ti, length))))
    
    Eq6 <- mean(regBetween$res^3) + p * sqrt(2 / pi) * (4/pi - 1) * 
      (sigu^3 / mean(sapply(Ti, function(x) length(x)^2)) + sigeta^3)
    
    Err <- Eq1^2 + Eq2^2 + Eq3^2 + Eq5^2 + Eq6^2
    
    return(Err)
  }
}


# Define the model function
get_SFpanelEasyABC <- function(Ti, K, model = 'cost') {
  function(theta) {
    N <- length(Ti)
    p <- ifelse(model == 'cost', 1, -1)
    N.T <- length(unlist(Ti))
    v <- rnorm(N.T, 0, theta[2])
    a <- rnorm(N, 0, theta[3])
    eta <- fdrtool::rhalfnorm(N, theta = sqrt(pi / 2) / theta[4])
    u <- fdrtool::rhalfnorm(N.T, theta = sqrt(pi / 2) / theta[5])
    rep_a <- sapply(1:N, function(i) rep(a[i], each = length(Ti[[i]]))) %>%
      unlist() %>% as.vector()
    rep_eta <- sapply(1:N, function(i) rep(eta[i], each = length(Ti[[i]]))) %>% 
      unlist() %>% as.vector()
    Resz <- theta[1] + v + rep_a + p * (u - mean(u)) + p * (rep_eta - mean(eta))
    SumStatZ <- SumStatOLS(Resz, Ti, K)
    return(SumStatZ)
  }
}

# Function to run ABC for SFA (AR-ABC method) in chunks
ABCpanelSFA_AR <- function(S, Hyp, a, Ti, K, OLSreg, regBetween, num_chunks = 5, model = 'cost') {
  prior <- matrix(runif(S * 5, rep(Hyp[c(TRUE, FALSE)], each = S), rep(Hyp[c(FALSE, TRUE)], each = S)), ncol = 5)
  plan(multisession, workers = availableCores() - 1)
  
  sum_stat_obs <- SumStatOLS(OLSreg$residuals, Ti, K)
  chunk_size <- S / num_chunks
  total_elapsed_time <- 0
  final_SelPrior <- NULL
  
  for (i in 1:num_chunks) {
    cat("Running AR-ABC chunk", i, "of", num_chunks, "\n")
    progressr::with_progress({
      start_time <- Sys.time()
      
      chunk_prior <- prior[((i-1) * chunk_size + 1):(i * chunk_size), ]
      z <- future_map(seq_len(nrow(chunk_prior)), ~ SFpanel(chunk_prior[.x, ], dist = "halfnormal", Ti = Ti, model = model), .options = furrr_options(seed = TRUE))
      SumStatsZs <- future_map(z, ~ SumStatOLS(.x, Ti = Ti, K = K))
      EucDist <- future_map_dbl(SumStatsZs, ~ dist(rbind(sum_stat_obs, .x)))
      
      chunk_OrdPrior <- cbind(chunk_prior[order(EucDist), ], EucDist[order(EucDist)])
      chunk_SelPrior <- chunk_OrdPrior[1:round(chunk_size * a), ]
      chunk_OrdDistSel <- sort(EucDist)[1:round(chunk_size * a)]
      
      elapsed_time <- Sys.time() - start_time
      total_elapsed_time <- total_elapsed_time + elapsed_time
      
      if (is.null(final_SelPrior)) {
        final_SelPrior <- chunk_SelPrior
      } else {
        final_SelPrior <- rbind(final_SelPrior, chunk_SelPrior)
      }
      
      # Print diagnostics
      cat("Chunk", i, "completed in", elapsed_time, "seconds\n")
      cat("Total elapsed time so far:", total_elapsed_time, "seconds\n")
      cat("Memory used:", pryr::mem_used(), "\n")
    })
  }
  
  ResultARABC <- list(SelPrior = final_SelPrior, total_elapsed_time = total_elapsed_time)
  return(ResultARABC)
}

# Define a function that runs the entire ABC pipeline for one application.
run_ABC_pipeline <- function(app_name, data_path, formula_list, prior_list, model_type_list,
                             time_limit_minutes = 60, chunk_size_AR = 500, 
                             chunk_size_MCMC = 500, chunk_size_SMC = 1000, 
                             tol = 0., folder_output = '.') {
  
  cat("===========================================\n")
  cat("Running ABC pipeline for application:", app_name, "\n")
  
  # 1. Read data
  data <- readRDS(data_path)[[app_name]]
  data <- cbind(y = data$y, data$X, id = data$loc$id, time = data$loc$time)
  
  # 2. Run OLS regression using the provided formula specification
  formula_spec <- formula_list[[app_name]]
  OLSreg <- lm(as.formula(formula_spec), data = data)
  cat("OLS regression summary:\n")
  print(summary(OLSreg))
  
  # 3. Extract time observations per individual using get_Ti
  Ti <- get_Ti(data)
  
  # 4. Compute observed summary statistics from the OLS residuals
  K <- length(OLSreg$coefficients)  # Number of parameters
  sum_stat_obs <- SumStatOLS(OLSreg$residuals, Ti, K)
  cat("Observed summary statistics:\n")
  print(sum_stat_obs)
  
  # 5. Define the ABC simulation model function
  #    get_SFpanelEasyABC returns a function that, given a theta, simulates residuals
  #    and computes the summary statistics (using SumStatOLS).
  model_type <- model_type_list[[app_name]]
  sim_model <- get_SFpanelEasyABC(Ti, K, model = model_type)
  
  # 6. Define the prior for theta (here a 5-dimensional parameter vector)
  toy_prior <- prior_list[[app_name]]
  
  # 7. Set time limit in seconds for each ABC algorithm run
  time_limit_sec <- time_limit_minutes * 60
  
  # 8. Initialize storage for results
  results_SMC  <- list()   # For SMC ABC
  
  # 11. Run SMC ABC in chunks using EasyABC's ABC_sequential
  cat("Starting SMC ABC simulation...\n")
  start_time_SMC <- Sys.time()
  while(as.numeric(difftime(Sys.time(), start_time_SMC, units = "secs")) < time_limit_sec) {
    chunk_SMC <- ABC_sequential(method = "Lenormand",
                                model = sim_model,
                                prior = toy_prior,
                                summary_stat_target = sum_stat_obs,
                                nb_simul = chunk_size_SMC,
                                alpha = 0.5,
                                p_acc_min = 0.03,
                                progress_bar = FALSE)
    results_SMC[[length(results_SMC) + 1]] <- chunk_SMC$param
    elapsed_SMC <- as.numeric(difftime(Sys.time(), start_time_SMC, units = "secs"))
    cat("SMC chunk completed. Elapsed time:", round(elapsed_SMC, 1), "seconds\n")
  }
  
  combined_SMC <- do.call(rbind, results_SMC)
  PostSMC <- coda::mcmc(combined_SMC)
  
  # 12. Save results for the application (optional)
  save(
    PostSMC, file = file.path(
      folder_output,
      sprintf("ABCpostChain_SMC_%s_time=%d_chunk_size=%d_%s_sa.RData", 
      app_name, time_limit_minutes, chunk_size_SMC, Sys.Date())
    )
  )
  
  cat("ABC pipeline completed for application:", app_name, "\n")
  return(list(PostSMC = PostSMC))
}