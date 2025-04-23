setwd('auto_apolo/')
rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))
library(mvtnorm)
library(stringr)
library(dplyr)


# Results halfnormal-normal GTRE ---------------------------------------------
FOLDER <- 'abc4sfa/1_simulation_resources'
ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_Ystats.RData'))
inputData <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_simData.RData'))
files <- list.files(
  file.path(FOLDER, 'Data/Outputs'), full.names = T, pattern = '2025-04-06'
)
output_file <- sprintf(
  "abc4sfa/4_error_measurements/data/output/simulations/hf/simulation_hn_misspecification_exp_%s.xlsx",
  Sys.Date()
)


scenarios <- str_extract(files, 's\\d+')

df <- NULL
# For each scenario.
for (scenario in scenarios){
  # Read results.
  outputData <- readRDS(files[grep(sprintf("%s_", scenario), files)])

  sims <- names(outputData[[scenario]])
  sigmas <- get_real_sigmas(inputData, scenario)

  est_sigmas <- NULL
  for (sim in sims){
    est_sigmas <- rbind(
      est_sigmas, 
      get_abc_sigmas(outputData[[scenario]][[sim]])
      )
  }
  df0 <- rbind(
    sigmas,
    colMeans(est_sigmas),
    colMeans(abs(est_sigmas / sigmas - 1)),
    colMeans((est_sigmas - sigmas)^2) %>% sqrt
  ) %>% as_tibble() %>% 
    rename(sigma_u = 1, sigma_eta = 2, sigma_v = 3, sigma_alpha = 4) %>% 
    mutate(
      scenario = scenario,
      type = c('Population', 'Estimate', 'Relative bias', 'RootMSE')
    ) %>% 
    relocate(scenario, type, sigma_v, sigma_alpha)
  
  df <- bind_rows(df, df0)
}

df %>% 
  fix_inverted_ids %>% 
  writexl::write_xlsx(output_file)


