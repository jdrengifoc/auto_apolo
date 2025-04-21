rm(list = ls())
set.seed(010101)

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))
library(mvtnorm)
library(stringr)


# Results halfnormal-normal GTRE ---------------------------------------------
FOLDER <- 'abc4sfa/1_simulation_resources'
ABC_inputs <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_Ystats.RData'))
inputData <- readRDS(file.path(FOLDER, 'Data/Inputs/BK_simData.RData'))
files <- list.files(
  file.path(FOLDER, 'Data/Outputs'), full.names = T, pattern = '2025-04-06'
)
output_file <- sprintf(
  "abc4sfa/4_error_measurements/data/output/simulation_hn_misspecification_exp_%s.xlsx",
  Sys.Date()
)
#' It is also a common practice in ABC to perform a regression adjustment after 
#' retained draws (Beaumont et al., 2002; Leuenberger and Wegmann, 2010;
#' Sisson et al., 2018)
# Preallocate.
# Without adjustment.
ID <- str_extract(files, 's\\d+')
n_ids <- length(ID)
RelativeBias <- matrix(0, n_ids, 3)
UpperBias <- matrix(0, n_ids, 3)
PearsonCorrCoef <- matrix(0, n_ids, 3)
RelativeMSE <- matrix(0, n_ids, 3)
# With adjustment.
RelativeBiasAdj <- matrix(0, n_ids, 3)
UpperBiasAdj <- matrix(0, n_ids, 3)
PearsonCorrCoefAdj <- matrix(0, n_ids, 3)
RelativeMSEAdj <- matrix(0, n_ids, 3)

tic <- Sys.time()
df <- NULL
# For each scenario.
for (id in ID){
  # Read results.
  outputData <- readRDS(files[grep(sprintf("%s_", id), files)])

  sims <- names(outputData[[id]])
  nsim <- length(sims)
  # Fix inverted sigmas in experiments.
  sigmas <- rev(inputData[[id]]$params$sigma)

  est_sigmas <- NULL
  for (sim in sims){
    est_ABCparams <- apply(outputData[[id]][[sim]]$ABCpostChain, 2, mean)
    est_sigmas <- bind_rows(est_sigmas, est_ABCparams[2:5])
  }
  df0 <- rbind(
    sigmas,
    colMeans(est_sigmas),
    colMeans(abs(est_sigmas / sigmas - 1)),
    colMeans(est_sigmas > sigmas)
  ) %>% as_tibble() %>% 
    rename(sigma_u = 1, sigma_eta = 2, sigma_v = 3, sigma_alpha = 4) %>% 
    mutate(
      scenario = readxl::read_excel(
        file.path(FOLDER, 'Data/Inputs/Badunenko&Kumbhakar.xlsx')
        ) %>% 
        filter(ID == id) %>% pull(ID_inverted),
      type = c('Population', 'Estimate', 'Relative bias', 'Upper bias')
    ) %>% 
    relocate(scenario, type, sigma_v, sigma_alpha)
  
  df <- bind_rows(df, df0)
}

df %>% mutate(
  across(
    starts_with('sigma'), ~round(., 2L)
    )
  ) %>% 
  writexl::write_xlsx(output_file)


