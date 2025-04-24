# Pipeline for Running ABC on Multiple Applications
# MainIntegration_EasyABC_Pipeline.R
# Set your working directory (adjust the path as needed)
FOLDER <- "abc4sfa/3_gibbs_resources"
FOLDER_OUTPUT <- file.path(FOLDER, "output")
FOLDER_INPUT <- file.path(FOLDER, "input")
data_path <- file.path(FOLDER_INPUT, "BK_empiricalData.RData")

source(file.path(FOLDER_INPUT, "Functions2.R"))

scenario <- "__scenario__"
distribution <- "hn"

output_file <- file.path(
  FOLDER_OUTPUT,
  sprintf("Gibbs_experiment_%s_%s_%s.RData", distribution, scenario, Sys.Date())
)

# run ---------------------------------------------------------------------

FOLDER_fixed <- "abc4sfa/1_simulation_resources"
# Read data.
ABC_sample <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_simData.RData'))
ABC_inputs <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_Ystats.RData'))[[scenario]][['params']]
# if (distribution == 'hn') {
#   ABC_sample <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_simData.RData'))
#   ABC_inputs <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_Ystats.RData'))[[scenario]][['params']]
# } else {
#   ABC_sample <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_exp_simData.RData'))
#   ABC_inputs <- readRDS(file.path(FOLDER_fixed, 'Data/Inputs/BK_exp_Ystats.RData'))[[scenario]][['params']]
# }

# Process data for `get_posterior()`
data_base <- list(
  X = ABC_sample$X,
  loc = tibble(
    id = rep(1:ABC_inputs$n, each = ABC_inputs$t),
    time = rep(1:ABC_inputs$t, times = ABC_inputs$n)
    ))

# Initialize.
results <- list()
sims <- names(ABC_sample[[scenario]])
sims <- sims[grepl('sim', sims)]

for (sim in sims) {
  data <- c(data_base, list(y = ABC_sample[[scenario]][[sim]]$y))
  postChain <- get_posterior(
    data, model = "production", n_samples = 25e3, burnin = 10e3, thinning = 10
    )
  results[[sim]][['postChain']] <- postChain
}

# Save results.
saveRDS(results, output_file)