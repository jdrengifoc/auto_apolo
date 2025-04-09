# Pipeline for Running ABC on Multiple Applications
# MainIntegration_EasyABC_Pipeline.R
# Set your working directory (adjust the path as needed)
FOLDER <- "abc4sfa/4_enrealidad_3_gibbs_resources"
FOLDER_OUTPUT <- file.path(FOLDER, "output")
FOLDER_INPUT <- file.path(FOLDER, "input")
data_path <- file.path("abc4sfa/2_empirical_applications_resources", "BK_empiricalData.RData")

# Source your helper functions.
source(file.path(FOLDER_INPUT, "Functions2.R"))

distribution <- "hfn"
app <- "__app__"

output_file <- file.path(
  FOLDER_OUTPUT,
  sprintf("Gibbs_%s_%s_%s.RData", distribution, app, Sys.Date())
)

# run ---------------------------------------------------------------------

results <- list()
results[[app]] <- list()
data <- readRDS(data_path)[[app]]
postChain <- get_posterior_time_limit(
  data, model = "cost", burnin_rate=0.3, time_limit_seconds = __time_limit_seconds__,
  max_na_iterations = 5, fixed_beta = TRUE, thinning = 10
  )
results[[app]][['postChain']] <- postChain
# Save results.
saveRDS(results, output_file)

