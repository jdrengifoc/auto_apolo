# Pipeline for Running ABC on Multiple Applications
# MainIntegration_EasyABC_Pipeline.R
# Set your working directory (adjust the path as needed)
library(dplyr)
FOLDER <- "abc4sfa/3_gibbs_resources"
FOLDER_OUTPUT <- file.path(FOLDER, "output")
FOLDER_INPUT <- file.path(FOLDER, "input")
data_path <- file.path("abc4sfa/2_empirical_applications_resources/input", "BK_empiricalData.RData")

# Source your helper functions.
source(file.path(FOLDER_INPUT, "Functions2.R"))

distribution <- "hfn"
app <- "usBanks"

output_file <- file.path(
  FOLDER_OUTPUT,
  sprintf("Gibbs_%s_%s_%s.RData", distribution, app, Sys.Date())
)

# run ---------------------------------------------------------------------

results <- list()
results[[app]] <- list()
data <- readRDS(data_path)[[app]]

fmla <- formula(
   ~ log(w1/w3)+log(w2/w3)+I((log(w1/w3))^2)+I((log(w2/w3))^2)+
    I(log(w1/w3)*log(w2/w3))+log(y1)+ log(y2)+log(y3)+log(y4)+log(y5)+
    I((log(y1))^2)+ I((log(y2))^2)+I((log(y3))^2)+I((log(y4))^2)+
    I((log(y5))^2)+I(log(y1)*log(y2))+I(log(y1)*log(y3))+I(log(y1)*log(y4))+
    I(log(y1)*log(y5))+I(log(y2)*log(y3))+I(log(y2)*log(y4))+I(log(y2)*log(y5))+
    I(log(y3)*log(y4))+I(log(y3)*log(y5))+I(log(y4)*log(y5))+I(log(y1)*log(w1/w3))+
    I(log(y1)*log(w2/w3))+I(log(y2)*log(w1/w3))+I(log(y2)*log(w2/w3))+
    I(log(y3)*log(w1/w3))+I(log(y3)*log(w2/w3))+I(log(y4)*log(w1/w3))+
    I(log(y4)*log(w2/w3))+I(log(y5)*log(w1/w3))+I(log(y5)*log(w2/w3))+trend
)
data$X <- model.matrix(fmla, data = data$X) %>% as_tibble()

postChain <- get_posterior(data, model = "cost", burnin=3.9e5, n_samples=1.3e6, thinning=20)
results[[app]][['postChain']] <- postChain
# Save results.
saveRDS(results, output_file)

