# Pipeline for Running ABC on Multiple Applications
# MainIntegration_EasyABC_Pipeline.R
print(getwd())
# Set your working directory (adjust the path as needed)
FOLDER_OUTPUT <- "tristeza_empirica_aux/output"
FOLDER_INTPUT <- "tristeza_empirica_aux/input"
data_path <- file.path(FOLDER_INTPUT, "BK_empiricalData_new.RData")


# Source your helper functions.
# functions_abc.R should define get_Ti, SumStatOLS, get_SFpanelEasyABC, etc.
source(file.path(FOLDER_INTPUT, "functions_abc.R"))

app_name <- "__app_name__"
formula_list <- list(
  swissRailWays = formula(
    y ~ LNQ2 + LNQ3 + LNNET + LNSTOP + LNPK + LNPL +
      YEAR_86 + YEAR_87 + YEAR_88 + YEAR_89 + YEAR_90 +
      YEAR_91 + YEAR_92 + YEAR_93 + YEAR_94 + YEAR_95 +
      YEAR_96 + YEAR_97
  ),
  spainDairy = formula(
    y ~ X1 + X2 + X3 + X4 + X11 + X12 + X13 + X14 +
      X22 + X23 + X24 + X33 + X34 + X44 +
      YEAR_94 + YEAR_95 + YEAR_96 + YEAR_97 + YEAR_98
  ),
  usElectricity = formula(
      y ~ l_fuel + l_labor + l_k +
      I(0.5 * l_fuel^2) +
      I(l_fuel * l_labor) +
      I(l_fuel * l_k) +
      I(0.5 * l_labor^2) +
      I(l_labor * l_k) +
      I(0.5 * l_k^2) +
      year_87 + year_88 + year_89 + year_90 +
      year_91 + year_92 + year_93 + year_94 +
      year_95 + year_96
  ),
  indonesianRice = formula(
    y ~ l_seed + l_urea + l_labour + l_land + DP + DV1 + DV2 + DSS
  ),
  usBanks = formula(
    y ~ log(w1/w3)+log(w2/w3)+I((log(w1/w3))^2)+I((log(w2/w3))^2)+
    I(log(w1/w3)*log(w2/w3))+log(y1)+ log(y2)+log(y3)+log(y4)+log(y5)+
    I((log(y1))^2)+ I((log(y2))^2)+I((log(y3))^2)+I((log(y4))^2)+
    I((log(y5))^2)+I(log(y1)*log(y2))+I(log(y1)*log(y3))+I(log(y1)*log(y4))+
    I(log(y1)*log(y5))+I(log(y2)*log(y3))+I(log(y2)*log(y4))+I(log(y2)*log(y5))+
    I(log(y3)*log(y4))+I(log(y3)*log(y5))+I(log(y4)*log(y5))+I(log(y1)*log(w1/w3))+
    I(log(y1)*log(w2/w3))+I(log(y2)*log(w1/w3))+I(log(y2)*log(w2/w3))+
    I(log(y3)*log(w1/w3))+I(log(y3)*log(w2/w3))+I(log(y4)*log(w1/w3))+
    I(log(y4)*log(w2/w3))+I(log(y5)*log(w1/w3))+I(log(y5)*log(w2/w3))+trend
    )
)

prior_list <- list(
  swissRailWays = list(
    c('unif', -12, 0),
    c('unif', 0.02, 0.6),
    c('unif', 0.05, 0.45),
    c('unif', 0.005, 0.12),
    c('unif', 0.08, 0.2)
  ),
  spainDairy = list(
    c('unif', -20, -2),
    c('unif', 0.05, 0.2),
    c('unif', 0.005, 0.15),
    c('unif', 0.01, 0.25),
    c('unif', 0.001, 0.15)
  ),
  indonesianRice = list(
    c('unif', -10, -1),
    c('unif', 0.02, 0.35),
    c('unif', 0.005, 0.22),
    c('unif', 0.05, 0.48),
    c('unif', 0.1, 0.75)
  ),
  usElectricity = list(
    c('unif', -7, 0),
    c('unif', 0.04, 0.45),
    c('unif', 0.06, 0.35),
    c('unif', 0.005, 0.12),
    c('unif', 0.02, 0.36)
  ),
  usBanks = list(
    c('unif', -1.71, -0.71),
    c('unif', 0.01, 0.1),
    c('unif', 0.06, 0.26),
    c('unif', 0.01, 0.18),
    c('unif', 0.15, 0.35)
  )
)


model_type_list <- list(swissRailWays = 'cost', 
                        spainDairy = 'production',
                        usElectricity = 'production',
                        indonesianRice = 'production',
                        usBanks = 'cost')

run_ABC_pipeline(
  app_name, data_path, formula_list, prior_list, model_type_list,
  time_limit_minutes = __time_limit_minutes__, tol = 0.01, 
  chunk_size_AR = 1e5, chunk_size_MCMC = 1e4, chunk_size_SMC = 5e3
  ) 

