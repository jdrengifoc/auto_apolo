
# Requirements ------------------------------------------------------------

SETUP_FOLDER <- 'abc4sfa/setup'
source(file.path(SETUP_FOLDER, 'requirements.R'))
source(file.path(SETUP_FOLDER, 'functions.R'))

FOLDER_OUTPUT <- "abc4sfa/4_error_measurements/data/output/simulations/exp"
dir.create(FOLDER_OUTPUT, recursive = T)


scenarios <- paste0('s', 1:16)
FOLDER_EXP_EFFICIENCIES <- '../../ABC4SFA/After/ABC4SFAv2/Data/Outputs/MH_final_este_si'
files <- list.files(FOLDER_EXP_EFFICIENCIES, full.names = T)

# Metrics efficiency simulations exp --------------------------------------

df <- NULL
for (scenario in scenarios) {
  file <- files[str_extract(files, 's\\d+') == scenario]
  mh_exp_results <- readRDS(file)[[scenario]]
  mh_exp_results$variances <- NULL
  mh_exp_results$betas <- NULL
  metrics <- get_metrics_efficiency(
    mh_exp_results, 
    name_real = 'EF_reales', 
    name_est = 'EF_abc'
    )
  
  df0 <- tibble(
    scenario = scenario, metric_name = names(metrics), metric_value = metrics
  )
  df <- bind_rows(df, df0)
}

df %>% 
  mutate(
    efficiency_type = str_split_i(metric_name, '_', -1L),
    metric_name = str_remove(metric_name, "_[^_]*$")
  ) %>% relocate(scenario, efficiency_type, metric_name, metric_value) %>% 
  fix_inverted_ids %>% 
  mutate(pp = str_extract(scenario, '\\d+') %>% as.integer) %>% 
  arrange(pp) %>% select(-pp) %>% 
  pivot_wider(names_from = 'scenario', values_from = 'metric_value') %>% 
  writexl::write_xlsx(
    sprintf(
      "%s/efficiency_metrics_abc_tun_exp_%s.xlsx",
      FOLDER_OUTPUT, Sys.Date()
    )
  )


# Count NAs ---------------------------------------------------------------

mh_exp_TE <- list()
for (scenario in scenarios) {
  file <- files[str_extract(files, 's\\d+') == scenario]
  mh_exp_results <- readRDS(file)[[scenario]]
  mh_exp_results$variances <- NULL
  mh_exp_results$betas <- NULL
  mh_exp_TE[[scenario]] <- mh_exp_results
}
count_nas(mh_exp_TE) %>% apply(2, sum) / 300
