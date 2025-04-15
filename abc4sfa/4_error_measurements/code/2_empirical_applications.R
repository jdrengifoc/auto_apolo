library(dplyr)
library(coda)
library(stringr)


theme4policies <- function(base_size = 16, base_family = "serif") {
  theme(
    axis.line = element_line(colour = 'black'),
    axis.ticks = element_line(colour = 'black'),
    axis.ticks.length = unit(0.5, 'lines'),
    
    text = element_text(
      family = base_family, face = "plain", colour = "black", size = base_size,
      hjust = 0.5, vjust = 1, angle = 0, lineheight = 1, margin = margin(),
      debug = FALSE),
    axis.title = element_text(colour = 'black', face = 'plain'),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(colour = 'black'),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
    
    plot.background = element_rect(fill = 'white'),
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(linetype = 'dotted', linewidth = 0.3,
                                    color = 'grey60'),
    panel.grid.minor.x = element_line(linetype = 'dotted', linewidth = 0.3, color =
                                        'grey60'),
    panel.border = element_blank(),
    plot.margin = margin(t = 1, b = 1, r = 2, l = 1, unit = 'lines'),
    
    legend.background = element_rect(linetype = 'solid', linewidth = 0.4,
                                     colour = NA, fill = NA),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = 0.5,
    legend.key = element_rect(fill = NA, colour = NA)
  )
}


setwd('auto_apolo/')
folder <- "abc4sfa/2_empirical_applications_resources/output"
files <- list.files(folder, "RData", full.names = T)
output_file <- sprintf(
  "abc4sfa/3_error_measurements/data/output/2_empirical_applications_%s.rds",
  Sys.Date()
)

testthat::test_that('Check missings', {
  for (file in files) {
    results <- load(file) %>% get
    n_nas <- sum(apply(results, 2, function(x) sum(is.na(x))))
    testthat::expect_equal(n_nas, 0)
  }
})



# Compare the results of the empirical applications for each pair of methods ----
apps <- str_split_i(basename(files), '_', -1) %>% 
  str_remove('\\..*') %>% unique
  
to_save <- list()
for (app in apps) {
  app_files <- files[str_detect(files, app)]
  to_save[[app]] <- list()
  
  # For each pair of files.
  for (file1 in app_files) {
    for (file2 in app_files) {
      i <- which(app_files == file1)
      j <- which(app_files == file2)
      
      if (i >= j) next
      method1 <- str_split_i(basename(file1), '_', -2)
      method2 <- str_split_i(basename(file2), '_', -2)
      sprintf("%s: %s Vs %s", app, method1, method2) %>% print
      summary1 <- (load(file1) %>% get %>% mcmc %>% summary)$statistics
      summary2 <- (load(file2) %>% get %>% mcmc %>% summary)$statistics
      
      absolute_dist <- summary1 - summary2
      relative_dist <- (summary1 / summary2 - 1)
      to_save[[app]][[sprintf("%s-%s", method1, method2)]] <- list(
        relative_dist = relative_dist,
        absolute_dist = absolute_dist
      )
    }
  }
}

saveRDS(to_save, output_file)



# Get data to plot  -------------------------------------------------------


library(ggplot2)
library(tidyr)

sig_level <- 0.05
df <- NULL
df_ci <- NULL 

apps <- str_split_i(basename(files), '_', 3) %>% 
  str_remove('\\..*') %>% unique

to_save <- list()
for (app in apps) {
  app_files <- files[str_detect(files, app)]
  to_save[[app]] <- list()
  
  # For each pair of files.
  for (file in app_files) {
      method <- str_split_i(basename(file), '_', 2)
      running_time <- str_split_i(basename(file), '_', 4) %>% 
        str_remove('.*=') %>% as.integer()
      chunk_size <- str_split_i(basename(file), '_', 6) %>% 
        str_remove('.*=') %>% as.integer()
      saving_date <- basename(file) %>% str_extract('\\d{4}-\\d{2}-\\d{2}')
      
      chains <- load(file) %>% get 
      point_estimate <- (chains %>% mcmc() %>% summary)$statistic[, 'Mean']
      cis <- chains %>% mcmc() %>% HPDinterval(prob = 1 - sig_level) %>% 
        as_tibble() %>% 
        mutate(
          point = point_estimate,
          parameter = c('beta0', 'sigma_v', 'sigma_alpha', 'sigma_eta', 'sigma_u')
          )
      
  
      cis <- cis %>% mutate(app = app, method = method) %>% 
        relocate(app, method, parameter)
      
      chains <- chains %>% as_tibble()
      names(chains) <- c('beta0', 'sigma_v', 'sigma_alpha', 
                         'sigma_eta', 'sigma_u')
      chains <- chains %>% 
        mutate(
          date = saving_date, app = app, 
          method = method, running_time = running_time, chunk_size = chunk_size,
          lambda = sigma_u / sigma_v,
          lambda_delta = sigma_eta / sigma_alpha,
          Lambda = sigma_eta / sigma_u
          ) %>% relocate(date, app, method, running_time, chunk_size)
      df_ci <- bind_rows(df_ci, cis)
      df <- bind_rows(df, chains)
  }
}

df_priors <- list(
  swissRailWays = list(
    beta0       = c(-12, 0),
    sigma_v     = c(0.02, 0.6),
    sigma_alpha = c(0.05, 0.45),
    sigma_eta   = c(0.005, 0.12),
    sigma_u     = c(0.08, 0.2)
  ),
  spainDairy = list(
    beta0       = c(-20, -2),
    sigma_v     = c(0.05, 0.2),
    sigma_alpha = c(0.005, 0.15),
    sigma_eta   = c(0.01, 0.25),
    sigma_u     = c(0.001, 0.15)
  ),
  indonesianRice = list(
    beta0       = c(-10, -1),
    sigma_v     = c(0.02, 0.35),
    sigma_alpha = c(0.005, 0.22),
    sigma_eta   = c(0.05, 0.48),
    sigma_u     = c(0.1, 0.75)
  ),
  usElectricity = list(
    beta0       = c(-7, 0),
    sigma_v     = c(0.04, 0.45),
    sigma_alpha = c(0.06, 0.35),
    sigma_eta   = c(0.005, 0.12),
    sigma_u     = c(0.02, 0.36)
  ),
  usBanks = list(
    beta0       = c(-1.71, -0.71),
    sigma_v     = c(0.01, 0.1),
    sigma_alpha = c(0.06, 0.26),
    sigma_eta   = c(0.01, 0.18),
    sigma_u     = c(0.15, 0.35)
  )
) %>% purrr::map_df(
  ~ tibble::enframe(.x, name = "parameter", value = "interval"), 
  .id = "app"
  ) %>%
  unnest_wider(interval, names_sep = "_") %>%
  rename(lower = interval_1, upper = interval_2)


# Plot posteriors density per app-method ----------------------------------

library(lubridate)
library(scales)
library(latex2exp)
library(gridExtra)

parameter_labels <- c(
  "sigma_v"      = TeX("$\\sigma_v$"),
  "sigma_alpha"  = TeX("$\\sigma_{\\alpha}$"),
  "sigma_eta"    = TeX("$\\sigma_{\\eta}$"),
  "sigma_u"      = TeX("$\\sigma_u$"),
  "lambda"       = TeX("$\\lambda$"),
  "lambda_delta" = TeX("$\\lambda_{\\delta}$"),
  "Lambda"       = TeX("$\\Lambda$")
)

df_plot_base <- df %>% 
  mutate(date = ymd(date)) %>% 
  filter(between(date, ymd('2025-03-20'), ymd('2025-03-31'))) %>% 
  select(-beta0) %>% 
  # distinct(app, date, method, chunk_size, running_time)
  pivot_longer(
    cols = matches('beta|sigma|ambda'), 
    names_to = 'parameter', values_to = 'value'
    ) %>% 
  mutate(
    # parameter = paste(app, parameter, sep = '_'),
    parameter = factor(
      parameter, 
      levels = names(parameter_labels),
      labels = parameter_labels
      ),
    app = case_when(
        app == "indonesianRice" ~ "IndonesianRice",
        app == "spainDairy" ~ "SpainDairy",
        app == "swissRailWays" ~ "SwissRailways",
        app == "usBanks" ~ "USBanks",
        app == "usElectricity" ~ "USElectricity"
    )
  )


# Sigmas ------------------------------------------------------------------

df_plot <- df_plot_base %>% filter(str_detect(parameter, 'sigma'))
fig_ir <- df_plot %>% 
  filter(app == 'IndonesianRice') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = "Indonesian Rice")

fig_sp <- df_plot %>% 
  filter(app == 'SpainDairy') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = "Spain Dairy")

fig_sr <- df_plot %>% 
  filter(app == 'SwissRailways') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = 'Density', title = 'Swiss Railways')

fig_ub <- df_plot %>% 
  filter(app == 'USBanks') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = 'US Banks')

fig_ue <- df_plot %>% 
  filter(app == 'USElectricity') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.title = element_blank(),
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = "Parameter Value", y = "", title = 'US Electricity')



ggsave(
  filename = 'abc4sfa/4_error_measurements/data/output/empirical_applications_grid_sigmas.png',
  plot = grid.arrange(fig_ir, fig_sp, fig_sr, fig_ub, fig_ue, nrow = 5),
  width = 18, height = 10, dpi = 300
  )


# Lambdas -----------------------------------------------------------------

df_plot <- df_plot_base %>% filter(!str_detect(parameter, 'sigma'))

fig_ir <- df_plot %>% 
  filter(app == 'IndonesianRice') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = "Indonesian Rice")

fig_sp <- df_plot %>% 
  filter(app == 'SpainDairy') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = "Spain Dairy")

fig_sr <- df_plot %>% 
  filter(app == 'SwissRailways') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = 'Density', title = 'Swiss Railways')

fig_ub <- df_plot %>% 
  filter(app == 'USBanks') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = NULL, y = "", title = 'US Banks')

fig_ue <- df_plot %>% 
  filter(app == 'USElectricity') %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.title = element_blank(),
    legend.position = 'none',
    strip.text = element_text(size = 10)
  ) +
  labs(x = "Parameter Value", y = "", title = 'US Electricity')



ggsave(
  filename = 'abc4sfa/4_error_measurements/data/output/empirical_applications_grid_lambdas.png',
  plot = grid.arrange(fig_ir, fig_sp, fig_sr, fig_ub, fig_ue, nrow = 5),
  width = 18, height = 10, dpi = 300
)


# Just SMC ----------------------------------------------------------------

df_plot %>% 
  ggplot(aes(x = value, color = method)) +
  facet_wrap(
    ~ app + parameter, 
    scales = "free", ncol = 7, 
    labeller = label_parsed
  ) + 
  geom_density(linewidth = 0.35) +
  scale_color_manual(
    values = c(
      "AR" = "#91B654", 
      "MCMC" = "#9B006E", 
      "SMC" = "#00929B"
    )
  ) +
  theme_minimal() + 
  theme(
    legend.title = element_blank(),
    legend.position = 'top',
    strip.text = element_text(size = 10)
  ) + 
  labs(
    y = 'Density', x = 'Parameter Value'
  )

# Plot credible intervals -------------------------------------------------

df_ci %>% 
  bind_rows(df_priors %>% mutate(method = 'Prior')) %>% 
  mutate(
    method = factor(
      method, levels = c("Prior", "AR", "MCMC", "SMC")
      )
    ) %>% 
  ggplot(aes(y = method, xmin = lower, xmax = upper, x = point, color = method)) +
  geom_errorbarh(height = 0.3) +  # Horizontal bars for confidence intervals
  geom_point(size = 2) +  # Point estimates
  facet_wrap(~app + parameter, scales = "free_x") +  # Separate plots for each parameter and app
  scale_color_manual(values = c("Prior" = "black", "AR" = '#9459C9', "MCMC" = '#2BA2CF', "SMC" = "#469500")) +  # Colores manuales
  theme_minimal() +
  theme_minimal() +
  theme(legend.position = "top") +  # Remove legend if method is already on the y-axis
  labs(x = "Estimate", y = "Method")


