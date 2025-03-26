library(dplyr)
library(coda)
library(stringr)

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
# For each empirical application.
apps <- str_split_i(basename(files), '_', -1) %>% 
  str_remove('\\..*') %>% unique
  
to_save <- list()
for (app in apps) {
  app_files <- files[str_detect(files, app)]
  to_save[[app]] <- list()
  
  # Initialize matrices.
  n <- length(app_files)
  relative_matrix <- matrix(
    0, n, n, 
    dimnames = list(
      str_split_i(basename(app_files), '_', -2), 
      str_split_i(basename(app_files), '_', -2)
      )
    )
  absolute_matrix <- relative_matrix
  
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
      
      relative_matrix[i, j] <- mean(abs(relative_dist))
      relative_matrix[j, i] <- relative_matrix[i, j]
      
      absolute_matrix[i, j] <- mean(abs(absolute_dist))
      absolute_matrix[j, i] <- absolute_matrix[i, j]
      
      
      #testthat::expect_lte(relative_matrix[i, j], 0.05)
      #testthat::expect_lte(absolute_matrix[i, j], 0.05)
    }
  }
  
  to_save[[app]][['relative_matrix']] <- relative_matrix
  to_save[[app]][['absolute_matrix']] <- absolute_matrix
}

saveRDS(to_save, output_file)
