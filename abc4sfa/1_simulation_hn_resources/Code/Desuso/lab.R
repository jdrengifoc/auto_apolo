rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
# Results -----------------------------------------------------------------
permanentTE <- function(sigma_u0){
  return(2*(1 - pnorm(sigma_u0)) * exp(sigma_u0^2 / 2))
}

mse <- function(y, yhat){
  return(mean((y-yhat)^2))
}

# Read data experiment.
# wrong sigma
data <- readRDS('Data/Data_apolo/ABC_par100k.RData')
data <- data$`1.1.p`
data <- readRDS('Data/Outputs/ABC_500k_s1.RData')
data <- data$s1
data <- readRDS('Data/Outputs/ABC_batch1M_tests12.RData')
data <- readRDS('Data/Outputs/ABC_par1M_s12.RData')
data <- data$`1.2.p`
data <- readRDS('Data/Outputs/ABC_par1M_s12.RData')
data <- data$`1.2.p`

#
data <- readRDS('Data/Outputs/ABC_batch10M_fixed_s12.RData')
data <- data$`1.2.p`
data <- readRDS('Data/Outputs/ABC_batch1M_fixed_s12.RData')
data <- data$`1.2.p`
data <- readRDS('Data/Outputs/ABC_batch500k_fixed_s12.RData')
data <- data$`1.2.p`


#MSE.
# Save params and delete from list.
params <- data$params
data$params <- NULL
data$X <- NULL

(nsim <- length(data))


# MSE ---------------------------------------------------------------------

real <- c(0.3, rep(0.04, 4))

MSE_min <- rep(0, nsim)
MSE_mean <- rep(0, nsim)
for (i in 1:nsim){
  name <- names(data)[i]
  media <- data[[name]]$sumary[1:5,1]
  minimo <- data[[name]]$ABCpostChain[1,1:5]
  MSE_min[i] <- mse(media, real)
  MSE_mean[i] <- mse(minimo, real)
}
data.frame(mean = mean(MSE_mean), min =  mean(MSE_min))
# Coverage ----------------------------------------------------------------

# Parameters.
opt_fun <- function(first_n){
  real_PTE <- permanentTE(params$sigma[4])
  sigs <- c(0.1, 0.05, 0.01)
  cover_per_alpha <- rep(0, length(sigs))
  
  for (j in 1:length(sigs)) {
    coverage_permEff <- rep(0,  nsim)
    sig <- sigs[j]
    for (i in 1:nsim){
      name <- names(data)[i]
      
      PTE <- permanentTE(data[[name]]$ABCpostChain[1:first_n, 5])
      lb <- quantile(PTE, sig/2)
      ub <- quantile(PTE, 1-sig/2)
      cover <- between(real_PTE, lb, ub)
      coverage_permEff[i] <- cover
    }
    cover_per_alpha[j] <- mean(coverage_permEff)
  }
  dist <- sum((1-sigs-cover_per_alpha)^2)
  return(dist)
}

coverage <- function(first_n){
  real_PTE <- permanentTE(params$sigma[4])
  sigs <- c(0.1, 0.05, 0.01)
  cover_per_alpha <- rep(0, length(sigs))
  length_per_alpha <- rep(0, length(sigs))
  
  for (j in 1:length(sigs)) {
    coverage_permEff <- rep(0,  nsim)
    length_permEff <- rep(0,  nsim)
    sig <- sigs[j]
    for (i in 1:nsim){
      name <- names(data)[i]
      
      PTE <- permanentTE(data[[name]]$ABCpostChain[1:first_n, 5])
      lb <- quantile(PTE, sig/2)
      ub <- quantile(PTE, 1-sig/2)
      cover <- between(real_PTE, lb, ub)
      coverage_permEff[i] <- cover
      length_permEff[i] <- ub - lb
    }
    cover_per_alpha[j] <- mean(coverage_permEff)
    length_per_alpha[j] <- mean(length_permEff)
  }
  
  df <- data.frame(confidence = 1-sigs, coverage = cover_per_alpha, length = length_per_alpha)
  return(df)
}




optimize(opt_fun, c(1, 100))
# Se cogen 
coverage(5)
coverage(100)

# Estimation --------------------------------------------------------------

# Preallocation.
MEDIA <- matrix(nrow = nsim, ncol = 6)
MEDIANA <- matrix(nrow = nsim, ncol = 6)
MINIMO <- matrix(nrow = nsim, ncol = 6)
# For each simulation gather data.
for (i in 1:nsim){
  name <- names(data)[i]
  MEDIA[i,] <- data[[name]]$sumary$mean_ABC
  MEDIANA[i,] <- data[[name]]$sumary[,2]
  MINIMO[i,] <- data[[name]]$ABCpostChain[100,]
}


# Visualize estimation distribution ---------------------------------------

df <- as.data.frame(MINIMO)
names_est <- c('B0', 'v', 'v0', 'u', 'u0', 'dist')
real <- c(params$beta[1], params$sigma, 0)
names(df) <- names_est

df_long <- pivot_longer(df, cols = everything(), values_to = 'Estimate') %>% 
  group_by(name) %>%
  mutate(Mean.est = mean(Estimate)) %>% 
  mutate(real = real[which(names_est%in%name)]) 

ggplot(df_long, aes(x = Estimate, fill = name)) +
  geom_histogram(alpha=0.4) +
  geom_vline(aes(xintercept = Mean.est, col = name), size = 0.5) +
  geom_vline(aes(xintercept = real), linetype="dashed",
             col = 'black') +
  facet_wrap(. ~ name, ncol = 2, scales = 'free_x') +
  labs(title = 'Estimation for 100 simulation',
       y = 'Count', x = 'Estimate') +
  theme(
    panel.spacing.x = unit(8, 'mm'),
    panel.spacing.y = unit(5, 'mm'),
    text = element_text(family = "sans", color = "grey20"),
    #panel.background = element_blank(),
    strip.text = element_text(size=10, face = 'bold', hjust = 0.5),
    #panel.grid.major = element_line(colour="#C5C6C6", size=0.15),
    #panel.grid.minor =  element_line(colour = "#E7E7E7"),
    plot.margin = unit(c(1,1,1,1),"cm"),
    plot.caption = element_text(size = 9,colour = "grey30", hjust = 1),
    plot.subtitle = element_text(size = 12, face = "italic",colour = "grey40"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_blank(),#element_text(size=14, face="bold"),
    legend.text = element_text(size=10),
    legend.position = "none",
    axis.title = element_text(size=11, face='bold'),
    axis.text = element_text(size=10)
  )

#' Por cada simulación calcular el intervalo y veo si cubre al poblacional.
#' Luego calcular la propación de cubiertos.
#' 
#' S = 1e7, S* = 100, n = 50, t = 6, nsim = 1.
#' #' S = 1e6, S* = 100, n = 50, t = 6, nsim = 1.

#' 



n <- 300
d <- 5 
1 / (n^(-d/2) * (1/log(n)))
1e8 * (n^(-d/2) * (1/log(n)))

# Lab ---------------------------------------------------------------------

N <- c(1e3, 1e4, 1e5, 1e6, 5e6)
time <- 1:length(N)

for (j in time){
  tic <- Sys.time()
  for (i in 1:10){
    pp <- matrix(rnorm(N[j]*3), ncol = 3)
  }
  time[j] <- Sys.time() - tic
  
}
plot(N[-1], time[-1], type = 'l')


# -------------------------------------------------------------------------


# s -----------------------------------------------------------------------
#' 1. Corregir B0.
#' 2. Calcular las ineficiencias individuales
#' 3. Pensarlo para datos de panel (between).
#' 4. Buscar datos de Kumba.
#' 5. Añadir random effects (pooling).
B0 <- 0.19
(B0fix <- B0 + 0.075 * (2 / pi)^0.5)
