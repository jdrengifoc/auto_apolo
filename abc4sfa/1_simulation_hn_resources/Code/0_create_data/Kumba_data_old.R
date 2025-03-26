# Create the empirical data employed by Badunenko & Kumbhakar 2016.
install.packages("micEcon", repos="http://R-Forge.R-project.org")
install.packages('sfaR')

library(fastDummies)
library(micEcon)
data("utility")
library(sfaR, include.only = 'swissrailways')
library(sfaR)
swissrailways
data("dairyspain")

# Initialize.
empiricalData <- list()
# Swiss Railway.
attach(swissrailways)
empiricalData$swissRailWays <- list()
empiricalData$swissRailWays$y <- LNCT
X <- cbind(YEAR, LNQ2, LNQ3, LNNET, LNSTOP, LNPL, LNPK)
empiricalData$swissRailWays$X <- dummy_cols(X, select_columns = 'YEAR', remove_first_dummy = TRUE,
                remove_selected_columns = TRUE)
empiricalData$swissRailWays$loc <- data.frame(id = ID, time = YEAR)
detach(swissrailways)

# Dairy Spain.
attach(dairyspain)
empiricalData$spainDairy <- list()
empiricalData$spainDairy$y <- MILK
X <- cbind(YEAR, COWS, LAND, LABOR, FEED,
           COWS * cbind(COWS, LAND, LABOR, FEED),
           LAND * cbind(LAND, LABOR, FEED),
           LABOR * cbind(LABOR, FEED),
           FEED * FEED
)
empiricalData$spainDairy$X <- dummy_cols(X, select_columns = 'YEAR', remove_first_dummy = TRUE,
           remove_selected_columns = TRUE)
empiricalData$spainDairy$loc <- data.frame(id = FARM, time = YEAR)
detach(dairyspain)

# Utility US electricity.
attach(utility)
empiricalData$usElectricity <- list()
empiricalData$usElectricity$y <- y
X <- cbind(year, fuel, labor, k, wf, wl, wk, tc)
empiricalData$usElectricity$X <- dummy_cols(X, select_columns = 'year', remove_first_dummy = TRUE,
                                            remove_selected_columns = TRUE)
empiricalData$usElectricity$loc <- data.frame(id = firm, time = year)
detach(utility)

# Indonesian rice farms.
library(plm)

data("RiceFarms")
attach(RiceFarms)
empiricalData$indonesianRice <- list()
empiricalData$indonesianRice$y <- goutput
# FALTA WET SEASON.
empiricalData$indonesianRice$X <- cbind(seed, urea, phosphate, labour = totlabor,
                                        land = size, DP = pesticide > 0,
                                        DV1 = varieties == 'high',
                                        DV2 = varieties == 'mixed',
                                        DSS = c(1, 0))
#' file:///Users/juanrengifo101/Library/CloudStorage/OneDrive-UniversidadEAFIT/_Publicaciones/ABC4SFA/Material/Bibliografi%CC%81a/2016_Badunenko_Kumbhakar/4.a.%20Erwido1990.pdf
#' pp 195 or 221.
#' years = W75/76, D76, W76/77, D77, W82/83, D83
# wet = 1, 0, 1, 0, 1, 0
empiricalData$indonesianRice$loc <- data.frame(
  id = rep(1:length(unique(id)), each = 6), # farm
  time = rep(1:6, times = length(unique(id))) # growingseason
  )

detach(RiceFarms)

saveRDS(empiricalData, 'After/ABC4SFAv2/Data/Inputs/BK_empiricalData.RData')

empiricalData <- readRDS('After/ABC4SFAv2/Data/Inputs/BK_empiricalData.RData')



