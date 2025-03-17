# Create the empirical data employed by Badunenko & Kumbhakar 2016.
install.packages("micEcon", repos="http://R-Forge.R-project.org")
install.packages('sfaR')

library(fastDummies)
library(micEcon)

library(sfaR, include.only = 'swissrailways')
library(sfaR)

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
empiricalData$spainDairy$y <- YIT
X <- cbind(YEAR, X1, X2, X3, X4, X11, X22, X33, X44,
           X12, X13, X14, X23, X24, X34)



empiricalData$spainDairy$X <- dummy_cols(X, select_columns = 'YEAR', remove_first_dummy = TRUE,
           remove_selected_columns = TRUE)
empiricalData$spainDairy$loc <- data.frame(id = FARM, time = YEAR)
detach(dairyspain)

# Utility US electricity.
data("utility")
attach(utility)
# NO CUADRAN LAS UNIDADES PARA NADA, NI SÉ QUE VARIABLES SON.
data.frame(Mean = apply(utility, 2, mean), Sd = apply(utility, 2, sd))
empiricalData$usElectricity <- list()
empiricalData$usElectricity$y <- log(y)
#X <- cbind(year, fuel, labor, k, wf, wl, wk, tc)
X <- cbind(year, l_fuel = log(fuel), l_labor = log(labor), l_k = log(k),
           l_fuel2 = I(0.5*log(fuel)^2), l_labor2 = I(0.5*log(labor)^2),
           l_k2 = I(0.5*log(k)^2),
           l_fuel_labor = I(log(fuel)*log(labor)), l_fuel_kI = (log(fuel)*log(k)),
           l_labor_k = I(log(labor)*log(k)))
empiricalData$usElectricity$X <- dummy_cols(X, select_columns = 'year', remove_first_dummy = TRUE,
                                            remove_selected_columns = TRUE)
empiricalData$usElectricity$loc <- data.frame(id = firm, time = year)
detach(utility)

# Indonesian rice farms.
library(plm)

data("RiceFarms")
attach(RiceFarms)
empiricalData$indonesianRice <- list()
empiricalData$indonesianRice$y <- log(goutput)
#' FALTA WET SEASON.
#' No se pone fosfato.
# empiricalData$indonesianRice$X <- cbind(seed, urea, phosphate, labour = totlabor,
#                                         land = size, DP = pesticide > 0,
#                                         DV1 = varieties == 'high',
#                                         DV2 = varieties == 'mixed',
#                                         DSS = c(1, 0))
empiricalData$indonesianRice$X <- data.frame(l_seed = log(seed), l_urea = log(urea),
                                        l_labour = log(totlabor),
                                        l_land = log(size), DP = pesticide > 0,
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

saveRDS(empiricalData, 'After/ABC4SFAv2/Data/Inputs/BK_empiricalData_new.RData')

empiricalData <- readRDS('After/ABC4SFAv2/Data/Inputs/BK_empiricalData_new.RData')

for (name in names(empiricalData)) {
  print(names(empiricalData[[name]]$X))
}
