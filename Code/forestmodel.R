# Forest model from Sierra et al. 2012 without the soil compartment

TCS_model <- function(t) {
  TAGB <- 111.51 * (1 - exp(-0.064 * t))^1.964
  TBB <- 37.665 * (1 - exp(-0.022 * t))^0.897
  TDM <- 6.615 / (1 + 2.363 * exp(-0.062 * t))
  TCS <- TAGB + TBB + TDM
  return(data.frame(
    t = t,
    TAGB = TAGB,
    TBB = TBB,
    TDM = TDM,
    TCS = TCS
  ))
}

# Create a time series (e.g., from 0 to 100 in steps of 1)
time_series <- seq(0, 135, by = 1)

# Simulate the model for the time series
TCS_values <- TCS_model(time_series)

# Plot the results
matplot(TCS_values$t,TCS_values[,2:5], type="l", lty=1,lwd=3, col=c(2:4,1),ylim=c(0,max(TCS_values)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topleft",c("TAGB", "TBB","TDM","TCS"),lty=1,col=c(2:4,1), lwd=3, bty="n")

#################################################################
# Let's try now with values for Eucalyptus plantations in Hawai'i
AGB=337.6 # use the lidar data obtained by Hayley
BGB=mean(c(22.4, 66)) # Reeves 2012 thesis and Selments et al. 2017
DM=12 # Reeves 2012 thesis

TCS_model2 <- function(t) {
  TAGB <- AGB * (1 - exp(-0.064 * t))^1.964
  TBB <- BGB * (1 - exp(-0.022 * t))^0.897
  TDM <- DM / (1 + 2.363 * exp(-0.062 * t))
  TCS <- TAGB + TBB + TDM
  return(data.frame(
    t = t,
    TAGB = TAGB,
    TBB = TBB,
    TDM = TDM,
    TCS = TCS
  ))
}

TCS_values2 <- TCS_model2(time_series)

# Plot the results
matplot(TCS_values2$t,TCS_values2[,2:5], type="l", lty=1,lwd=3, col=c(2:4,1),ylim=c(0,max(TCS_values2)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topleft",c("TAGB", "TBB","TDM","TCS"),lty=1,col=c(2:4,1), lwd=3, bty="n")

##################################################################
# Pakukui growth model 
kukui_model <- function(t, d) {
  DBH <- -9.021135 + 16.757157 * log(t) # from Lincoln 2023
  HEIGHT <- -3.635323 + 8.5568025 * log(t) # from Lincoln 2023
  AGB <- 0.0673 * (d*DBH^2*HEIGHT)^0.976  # from Chave et al. 2014
  return(data.frame(
    t = t,
    DBH = DBH, 
    HEIGHT = HEIGHT,
    AGB = AGB
  ))
}

wooddensities=c(230, 310, 330, 440, 490)/1000 # wood density from Krisnawati et al. 2011

kukui_values = lapply(wooddensities, function (density) kukui_model(time_series, density)) 

# Plot the results
matplot(time_series,kukui_values[[5]]$AGB, type="l", lty=1,lwd=3, 
        ylab="Aboveground biomass (kg/ha)", xlab="Time (years)")
lines(time_series,kukui_values[[1]]$AGB, type="l", lty=1, col = 2, lwd=3)
lines(time_series,kukui_values[[2]]$AGB, type="l", lty=1, col = 3, lwd=3)
lines(time_series,kukui_values[[3]]$AGB, type="l", lty=1, col = 4, lwd=3)
lines(time_series,kukui_values[[4]]$AGB, type="l", lty=1, col = 5, lwd=3)
legend("topleft",c("230 kg/m3", "310 kg/m3", "330 kg/m3", "440 kg/m3", "490 kg/m3"),lty=1,col=c(2:5,1), lwd=3, bty="n")

# NEED BELOWGROUND BIOMASS AND KUKUI CARBON CONTENT