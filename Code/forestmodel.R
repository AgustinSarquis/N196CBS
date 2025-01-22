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
# Pakukui growth model from Lincoln 2023

# time series of kukui biomass (ton/ha)
TB=c( 0.0003, 0.0100, 0.0317, 0.0613, 0.0959, 0.1334, 0.1729, 0.2135, 0.2547, 0.2962, 0.3378,
0.3793, 0.4206, 0.4616, 0.5022, 0.5425, 0.5824, 0.6218, 0.6609, 0.6994, 0.7376, 0.7752,
0.8125 ,0.8493, 0.8857, 0.9217, 0.9573 ,0.9924 ,1.0272)

# Plot the results
matplot(c(2:30), TB, type="l", lty=1,lwd=3, 
        ylab="Total biomass (ton/ha)", xlab="Time (years)")
