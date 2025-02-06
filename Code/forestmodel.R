# BUSINESS AS USUAL SCENARIO
# The simulated forest is 60 years old, so we run the simulation starting 60 years from 2025, until 2100
t=seq(1965:2100)
# We use data from Eucalyptus plantations in Hawai'i
AGB=337.6 # aboveground biomass carbon (tons ha-1) lidar data obtained by Hayley
BGB=mean(c(22.4, 66)) # belowground biomass carbon (tons ha-1) averaged from Reeves 2012 thesis and Selments et al. 2017
DM=12 # litter carbon (tons ha-1) from Selments et al. 2017
# Forest model from Sierra et al. 2012 without the soil compartment adjusted to Hawai'i data
# CAN BE IMPROVED: NEED A EUCALYPTUS GROWTH MODEL
forestmodel <- function(t) {
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
forestvalues <- forestmodel(t)
# Plot the results
matplot(forestvalues$t,forestvalues[,2:5], type="l", lty=1,lwd=3, col=c(2:4,1),ylim=c(0,max(TCS_values2)), 
        ylab=" Eucalyptus Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("right",c("TAGB", "TBB","TDM","TCS"),lty=1,col=c(2:4,1), lwd=3, bty="n")

##################################################################

# UTOPIAN SCENARIO

# Imagine Eucalyptus are cleared and turned into biochar (see script soilmodel.R)
# In that place, Kukui (Aleurites moluccanus) are planted as an agroforestry management 
# growth model (only aboveground biomass)
# CAN BE IMPROVED: could try different belowground models. Also, lacking dead biomass model. Also, Kukui C content
# transform to Mg from kg dividing by 1000
kukuimodel <- function(t, d) {
  DBH <- -9.021135 + 16.757157 * log(t) # from Lincoln 2023 (in cm)
  HEIGHT <- -3.635323 + 8.5568025 * log(t) # from Lincoln 2023 (in m)
  AGB <- 0.0673 * (d*DBH^2*HEIGHT)^0.976  # from Chave et al. 2014 (in kg/ha)
  BGB =   0.016*DBH^2.597  # From Waring & Powers 2017 (in kg/ha) CAN BE IMPROVED
  TCS = (AGB + BGB)*0.45 # convert to carbon (standard placeholder for now)
  return(data.frame(
    t = t,
    DBH = DBH, 
    HEIGHT = HEIGHT,
    AGB = AGB/1000, # convert to Mg from kg
    BGB = BGB/1000,
    TCS = TCS/1000
  ))
}
# wood density from Krisnawati et al. 2011 
wooddensities=c(230, 310, 330, 440, 490) # (kg/m3) 
wooddensities=c(wooddensities, mean(wooddensities))/1000 # divide by 1000 because in the model it is in g/cm3

# this simulation starts in 2025, so it's only 76 years
tb=seq(2025:2100)
kukuivalues = lapply(wooddensities, function (density) kukuimodel(tb, density)) 
# Plot the results for different values of wood density (as a sensitivity analysis)
matplot(tb,kukuivalues[[6]]$AGB, type="l", lty=1,lwd=3, 
        ylab="Aboveground biomass carbon (Mg/ha)", xlab="Time (years)")
lines(tb,kukuivalues[[1]]$AGB, type="l", lty=1, col = 2, lwd=3)
lines(tb,kukuivalues[[2]]$AGB, type="l", lty=1, col = 3, lwd=3)
lines(tb,kukuivalues[[3]]$AGB, type="l", lty=1, col = 4, lwd=3)
lines(tb,kukuivalues[[4]]$AGB, type="l", lty=1, col = 5, lwd=3)
lines(tb,kukuivalues[[5]]$AGB, type="l", lty=1, col = 6, lwd=3)
legend("topleft",c("average","230 kg/m3", "310 kg/m3", "330 kg/m3", "440 kg/m3", "490 kg/m3"),lty=1,col=c(1:6), lwd=3, bty="n")
# Plot the results for average wood density
matplot(tb,kukuivalues[[6]]$TCS, type="l", lty=1,lwd=3, 
        ylab="Total Kukui carbon (Mg/ha)", xlab="Time (years)")
lines(tb,kukuivalues[[6]]$AGB*0.45, type="l", lty=1, col = 2, lwd=3)
lines(tb,kukuivalues[[6]]$BGB*0.45, type="l", lty=1, col = 3, lwd=3)    
legend("topleft",c("Total", "Aboveground", "Belowground"),lty=1,col=c(1:3), lwd=3, bty="n")

# Kukui values are smaller than Eucalyptus, but generally are half their height, so it could be right.