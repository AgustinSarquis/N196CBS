library(ggplot2)
library(plyr)
library(dplyr)
# BUSINESS AS USUAL SCENARIO
# The simulated forest is 60 years old, so we run the simulation starting 60 years from 2025, until 2100
t=seq(1965:2100)
# We use data from Eucalyptus plantations in Hawai'i
# CAN BE IMPROVED: instead of starting with this amount of C, how about constraining the model
# to reach this amount whenever the data were taken
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
matplot(forestvalues$t,forestvalues[,2:5], type="l", lty=1,lwd=3, col=c(2:4,1), 
        ylab=" Eucalyptus Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("right",c("TAGB", "TBB","TDM","TCS"),lty=1,col=c(2:4,1), lwd=3, bty="n")

# A SECOND VERSION:
# Using repeated growth plot data from the Division of Forestry and Wildlife (DOFAW)
# For two locations near HuiMAU 
# For four Eucalyptus species:
# ES: Eucalyptus saligna
# EM: Eucalyptus microcorys
# ER: Eucalyptus robusta
# EF: Eucalyptus resinifera
# Compare data from DOFAW with measuremnts made at huiMAU on 04/08/2025
matplot(DOFAWfiltered$age,DOFAWfiltered$DBH, pch=16, col=2, 
        ylab="Eucalyptus DBH (cm)", xlab="Years")
eucs=huiMAU04082025[1:10,]
points(eucs$age, eucs$DBHcm, pch=16, col=3)
legend("bottomright",c("DOFAW", "huiMAU"),pch=16,col=2:3, cex = 0.9)








# Use allometric equations in Kaye et al. 2000 developed for Eucalyptus in Hawai'i
biomass= function(DBH) {
  wood=0.0062*DBH^3.3178
  crown=0.0082*DBH^2.2095
  root=0.029*DBH^2.432 # in Waring and Powers 2017 for Eucalyptus in Kenya
  biomass = wood + crown + root
  return(data.frame(
    wood=wood,
    crown=crown,
    root=root,
    total=biomass
  ))
}
# data=na.omit(DOFAWfiltered)
# eucbiomass=biomass(data$'DBH..cm.') # kg
# data=cbind(data, eucbiomass)
# write.csv(data, 'C:/Users/asarq/Documents/GitHub/N196CBS/Data/DOFAWfiltered.csv')
byplot = ddply(DOFAWfiltered, .(Plot, age),summarise, 
               woodarea = sum(wood)/0.04,
               crownarea = sum(crown)/0.04,
               rootarea = sum(root)/0.04, 
               totalarea=sum(total)/0.04
               ) 
ggplot(data=byplot, aes(age, totalarea/1000,color=as.factor(Plot))) +
  geom_point(size=5) +
  ylab('Total biomass (Mg/ha)') +
  theme_bw()
ggplot(data=byplot, aes(age, (woodarea+crownarea)/1000,color=as.factor(Plot))) +
  geom_point(size=5) +
  ylab('AGB (Mg/ha)')  +
  theme_bw()+ 
  geom_point(data = data.frame(age = 31, agb = 689.1), # add data point of lidar measured by HPC
                         aes(x = age, y = agb), 
                         color = "black", 
                         size = 5)
# try with Chave et al. 2014 equation instead of Kaye's
biomass2= function(DBH, Height) {
  AGB <- 0.0673 * (d*DBH^2*Height)^0.976
  BGB=0.029*DBH^2.432 # in Waring and Powers 2017 for Eucalyptus in Kenya
  total = AGB + BGB
  return(data.frame(
    AGB=AGB,
    BGB=BGB,
    total=total
  ))
}
d=mean(0.803, 0.679, 0.87) #specific gravity (g/cm3) from moisturemeters.com.au/pages/australian-species-spacific-gravity-table
eucbiomass2=biomass2(DOFAWfiltered$DBH, DOFAWfiltered$Height)
#compare predictions
plot(eucbiomass$total, eucbiomass2$total)
cor(eucbiomass$total, eucbiomass2$total) # r2: 0.968
# both methods work similarly, so I'm sticking with the Hawaiian equations in Kaye et al. 2000
# Fitting a Gompertz logistic model to total biomass data.
time=c(0, byplot$age)
biomass=c(0,byplot$totalarea/1000) # transform to Mg
gompertz <- nls(biomass ~ A * exp(-B * exp(-k * time)),
                      start = list(A = 2084, B = 3, k = 0.022), # initial values for iterative optimization
                      data = data.frame(time, biomass))
summary(gompertz)
# residual analysis
plot(time, residuals(gompertz), main="Residuals vs Time", ylab="Residuals")
hist(residuals(gompertz), main="Histogram of Residuals", xlab="Residuals")
# pseudo R2
rss <- sum(residuals(gompertz)^2)
tss <- sum((biomass - mean(biomass))^2)
pseudo_R2 <- 1 - (rss/tss)
# fitted curve vs data
plot(time, biomass, main="Gompertz Model Fit", xlab="Time", ylab="Biomass")
lines(time, predict(gompertz), col="red", lwd=2)
# final equation
forestmodel2=function(t) {
  biomass=2223 * exp(-2.946 * exp(-0.04306 * t)) # parameters from previous optimization
  return(data.frame(
    time=t,
    biomass=biomass
  ))
}
# apparently the forest was planted in 1994
t=seq(1994:2100)
forestvalues2=forestmodel2(t)
years=seq(from=1994, to=2100, by=1)
matplot(t,forestvalues2$biomass, type="l", lty=1,lwd=3, col=c(2), 
        ylab=" Eucalyptus biomass (Mg C/ha)", xlab="Years")
points(time, biomass)
# fit a model for BGB
BGB=c(0,byplot$rootarea/1000) # transform to Mg
BGBgompertz <- nls(BGB ~ A * exp(-B * exp(-k * time)),
                start = list(A = 222, B = 3, k = 0.022), # initial values for iterative optimization
                data = data.frame(time, BGB))
summary(BGBgompertz)
# residual analysis
plot(time, residuals(BGBgompertz), main="Residuals vs Time", ylab="Residuals")
hist(residuals(BGBgompertz), main="Histogram of Residuals", xlab="Residuals")
# pseudo R2
rss <- sum(residuals(BGBgompertz)^2)
tss <- sum((BGB - mean(BGB))^2)
pseudo_R2 <- 1 - (rss/tss)
# fitted curve vs data
plot(time, BGB, main="Gompertz Model Fit", xlab="Time", ylab="BGB")
lines(time, predict(BGBgompertz), col="red", lwd=2)
# final equation
BGBmodel2=function(t) {
  BGB=198.9191 * exp(-2.76846 * exp(-0.06354 * t)) # parameters from previous optimization
  return(data.frame(
    time=t,
    BGB=BGB
  ))
}
# apparently the forest was planted in 1994
t=seq(1994:2100)
BGBvalues2=BGBmodel2(t)
years=seq(from=1994, to=2100, by=1)
matplot(t,BGBvalues2$BGB, type="l", lty=1,lwd=3, col=2, ylim=c(0, max(BGB)),
        ylab=" Eucalyptus root biomass (Mg C/ha)", xlab="Years")
points(time, BGB)



##################################################################

# UTOPIAN SCENARIO

# Imagine Eucalyptus are cleared and turned into biochar (see script soilmodel.R)
# In that place, Kukui (Aleurites moluccanus) are planted as an agroforestry management 
# CAN BE IMPROVED: could try different belowground models. Also, lacking dead biomass model. Also, Kukui C content
# OJO: the formula from Chave is not in kg/ha. Maybe what Lincoln did is incorrect.
kukuimodel <- function(t, d) {
  DBH <- -9.021135 + 16.757157 * log(t) # from Lincoln 2023 (in cm)
  HEIGHT <- -3.635323 + 8.5568025 * log(t) # from Lincoln 2023 (in m)
  AGB <- 0.0673 * (d*DBH^2*HEIGHT)^0.976  # formula 4 in Chave et al. 2014 (in kg/ha) !!!!!!!!
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
# CAN BE IMPROVED:
# Kukui values are smaller than Eucalyptus, but generally are half their height, so it could be right.