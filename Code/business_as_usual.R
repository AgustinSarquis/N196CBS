#################################
# SCENARIO 1: BUSINESS AS USUAL #
#################################

# load libraries
library(SoilR)
library(plyr)

# This is the simulation of 1 ha of an unmanaged Eucalyptus plantation in the Hamakua Coast, Hawai'i, planted in 1994. 

# Tree biomass

# Using repeated growth plot data from the Division of Forestry and Wildlife (DOFAW)
# Use aboveground biomass allometric equations in Kaye et al. 2000 developed for Eucalyptus in Hawai'i and
# root allometric equation in Kuyah et al. 2013 for Eucalyptus in Kenya
euc_allometry= function(DBH) {
  wood=0.0062*DBH^3.3178
  crown=0.0082*DBH^2.2095
  root=0.029*DBH^2.432 
  biomass = wood + crown + root
  return(data.frame(
    wood=wood,
    crown=crown,
    root=root,
    total=biomass
  ))
}

# Calculate tree biomass (kg) for DOFAW samples and store in data frame
# data=na.omit(DOFAWfiltered)
# eucbiomass=euc_allometry(data$'DBH..cm.') # kg
# data=cbind(data, eucbiomass)
# write.csv(data, 'C:/Users/asarq/Documents/GitHub/N196CBS/Data/DOFAWfiltered.csv')

# With results, get total biomass per area by plot (0.04 hectares)
byplot = ddply(DOFAWfiltered, .(Plot, age),summarise, 
               woodarea = sum(wood)/0.04,
               crownarea = sum(crown)/0.04,
               rootarea = sum(root)/0.04, 
               totalarea=sum(total)/0.04
) 

# Fitting a Gompertz logistic growth model to total biomass data.
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

# final Ecualyptus growth equation (Mg/ha)
euc_model=function(t) {
  biomass=2223 * exp(-2.946 * exp(-0.04306 * t)) # parameters from previous optimization
  return(data.frame(
    time=t,
    biomass=biomass
  ))
}

# Start simulation from 1994 until 2100
t=seq(1994:2100)
euc_values=euc_model(t) # in Mg/ha, multiply by 0.482 to convert to C (Kaye et al. 2000)
years=seq(from=1994, to=2100, by=1)
matplot(t,euc_values$biomass, type="l", lty=1,lwd=3, col=c(2), 
        ylab=" Eucalyptus biomass (Mg/ha)", xlab="Years")
points(time, biomass) # add DOFAW plot points for reference

# Now fit a similar model only for aboveground biomass
AGB=c(0,(byplot$crownarea+byplot$woodarea)/1000) # transform to Mg
AGBgompertz <- nls(AGB ~ A * exp(-B * exp(-k * time)),
                   start = list(A = 1862, B = 3, k = 0.022), # initial values for iterative optimization
                   data = data.frame(time, AGB))
summary(AGBgompertz)
# residual analysis
plot(time, residuals(AGBgompertz), main="Residuals vs Time", ylab="Residuals")
hist(residuals(AGBgompertz), main="Histogram of Residuals", xlab="Residuals")
# pseudo R2
rss <- sum(residuals(AGBgompertz)^2)
tss <- sum((AGB - mean(AGB))^2)
pseudo_R2 <- 1 - (rss/tss)
# final equation
euc_AGBmodel=function(t) {
  AGB=2043 * exp(-3.01 * exp(-0.0415 * t)) # parameters from previous optimization
  return(data.frame(
    time=t,
    AGB=AGB
  ))
}
# Start simulation
euc_AGBvalues=euc_AGBmodel(t) # in Mg/ha, multiply by 0.482 to convert to C (Kaye et al. 2000)
matplot(t,euc_AGBvalues$AGB, type="l", lty=1,lwd=3, col=2, ylim=c(0, max(AGB)),
        ylab=" Eucalyptus aboveground biomass (Mg/ha)", xlab="Years")
points(time, AGB)

# Now fit a similar model only for roots to compute C inputs in the soil model
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
# final equation
euc_BGBmodel=function(t) {
  BGB=198.9191 * exp(-2.76846 * exp(-0.06354 * t)) # parameters from previous optimization
  return(data.frame(
    time=t,
    BGB=BGB
  ))
}
# Start simulation
euc_BGBvalues=euc_BGBmodel(t) # in Mg/ha, multiply by 0.482 to convert to C (Kaye et al. 2000)
matplot(t,euc_BGBvalues$BGB, type="l", lty=1,lwd=3, col=2, ylim=c(0, max(BGB)),
        ylab=" Eucalyptus root biomass (Mg/ha)", xlab="Years")
points(time, BGB)

# Soil 

# the next model was fit with SIDb data for non native forests in Andisols in the Crow2019a entry.
# sites 46, 47, 48, 70, 71 and 72
# parameters for a one pool model (in Azizi-Rad et al. this is the best model)
ks=c(0.00138, 0.00169, 0.00142, 0.00214, 0.00327, 0.0026023)
# inputs from roots (C Mg / ha)
# Giardina et al. 2014 proposed a flux of 0.4 Mg C ha-1 y-1 up to 91.5 cm from roots to soils at steady state
# This represents 2% of total below-ground carbon flux 
# We define inputs using 2% of the BGB model
inputs=function (t) {0.02*euc_BGBmodel(t)[,2]}
In=as.data.frame(cbind(t, inputs(t)))
# This is the SOC stock in Eucalyptus plantations in Mg ha-1 up to 15 cm (from CIG project on site)
SOC=72.38
# Start simulations
soilmodel=OnepModel(t, mean(ks), In, C0=SOC)
Ct=getC(soilmodel)
matplot(t,Ct, type="l", lty=1,lwd=3, col=4, 
        ylab="Soil C stocks (Mg C/ha)", xlab="Time (years)")

# finally add emissions from the WARM model
