#####################
# BUSINESS AS USUAL #
#####################

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
euc_values=euc_model(t)
years=seq(from=1994, to=2100, by=1)
matplot(t,euc_values$biomass, type="l", lty=1,lwd=3, col=c(2), 
        ylab=" Eucalyptus biomass (Mg C/ha)", xlab="Years")
points(time, biomass) # add DOFAW plot points for reference

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
euc_BGBvalues=euc_BGBmodel(t)
matplot(t,euc_BGBvalues$BGB, type="l", lty=1,lwd=3, col=2, ylim=c(0, max(BGB)),
        ylab=" Eucalyptus root biomass (Mg C/ha)", xlab="Years")
points(time, BGB)

# Soil 

# Using the model proposed in Crow et al. 2015 for a native forest near the area of our study site
# We define the matrix of parameters where the diagonal parameters are the decomposition constants for each pool,
# and the rest are coefficients of transference from one pool to another.
# This model has 5 pools defined from fractionation at a depth of 0-15 cm.
A=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433,
           0, 0, -0.008, 0, 0.008*0.0277,
           0, 0, 0, -0.0032, 0,
           0, 0, 0, 0, -0.001), nrow=5)
# This is the SOC stock in Eucalyptus plantations in Mg ha-1 up to 15 cm (from CIG project on site)
SOC=72.38
# This is the total SOC % in each fraction (Crow et al. 2015)
OC=sum(0.52, 3.07, 1.6, 8.36, 2.54) 
# This is the proportion of organic C per pool
pools=c(0.52/OC, 3.07/OC, 1.6/OC, 8.36/OC, 2.54/OC) 
# These are the stocks of SOC per pool
ivList=pools*SOC
# These are the C input fluxes in Mg C ha-1 y-1
# Giardina et al. 2014 proposed a flux of 0.4 Mg C ha-1 y-1 up to 91.5 cm from roots to soils at steady state
# This represents 2% of total below-ground carbon flux 
# We define inputs using 2% of the euc_BGBmodel
inputs=function (t) {0.02*euc_BGBmodel(t)[,2]}
# Multiplied by the proportion that enters each pool (from Crow et al. 2015) we get inputs per pool
inputFluxes=as.data.frame(cbind(t, 0.9184*inputs(t), 0.0045*inputs(t), 0.0755*inputs(t), 0*inputs(t), 0*inputs(t)))

# Start simulation
soilmodel=Model(t, A, ivList, inputFluxes )
# Get stocks dynamics
Ct=cbind(rowSums(getC(soilmodel)),getC(soilmodel))
matplot(years,Ct, type="l", lty=1,lwd=3, col=1:6, 
        ylab="Carbon stocks (Mg C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5"),lty=1,col=1:6, lwd=3, bty="n", cex = 0.7)
abline(v = 2025, col = "black", lty = 2)

# alternative. the previous model overestimates final C
# the next model was fit with SIDb data for non native forests in Andisols in the Crow2019a entry.
# sites 46, 47, 48, 70, 71 and 72
# parameters for a two pool series model
k1=c(0.03907, 0.02842, 0.02411, 0.01062, 0.00331, 0.0068647)
k2=c(0.00102, 0.00125, 0.00108, 0.00162, 0.000965, 0.0020263)
gamma=c(0.141, 0.093, 0.058, 0.117, 0.987, 0.1807374)
a21=c(0.688883772, 0.424837969, 0.205219313, 0.170861778, 0.0000944,0.1129021)
# inputs from roots
In=as.data.frame(cbind(t, inputs(t), 0*inputs(t)))
# Start simulation
soilmodel46=TwopSeriesModel(t, ks=c(k1[1], k2[1]), a21[1], In, C0=c(SOC*gamma[1], SOC*(1-gamma[1])))
soilmodel47=TwopSeriesModel(t, ks=c(k1[2], k2[2]), a21[2], In, C0=c(SOC*gamma[2], SOC*(1-gamma[2])))
soilmodel48=TwopSeriesModel(t, ks=c(k1[3], k2[3]), a21[3], In, C0=c(SOC*gamma[3], SOC*(1-gamma[3])))
soilmodel70=TwopSeriesModel(t, ks=c(k1[4], k2[4]), a21[4], In, C0=c(SOC*gamma[4], SOC*(1-gamma[4])))
soilmodel71=TwopSeriesModel(t, ks=c(k1[5], k2[5]), a21[5], In, C0=c(SOC*gamma[5], SOC*(1-gamma[5])))
soilmodel72=TwopSeriesModel(t, ks=c(k1[6], k2[6]), a21[6], In, C0=c(SOC*gamma[6], SOC*(1-gamma[6])))
# Get stocks dynamics
Ct=cbind(rowSums(getC(soilmodel71)),getC(soilmodel71))
matplot(t,Ct, type="l", lty=1,lwd=3, col=1:6, 
        ylab="Carbon stocks (Mg C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2"),lty=1,col=1:3, lwd=3, bty="n", cex = 0.7)
