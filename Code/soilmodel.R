library(SoilR)

# This is a simulation of soil organic C dynamics in a Eucalyptus plantation from the model proposed in Crow et al. 2015.
# We define a time frame of 100 years
t=seq(from=1, to=100)
# We define the matrix of parameters where the diagonal parameters are the decomposition constants for each pool,
# and the rest are coefficients of transference from one pool to another.
# This model has 5 pools defined from fractionation techniques in the publication.
A=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433,
           0, 0, -0.008, 0, 0.008*0.0277,
           0, 0, 0, -0.0032, 0,
           0, 0, 0, 0, -0.001), nrow=5)
# This is the SOC stock in Eucalyptus plantations in tons ha-1 up to 1 m (from Reeves 2012 thesis)
SOC=593 
# This is the total SOC % up to 15 cm (Crow et al. 2015)
OC=sum(0.52, 3.07, 1.6, 8.36, 2.54) 
# This is the proportion of organic C per pool
pools=c(0.52/OC, 3.07/OC, 1.6/OC, 8.36/OC, 2.54/OC) 
# These are the stocks of SOC per pool
ivList=pools*SOC
# These are the C input fluxes in tons C ha-1 y-1
# Giardina et al. 2014 proposed a flux of 4.8 Mg C ha-1 y-1 up to 91.5 cm from plants to soils
# Multiplied by the proportion that enters each pool (from Crow et al. 2015) we get inputs per pool
inputFluxes=c(0.4*0.9184, 0.4*0.0045, 0.4*0.0755, 0, 0) 
# The model
model=Model(t, A, ivList, inputFluxes)
# Get stocks dynamics
Ct=cbind(rowSums(getC(model)),getC(model))
matplot(t,Ct, type="l", lty=1,lwd=3, col=1:6,ylim=c(0,max(Ct)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topright",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5"),lty=1,col=1:6, lwd=3, bty="n")

#############################################################################################
# Imagine a scenario where the Eucalyptus plantation is harvested to make biochar which is applied to the soil on site
# tree harvest rate in Mg ha-1 y-1 (Gonzalez Garcia et al. 2009)
harvestrate=11.4*549/1000
# biochar yield % from wood (Sahoo et al. 2021)
yield=mean(c(0.20, 0.065, 0.13, 0.21))
# biochar fixed carbon % (Sahoo et al. 2021)
biocharC=mean(c(0.76, 0.89, 0.58, 0.83))
# how much C in biochar can be obtained with this tree harvesting rate?
# biochar C produced from tree harvest in Mg C ha-1 y-1
biocharrate=harvestrate*yield*biocharC # this value is low compared to forest inputs (0.72 VS. 4.8)

# Now imagine all trees are harvested at once and turned into biochar
# We use Eucalyptus grandis AGB of 111.9 Mg ha-1 (Reeves 2012)
biochartotal=111.9*yield*biocharC # this is more (12.95 VS. 0.72)
# this amount of biochar is below the maximum application threshold (20 Mg ha-1, Wang et al. 2016),
# so it could be applied at once to the soil

# new soil model without inputs from the forest, and a new pool that represents biochar (decay rate from Wang et al. 2016)
# no input fluxes because the forest is gone and biochar is applied once
Ab=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433, 0,
           0, 0, -0.008, 0, 0.008*0.0277, 0,
           0, 0, 0, -0.0032, 0, 0,
           0, 0, 0, 0, -0.001, 0,
           0, 0, 0, 0, 0, -0.004), nrow=6)
ivListb=c(ivList,biochartotal)
inputFluxesb=rep(0,6) 

modelb=Model(t, Ab, ivListb, inputFluxesb)

Ctb=cbind(rowSums(getC(modelb)),getC(modelb))
matplot(t,Ctb, type="l", lty=1,lwd=3, col=1:7,ylim=c(0,max(Ctb)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topright",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5", "Biochar"),lty=1,col=1:7, lwd=3, bty="n")
