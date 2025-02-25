library(SoilR)

# BUSINESS AS USUAL SCENARIO

# This is a simulation of soil organic C dynamics in a Eucalyptus plantation from the model proposed in Crow et al. 2015.
# The simulated forest is 60 years old, so we run the simulation starting 60 years from 2025, until 2100
t=seq(1965:2100)
# We define the matrix of parameters where the diagonal parameters are the decomposition constants for each pool,
# and the rest are coefficients of transference from one pool to another.
# This model has 5 pools defined from fractionation techniques in the publication.
A=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433,
           0, 0, -0.008, 0, 0.008*0.0277,
           0, 0, 0, -0.0032, 0,
           0, 0, 0, 0, -0.001), nrow=5)
# This is the SOC stock in Eucalyptus plantations in tons ha-1 up to 1 m (from Reeves 2012 thesis)
# CAN BE IMPROVED: instead of starting with this amount of SOC, how about constraining the model
# to reach this amount in 2012 (or whenever the data were taken)
SOC=593 
# This is the total SOC % up to 15 cm (Crow et al. 2015)
OC=sum(0.52, 3.07, 1.6, 8.36, 2.54) 
# This is the proportion of organic C per pool
pools=c(0.52/OC, 3.07/OC, 1.6/OC, 8.36/OC, 2.54/OC) 
# These are the stocks of SOC per pool
ivList=pools*SOC
# These are the C input fluxes in tons C ha-1 y-1
# Giardina et al. 2014 proposed a flux of 15.5 Mg C ha-1 y-1 up to 91.5 cm from plants to soils
inputs=2.7+12.4+0.4
# Multiplied by the proportion that enters each pool (from Crow et al. 2015) we get inputs per pool
inputFluxes=c(inputs*0.9184, inputs*0.0045, inputs*0.0755, 0, 0) 
# The model
soilmodel=Model(t, A, ivList, inputFluxes)
# Get stocks dynamics
Ct=cbind(rowSums(getC(soilmodel)),getC(soilmodel))
matplot(t,Ct, type="l", lty=1,lwd=3, col=1:6,ylim=c(0,max(Ct)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5"),lty=1,col=1:6, lwd=3, bty="n")

#############################################################################################

# UTOPIAN SCENARIO

# Imagine a scenario where the Eucalyptus plantation is harvested to make biochar which is applied to the soil on site
# The Eucalytpus are then replaced by Pakukui agroforestry
# this simulation starts in 2025, so it's only 76 years
tb=seq(2025:2100)
# biochar yield % from wood (Sahoo et al. 2021)
yield=mean(c(0.20, 0.065, 0.13, 0.21))
# biochar fixed carbon % (Sahoo et al. 2021)
biocharC=mean(c(0.76, 0.89, 0.58, 0.83))
# Now imagine all trees are harvested at once and turned into biochar
# We use Eucalyptus grandis aboveground carbon from the forest model at year 60
# and convert to biomass using 48.2 % Eucalyptus wood C content (Ryan et al. 2004)
AGB60=forestvalues$TAGB[60]/0.482
biochartotal=AGB60*yield*biocharC
# the maximum biochar application threshold in Wang et al. (2016) is 20 % of SOC (CAN BE IMPROVED: maybe more?)
# considering soil C at year 60 in the BaU model, that value would be:
sum(getC(soilmodel)[60,])*0.2
# we can then apply all the biochar since it's below the threshold
# new soil model with a new pool that represents a single application of biochar (decay rate from Wang et al. 2016)
# this assumes that changing the vegetation and applying biochar does not change the behavior of the soil. 
# CAN BE IMPROVED!
Ab=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433, 0,
           0, 0, -0.008, 0, 0.008*0.0277, 0,
           0, 0, 0, -0.0032, 0, 0,
           0, 0, 0, 0, -0.001, 0,
           0, 0, 0, 0, 0, -0.004), nrow=6)
# since this simulation starts 60 years after the previous one, we use those stocks at year 60 as initial conditions for this, plus biochar
ivListb=c(getC(soilmodel)[60,],biochartotal)
# FOR NOW KUKUI INPUTS ARE EQUAL TO EUCALYPTUS', BUT CAN BE IMPROVED
inputFluxesb=c(inputFluxes,0)

soilmodelb=Model(tb, Ab, ivListb, inputFluxesb)

Ctb=cbind(rowSums(getC(soilmodelb)),getC(soilmodelb))
matplot(tb,Ctb, type="l", lty=1,lwd=3, col=1:7,ylim=c(0,max(Ctb)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5", "Biochar"),lty=1,col=1:7, lwd=3, bty="n")
