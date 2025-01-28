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
SOC=593 
# This is the total SOC % up to 15 cm (Crow et al. 2015)
OC=sum(0.52, 3.07, 1.6, 8.36, 2.54) 
# This is the proportion of organic C per pool
pools=c(0.52/OC, 3.07/OC, 1.6/OC, 8.36/OC, 2.54/OC) 
# These are the stocks of SOC per pool
ivList=pools*SOC
# These are the C input fluxes in tons C ha-1 y-1
# Giardina et al. 2014 proposed a flux of 23.6 Mg C ha-1 y-1 up to 91.5 cm from plants to soils
inputs=4.4+18.8+0.4
# Multiplied by the proportion that enters each pool (from Crow et al. 2015) we get inputs per pool
inputFluxes=c(inputs*0.9184, inputs*0.0045, inputs*0.0755, 0, 0) 
# The model
model=Model(t, A, ivList, inputFluxes)
# Get stocks dynamics
Ct=cbind(rowSums(getC(model)),getC(model))
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
# We use Eucalyptus grandis AGB of 111.9 Mg ha-1 (Reeves 2012)
# THIS CAN BE FINETUNED USING THE AMOUNT OF BIOMASS AVAILABE IN 2025 FROM THE NEW MODEL!!!
biochartotal=111.9*yield*biocharC
# this amount of biochar is below the maximum application threshold (20 Mg ha-1, Wang et al. 2016),
# so it could technically be applied at once to the soil

# new soil model with a new pool that represents a single application of biochar (decay rate from Wang et al. 2016)
# this asumes that changing the vegetation and applying biochar does not change the behavior of the soil. CAN BE IMPROVED!
Ab=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433, 0,
           0, 0, -0.008, 0, 0.008*0.0277, 0,
           0, 0, 0, -0.0032, 0, 0,
           0, 0, 0, 0, -0.001, 0,
           0, 0, 0, 0, 0, -0.004), nrow=6)
# since this simulation starts 60 years after the previous one, we use those stocks at year 60 as initial conditions for this, plus biochar
ivListb=c(getC(model)[60,],biochartotal)
# FOR NOW KUKUI INPUTS ARE EQUAL TO EUCALYPTUS', BUT CAN BE IMPROVED
inputFluxesb=c(inputFluxes,0)

modelb=Model(tb, Ab, ivListb, inputFluxesb)

Ctb=cbind(rowSums(getC(modelb)),getC(modelb))
matplot(tb,Ctb, type="l", lty=1,lwd=3, col=1:7,ylim=c(0,max(Ctb)), 
        ylab="Carbon stocks (Ton C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5", "Biochar"),lty=1,col=1:7, lwd=3, bty="n")
