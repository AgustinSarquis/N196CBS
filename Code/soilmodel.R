library(SoilR)

# BUSINESS AS USUAL SCENARIO

# This is a simulation of soil organic C dynamics in a Eucalyptus plantation from the model proposed in Crow et al. 2015.
# We define the matrix of parameters where the diagonal parameters are the decomposition constants for each pool,
# and the rest are coefficients of transference from one pool to another.
# This model has 5 pools defined from fractionation techniques in the publication.
A=matrix(c(-0.135, 0, 0, 0.135*0.9958, 0,
           0, -0.0033, 0, 0, 0.0033*0.8433,
           0, 0, -0.008, 0, 0.008*0.0277,
           0, 0, 0, -0.0032, 0,
           0, 0, 0, 0, -0.001), nrow=5)
# This is the SOC stock in Eucalyptus plantations in tons ha-1 up to 1 m (from Crow et al. 2016)
SOC=543.4 
# FOR LATER CONSIDERATION: this amount of SOC was measured in 2011, when the plantation was between
# 7-10 years old (say 8.5 years on average). When accounting for time in the BaU scenario, 
# deduct this amount of time. That is, if HuiMAU's plantations are 31 years old, start the BaU simulation
# at 23 years prior to account for the first 8 years of SOC accrual.
t=seq(2002:2100)
# This is the total SOC % up to 15 cm (Crow et al. 2015)
OC=sum(0.52, 3.07, 1.6, 8.36, 2.54) 
# This is the proportion of organic C per pool
pools=c(0.52/OC, 3.07/OC, 1.6/OC, 8.36/OC, 2.54/OC) 
# These are the stocks of SOC per pool
ivList=pools*SOC
# These are the C input fluxes in tons C ha-1 y-1
# Giardina et al. 2014 proposed a flux of 0.4 Mg C ha-1 y-1 up to 91.5 cm from roots to soils at steady state
# This represents 2% of total below-ground carbon flux 
# We define inputs using 2% of the BGB model
inputs=function (t) {0.02*BGBmodel2(t)[,2]}
# Multiplied by the proportion that enters each pool (from Crow et al. 2015) we get inputs per pool
inputFluxes=as.data.frame(cbind(t, 0.9184*inputs(t), 0.0045*inputs(t), 0.0755*inputs(t), 0*inputs(t), 0*inputs(t)))
# The model
soilmodel=Model(t, A, ivList, inputFluxes )
# Get stocks dynamics
Ct=cbind(rowSums(getC(soilmodel)),getC(soilmodel))
years=seq(from=2002, to=2100)
matplot(years,Ct, type="l", lty=1,lwd=3, col=1:6, 
        ylab="Carbon stocks (Mg C/ha)", xlab="Time (years)")
legend("topleft",c("Total SOC", "Pool 1","Pool 2","Pool 3", "Pool 4", "Pool 5"),lty=1,col=1:6, lwd=3, bty="n", cex = 0.7)
abline(v = 2025, col = "black", lty = 2)

inputFluxes2=c(0.9184*4, 0.0045*4, 0.0755*4, 0, 0)
soilmodel=Model(t, A, ivList, inputFluxes2)

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
