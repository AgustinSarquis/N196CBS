#############################
# SCENARIO 2: TREE CLEARING #
#############################

# load libraries
library(SoilR)
library(plyr)

# This is the simulation of 1 ha of the Eucalyptus plantation in scenario 1 after it being cleared, chipped and the chips then 
# layed on top of the soil in 2025 until the end of the century
t=seq(2025:2100)

# Chip mass

# We assume total aboveground biomass is fully converted to chips
# So we take total aboveground biomass Mg / ha at year 2025 from scenario 1
# We multiply by 0.482 to convert to C (Kaye et al. 2000)
chipC=euc_AGBvalues[32,2]*0.482 # 443.5 C Mg / ha
# In Mackensen et al. (2003) they present mass loss for 3 of the Eucalyptus species in the DOFAW data (E. grandis, microcorys and saligna)
# they show mass loss (x) over a period of 0.42 year, so in this way we compute decomposition rate (k)
k=function(x) log((100-x)/100)/0.42
massloss=c(8.5, 8.7, 25, 3.3, 18) 
ks=sapply(massloss, k)
chip_model=OnepModel(t, C0=chipC, mean(ks), In=0)
chip_Ct=getC(chip_model)
matplot(t,chip_Ct, type="l", lty=1,lwd=3, col=3, 
        ylab="Chip C stocks (Mg C/ha)", xlab="Time (years)")

# Roots and stumps

# after clearing, roots and stumps are left in place
# we simulate the decomposition of below ground biomass in 2025 from the previous scenario
rootC=euc_BGBvalues[32,2]*0.482 # 66.7 C Mg / ha
# we use the decomposition rate for broad leaf trees of 0.45 y-1 in Silver & Miya 2001
kr=0.45
root_model=OnepModel(t, C0=rootC, k=kr, In=0)
root_Ct=getC(root_model)
matplot(t,root_Ct, type="l", lty=1,lwd=3, col=5, 
        ylab="Root C stocks (Mg C/ha)", xlab="Time (years)")

# ALTERNATIVELY, we can consider roots decompose and a portion is stabilized in SOM
# we can fit a two-pool series model with one pool being roots that have a transfer coefficient of 0.02
# this is an approximation from Giardina et al. (2014), a placeholder for now

# Soil and roots

ksoil=mean(ks)
a21=0.02
# This is the SOC stock in 2025
SOC25=Ct[32]
# Start simulations
soil_model2=TwopSeriesModel(t, ks=c(kr, ksoil), In=0, C0=c(rootC,SOC25), a21=a21)
Ct2=getC(soil_model2)
matplot(t,Ct2, type="l", lty=1,lwd=3, col=c(2,3), 
        ylab="Soil C stocks (Mg C/ha)", xlab="Time (years)")

# finally add emissions from the WARM model