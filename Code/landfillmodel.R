################
# North 196 CBS#
################

library(SoilR)
library(FME)

# LANDFILLS

# Kona
# BaU: bovine cadavers are taken from Kulana Foods to the landfill and buried.
# Proportion of mass loss from buried bovine skeletal muscle tissue (Stokes et al. 2013)
# Multiplied by 1.8 tons of C in carcasses taken to the landfill each week (9 tons from EPA * 0.2 C content in DeBruyn et al. 2025)
df=na.omit(beefDecayStokes2013[,c(1,7)]) # beef mass loss in tons C/d
# Load function for a single pool model (from Sarquis & Sierra 2023)
onepFit=function(timeSeries, initialCarbon, inipars){
  complete=data.frame(time=timeSeries[,1],Ct=timeSeries[,2])
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
    Func=function(pars){
    mod=SoilR::OnepModel(t=tt,k=pars[1], C0=initialCarbon, In=0)
    Ct=SoilR::getC(mod)
    Rt=SoilR::getReleaseFlux(mod)
    return(data.frame(time=tt, Ct=Ct))
  }
    costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
    Fit=modFit(f=costFunc, p=inipars, method="Marq", lower= 0, upper=Inf)
  bestMod=Func(pars=Fit$par)
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
    SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=initialCarbon, In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TransitTime=1/(Fit$par)))
}
beefdecay=onepFit(df,1.8, inipars=0.5) # fit the model
# beef decomposition rate 
k=beefdecay$FMEmodel$par

# Use k rate from previous model to simulate decomposition at landfill level over 100 years 
weeks=365/7 # weeks in a year
# inputs according to Kulana Foods are around 9 Mg of carcasses per week (1.8 Mg of C)
In=1.8*weeks # inputs of carcass C Mg per year
# start in a landfill without carcasses (C0=0)
KonaBaU=OnepModel(1:100, k, 0, In)
plot(KonaBaU)
plot(getReleaseFlux(KonaBaU))
getC(KonaBaU)

# this is just to corroborate that the model is accurate
m100=OnepModel(1:100, k, In, 0)
m99=OnepModel(1:99, k, In, 0)
m98=OnepModel(1:98, k, In, 0)
m97=OnepModel(1:97, k, In, 0)
m96=OnepModel(1:96, k, In, 0)
m95=OnepModel(1:95, k, In, 0)
m94=OnepModel(1:94, k, In, 0)
m93=OnepModel(1:93, k, In, 0)
m92=OnepModel(1:92, k, In, 0)
m91=OnepModel(1:91, k, In, 0)
m90=OnepModel(1:90, k, In, 0)

C100=getC(m100)
C99=c(0, getC(m99))
C98=c(0, 0, getC(m98))
C97=c(0, 0, 0, getC(m97))
C96=c(0, 0, 0, 0, getC(m96))
C95=c(0, 0, 0, 0, 0, getC(m95))
C94=c(0, 0, 0, 0, 0, 0, getC(m94))
C93=c(0, 0, 0, 0, 0, 0, 0, getC(m93))
C92=c(0, 0, 0, 0, 0, 0, 0, 0, getC(m92))
C91=c(0, 0, 0, 0, 0, 0, 0, 0, 0, getC(m91))
C90=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, getC(m90))

all=rowSums(cbind(C100, C99, C98, C97, C96, C95, C94, C93, C92, C91, C90))
head(all)
head(getC(KonaBaU))
