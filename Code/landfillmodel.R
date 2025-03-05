################
# North 196 CBS#
################

# LANDFILLS

# Kona
# BaU: bovine cadavers are taken from Kulana Foods to the landfill and buried.
# Respiration from buried bovine skeletal muscle tissue (Stokes et al. 2013)
library(SoilR)
library(FME)
df=na.omit(beefDecayStokes2013[,c(1,4)]) # beef mass loss in C %/d
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
beefdecay=onepFit(df,100, inipars=2.5) # fit the model
# beef decomposition rate (C g/d)
k=beefdecay$FMEmodel$par

# Use k rate from previous model to simulate decomposition at landfill level over 100 years (36500 days)
# inputs according to Kulana Foods are around 9 Mg of carcasses per week
# converted to daily rates (1.28), to C (0.2) and to g (1e6)
In=1.28*0.2*1e6
# start in a landfill without carcasses (C0=0)
KonaBaU=OnepModel(1:36500, k, 0, In)
plot(KonaBaU)
plot(getReleaseFlux(KonaBaU))
# another way
u=onePoolFlux(df, 100, 0.1)

df=na.omit(beefDecayStokes2013[,c(1,5)]) # beef mass loss in C %/d
# another way
lin = lm(df$Mt ~ df$days)
summary(lin)
coef(lin)[2]
#another way
df=na.omit(beefDecayStokes2013[,c(1,2)]) # beef mass loss in C %/d

chu=onePoolFlux(df, 80000, 0.1)
getReleaseFlux(chu$SoilRmodel)
