# Code for parameter estimation of fractionation data

#First load these packages into your R session. If not installed yet, simply type install.packages("name"), with the name of the package inside quotations
library(SoilR)
library(FME)
library(forecast)

#Load the dataset. Change address to a specific directory in your computer
and=read.csv("~/Dropbox/SoilR_Crow/AND.csv")

#Load the model. Change address to a specific directory in your computer
source("~/Dropbox/SoilR_Crow/CrowModel.R")

#Create objects with the year of sampling and sequence of years for simulation. 
yrsample=2010
years=seq(-20000,2010)

#Now the C stock for the different land uses
Cnat=33.443 # C stock native forest (Mg C ha-1)
Ceuc=120.467 # C stock Eucalyptus plantation (Mg C ha-1)
Cpas=86.963 # C stock pasture (Mg C ha-1)

#C input numbers
TBCF=12.4 #Total Belowground C Flux (Mg C ha-1 yr-1) from Fig.4 Nat Clim Ch paper.
LF=2.7 #Litter Fall C (Mg C ha-1 yr-1) from Fig.4 Nat Clim Ch paper.
RI=0.4 #Root inputs to mineral soil (Mg C ha-1 yr-1) from Fig.4 Nat Clim Ch paper.

################################################################################
# Preparation of bomb-curve dataset

#Bind pre- and postbomb curves
calcurves=bind.C14curves(prebomb=IntCal13,postbomb=Hua2013$NHZone2,time.scale="AD")
bombcurve=calcurves[calcurves[,1]>=-10000,1:2]

#Forecast atmospheric Delta14C values up to 2014
yrs=seq(1966,2009.5,by=1/4) # A series of years by quarters
nz2=spline(Hua2013$NHZone2[,c(1,4)],xout=yrs) #Spline interpolation of the NH_Zone 2 dataset at a quaterly basis
nhz2=ts((nz2$y-1)*1000,start=1966,freq=4) #Transformation into a time-series object
m=ets(nhz2) #Fits an exponential smoothing state space model to the time series
f2=forecast(m,h=2*4) #Uses the fitted model to forecast 2 years into the future

#Bind observed and forecasted curves
bc=data.frame(Year=c(bombcurve[-dim(bombcurve)[1],1],seq(2009.75,2011.5,by=0.25)),
              Delta14C=c(bombcurve[-dim(bombcurve)[1],2],as.numeric(f2$mean)))
################################################################################

### One pool estimation

vtF=Vectorize(turnoverFit, vectorize.args=c("obsC14","obsyr","yr0","Fatm"))

#Subset the data for the different land uses
Nat=subset(and,LandUse=="NAT")
Pas=subset(and,LandUse=="PAS")
Euc=subset(and,LandUse=="EUC")

#-------------------------------------------------------------------------------------------------
# WARNING: This part if the code doesn't run. The vectirize function can't run for all data points in the dataset. It needs to be run for each individual data point.
#Calculate the turnover times
#tNat=vtF(obsC14=Nat[,5],obsyr=2009,In=(LF+RI)*Nat[,9],C0=Cnat*Nat[,9],yr0=-10000,by=1,plot=FALSE)
tPas=vtF(obsC14=Pas[,5],obsyr=2009,In=(LF+RI)*Pas[,9],C0=Cpas*Pas[,9],yr0=-10000,by=1,plot=FALSE)
tEuc=vtF(obsC14=Euc[,5],obsyr=2009,In=(LF+RI)*Euc[,9],C0=Ceuc*Euc[,9],yr0=-10000,by=1,plot=FALSE)
tNat=vtF(obsC14=Nat[,5],obsyr=2009,yr0=-10000,Fatm=bc,by=1,plot=FALSE)

#Put the data together
Nat=cbind(Nat,tau1=as.numeric(tNat[1,]),tau2=as.numeric(tNat[2,]))
Pas=cbind(Pas,tau1=as.numeric(tPas[1,]),tau2=as.numeric(tPas[2,]))
Euc=cbind(Euc,tau1=as.numeric(tEuc[1,]),tau2=as.numeric(tEuc[2,]))

Out=rbind(Nat,Pas,Euc)

#Take the output and save it as csv file. Change the location to a specific directory in your computer.
write.csv(Out,"~/Documents/Teach/Radiocarbon/Crow/AndOut.csv",row.names=FALSE)
#--------------------------------------------------------------------------------------------------


############################
# Multiple-pool parameter estimation. Method 4

#Subset the data for method 4
M4N=subset(and,Fract_Method==4 & LandUse=="NAT")

#Plot the data for the different fractions
plot(Hua2013$NHZone2[,1:2],type="l",ylim=c(-200,1000), ylab="Delta14C (per mil)")
points(rep(yrsample,5),M4N$Delta14C[-1],pch=20,col=1:5)
legend("topright",as.character(M4N$Fraction[-1]),pch=20,col=1:5,bty="n")

#Run the model with some first-guess parameter values
Natmod=CrowModel(
  t=years,
  ks=c(k1=1/2, k2=1/300, k3=1/150, k4=1/500, k5=1/1000),
  In=RI,
  gamma=c(gamma1=0.7, gamma2=0.2),
  alpha=c(alpha41=0.2, alpha52=0.05,alpha53=0.05),
  inputFc=bc
)
# R14m=getF14R(Natmod)
# C14m=getF14C(Natmod)
C14t=getF14(Natmod)

#Plot model and data
plotC14Pool(years,C14t,bc,col=2:6,ylim=c(-200,1000),xlim=c(1900,2010))
points(rep(yrsample,5),M4N$Delta14C[-1],pch=20,col=2:6)
arrows(rep(yrsample,5),M4N$Delta14C[-1]-M4N$D14Cerror[-1],rep(yrsample,5),M4N$Delta14C[-1]+M4N$D14Cerror[-1],
       angle=90,length=0.1,col=2:6)
legend("topright",as.character(M4N$Fraction[-1]),pch=20,col=2:6,bty="n")

#Create a function for optimization
M4func=function(pars){
  mod=CrowModel(t=years,ks=pars[1:5],In=RI,gamma=c(gamma1=pars[9], gamma2=pars[10]),
      alpha=c(pars[1]*pars[6],pars[2]*pars[7],pars[3]*pars[8]),inputFc=bc,pass=TRUE)
  C14t=getF14(mod)
  Ct=getC(mod)
  return(data.frame(Fractions=c(1:5),C14t=C14t[years==yrsample], Ct=Ct[years==yrsample]))
}

#Prepare the data for optimization
obsM4NCt14=data.frame(Fractions=c(1:5),C14t=M4N$Delta14C[-1],sdC14=M4N$D14Cerror[-1])
obsM4NCt=data.frame(Fractions=c(1:5),Ct=Cnat*M4N$Fractions[-1],sdC=Cnat*M4N$Frac_SE[-1])

#Cost fuction
costM4N=function(pars){
  modOut=M4func(pars)
  cost=modCost(model=modOut,obs=obsM4NCt14,err="sdC14",x="Fractions")
  return(modCost(model=modOut,obs=obsM4NCt,err="sdC",x="Fractions",cost=cost))
}

#Initial guess of parameter values
inipars=c(k1=1/2, k2=1/300, k3=1/150, k4=1/500, k5=1/1000, alpha41=0.2, 
          alpha52=0.05,alpha53=0.05, gamma1=0.7, gamma2=0.2)

#Fit model to data. Takes some time to run.
M4Nfit=modFit(f=costM4N,p=inipars,method="Newton", lower=rep(0,10),upper=c(3,rep(1,9)))

#Best parameter set and model
bestpar=M4Nfit$par
bestmod=CrowModel(t=years,ks=bestpar[1:5],In=RI,gamma=c(gamma1=bestpar[9], gamma2=bestpar[10]),
                  alpha=c(bestpar[1]*bestpar[6],bestpar[2]*bestpar[7],bestpar[3]*bestpar[8]),inputFc=bc)

bestC14t=getF14(bestmod)

#Plot with best parameter set5
par(mar=c(5,5,4,2))
plotC14Pool(years,bestC14t,bc,col=2:6,ylim=c(-200,1000),xlim=c(1900,2010),ylab=expression(paste(Delta^14,"C ","(\u2030)")))
points(rep(yrsample,5),M4N$Delta14C[-1],pch=20,col=2:6)
arrows(rep(yrsample,5),M4N$Delta14C[-1]-M4N$D14Cerror[-1],rep(yrsample,5),M4N$Delta14C[-1]+M4N$D14Cerror[-1],
       angle=90,length=0.1,col=2:6)
legend("topright",as.character(M4N$Fraction[-1]),pch=20,col=2:6,bty="n")

bestCt=getC(bestmod)
matplot(years,bestCt,type="l",lty=1,ylim=c(0,40),ylab="Carbon stores (Mg C/ha)")
points(rep(2010,5),Cnat*M4N$Fractions[-1],col=1:5,pch=20)
legend("topright",as.character(M4N$Fraction[-1]),pch=20,col=2:6,bty="n")

bestout=data.frame(Year=years,AtmC14=bc)
bestout=cbind(years,bestCt,bestC14t)
colnames(bestout)<-c("Years","C_FLF","C_OLF","C_1.8-2.0","C_2-2.4","C_>2.4",
                     "C14_FLF","C14_OLF","C14_1.8-2.0","C14_2-2.4","C14_>2.4")
write.csv(bestout,"~/Documents/Teach/Radiocarbon/Crow/bestout.csv")

#Bayesian optimization
var0 <- M4Nfit$var_ms_unweighted
cov0 <- summary(M4Nfit)$cov.scaled # The covariance matrix can be used for the jump, but wasn't used in this example.
MCMC <- modMCMC(f=costM4N, p = M4Nfit$par, niter = 2500, jump = NULL, var0 = var0, wvar0 = 0, #updatecov = 50,
                lower=rep(0,10),upper=c(4,rep(1,9)))

#Uncertainty values for the parameters
load("~/Dropbox/SoilR_Crow/MCMC.RData")
summary(MCMC)
sm=summary(MCMC)
#Correlation between posterior distributions
pairs(MCMC,nsample=500)

save(MCMC,file="~/Documents/Teach/Radiocarbon/Crow/MCMC.RData")

#build a matrix with the mean values
A=-1*diag(sm[1,1:5])
A[4,1]<-sm[1,6]*sm[1,1]
A[5,2]<-sm[1,7]*sm[1,2]
A[5,3]<-sm[1,8]*sm[1,3]
IN=matrix(c(sm[1,9:10],1-sum(sm[1,9:10]),0,0),nrow=5, ncol=1)

#compute the age and transit time distributions
age=seq(0,1000)
AD=systemAge(A=A,u=as.numeric(IN),a=age)
TT=transitTime(A=A,u=as.numeric(IN),a=age)

#plot
plot(age,AD$systemAgeDensity,type="l")
abline(v=AD$meanSystemAge,col=2)
abline(v=AD$quantilesSystemAge[2],col=4)
plot(age,TT$transitTimeDensity,type="l")
abline(v=TT$meanTransitTime,col=2)
abline(v=TT$quantiles[2],col=4)

#AD of indiv pools
matplot(age,AD$poolAgeDensity,type="l",lty=1,col=1)

par(mfrow=c(3,2))
plot(age,AD$poolAgeDensity[,1],type="l")
plot(age,AD$poolAgeDensity[,2],type="l")
plot(age,AD$poolAgeDensity[,3],type="l")
plot(age,AD$poolAgeDensity[,4],type="l")
plot(age,AD$poolAgeDensity[,5],type="l")

