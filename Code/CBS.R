# Let's calculate CS and CBS for a specific case: Hawaii case study from Crow & Sierra 2022

library(expm)
library(SoilR)

# Construct the A matrix with decay constants and transference coefficients
ks= c(k1=0.299, k2=0.306, k3=0.058)
alpha21=0.005; alpha12=0.008;alpha32=0.501; alpha23=0.036
A=diag(-ks)
A[2,1]=alpha21*ks[1]; A[1,2]=alpha12*ks[2]; A[3,2]=alpha32*ks[2]; A[2,3]=alpha23*ks[3]

####################### TRANSIT TIME AND SYSTEM AGE

TT=transitTime(A=A, u=c(2.8+0.72,0,0)) 
Age=systemAge(A=A, u=c(2.8+0.72,0,0))
# How long do new C inputs stay in the soil on average? (mean TT)
TT$meanTransitTime
# How long does half of the C in new inputs stay in the soil? (median TT)
TT$quantiles[2]

################################ CARBON SEQUESTRATION

# STEADY STATE SYSTEM
# This is for determining bulk CS 
# C remaining at steady state 
M=function(t,B,u){expm(t*B)%*%u} # Eq. 13 from Sierra et al. 2021
yr=seq(0,100, by=0.1) # Set time
Mt=t(sapply(yr,FUN=M, B=A, u=c(2.8+0.72,0,0))) # u from Crow (2016), root inputs + turnover and loss
totalMt=rowSums(Mt)
# C sequestration (at 5, 10 and 20 years)
# This answers the question: What is the amount of ecosystem C inputs stored over x amount of years?
CS=function(B,H, u){sum((solve(B)%*%(expm(H*B) - diag(1, dim(B)[1])))%*%(u))} # Eq. 19 from Sierra et al. 2021
CSt=sapply(c(5,10,20), FUN=CS, B=A, u=c(2.8+0.72,0,0))
# plot
plot(yr, totalMt, type="l", xlim=c(0,20), ylim=c(0,3.5), lwd=2,xlab="Time since C addition to soil (years)",  
     ylab=expression(paste("Mass remaining (Mg C ",ha^-1, ")" )), bty="n")
polygon(x=c(yr[1:i20], rev(yr[1:i20])), y=c(totalMt[1:i20], rep(0,i20)),col='hotpink')
abline(v=c(5,10,20),lty=2)

# This is for determining proportional CS
# C remaining at steady state as the proportion of C remaining after the t of ﬁxation (eq. 14)
M1=function(t,B,u){expm(t*B)%*%(u/sum(u))} 
M1t=t(sapply(yr,FUN=M1, B=A, u=c(1,0,0)))
totalM1t=rowSums(M1t)
# This answers the question: What is the amount of one unit of C stored over x amount of years?
CS1=sapply(c(5,10,20), FUN=CS, A=A, u=c(1,0,0))
# plot
plot(yr,totalM1t, type="l", xlim=c(0,20), ylim=c(0,1), lwd=2, xlab="Time since C addition to soil (years)", 
     ylab=expression(paste("Mass remaining (Mg C ",ha^-1, ")" )), bty="n")
i20=which(yr==20.0); i5=which(yr==5.0)
polygon(x=c(yr[1:i20], rev(yr[1:i20])), y=c(totalMt[1:i20], rep(0,i20)),col="hotpink")
abline(v=c(5,10,20),lty=2)

# NON STEADY STATE SYSTEM WITH TRANSIENT ACCUMULATION OF C
# For a trajectory of accumulated annual pulses
yr=seq(0,100, by=1)
Mt=t(sapply(yr,FUN=M, B=A, u=c(2.8+0.72,0,0)))
totalMt=rowSums(Mt)
cm=matrix(NA,nrow=length(totalMt),ncol=100)
for(i in 1:100){
  cm[,i]<-c(rep(NA,i-1),totalMt[1:(length(totalMt)-i+1)])
}
# CS cummulative pulses
CSpulses<-cumsum(sapply(yr, FUN=CS, B=A, u=c(2.8+0.72,0,0)))
CSpulses[which(yr==20 | yr==40 | yr==100)]
# plot
# The areas under the curve of each pulse accumulate the amount of C and the time it is retained in an ecosystem
# The areas under the curve of individual pulses can be summed over the time period of interest to obtain CS
plot(yr,totalMt,type="l", ylim=c(0,15),xlim=c(0,50), xlab="Time (yr)", ylab=expression(paste("Cummulative pulses (Mg C h", a^-1,")")),bty="n")
for(i in 2:50){
  lines(yr,rowSums(cm[,1:i]))
}
lines(yr,rowSums(cm,na.rm = TRUE),col=2, lty=2)

############################# CLIMATE BENEFIT OF SEQUESTRATION

# STEADY STATE SYSTEM
# convolution between the impulse response function of atmospheric CO2 and the carbon returning from ecosystems to the atmosphere
# convolution: a mathematical operation of 2 functions producing a 3rd function
convolutionfun=function(t, t0=0, f, g){
  function(t){integrate(function(u,t){f(t-u)*g(u)}, lower=t0, upper=t, t)$value}
} 
# CBS, equation 11 in Sierra et al. 2021
CBSfun=function(yr, t0, kCO2, ha, B, u, xi, gamma){ 
  rfun=function(t){sum((-colSums(xi*B)%*%expm(t*(xi*B))%*%(gamma*u)))}
  rfun=Vectorize(rfun)
  cv=convolutionfun(t, f=ha, g=rfun)
  integrand=function(t){-sum(gamma*u)*ha(t)+cv(t)}
  integrand=Vectorize(integrand)
  function(yr){kCO2*(integrate(integrand,lower=t0,upper=yr)$value)}
}
yr=seq(0,100, by=1)
# Radiative efficiency of one MgC in W m-2 (from Joos et al. 2013)
RE1Mg=5.35*(1/389)*(1/2.123e9)
# Impulse response function (eq. 2 in Sierra et al. 2021), fate of C entering the atmosphere
IRF_func=function(t, a0, as, tau){ 
  a0+sum(as*exp(-t/tau))
}
# Impulse response function modified after Millar et al. 2017
# using this IRF is preferred to avoid mathematical problems (see Sierra et al. 2021)
IRF_PD100=Vectorize(function(x){IRF_func(x, a0=0, as=c(0.2173,0.2240, 0.2824, 0.2763), tau=c(1000000,394.4, 36.54, 4.304))})
CBS=sapply(yr, CBSfun(t0=0,kCO2=RE1Mg, ha=IRF_PD100, B=A, u=c(2.8+0.72,0,0),xi=1,gamma=1))
# What is the amount of warming mitigated by soil C storage in the ecosystem after x amount of years (W m-2 year)?
CBS[which(yr==20 | yr==40 | yr==100)]
# absolute global warming potential (eq. 4 in Sierra et al. 2021)
AGWP=function(x){integrate(IRF_PD100, 0, x)$value * RE1Mg} 
# Computation of AGWP for different Time Horizons
agwp_yr=sapply(yr, AGWP)
# plot: CBS of one year of productivity over time
plot(yr, CBS*1e10, type="l", xlab="Time horizon (yr)",lwd=2,cex.axis=1.3, cex.lab=1.3,
     ylab=expression(paste("CBS (", 10^-10," W ",m^-2, " yr)")), ylim=c(-0.7, 0), col='hotpink', bty="n")
abline(v=c(10,40,100),lty=2, col="gray")
# plot: absolute CBS compared to the AGWP of 1 Mg C of C02 emission
plot(yr,agwp_yr*1e10,type="l",cex.axis=1.3, cex.lab=1.3,
     ylab=expression(paste("| CBS | or AGWP (", 10^-10," W ", m^-2, " yr)")),
     xlab="Time horizon (yr)", col=gray(0.2),lwd=2, bty="n")
lines(yr,abs(CBS)*1e10,col='hotpink',lwd=2)
abline(v=c(10,40,100),lty=2, col="gray")
legend("topleft", c("AGWP of an emission","Hawaiian Mollisol"), 
       lty=c(1,1,2), lwd=2, col=c(gray(0.2),'hotpink'), bty="n")

# NON STEADY STATE SYSTEM WITH TRANSIENT ACCUMULATION OF C
# Cummulative CBS of pulses
# Run transient simulations with zero initial conditions
times<-seq(0,100,by=0.1)
tr<-Model(t=times,A=A, ivList = rep(0,3), inputFluxes=c(2.8+0.72,0,0))
rTr<-rowSums(getReleaseFlux(tr))
smr<-splinefun(x=times, y=(2.8+0.72)-rTr)
CBSfun3=function(TH, t0, kCO2, ha, smrfun){
  smrfun=Vectorize(smrfun)
  cv=convolutionfun(t, f=ha, g=smrfun,t0=t0)
  cv=Vectorize(cv)
  function(TH){-kCO2*(integrate(cv,lower=t0,upper=TH)$value)}
}
CBS_Tr_fun=CBSfun3(TH=times,t0=0,kCO2 = RE1Mg, ha=IRF_PD100, smrfun = smr)
CBS_Tr<-sapply(times, FUN=CBS_Tr_fun)
# plot
plot(times,CBS_Tr,type="l")

# What is the amount of warming mitigated by soil C storage in the ecosystem after x amount of years (W m-2 year)?
CBS_Tr[which(times==20 | times==40 | times==100)]
