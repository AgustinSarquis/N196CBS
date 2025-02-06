# Let's get the entire ecosystem C dynamics from the soil and forest models to compute CS and CBS

TEC=as.data.frame(cbind(t=(2025:2100),
  BAU=forestvalues$TCS[60:136], 
  TEC2=Ctb[25:100,1], # total ecosystem C after tree harvesting is only soil C with biochar ammendments
  forest=TCS_values2[60:135, 5]
  ))
TEC=cbind(TEC, TEC1=rowSums(TEC[,c(2,4)])) #total ecosystem C with an intact forest are the soil and forest models together
#write.csv(TEC, 'C:/Users/asarq/Documents/TEC.csv')

TEC<-read.csv("~/tmp/TEC.csv")

# Plot the results
matplot(TEC$t,TEC[,c(4,6)], type="l", lty=1,lwd=3, col=c(2,3),ylim=c(0,1000), 
        ylab="Carbon stocks (Mg C/ha)", xlab="Year")
legend("topleft",c("After tree harvesting", "Intact plantation"),lty=1,col=c(2,3), lwd=3, bty="n")
polygon(c(TEC$t, rev(TEC$t)), c(TEC[,6], rep(0, length(TEC[,5]))),
        density = 20, angle = 45, border = NA, col = "green")
polygon(c(TEC$t, rev(TEC$t)), c(TEC[,4], rep(0, length(TEC[,3]))),
        density = 20, angle = 135, border = NA, col = "red")

# Compute Carbon Sequestration (CS) as the area under the curves at the end of the simulation (year 2100)
# For the intact plantation
t0<-head(TEC$t, 1)
tf<-tail(TEC$t, 1)
CS1 <- integrate(splinefun(TEC$t, TEC$TEC1), lower=t0, upper=tf)$value # 69309.1 Mg C ha-1 * y 
CS1 # correct value
square1<-((tf-t0)* max(TEC$TEC1)) # test accuracy approximating to the area of a rectangle
CS1-square1 # Loss of mass*time in comparison to keeping the system completely intact since t0

# For the plantation after being entirely harvested
CS2 <- integrate(splinefun(TEC$t, TEC$TEC2), lower=t0, upper=tf)$value # 39155.2 Mg C ha-1 * y 
square2<-((tf-t0) * max(TEC$TEC2)) # test accuracy approximating to the area of a rectangle
CS2-square2


# function to get CS for specific time horizons
CS=function(x,y, t){integrate(splinefun(x, y), lower=t0, upper=t)$value} 
# CS for year 1 to 76 (years 2025 to 2100)
CSt1=unlist(sapply(TEC$t, FUN=CS, x=TEC$t, y=TEC$TEC1)) # intact plantation
CSt2=unlist(sapply(TEC$t, FUN=CS, x=TEC$t, y=TEC$TEC2)) # harvested plantation
CSdf=as.data.frame(cbind(year=TEC$t, CSt1, CSt2))
matplot(CSdf$year,CSdf[,c(2,3)], type="l", lty=1,lwd=3, col=c(3,2), 
        ylab="Carbon Sequestration (Mg C/ha y)", xlab="Year")
legend("topleft",c("Intact plantation", "After tree harvesting"),lty=1,col=c(3,2), lwd=3, bty="n")

# Now calculate CBS (Climate Benefit of Sequestration)
# Radiative efficiency of one MgC in W m-2 (from Joos et al. 2013)
RE1Mg=5.35*(1/389)*(1/2.123e9)
# Impulse response function (eq. 2 in Sierra et al. 2021), fate of C entering the atmosphere
IRF_func=function(t, a0, as, tau){ 
  a0+sum(as*exp(-t/tau))
}
# Impulse response function modified after Millar et al. 2017
IRF_PD100=Vectorize(function(x){IRF_func(x, a0=0, as=c(0.2173,0.2240, 0.2824, 0.2763), tau=c(1000000,394.4, 36.54, 4.304))})
# convolution between the impulse response function of atmospheric CO2 and the carbon returning from ecosystems to the atmosphere
convolutionfun=function(t, t0=0, f, g){
  function(t){integrate(function(u,t){f(t-u)*g(u)}, lower=t0, upper=t, t)$value}
} 
# CBS function
CBSfun=function(TH, t0, kCO2, ha, smrfun){
  smrfun=Vectorize(smrfun)
  cv=convolutionfun(t, f=ha, g=smrfun,t0=t0)
  cv=Vectorize(cv)
  function(TH){-kCO2*(integrate(cv,lower=t0,upper=TH)$value)}
}

# For the intact plantation
smr1=splinefun(TEC$t, c(0,diff(TEC$TEC1)))  
CBS_Tr_fun=CBSfun(TH=TEC$t,t0=2025,kCO2 = RE1Mg, ha=IRF_PD100, smrfun = smr1)
CBS_Tr1<-sapply(TEC$t, FUN=CBS_Tr_fun)
# For the harvested plantation (emissions from fuels are subtracted in the dataframe (0.0297 C Mg ha-1 y-1))
TEC[1,4]=TEC[1,4]-0.0297
smr2=splinefun(TEC$t, c(0,diff(TEC$TEC2)))  
CBS_Tr_fun=CBSfun(TH=TEC$t,t0=2025,kCO2 = RE1Mg, ha=IRF_PD100, smrfun = smr2)
CBS_Tr2<-sapply(TEC$t, FUN=CBS_Tr_fun)

# plot
plot(TEC$X,CBS_Tr1,type="l",col = "green", lwd=3,
     ylab=expression(paste("CBS ("," W ", ha^-1, " yr)")), 
     xlab="Time horizon (yr)")
lines(TEC$X,CBS_Tr2,col='red',lwd=3)
abline(0, 0, lty='dashed')
legend("topleft",c("After tree harvesting", "Intact plantation"),lty=1,col=c('red','green'), lwd=3, bty="n")



