# Let's get the entire ecosystem C dynamics from the soil and forest models to compute CS and CBS
TEC=as.data.frame(cbind(t=2025:2100,
  BAU=forestvalues$TCS[61:136]+Ct[61:136,1], 
  Utopy= c(0,na.omit(kukuivalues[[6]]$TCS))+Ctb[,1] 
  ))
#write.csv(TEC, 'C:/Users/asarq/Documents/GitHub/N196CBS/Data/TEC.csv')
#TEC<-read.csv("~/tmp/TEC.csv")

# Plot the results
matplot(TEC$t,TEC[,c(2,3)], type="l", lty=1,lwd=3, col=c('#ff6b35','#53c98b'), ylim = c(0,2600),
        ylab="Carbon stocks (Mg C/ha)", xlab="Year")
legend("topleft",c("Eucalyptus plantation", "Pakukui agroforestry"),lty=1,col=c('#ff6b35','#53c98b'), lwd=3, bty="n")
polygon(c(TEC$t, rev(TEC$t)), c(TEC[,2], rep(0, length(TEC[,2]))),
        density = 20, angle = 45, border = NA, col = 4)
polygon(c(TEC$t, rev(TEC$t)), c(TEC[,3], rep(0, length(TEC[,3]))),
        density = 20, angle = 135, border = NA, col = 7)

# Compute Carbon Sequestration (CS) as the area under the curves at the end of the simulation (year 2100)
# For the business as usual scenario
t0=2025
tf=2100
CS1 <- integrate(splinefun(TEC$t, TEC$BAU), lower=t0, upper=tf)$value # 210992 Mg C ha-1 * y 
# For the utopian scenario
CS2 <- integrate(splinefun(TEC$t, TEC$Utopy), lower=t0, upper=tf)$value # 182208 Mg C ha-1 * y 

# function to get CS for specific time horizons
CS=function(x,y, t){integrate(splinefun(x, y), lower=t0, upper=t)$value} 
# CS for year 1 to 76 (years 2025 to 2100)
CSt1=unlist(sapply(TEC$t, FUN=CS, x=TEC$t, y=TEC$BAU)) 
CSt2=unlist(sapply(TEC$t, FUN=CS, x=TEC$t, y=TEC$Utopy)) 
CSdf=as.data.frame(cbind(year=TEC$t, CSt1, CSt2))
matplot(CSdf$year,CSdf[,c(2,3)], type="l", lty=1,lwd=3, col=c('#ff6b35','#53c98b'), 
        ylab="Carbon Sequestration (Mg C/ha y)", xlab="Year")
legend("topleft",c("Eucalyptus plantation", "Pakukui Agroforestry"),lty=1,col=c('#ff6b35','#53c98b'), lwd=3, bty="n")

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

# For the business as usual scenario
smr1=splinefun(TEC$t, c(0,diff(TEC$BAU)))  
CBS_Tr_fun=CBSfun(TH=TEC$t,t0,kCO2 = RE1Mg, ha=IRF_PD100, smrfun = smr1)
CBS_Tr1<-sapply(TEC$t, FUN=CBS_Tr_fun)
# For the utopian scenario 
# t=0 represents emissions from clearing eucalyptus and biochar production
# this includes emissions from fuel calculated by Darshi 
# COULD BE IMPROVED: it doesn't include CO2 equivalents
# COULD BE IMPROVED: I ended up using values for a density of 1598 m3 which matched AGB60 (wood density of 0.42 from Turn et al., 2014)
fuelC=-2.5 # C Mg ha-1 (one time emission)
# and also includes the C lost during biochar production calculated from biochar efficiency (eucClost)
eucClost=-(1-yield)*AGB60*0.482 # convert Eucalyptus biomass back to C
smr2=splinefun(TEC$t, c(fuelC+eucClost,diff(TEC$Utopy)))  
CBS_Tr_fun=CBSfun(TH=TEC$t,t0,kCO2 = RE1Mg, ha=IRF_PD100, smrfun = smr2)
CBS_Tr2<-sapply(TEC$t, FUN=CBS_Tr_fun)

# plot
plot(TEC$t,CBS_Tr1,type="l",col = '#ff6b35', lwd=3,
     ylab=expression(paste("CBS ("," W ", ha^-1, " yr)")), 
     xlab="Time horizon (yr)")
lines(TEC$t,CBS_Tr2,col='#53c98b',lwd=3)
abline(0, 0, lty='dashed')
legend("bottomleft",c("Eucalyptus plantation", "Pakukui agroforestry"),lty=1,col=c('#ff6b35','#53c98b'), lwd=3, bty="n")

# plot absolute difference between the 2, using BaU as the baseline
plot(TEC$t,abs(CBS_Tr1-CBS_Tr2),type="l",col = '#cb71ff', lwd=3,
     ylab=expression(paste("|CBS| ("," W ", ha^-1, " yr)")), 
     xlab="Time horizon (yr)")
abline(0, 0, lty='dashed')
