# Total ecosystem C per scenario starting in 2025
years=seq(from=2025, to=2100, by=1)
TEC = data.frame(years, Ct1[32:107,3], Ct2[,4])
colnames(TEC) <- c("years", "scenario1", "scenario2")

matplot(TEC$years,TEC[,2:3], type="l", lty=1,lwd=3, col=c(2,3), 
        ylab="Total Ecosystem C (Mg C/ha)", xlab="Years", ylim=c(0,1400))
legend("topleft",c("Business as usual", "Clearing and chipping"),lty=1,col=c(2,3), lwd=3, bty="n", cex = 0.7)
