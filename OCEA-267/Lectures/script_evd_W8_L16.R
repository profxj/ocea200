setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ESCI167/Rscripts_lectures")

library(evd)

####
data=read.csv("usgs_03335500_peakflow.csv",header=T)
names(data)=c("Date","Pflow")

mflow=data$Pflow
year=seq(1981,2015)

options(scipen = 999)# just doing this to avoid scientific notation on my plots

# Plot the time series
plot(year,mflow,type="h",xlab="",ylab="Annual peak streamflow (cfs)")

# Plot the empirical CDF
x.ord=sort(mflow)
n=length(x.ord)
F_x=(1:n)/(n+1)
plot(x.ord, F_x, type = "s", ylim = c(0, 1), xlab="Annual peak streamflow (cfs)", 
     ylab = "F(x)", main = "",lwd=2)

# Plot the return periods
T=1/(1-F_x)
plot(T,x.ord,ylim=c(10000,90000),xlim=c(0,60),col="black",xlab="Return period (yrs)",
     ylab="Annual peak streamflow (cfs)",lwd=2)

#estimate EVD parameters
m1 = fgev(mflow)

#calculate the GEV CDF
x = seq(10000,90000,1)
F1_x = pgev(x,loc=44574.4700,scale=13492.2966,shape=-0.1604)

#Gumbel CDF
m2 = fgumbel(mflow)
F2_x = pgumbel(x,loc=44574,scale=13492)

# superpose GEV fit to empirical CDF
plot(x.ord, F_x, type = "s", ylim = c(0, 1), xlab="Annual peak streamflow (cfs)", 
     ylab = "F(x)", main = "",lwd=2,cex.lab=1.5,cex.axis=1.5)
lines(x,F1_x,col="lightblue",lwd=2)
lines(x,F2_x,col="gold",lwd=2)
labels=c('Empirical CDF',"GEV","Gumbel")
colors=c("black","lightblue","gold")
legend("topleft", inset=.05, labels,lty=c(1,1,1),lwd=c(2,2,1),col=colors)

#probability of exceedance and return period of the largest event
pe = 1-pgumbel(85700,loc=44574,scale=13492)
1/(pe)
