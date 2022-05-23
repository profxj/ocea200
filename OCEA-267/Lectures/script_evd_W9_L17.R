setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ESCI167/Rscripts_lectures")
#####
# This script provides an example of extreme value change analysis
# for the sea surface temperature time series collected at
# the Scripps Pier.
#####

# load packages
library(extRemes)
library(distillery)

##### 
# Data preparation
#####

# Reads the ENSO data
nino = read.table("enso_soi_ts.txt")

# Extract annual means 
nino_am = rowMeans(nino[48:149,2:13])

#load the Scripps data 
scripps = read.csv("SIO_TEMP_1916_201905.csv", skip = 26, header = T)

# get rid of the extra columns and obs taken in 1916 and 2019 (those years are incomplete)
scripps = subset(scripps, YEAR > 1916 & YEAR < 2019, select = -c(10:14))

# histogram of the daily values
hist(scripps$SURF_TEMP_C,xlab="Daily SST (°C)",main=" ")

# get annual maximums and annual means - initialise
year = min(scripps$YEAR):max(scripps$YEAR) 
ny = length(year)
max_sst = vector(mode="numeric",length=ny)
mean_sst = vector(mode="numeric",length=ny)

# loop through years
k = 1
for (i in 1:ny){
  cy = scripps$YEAR == year[i] 
  max_sst[k] = max(scripps$SURF_TEMP_C[cy],na.rm=T) # mt is the SSt maximum for each year
  mean_sst[k] = mean(scripps$SURF_TEMP_C[cy],na.rm=T) # meant is the SST mean for each year
  k = k+1
}

hist(max_sst,xlab="Max SST (°C)",main=" ")

# create a data frame, that has year, a trend, max temp and Nino 3.4SST Index
lintr = 1:ny
scripps_nino = data.frame(year, max_sst, nino_am, lintr)
names(scripps_nino) = c("Year", "Max.SST", "Nino3.4", "Trend")

#####
# GEV analysis
#####

# Fits a GEVD to the data of the scripps dataset - this is the "H0" fit with no change
fit_h0 = fevd(Max.SST, scripps_nino, units = "degC")
ci_h0 = ci(fit_h0, type="parameter") # confidence interval on the three GEV params

# Fits a GEV model with linear trend in location param
fit_loctr = fevd(Max.SST, scripps_nino, location.fun = ~Trend,
             units = "degC")
ci_loctr = ci(fit_loctr,type="parameter")
# Trend (mu1) significantly different from zero

# Likelihood ratio test comparing the two fits
lr.test(fit_loctr,fit_h0)
# consistent with ci intervals above

#return periods
T=return.level(fit_loctr)
T=return.level(fit_loctr,return.period = c(2,20,100,500))

#visualize return periods
plot(scripps_nino$Year,scripps_nino$Max.SST,ylim=c(20,27),type="l",xlab="Time",ylab=" Max SST (°C)")
lines(scripps_nino$Year,T[,1],col="blue",lwd=2)
lines(scripps_nino$Year,T[,2],col="gold",lwd=2)
lines(scripps_nino$Year,T[,3],col="red",lwd=2)
labels = c("T=2","T=20","T=100")
colors = c("blue","gold","red")
legend("bottomright",labels,lwd=c(3,3,3),col=colors)

# We can do a likelihood ratio test here because the models are nested
# but we could also rely on AIC or BIC
fit_h0  #AIC is 309.1 (and -loglik is 151.6)
fit_loctr#AIC is 288.3 (and -loglik is 140.1)

# LR test by hand
lambda = 2*(-140.1 - -151.6)
#lambda follows a chisq dist with 1 df

# plot the chisq(1) distribution
x = seq(0,30,0.1)
dx = dchisq(x,df=1)
plot(x,dx,type="l")
abline(v=lambda,col="red")
pval = pchisq(lambda,df=1,lower.tail=F)
####

# Fit a GEV model with ENSO as a covariate
fit_enso = fevd(Max.SST, scripps_nino, location.fun = ~Nino3.4,
             units = "deg C")
ci_enso = ci(fit_enso,type="parameter")
lr.test(fit_enso,fit_h0)
# ENSO not a good covariate here

##### 
# Plots
#####

# Plot of daily, annual mean and annual max SST
par(mfrow=c(1,2))
plot(scripps$YEAR,scripps$SURF_TEMP_C,type="l",xlab="Year",ylab="SST (°C)",ylim=c(10,33),main="a)")
lines(year,max_sst,type="l",col="red",lwd=3)
lines(year,mean_sst,type="l",lwd=3)
labels=c('Daily','Annual mean',"Annual maximum")
colors=c("black","black","red")
legend("topleft",labels,lwd=c(1,3,3),col=colors)

# visualize the change in location parameter by plotting the distribution start vs end of the record
xt = seq(from=20,to=28,by=0.1)
xd1 = devd(xt,loc=ci_loctr[1,2],scale=ci_loctr[3,2],shape=ci_loctr[4,2],type="GEV")
xd2 = devd(xt,loc=ci_loctr[1,2]+ci_loctr[2,2]*102,scale=ci_loctr[3,2],shape=ci_loctr[4,2],type="GEV")
plot(xt,xd1,type="l",lwd=3,xlab=" Max SST (°C)",ylab="Density",col="orange",ylim=c(0,0.45),main="b)")
lines(xt,xd2,type="l",lwd=3,col="red")
labels = c('1917',"2018")
colors = c("orange","red")
legend("topright",labels,lwd=c(3,3),col=colors)
arrows(22.2, 0.42, 23.8, 0.42,lwd=3,lty=6)
