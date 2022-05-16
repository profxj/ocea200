setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ESCI167/Rscripts_lectures")

library(lmtest)
library(nlme)

##### Differencing
# Read the MLO CO2 data 
mlo = read.csv("co2_mm_mlo.csv")
names(mlo) = c("year","month","dd","average")

# Create a time series (ts) object from the CO2 data
co2 = ts(data = mlo$average, frequency = 12, start = c(1958, 3))

co2_d1 = diff(co2,differences=1)
plot(co2_d1, ylab = expression("CO"[2]~" first differences"))
#still shows a trend

co2_d2 = diff(co2,differences=2)
plot(co2_d2, ylab = expression("CO"[2]~" second differences"))
#trend gone, but seasonal cycle is there

co2_d12 = diff(co2_d2,differences=12)
plot(co2_d12, ylab = expression("CO"[2]~" seasonal differences"))
#now appears stationary


##### GLS
#simulating a time series
set.seed(6000)
x = arima.sim(n = 50, model = list(order = c(1, 0, 0), ar=0.8, sd = 0.1)) #+ seq(1,100,1)*0.1

#export it as a .txt file
write.table(x,file="example_gls_W6_L12.txt",quote=F,col.names=F,row.names=F)
plot(x)
# Say x is a climate anomaly that we monitor, we want to determine whether a trend is present.

#Start with ordinary least squares
time=seq(1,50,1)
ols.m = lm(x~time)
abline(ols.m,col="red")
summary(ols.m)
dwtest(ols.m)
#We find a significant trend, but also the resioduals are not independent

#Let's look at their structure
acf(ols.m$residuals)
pacf(ols.m$residuals)
# Seems AR(1)

gls.m = gls(x~time,correlation = corARMA(p=1))
summary(gls.m)
#The trend is no longer significant

#Visualize the difference
plot(x)
abline(ols.m,col="red")
abline(gls.m,col="blue")