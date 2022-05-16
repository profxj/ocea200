setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ESCI167/Rscripts_lectures")

library(lmtest)
library(nlme)

# Firt we will download and plot the data
data = read.csv("PM25West.csv")
plot(data$Year,data$Mean,type="l",xlab="Time",ylab="PM2.5")

#Fit a linear model and visualize
m1 = lm(data$Mean~data$Year)
summary(m1)
abline(m1, col="red")

# Durbin-Watson test
dwtest(m1)
# Independence is rejected

# Let's look at the acf and pacf
acf(m1$residuals)
pacf(m1$residuals)
# but the acf and pacf do not suggest autocorrelation

# Need to also check normality and constant variance
shapiro.test(m1$residuals)
bartlett.test(m1$residuals,c(rep(1,10),rep(2,10)))
#these assumptions are respected

# What do we do?

#Option 1: Stick with the OLS model fitted above, Report the results of the Durbin-Watson test
#and report the acf and pacf plots to support sticking with this approach despite 
#the DW test results.

#Option 2: Use a GLS model:
x = data$Mean
time = data$Year
m2 = gls(x~time,correlation = corARMA(p=1))
summary(m2)
# The trend estimate is essentially the same and still significant with GLS.

plot(data$Year,data$Mean,type="l",xlab="Time",ylab="PM2.5")
abline(m1, col="red")
abline(m2,col="blue")
legend("topright",c("OLS","GLS"),col=c("red","blue"),lty=c(1,1))

# Using either option 1 or 2, doesn't affect the conclusions here.