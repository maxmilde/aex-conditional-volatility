library(tidyverse)
library(tseries) #adf.test, kpss.test, jarque.bera.test
library(xts)
library(forecast) #auto.arima, Acf, Pacf, checkresiduals
library(moments) #kurtosis, skewness
library(aTSA) #arch.test
library(rugarch)

rm(list=ls())


#Load/clean data
data_raw <- read.csv("^AEX.csv")
head(data_raw)
nrow(data_raw)
summary(data_raw)

data_raw$Index <- as.Date(data_raw$Index)
data_raw <- data_raw %>% rename(price=AEX.Close)

data_raw[is.na(data_raw$price),] #2nd Jan was a public holiday - safe to remove NA
data <- na.omit(data_raw)

ggplot(data, aes(x=Index, y=price)) + 
  geom_line() + 
  labs(title="AEX Close Price", x="Date", y="Close") +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))


#Returns
price <- xts(data$price, order.by=data$Index)
log_ret <- na.omit(diff(log(price)))

plot(log_ret, ylab="Log ret", xlab="Date", type="l")
title(main = "Log returns AEX", line = 2)
#shows volatility clustering

#Stationarity test
tseries::adf.test(price) #ADF rejects a unit root, but conflicts with KPSS; likely a near-integrated process
tseries::kpss.test(price) #strongly rejects stationarity (contradicts adf)

tseries::adf.test(log_ret) #non-stationarity rejected
tseries::kpss.test(log_ret) #stationarity not rejected

Acf(log_ret, lag.max=50, main="ACF")
Pacf(log_ret, lag.max=50, main="PACF") #there is autocorrelation present (but not systematic)


#Modelling mean
(fit_auto <- auto.arima(as.numeric(log_ret), ic=c("bic"), stationary=TRUE))
(fit_manual11 <- arima(as.numeric(log_ret), order=c(1,0,1)))
(fit_manual22 <- arima(as.numeric(log_ret), order=c(2,0,2)))
(fit_manual55 <- arima(as.numeric(log_ret), order=c(5,0,5))) #no manually fitted model provides stable or meaningful results that are not redundant
#coefficients cancel each other - overparametrization

#Acf(fit_manual11$residuals)
#Pacf(fit_manual11$residuals)
#Acf(fit_manual22$residuals)
#Pacf(fit_manual22$residuals)
#Acf(fit_manual55$residuals)
#Pacf(fit_manual55$residuals)
#no systematic reduction in autocorrelations

(fit_manual00 <- arima(as.numeric(log_ret), order=c(0,0,0))) 
#not plotting acf/pacf of residuals of the ARMA(0,0) fit as the graphs would not be any different
#residual autocorrelation likely due to heteroskedasticity


#Ljung-Box: test whether ACF at first 25 lags is jointly 0
Box.test(fit_manual00$residuals, lag=25, type=c("Ljung-Box"), fitdf=length(fit_manual00$coef)) #null of no resid. autocorr. rejected at 5%; lag=25
checkresiduals(fit_manual00) #no significant evidence of autocorrelation; lag=10

#Box.test(fit_manual11$residuals, lag=25, type=c("Ljung-Box"), fitdf=length(fit_manual11$coef))
#Box.test(fit_manual22$residuals, lag=25, type=c("Ljung-Box"), fitdf=length(fit_manual22$coef))
#Box.test(fit_manual55$residuals, lag=25, type=c("Ljung-Box"), fitdf=length(fit_manual55$coef))
#large order ARMA leads to stronger rejection in Box test an higher order leads to reduction in df -> larger model is worse


#normality of residuals
resid00 <- fit_manual00$residuals
resid00_norm <- resid00/sqrt(fit_manual00$sigma2)
qqnorm(resid00_norm)
qqline(resid00_norm) #far from normal

jarque.bera.test(resid00) #strong rejection of the null of normality
kurtosis(resid00) 
#kurtosis(log_ret) #same since we fit arma(0,0)
skewness(resid00)


#Summary statistics for the return series
ret <- as.numeric(log_ret)

summary_stats <- c(
  observations = length(ret),
  mean = mean(ret),
  sd = sd(ret),
  min = min(ret),
  max = max(ret),
  skewness = skewness(ret),
  kurtosis = kurtosis(ret),
  excess_kurtosis = kurtosis(ret) - 3)

summary_stats
jarque.bera.test(ret)



#Modelling volatility
resid_sqrd00 <- resid00^2
plot.ts(resid_sqrd00, ylab="Squared residuals") #definite volatility clustering

Acf(resid_sqrd00, lag.max=50, main="ACF")
Pacf(resid_sqrd00, lag.max=50, main="PACF")
#obvious serial dependence in the second moment -> GARCH

arch.test(fit_manual00) #heteroscedasticity, strong rejection of null of no resid^2 autocorrelation -> ARCH effects
Box.test(resid_sqrd00, lag=25, type=c("Ljung-Box")) #null of no autocorrelation rejected -> ARCH effects


#Estimating GARCH (normal distrib errors)
garch11 <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                      variance.model=list(garchOrder=c(1,1)))
fit11 <- ugarchfit(spec=garch11, data=log_ret)
fit11 #no ARCH effects remaining at lower lags; significant Sign Bias: leverage effect -> fit TARCH; Goodness-of-Fit: normality rejection
#high persistence
#report on statistics in paper 


#Persistence (alpha + beta)
coef11 <- coef(fit11)
coef11["alpha1"] + coef11["beta1"]


#Plotting estimated volatility
residAG <- as.numeric(residuals(fit11, standardize=FALSE)) #resid from mean eq.
plot(abs(residAG), type="l", col="grey85")
lines(array(sigma(fit11)), col="blue") #estimated sigma
abline(h=sqrt(fit_manual00$sigma2), col="black")

#adding 20-day rolling mean of absolute residuals 

w=20 #window


rolling_vol <- stats::filter(abs(residAG), rep(1/w, w), sides=1)
rolling_vol <- c(NA, head(as.vector(rolling_vol), -1))  #lag by one to avoid look-ahead bias


lines(rolling_vol, col="red")
legend("topleft", legend=c("Absolute residuals", "GARCH sigma", "Unconditional SD (ARMA)", "20-day rolling mean"),
       col=c("grey", "blue", "black", "red"), lty=1)


#Checking for remaining ARCH effects
residAG_stand <- residuals(fit11, standardize=TRUE)
sqr_residAG_stand <- residuals(fit11, standardize=TRUE)^2

Acf(residAG_stand, lag.max=30, main="ACF standardized residuals")
Pacf(residAG_stand, lag.max=30, main="PACF standardized residuals")

Acf(sqr_residAG_stand, lag.max=30, main="ACF squared standardized residuals")
Pacf(sqr_residAG_stand, lag.max=30, main="PACF squared standardized residuals")
#some ARCH effects appear to remain at high lags


#Checking for normality of standardised residuals
hist(residAG_stand, freq=FALSE, nclass=40)
curve(dnorm(x, mean=0, sd=1),
      col="red", lwd=2, add=TRUE)
jarque.bera.test(residAG_stand) #fat-tails
kurtosis(residAG_stand) #kurtosis dropped from fitting garch -> still above 4, student-t errors more appropriate
skewness(residAG_stand)
#normality of resid still rejected but distrib. of innovations much closer to normal


#GARCH with t-distrib
garch11_t <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                        variance.model=list(garchOrder=c(1,1)),
                        distribution.model="std")
fit11_t <- ugarchfit(spec=garch11_t, data=log_ret)
fit11_t
#alpha no longer signif under robust errors
#t-distrib wins over normal distrib under AIC and BIC and GoF
#signif. shape parameter - fat tails
#nyblom - parameters not stable over time (omega) - structural breaks
cbind(infocriteria(fit11), infocriteria(fit11_t))


#Persistence (alpha + beta)
coef11_t <- coef(fit11_t)
coef11_t["alpha1"] + coef11_t["beta1"]


plot(array(sigma(fit11)), col="blue", type="l")
lines(array(sigma(fit11_t)), col="lightblue")
#no visible change in estimates -> distribution of errors has little effect on vola estimates (QMLE)


#Checking if another GARCH specification would be more appropriate
p_max <- 5
q_max <- 5
ic_min <- Inf
best_p <- 0
best_q <- 0

for (i1 in 1:p_max){
  for (i2 in 1:q_max){
    spec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                       variance.model=list(garchOrder=c(i1,i2)),
                       distribution.model="std")
    fit_temp <- tryCatch(suppressWarnings(ugarchfit(spec=spec, 
                                                    data=log_ret, 
                                                    solver="hybrid")),
                         error=function(e) NULL)
    bic <- infocriteria(fit_temp)[2]
    if (bic < ic_min){
      ic_min <- bic
      best_p <- i1
      best_q <- i2
    }
  }
}
c(best_p, best_q) #we have selected the best GARCH specification according to BIC


#TARCH: GJR-GARCH
tarch_spec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                         variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                         distribution.model="std")
fit_tarch <- ugarchfit(spec=tarch_spec, data=log_ret)
fit_tarch


#gamma not signif. under robust errors
#alpha=0 -> only negative shocks impact next period's variance
#however, this is inconsistent with the significant positive bias and insignificant gamma coef.
#beta highly signif.
#Ljung-Box, ARCH-LM insignificant -> no remaining ARCH effects
cbind(infocriteria(fit11), infocriteria(fit11_t), infocriteria(fit_tarch))
#tarch best based on all info criteria


#Persistence (alpha + beta + gamma/2)
coef_tarch <- coef(fit_tarch)
coef_tarch["alpha1"] + coef_tarch["beta1"] + coef_tarch["gamma1"] / 2

#Likelihood ratio test: GARCH(1,1)-t vs GJR-GARCH(1,1)-t (nested, df=1)
lr_stat <- 2 * (likelihood(fit_tarch) - likelihood(fit11_t))
lr_pval <- pchisq(lr_stat, df=1, lower.tail=FALSE)
lr_stat
lr_pval



#News Impact Curve: comparing symmetric and asymmetric volatility response
ni_garch <- newsimpact(fit11_t)
ni_tarch <- newsimpact(fit_tarch)

plot(ni_tarch$zx, ni_tarch$zy, type = "l", col = "red", lwd = 2,
     main = "News Impact Curve", xlab = "z(t-1)", ylab = "sigma^2(t)")
lines(ni_garch$zx, ni_garch$zy, col = "blue", lwd = 2)

legend("top",
       legend = c("GJR-GARCH(1,1) Student-t", "GARCH(1,1) Student-t"),
       col = c("red", "blue"), lty = 1, lwd = 2)

plot(abs(residAG), type="l", col="grey85")
lines(array(sigma(fit_tarch)), col="red")
lines(array(sigma(fit11_t)), col="blue")
legend("topleft", legend=c("Absolute residuals", "TARCH Student-t", "GARCH Student-t"),
       col=c("grey85", "red", "blue"),
       lty=1)



#Diagnostics for the preferred model (GJR-GARCH with student-t errors)
std_tarch <- as.numeric(residuals(fit_tarch, standardize = TRUE))

#Serial dependence in standardized residuals
Acf(std_tarch, lag.max = 30, main = "ACF Std. Residuals (GJR-GARCH)")
Pacf(std_tarch, lag.max = 30, main = "PACF Std. Residuals (GJR-GARCH)")

#Remaining ARCH effects in squared standardized residuals
Acf(std_tarch^2, lag.max = 30, main = "ACF Squared Std. Residuals (GJR-GARCH)")
Pacf(std_tarch^2, lag.max = 30, main = "PACF Squared Std. Residuals (GJR-GARCH)")

#Ljung-Box tests on residuals and squared residuals
Box.test(std_tarch, lag = 25, type = "Ljung-Box")
Box.test(std_tarch^2, lag = 25, type = "Ljung-Box")

#Distributional diagnostics
jarque.bera.test(std_tarch)
skewness(std_tarch)
kurtosis(std_tarch)

qqnorm(std_tarch, main = "QQ-Plot: GJR-GARCH Std. Residuals")
qqline(std_tarch)

