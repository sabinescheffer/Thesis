library('survival')
library('survminer')
library('ggplot2')
library('ggfortify')
library('ranger')
library('tidyr')
library('rjags')
library('morse')
library('flexsurvcure')

#Clear all variables in workspace 
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#sink output to file
sink(file = "WEL_3_output.txt", split = TRUE, append = FALSE)

#Read data
WEL_3 <- read.csv("WEL_3.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
summary(WEL_3)

#Transform Click Flag into binary variable
WEL_3$FLG_CLK <- factor(ifelse(WEL_3$FLG_CLK == 'Y', 1, 0))
WEL_3$status <- as.numeric(as.character(WEL_3$FLG_CLK))
num_clicks <- table(WEL_3$status)
num_clicks

#Barplot of clicks
barplot(num_clicks)

#Transforming to dates
WEL_3$DAT_1ER_CLK <- as.POSIXct(WEL_3$DAT_1ER_CLK, format = "%d%b%Y:%H:%M:%S") 
WEL_3$DAT_1ER_CLK <- format(WEL_3$DAT_1ER_CLK, "%d-%m-%y %H:%M")
WEL_3$DAT_EXC <- as.POSIXct(WEL_3$DAT_EXC, format = "%d%b%Y")
WEL_3$DAT_EXC <- format(WEL_3$DAT_EXC, format="%d-%m-%y %H:%M")

#Calculating time interval between sending and clicking
times<- (difftime(as.POSIXct(WEL_3$DAT_1ER_CLK, format = "%d-%m-%y %H:%M", tz=Sys.timezone()), as.POSIXct(WEL_3$DAT_EXC, "%d-%m-%y %H:%M", tz=Sys.timezone()), units = 'mins')-365)

#for training data: create times using first 2 days censoring
min_train=60*24*2
WEL_3$times_train <- replace_na(times, min_train)
WEL_3$status_train <- WEL_3$status*as.numeric(WEL_3$times_train <= min_train)

table(WEL_3$status_train)/length(WEL_3$status_train)

#for testing: create times using threshold of 50 days for censoring
min_test=60*24*50
WEL_3$times_test<-replace_na(times, min_test)
table(WEL_3$status)/length(WEL_3$status)

WEL_3$times_train <- as.numeric(WEL_3$times_train)
WEL_3$times_test <- as.numeric(WEL_3$times_test)

#Survival object
sobj <- Surv(WEL_3$times_train, WEL_3$status_train) 

# Look at 10 first elements
sobj[0:10]

# Look at summary
summary(sobj)

# Look at structure
str(sobj)

#Kaplan-Meier estimate for training sample
km_train <- survfit(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, data = WEL_3)
summary(km_train, times = seq(0,min_train,200))

plot(km_train, xlab="Time (minutes)", ylab="Survival Probability", conf.int = TRUE)

#Kaplan Meier for test sample, you can use this as a benchmark to test the other models
km_test <- survfit(Surv(WEL_3$times_test, WEL_3$status) ~ 1, data = WEL_3)
summary<-summary(km_test, times = seq(0,72000,1440), digits = 8)
capture.output(summary, file = "summaryKMtestWEL_3.txt")

plot(km_test, xlab="Time (minutes)", ylab="Survival Probability", xlim=c(0,65000), conf.int = TRUE)

# Exponential, Weibull, log-logistic and log-normal parametric model coefficients
# Opposite signs from Stata results, Weibull results differ; same as SAS

exponential <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, dist="exponential")
summary(exponential)
AIC(exponential)
BIC(exponential)

weibull <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, dist="weibull")
summary(weibull)
AIC(weibull)
BIC(weibull)

loglogistic <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, dist="loglogistic")
summary(loglogistic)
AIC(loglogistic)
BIC(loglogistic)

lognormal<- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, dist="lognormal")
summary(lognormal)
AIC(lognormal)
BIC(lognormal)

# Exponential, Weibull, log-logistic and log-normal parametric model coefficients with covariates
# Opposite signs from Stata results, Weibull results differ; same as SAS

exponential_cov <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ WEL_3$FB_TIER + WEL_3$NBR_OUV, dist="exponential")
summary(exponential_cov)
AIC(exponential_cov)
BIC(exponential_cov)

weibull_cov <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ WEL_3$FB_TIER + WEL_3$NBR_OUV, dist="weibull")
summary(weibull_cov)
AIC(weibull_cov)
BIC(weibull_cov)

loglogistic_cov <- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~ WEL_3$FB_TIER + WEL_3$NBR_OUV, dist="loglogistic")
summary(loglogistic_cov)
AIC(loglogistic_cov)
BIC(loglogistic_cov)

lognormal_cov<- survreg(Surv(WEL_3$times_train, WEL_3$status_train) ~WEL_3$FB_TIER + WEL_3$NBR_OUV, dist="lognormal")
summary(lognormal_cov)
AIC(lognormal_cov)
BIC(lognormal_cov)

# Split hazard models: some fraction of the population will never click
split_exponential<- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="exp", mixture=T)
print(split_exponential)
BIC(split_exponential)

split_weibull <- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="weibullPH", mixture=T)
print(split_weibull)
BIC(split_weibull)

split_loglogistic<- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="llogis", mixture=T)
print(split_loglogistic)
BIC(split_loglogistic)

split_lognormal <- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="lnorm", mixture=T)
print(split_lognormal)
BIC(split_lognormal)

# Split hazard models: some fraction of the population will never click
split_exponential_cov<- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="exp", mixture=T, log=T, anc = list(shape = ~ WEL_3$FB_TIER + WEL_3$NBR_OUV))
print(split_exponential_cov)
BIC(split_exponential_cov)

split_weibull_cov <- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="weibullPH", mixture=T, log=T, anc = list(shape = ~ WEL_3$FB_TIER + WEL_3$NBR_OUV))
print(split_weibull_cov)
BIC(split_weibull_cov)

split_loglogistic_cov<- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="llogis", mixture=T, log=T, anc = list(shape = ~ WEL_3$FB_TIER + WEL_3$NBR_OUV))
print(split_loglogistic_cov)
BIC(split_loglogistic_cov)

split_lognormal_cov <- flexsurvcure(sobj~1, data=WEL_3, link="logistic", dist="lnorm", mixture=T, log=T, anc = list(shape = ~ WEL_3$FB_TIER + WEL_3$NBR_OUV))
print(split_lognormal_cov)
BIC(split_lognormal_cov)

# theta is the estimate of people who will never click

#Predict
pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(weibull, newdata=data.frame(ed=1), type='quantile',
                 p=pct, se=TRUE)
matplot(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
              ptime$fit - 2*ptime$se.fit), 1-pct,
        xlab="minutes", ylab="Survival", type='l', lty=c(1,2,2), col=1)
dev.off()
plot(km_test, xlab="Time (minutes)", ylab="Survival Probability", conf.int = TRUE)
lines(predict(exponential, newdata=data.frame(ed=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red")
lines(predict(weibull, newdata=data.frame(ed=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="blue")
lines(predict(loglogistic, newdata=data.frame(ed=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="green")
lines(predict(lognormal, newdata=data.frame(ed=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="orange")
legend("bottomleft", lty=c(1,1,1), lwd=c(2,2,2), cex=0.5,
       col=c("red", "blue", "green", "orange"),
       c("Exponential","Weibull","Loglogistic", "Lognormal"))

#Plot Best Fitting Models
plot(km_test, xlab="Time (minutes)", ylab="Survival Probability", main="Best Fitting models for WEL_3", conf.int = FALSE, ylim=c(0.4,1.0),xlim=c(0,65000), width=880, height=710, res=300, family="serif")
lines(split_lognormal_cov, col="green")
lines(split_loglogistic, col="blue")
lines(predict(lognormal, newdata=data.frame(ed=1),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col="red")
op <- par(family = "serif")
legend("bottomleft", lty=c(1,1,1), lwd=c(2,2,2), cex=0.8,
       col=c("black", "green", "blue", "red"),
       c("KM Test","Split-Log-Normal-Cov","Split-Log-Logistic","Log-Normal"))
par(op, cex=1.0)

sink()
