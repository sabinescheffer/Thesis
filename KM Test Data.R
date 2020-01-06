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
sink(file = "KMtest_output.txt", split = TRUE, append = FALSE)

#Read data
PLC <- read.csv("PLC.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
summary(PLC)
AMEX_NLA <- read.csv("AMEX_NLA.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
summary(AMEX_NLA)
WEL_3 <- read.csv("WEL_3.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)
summary(WEL_3)


#Transform Click Flag into binary variable
PLC$FLG_CLK <- factor(ifelse(PLC$FLG_CLK == 'Y', 1, 0))
PLC$status <- as.numeric(as.character(PLC$FLG_CLK))
PLC_num_clicks <- table(PLC$status)
PLC_num_clicks

AMEX_NLA$FLG_CLK <- factor(ifelse(AMEX_NLA$FLG_CLK == 'Y', 1, 0))
AMEX_NLA$status <- as.numeric(as.character(AMEX_NLA$FLG_CLK))
AMEX_NLA_num_clicks <- table(AMEX_NLA$status)
AMEX_NLA_num_clicks

WEL_3$FLG_CLK <- factor(ifelse(WEL_3$FLG_CLK == 'Y', 1, 0))
WEL_3$status <- as.numeric(as.character(WEL_3$FLG_CLK))
WEL_3_num_clicks <- table(WEL_3$status)
WEL_3_num_clicks

#Barplot of clicks
barplot(PLC_num_clicks)
barplot(AMEX_NLA_num_clicks)
barplot(WEL_3_num_clicks)

#Transforming to dates
PLC$DAT_1ER_CLK <- as.POSIXct(PLC$DAT_1ER_CLK, format = "%d%b%Y:%H:%M:%S") 
PLC$DAT_1ER_CLK <- format(PLC$DAT_1ER_CLK, "%d-%m-%y %H:%M")
PLC$DAT_ENV_CMU <- as.POSIXct(PLC$DAT_ENV_CMU, format = "%d%b%Y")
PLC$DAT_ENV_CMU <- format(PLC$DAT_ENV_CMU, "%d-%m-%y %H:%M")

AMEX_NLA$DAT_1ER_CLK <- as.POSIXct(AMEX_NLA$DAT_1ER_CLK, format = "%d%b%Y:%H:%M:%S") 
AMEX_NLA$DAT_1ER_CLK <- format(AMEX_NLA$DAT_1ER_CLK, "%d-%m-%y %H:%M")
AMEX_NLA$DAT_ENV_CMU <- as.POSIXct(AMEX_NLA$DAT_ENV_CMU, format = "%d%b%Y")
AMEX_NLA$DAT_ENV_CMU <- format(AMEX_NLA$DAT_ENV_CMU, "%d-%m-%y %H:%M")

WEL_3$DAT_1ER_CLK <- as.POSIXct(WEL_3$DAT_1ER_CLK, format = "%d%b%Y:%H:%M:%S") 
WEL_3$DAT_1ER_CLK <- format(WEL_3$DAT_1ER_CLK, "%d-%m-%y %H:%M")
WEL_3$DAT_EXC <- as.POSIXct(WEL_3$DAT_EXC, format = "%d%b%Y")
WEL_3$DAT_EXC <- format(WEL_3$DAT_EXC, format="%d-%m-%y %H:%M")


#Calculating time interval between sending and clicking
PLC_times<- (difftime(as.POSIXct(PLC$DAT_1ER_CLK, format = "%d-%m-%y %H:%M", tz=Sys.timezone()), as.POSIXct(PLC$DAT_ENV_CMU, "%d-%m-%y %H:%M", tz=Sys.timezone()), units = 'mins'))-582
AMEX_NLA_times<- (difftime(as.POSIXct(AMEX_NLA$DAT_1ER_CLK, format = "%d-%m-%y %H:%M", tz=Sys.timezone()), as.POSIXct(AMEX_NLA$DAT_ENV_CMU, "%d-%m-%y %H:%M", tz=Sys.timezone()), units = 'mins'))-362
WEL_3_times<- (difftime(as.POSIXct(WEL_3$DAT_1ER_CLK, format = "%d-%m-%y %H:%M", tz=Sys.timezone()), as.POSIXct(WEL_3$DAT_EXC, "%d-%m-%y %H:%M", tz=Sys.timezone()), units = 'mins')-365)

#for training data: create times using first 2 days censoring
min_train=60*24*2

PLC$times_train <- replace_na(PLC_times, min_train)
PLC$status_train <- PLC$status*as.numeric(PLC$times_train <= min_train)
table(PLC$status_train)/length(PLC$status_train)

AMEX_NLA$times_train <- replace_na(AMEX_NLA_times, min_train)
AMEX_NLA$status_train <- AMEX_NLA$status*as.numeric(AMEX_NLA$times_train <= min_train)
table(AMEX_NLA$status_train)/length(AMEX_NLA$status_train)

WEL_3$times_train <- replace_na(WEL_3_times, min_train)
WEL_3$status_train <- WEL_3$status*as.numeric(WEL_3$times_train <= min_train)
table(WEL_3$status_train)/length(WEL_3$status_train)

#Survival object
PLC_sobj<- Surv(PLC$times_train, PLC$status_train) 
AMEX_NLA_sobj <- Surv(AMEX_NLA$times_train, AMEX_NLA$status_train) 
WEL_3_sobj <- Surv(WEL_3$times_train, WEL_3$status_train) 

# Look at 10 first elements
PLC_sobj[0:10]
AMEX_NLA_sobj[0:10]
WEL_3_sobj[0:10]

# Look at summary
summary(PLC_sobj)
summary(AMEX_NLA_sobj)
summary(WEL_3_sobj)

#Kaplan-Meier estimate for training sample
km_PLC <- survfit(Surv(PLC$times_train, PLC$status_train) ~ 1, data = PLC)
summary(km_PLC, times = seq(0,72000,1440))

km_AMEX_NLA <- survfit(Surv(AMEX_NLA$times_train, AMEX_NLA$status_train) ~ 1, data = AMEX_NLA)
summary(km_AMEX_NLA, times = seq(0,72000,1440))

km_WEL_3 <- survfit(Surv(WEL_3$times_train, WEL_3$status_train) ~ 1, data = WEL_3)
summary(km_WEL_3, times = seq(0,72000,14400))

#PLot KM for all campaigns on training data
plot(km_PLC, xlab="Time (minutes)", ylab="Survival Probability", main="KM training estimate for FB Campaigns", conf.int = FALSE, ylim=c(0.8,1.0), width=880, height=710, res=300, family="serif")
lines(km_AMEX_NLA, col ='red', conf.int=FALSE)
lines(km_WEL_3, col = 'green', conf.int=FALSE)
op <- par(family = "serif")
legend("bottomleft", lty=c(1,1,1), lwd=c(2,2,2), cex=0.5,
       col=c("black", "red", "green"),
       c("PLC","AMEX_NLA","WEL_3"))
par(op, cex=1.0)

#for testing: create times using threshold of 50 days for censoring
min_test=60*24*50

PLC$times_test<-replace_na(PLC_times, min_test)
table(PLC$status)/length(PLC$status)
PLC$times_train <- as.numeric(PLC$times_train)
PLC$times_test <- as.numeric(PLC$times_test)

AMEX_NLA$times_test<-replace_na(AMEX_NLA_times, min_test)
table(AMEX_NLA$status)/length(AMEX_NLA$status)
AMEX_NLA$times_train <- as.numeric(AMEX_NLA$times_train)
AMEX_NLA$times_test <- as.numeric(AMEX_NLA$times_test)

WEL_3$times_test<-replace_na(WEL_3_times, min_test)
table(WEL_3$status)/length(WEL_3$status)
WEL_3$times_train <- as.numeric(WEL_3$times_train)
WEL_3$times_test <- as.numeric(WEL_3$times_test)

#Kaplan Meier for test sample, you can use this as a benchmark to test the other models
km_test_PLC <- survfit(Surv(PLC$times_test, PLC$status) ~ 1, data = PLC)
summary(km_test_PLC, times = seq(0,72000,1440))

km_test_AMEX <- survfit(Surv(AMEX_NLA$times_test, AMEX_NLA$status) ~ 1, data = AMEX_NLA)
summary(km_test_AMEX, times = seq(0,72000,1440))

km_test_WEL3 <- survfit(Surv(WEL_3$times_test, WEL_3$status) ~ 1, data = WEL_3)
summary(km_test_WEL3, times = seq(0,72000,1440))

  #PLot KM for all campaigns on test data
  plot(km_test_PLC, xlab="Time (minutes)", ylab="Survival Probability", main="KM test estimate for FB Campaigns test", conf.int = FALSE, ylim=c(0.8,1.0), xlim=c(0,65000), width=880, height=710, res=300, family="serif")
  lines(km_test_AMEX, col ='red', conf.int=FALSE)
  lines(km_test_WEL3, col = 'green', conf.int=FALSE)
  op <- par(family = "serif")
  legend("bottomleft", lty=c(1,1,1), lwd=c(2,2,2), cex=0.5,
         col=c("black", "red", "green"),
         c("PLC","AMEX_NLA","WEL_3"))
  par(op, cex=1.0)
  

