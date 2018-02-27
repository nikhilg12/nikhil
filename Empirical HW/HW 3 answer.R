#HW 3 Cohort 1: Bingjie Hu, Nikhil Guruji, Jialiang Le


library(data.table)
library(zoo)
library(xts)
library(lubridate)
#Process Data
ffactor <- fread("https://raw.githubusercontent.com/ErinHu95/HW/master/F-F_Research_Data_Factors.CSV")
indRtn <- fread("https://raw.githubusercontent.com/ErinHu95/HW/master/48_Industry_Portfolios.CSV")
ffactor$V1 <- paste(ffactor$V1,"01",sep = "")
indRtn$V1 <- paste(indRtn$V1,"01",sep = "")
ffactor$V1 <- ymd(ffactor$V1)
indRtn$V1 <- ymd(indRtn$V1)
class(ffactor$V1)
ffactor <- xts(ffactor[,2:ncol(ffactor)],order.by = ffactor$V1)
indRtn <- xts(indRtn[,2:ncol(indRtn)],order.by = indRtn$V1)
ffactor <- ffactor["1960-01-01::2016-01-01"]
indRtn <- indRtn["1960-01-01::2016-01-01"]
indRtn[indRtn==-99.99] <- NA
na_flag <- apply(is.na(indRtn),2,sum)
indRtn <- indRtn[,which(na_flag == 0)]
ind_excess_rtn <- data.frame()
for(i in 1:ncol(indRtn)){
  ind_excess_rtn <- cbind(ind_excess_rtn,indRtn[,i]-ffactor[,"RF"])
}
#Compute eigenvalues
varMatrix <- cov(ind_excess_rtn)
ev <- eigen(varMatrix)
evalue <- ev$values
frac <- evalue/sum(evalue)
barplot(frac,col = "lightblue")

#2(a)
f3frac <- sum(frac[1:3])
#2(b)
pca <- prcomp(ind_excess_rtn,scale = T)
lds1 <- t(pca$rotation)
a <- t(ind_excess_rtn)
f3 <- lds1%*%a
f3mean <- apply(f3[1:3,],1,mean)
f3std <- apply(f3[1:3,],1,sd)
f3cor <- cor(cbind(f3[1,],f3[2,],f3[3,]))
#2(c)
pca3 <- prcomp(ind_excess_rtn,scale = T,rank. = 3)
lds <- pca3$rotation
pred_rtn <- (lds)%*%f3mean
rlz_rtn <- apply(ind_excess_rtn,2,mean)
plot(x = t(pred_rtn),y = rlz_rtn,xlim = c(0,1),ylim = c(0,1),xlab = "predicted return",ylab = "realized return",col = "blue")
abline(a = 0,b = 1,col = "red")
#2(d)
csRsqr <- 1-(var(rlz_rtn-pred_rtn)/var(rlz_rtn))

#3
portRtn <- fread("https://raw.githubusercontent.com/ErinHu95/HW/master/25_Portfolios_5x5.CSV")
portRtn$V1 <- paste(portRtn$V1,"01",sep = "")
portRtn$V1 <- ymd(portRtn$V1)
portRtn <- xts(portRtn[,2:ncol(portRtn)],order.by = portRtn$V1)
portRtn <- portRtn["1960-01-01::2016-01-01"]
port_excess_rtn <- data.frame()
for(i in 1:ncol(portRtn)){
  port_excess_rtn <- cbind(port_excess_rtn,portRtn[,i]-ffactor[,"RF"])
}

#(a)
#Compute eigenvalues
varMatrix3a <- cov(port_excess_rtn)
ev3a <- eigen(varMatrix3a)
evalue3a <- ev3a$values
frac3a <- evalue3a/sum(evalue3a)
barplot(frac3a,col = "lightblue")

#(b)
pct <- sum(frac3a[1:5])
pct

#1
#(a)
#Construct a portfolio B using the portfolios that 'mimick' the factors:
#invest 0.5 in factor 1 mimick portfolio
#invest 0.75 in factor 2 mimick portfolio
#invest 1-0.5-0.75 in risk-free rate
#As a result, we copy a portfolio with the same factors as A. We could arbitrage 
#by shorting B and longing A to get risk-free alpha of A.
#The profit will be 1%

#(b)
#-1%
exp_rtn <- 0.5*0.06+0.75*(-0.02)
