library(data.table)
library(readxl)
library(dplyr)
library(xts)
library(zoo)

#1.1 Estimation of the tail from power law
gainsdata=read_excel("D:/MSFE/MFE/Q3/Financial Risk Measurement and Management/Week 3/homework3_data.xls")

#turn it into xts format to separate data by data
gainsdata=xts(gainsdata$gain,order.by = gainsdata$date)
gainsdata0708=gainsdata["2007::2008"]
gainsdata06=gainsdata["2006"]

plvar99=c()
plvar999=c()
for(i in 1:length(gainsdata0708)){
  #will be dealing with different set of data in each iteration
  currentdata=sort(c(as.vector(gainsdata06),as.vector(gainsdata0708[1:i])))
  
  #worst 5%
  worst5=ceiling((length(gainsdata06)+i)*0.05)
  u=currentdata[worst5]
  
  #data that is in the tail (in the worst 5%)
  v=currentdata[currentdata<=u]
  
  #likelihood function
  mlepl=function(x){sum(-(log(1/x[2]) - (1/x[1] + 1)*log(1 + x[1]*(v-u)/x[2])))}
  
  #optim function can minimize, not maximize, so we minimize the negative of log likelihood
  ksibeta=optim(c(0.5,1000),mlepl)$par
  ksi=ksibeta[1]
  beta=ksibeta[2]
  
  #var at different confidence levels
  plvar99[i]=u+(beta/ksi)*((0.99/0.95)^(-ksi)-1)
  plvar999[i]=u+(beta/ksi)*((0.999/0.95)^(-ksi)-1)
}
plot(plvar99,type='l',ylab="",xlab="No. of days after December 2006",main="VaR from Power Law")
lines(plvar999,col='blue')
legend("topright",legend=c("99%","99.9%"),col=c("black","blue"),lty=c(1,1))

#1.3 Calculate VaR through ARCH process for volatility

archvol=sqrt(rollapply(gainsdata0708^2,width=30,FUN=mean,fill=NA,align="right"))
archvar99=-qnorm(0.99)*archvol

plot(archvar99,type='l',main="99% VaR using ARCH(30) model")

#1.4 EWMA process for volatility
lambda=0.995
ewmavol=rep(NA,length(gainsdata0708))
ewmavol[1]=gainsdata0708[1]/100 #starting value
for(i in 2:length(gainsdata0708)){
  ewmavol[i]=sqrt(lambda*ewmavol[i-1]^2 + (1-lambda)*(gainsdata0708[i-1]/100)^2)
}

ewmavar99=-qnorm(0.99)*ewmavol
plot(ewmavar99,type='l',ylab="",xlab="No. of days after December 2006")

#1.5 GARCH Process for volatility estimated by MLE
longrunvol=sqrt(mean((gainsdata0708/100-mean(gainsdata0708/100))^2))

mlegarch=function(x){
  garchvol=c()
  garchvol[1]=sqrt(longrunvol/(1-x[2]-x[3]))
  for(i in 2:length(gainsdata0708)){
    garchvol[i]=sqrt(x[1] + x[2]*gainsdata0708[i-1]^2 + x[3]*garchvol[i-1]^2)
  }
  sum(log(garchvol^2) + (gainsdata0708/garchvol)^2)
}
omalbet=optim(c(0.2*longrunvol,0.4,0.4),mlegarch)$par
omega=omalbet[1]
alpha=omalbet[2]
betagarch=omalbet[3]

#calculating volatilites based on the parameters obtained above
garchvol=c()
garchvol[1]=sqrt(longrunvol/(1-alpha-betagarch))
for(i in 2:length(gainsdata0708)){
  garchvol[i]=sqrt(omega + alpha*gainsdata0708[i-1]^2 + betagarch*garchvol[i-1]^2)
}

garchvar99=-qnorm(0.99)*garchvol
plot(garchvar99,type='l',ylab="",xlab="No. of days after December 2006")


#1.6 Comparison of all the processes
vars=cbind(archvar99,ewmavar99*100,garchvar99,plvar99)

plot(vars[,1],type='l',grid.col = "white",ylim=c(-500,0),main="99% 1-day VaR for different models")
lines(vars[,2],type='l',col='red')
lines(vars[,3],type='l',col='green')
lines(vars[,4],type='l',col='blue')
legend("bottomleft",legend=c("ARCH(30)","EWMA","GARCH","Power Law"),
       col=c("black","red","green","blue"),lty=c(1,1,1,1),bty='n')
