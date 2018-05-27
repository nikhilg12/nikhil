#Nikhil Guruji

#Goal is to calculate Exponentially Weighed Historical VaR and the usual historical VaR and compare the 
# exceptions that occur in the two methods, then calculate the VaR through bootstrapping method

library(readxl)
library(xts)
library(data.table)
library(moments)
library(tseries)

gainsdata=read_excel("homework3_data.xls")
gainsdata=xts(gainsdata$gain,order.by = gainsdata$date)

gainsdata0708=gainsdata["2007::2008"]
gainsdata06=gainsdata["2006"]

gainsdata=as.data.table(gainsdata)

#function to calculate the Exponentially Weighted VaR
ewvar=function(retdata){
  n=length(retdata)
  
  #calculating the weights
  weights=seq(1:n)
  weights=sapply(weights,function(x) lambda^(n-x))
  weights=weights/sum(weights)
  
  retdata=as.data.table(cbind(retdata,weights))
  setnames(retdata,c("Returns","Weights"))
  retdata=retdata[order(Returns)] #order by the returns from lowest to highest
  retdata$cumweights=cumsum(retdata$Weights) #take the cumulative weights
  retdata[which(cumweights>0.01)[1],Returns] #return whose cum probability is around 0.01
}

histvar0708=c()
ewvar0708=c()
ci_hist_par=boot_ci_ew=boot_ci_hist=matrix(nrow=length(gainsdata0708),ncol=2)
lambda=0.995
pb <- txtProgressBar(min = 1, max = length(gainsdata0708), style = 3)
for(i in 1:length(gainsdata0708)){
  #will be dealing with different set of data in each iteration
  currentdata=c(as.vector(gainsdata06),as.vector(gainsdata0708[1:i]))
  
  #historical VaR
  worst=ceiling((length(gainsdata06)+i)*0.01) #calculate the index which will be the worst 1% return
  histvar0708[i]=sort(currentdata)[worst] # the worst 1% return
  n=length(currentdata)
  
  #Exponentially weighted VaR
  ewvar0708[i]=ewvar(currentdata)
  
  #Parametric confidence interval (assuming normal distribution)
  mu=mean(currentdata)
  sig=sd(currentdata)
  x=mu+sig*qnorm(0.01) #Because 99% 1 day VaR
  fx=dnorm(x,mean=mu,sd=sig)
  sdvar=(1/fx)*sqrt(0.99*0.01/n)
  ci_hist_par[i,]=c(histvar0708[i]+sdvar*qnorm(0.95),histvar0708[i]-sdvar*qnorm(0.95))
  
  #bootstrapping confidence interval for Exponentially Weighted
  ewsamples=matrix(sample(currentdata,size=n*1000,replace=TRUE),ncol=1000)
  ewvars=apply(ewsamples,2,ewvar)
  boot_ci_ew[i,]=c(sort(ewvars)[25],sort(ewvars)[975])
  
  #bootstrapping confidence interval for historical VaR
  histsamples=matrix(sample(currentdata,size=n*1000,replace=TRUE),ncol=1000)
  histvars=apply(histsamples,2,function(x) sort(x)[worst])
  boot_ci_hist[i,]=c(sort(histvars)[25],sort(histvars)[975])
  setTxtProgressBar(pb,i)
}
close(pb)

gainsdata0708=cbind(gainsdata0708,histvar0708,ewvar0708) #doing this to make sure index of histvar vector is in dates
plot(gainsdata0708[,1],type='l',main="Historical VaR")
lines(gainsdata0708[,2],col='red')

plot(gainsdata0708[,1],type='l',main="EW VaR")
lines(gainsdata0708[,3],col='red')

exceptionshist=length(which(gainsdata0708[,2]>gainsdata0708[,1]))
exceptionsew=length(which(gainsdata0708[,3]>gainsdata0708[,1]))
exceptionshist
exceptionsew

hist(gainsdata$V1,breaks=100)

#taking volatility & mean over past 21 trading days
gainsdata$Monthlyvol=rollapply(gainsdata$V1,21,sd,fill=NA,align="right")
gainsdata$Monthlymean=rollapply(gainsdata$V1,21,mean,fill=NA,align="right")
hist((gainsdata$V1-gainsdata$Monthlymean)/gainsdata$Monthlyvol,breaks = 100)

#Assessing the difference in the confidence intervals

plot(ci_hist_par[,2]-ci_hist_par[,1],type='l')
plot(boot_ci_ew[,2]-boot_ci_ew[,1],type='l')
plot(boot_ci_hist[,2]-boot_ci_hist[,1],type='l')

#Jarque Bera Test for normality

jb = function(test){
  if(test>qchisq(0.95,df=2)){
    sprintf("Null is rejected and distribution is not normal. statistic = %f",test)
  }
  else sprintf("Null cannot be rejected and the distribution is normal. Statistic = %f",test)
}
skew_test_norm = skewness((gainsdata$V1-gainsdata$Monthlymean)/
                            gainsdata$Monthlyvol,na.rm = TRUE)/
  sqrt(6/length(which(!is.na(gainsdata$Monthlyvol))))
kurt_test_norm = (kurtosis((gainsdata$V1-gainsdata$Monthlymean)/
                             gainsdata$Monthlyvol,na.rm = TRUE)-3)/
  sqrt(24/length(which(!is.na(gainsdata$Monthlyvol))))


jb_norm= skew_test_norm^2 + kurt_test_norm^2

jb(jb_norm)

skew_test_gains = skewness(gainsdata$V1)/sqrt(6/nrow(gainsdata))
kurt_test_gains = (kurtosis(gainsdata$V1)-3)/sqrt(24/nrow(gainsdata))


jb_gains= skew_test_gains^2 + kurt_test_gains^2

jb(jb_gains)

jarque.bera.test(gainsdata$V1)
jarque.bera.test(which(!is.na((gainsdata$V1-gainsdata$Monthlymean)/
                   gainsdata$Monthlyvol)))
