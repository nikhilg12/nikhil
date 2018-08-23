#DoubleLine ForwardRate 
library(readxl)
library(data.table)
library(lubridate)
library(zoo)
library(readr)
library(R.matlab)
library(matlib)
library(MASS)
library(optimx)
library(matrixcalc)
#Take Maturity & Fed as inputs.
Maturity<-c(3,6,12,24,60,84,120)
Fed <-fread("feds200628.xls.csv",skip=9)


#Compute Monthly Forward Rate
MonthlyForwardRate<-function(Fed,Maturity){
  Fed<-data.table(Fed)
  setnames(Fed,'V1','Date')
  Fed[,Year:=year(Date)][,Month:=month(Date)][,Day:=day(Date)]
  Fed<-Fed[Year>=1990]
  setorder(Fed,Date)
  Fed[,MaxDay:=max(Day),by=c("Year","Month")]  #Select the last line of each month.
  GSWdata<-Fed[which(Day==MaxDay),c("Year","Month","Day","BETA0","BETA1","BETA2","BETA3","TAU1","TAU2")]
  GSWdata[!(is.na(BETA0)|is.na(BETA1)|is.na(BETA2)|is.na(BETA3)|is.na(TAU1)|is.na(TAU2)), ] #Remove NA.
  
  if(GSWdata[nrow(GSWdata),3]<31) GSWdata[nrow(GSWdata),]=NULL  #Check if the last observation is the end of the month.
  GSWdata[,YearMon:=Year*100+Month]
  
  #Calculate the annualized forward rate.
  #NSS Function
  nelson=function(t,beta0=GSWdata$BETA0,beta1=GSWdata$BETA1,beta2=GSWdata$BETA2,beta3=GSWdata$BETA3,tau1=GSWdata$TAU1,tau2=GSWdata$TAU2){
    beta0+beta1*(1-exp(-t/tau1))/(t/tau1)+beta2*((1-exp(-t/tau1))/(t/tau1)-exp(-t/tau1))+beta3*((1-exp(-t/tau2))/(t/tau2)-exp(-t/tau2))
  }
  
  n=Maturity/12
  
  Yield1<-sapply(n,nelson)        #Apply the NNS model to get yield for each maturity.
  Yield1<-n*t(Yield1)
  
  n2=(Maturity+1)/12
  
  
  Yield2<-sapply(n2,nelson)
  Yield2<-n2*t(Yield2)
  
  ForwardRates=(Yield2-Yield1)*12
  ForwardRates<-t(ForwardRates)
  return(ForwardRates)
}


forwardrates<-MonthlyForwardRate(Fed,Maturity)
forwardrates<-t(forwardrates)

#Import parameters Wu and Xia used.
parameters<-unlist(readMat("parameters.mat"))  


#Kalman Filter (non-extended for GATSM) 
KF_GATSM<-function(parameters)
{
  
  #Parameterization
  T=ncol(forwardrates)
  J=Maturity
  rhoP = parameters[1:9]; 
  rhoP = matrix(rhoP,ncol=3) 
  muP = parameters[10:12] 
  rhoQ1 = parameters[13]
  rhoQ2 = parameters[14]
  sigma = matrix(c(abs(parameters[15]), 0, 0 , parameters[16], abs(parameters[18]) ,0, parameters[17], parameters[19], abs(parameters[20])),nrow=3,byrow=TRUE) 
  omega = sigma%*%t(sigma)
  delta0 = parameters[21];
  omegaM<-matrix(0,nrow=nrow(forwardrates),ncol=nrow(forwardrates))
  diag(omegaM)=parameters[22]^2
  I3=matrix(0,ncol=3,nrow=3) 
  diag(I3)=1
  
  
  sigma_test_mat<-matrix(c(abs(parameters[15]),  parameters[16],  parameters[17] , parameters[16], abs(parameters[18]) , parameters[19], parameters[17], parameters[19], abs(parameters[20])),nrow=3,byrow=TRUE) 
  #Restrictions on parameters.(return big values if restrictions not met)
  if(abs(parameters[1])>1 |abs(parameters[5])>1 |abs(parameters[9])>1) return(10^8) #Rho_P
  if(abs(parameters[13])>1 |abs(parameters[14])>1) return(10^8) #Rho_Q
  if(parameters[22]<0) return(10^8) #OmegaM
  
  #Initialization
  X1<-matrix(NA,nrow=3,ncol=T+1)
  X2<-matrix(NA,nrow=3,ncol=T+1)  #Here we assume X1 and X2 are one the SAME time horizon.
  X1[,1]=0  #X1 starts at 0.
  
  
  #Compute aJ and bJ
  JJ=seq(1,max(J))
  bn=cbind(rhoQ1^JJ,rhoQ2^JJ,JJ*rhoQ2^(JJ-1)) #rho^J in Appendix A.
  Bn=rbind(c(1,1,0),bn)       
  Bn<-apply(Bn,2,cumsum)      #Summation of rho^J
  Bn<-Bn[-nrow(Bn),]          #Summation up to (n-1)
  delta1=c(1,1,0)
  aJ=c()
  
  for(i in 1:length(J))
  {
    aJ[i]<-delta0-0.5*t(Bn[J[i],])%*%((sigma%*%t(sigma))%*%Bn[J[i],])/1200
  }
  
  
  bJ=cbind(rhoQ1^J,rhoQ2^J,J*rhoQ2^(J-1)) 
  
  V1<-array(NA,dim=c(3,3,T+1))  #V1 is a 3-D arry of covariance at each time
  V2<-array(NA,dim=c(3,3,T+1))  #So is V2.
  V1[,,1]=0
  diag(V1[,,1])=100
  V1[3,3,1]<-V1[3,3,1]/144
  loglikvec<-c()
  for(i in 1:T)
  {
    F_h=aJ+bJ%*%X1[,i]  #Prediction of forward rate based on Factors. 
    err=forwardrates[,i]-F_h #Error in prediction.
    
    H=bJ 
    S=H%*%V1[,,i]%*%t(H)+omegaM
    
    if (rcond(S) < 1e-8 | is.na(rcond(S)))   return(10^8)
    
    flag=det(S)
    if(!is.finite(flag) | (flag<=0)) return(10^9)
    
    
    InvS=solve(S)
    
    
    K=V1[,,i]%*%t(H)%*%InvS #Kalman Gain using guessed value of V1
    
    loglikvec[i] = length(J)*log(2*pi) + log(det(S)) + t(err)%*%InvS%*%err 
    loglikvec[i] = -1/2*loglikvec[i]
    
    X2[,i+1]=X1[,i]+K%*%err #Correct the guess with Kalman Gain.
    X1[,i+1]=muP+rhoP%*%X2[,i+1] #Update next value of state variable(X1).
    I3=matrix(0,ncol=3,nrow=3) 
    diag(I3)=1
    
    V2[,,i+1]=(I3-K%*%H)%*%V1[,,i] #Update the covariance matrix.
    V1[,,i+1]=rhoP%*%V2[,,i+1]%*%t(rhoP)+sigma%*%t(sigma)  #Predict the covariance for next step.
    
  }
  llf=sum(loglikvec)
  return(-llf)
}

KF_GATSM(parameters)

#Optimization
set.seed(20)
startvalue<-rnorm(22)
startvalue[21]=12
test1<-optim(startvalue,KF_GATSM,method="BFGS") 
parameters_test<-c(test1$par)
parameters_test

#Import parameters Wu and Xia used.
parameters_2<-unlist(readMat("parameters_rlb.mat"))  

#Extended Kalman Filter
EKF_SRTSM<-function(parameters)
{
  #Parameterization
  T=ncol(forwardrates)
  J=Maturity
  rhoP = parameters[1:9]; 
  rhoP = matrix(rhoP,ncol=3) 
  muP = parameters[10:12] 
  rhoQ1 = parameters[13]
  rhoQ2 = parameters[14]
  sigma = matrix(c(abs(parameters[15]), 0, 0 , parameters[16], abs(parameters[18]) ,0, parameters[17], parameters[19], abs(parameters[20])),nrow=3,byrow=TRUE); 
  omega = sigma%*%t(sigma)
  delta0 = parameters[21];
  omegaM<-matrix(0,nrow=nrow(forwardrates),ncol=nrow(forwardrates))
  diag(omegaM)=parameters[22]^2
  rlb=0.25
  
  #Restrictions on parameters.(return big values if restrictions not met)
  if(abs(parameters[1])>1 |abs(parameters[5])>1 |abs(parameters[9])>1) return(10^8) #Rho_P
  if(abs(parameters[13])>1 |abs(parameters[14])>1) return(10^8) #Rho_Q
  if(parameters[22]<0) return(10^8) #OmegaM
  
  #Initialization
  X1<-matrix(NA,nrow=3,ncol=T+1)
  X2<-matrix(NA,nrow=3,ncol=T+1)  #Here we assume X1 and X2 are one the SAME time horizon.
  X1[,1]=0  #X1 starts at 0.
  
  
  #Compute aJ and bJ
  JJ=seq(1,max(J))
  bn=cbind(rhoQ1^JJ,rhoQ2^JJ,JJ*rhoQ2^(JJ-1)) #rho^J in Appendix A.
  Bn=rbind(c(1,1,0),bn)       
  Bn<-apply(Bn,2,cumsum)      #Summation of rho^J
  Bn<-Bn[-nrow(Bn),]          #Summation up to (n-1)
  delta1=c(1,1,0)
  aJ=c()
  SigmaJ=c()
  for(i in 1:length(J))
  {
    aJ[i]<-delta0-0.5*t(Bn[J[i],])%*%((sigma%*%t(sigma))%*%Bn[J[i],])/1200
  }
  
  #Compute SigmaJ.
  cn=rbind(c(1,1,0),bn)
  cn=cn[-nrow(cn),]
  Sigma_step=c()
  SigmaJ=c()
  
  for(i in 1:length(J))
  {
    for(j in 1:Maturity[i])
    {
      Sigma_step[j]<-t(cn[j,])%*%omega%*%cn[j,]
    }
    SigmaJ[i]=sqrt(sum(Sigma_step))
  }
  
  bJ=cbind(rhoQ1^J,rhoQ2^J,J*rhoQ2^(J-1)) 
  V1<-array(NA,dim=c(3,3,T+1))  #V1 is a 3-D arry of covariance at each time
  V2<-array(NA,dim=c(3,3,T+1))  #So is V2.
  V1[,,1]=0
  diag(V1[,,1])=100
  V1[3,3,1]<-V1[3,3,1]/144
  loglikvec<-c()
  
  for(i in 1:T)
  {
    F_h=aJ+bJ%*%X1[,i]  #Prediction of forward rate based on Factors. 
    Z1_temp=(F_h-rlb)/SigmaJ
    Z2_temp=rlb+(F_h-rlb)*pnorm(Z1_temp)+SigmaJ*dnorm(Z1_temp)
    
    err=forwardrates[,i]-Z2_temp #Error in prediction.
    
    H=rep(pnorm(Z1_temp))*bJ
    
    
    S=H%*%V1[,,i]%*%t(H)+omegaM
    
    if (rcond(S) < 1e-8 | is.na(rcond(S)))   return(10^8)
    
    flag=det(S)
    if(!is.finite(flag) | (flag<=0)) return(10^9)
    
    InvS=solve(S)
    
    K=V1[,,i]%*%t(H)%*%InvS #Kalman Gain using guessed value of V1
    
    loglikvec[i] = length(J)*log(2*pi) + log(det(S)) + t(err)%*%InvS%*%err 
    loglikvec[i] = -1/2*loglikvec[i]
    
    X2[,i+1]=X1[,i]+K%*%err #Correct the guess with Kalman Gain.
    X1[,i+1]=muP+rhoP%*%X2[,i+1] #Update next value of state variable(X1).
    I3=matrix(0,ncol=3,nrow=3) 
    diag(I3)=1
    
    V2[,,i+1]=(I3-K%*%H)%*%V1[,,i] #Update the covariance matrix.
    V1[,,i+1]=rhoP%*%V2[,,i+1]%*%t(rhoP)+sigma%*%t(sigma)  #Predict the covariance for next step.
  }
  llf=sum(loglikvec)
  X2<-X2[,-1] #Remove NA column.
  SR=t(delta0+matrix(c(1,1,0),nrow=1)%*%X2) #1st column of X2 is NA. Third row is muliply by 0.
  Xf=X2 #Factors.
  
  return(-llf)
}

#Optimize SRTSM.
EKF_SRTSM(parameters_2)
set.seed(222)
start=rnorm(22)
start[21]=12
#Optimization
test3<-optim(parameters_test,EKF_SRTSM,method="BFGS") 


#WX method

EKF_Shadow_Rate=function(parameters,index=1){
  #Parameterization
  T=ncol(forwardrates)
  J=Maturity
  rhoP = parameters[1:9]; 
  rhoP = matrix(rhoP,ncol=3) 
  muP = parameters[10:12] 
  rhoQ1 = parameters[13]
  rhoQ2 = parameters[14]
  sigma = matrix(c(abs(parameters[15]), 0, 0 , parameters[16], abs(parameters[18]) ,0, parameters[17], parameters[19], abs(parameters[20])),nrow=3,byrow=TRUE); 
  omega = sigma%*%t(sigma)
  delta0 = parameters[21];
  omegaM<-matrix(0,nrow=nrow(forwardrates),ncol=nrow(forwardrates))
  diag(omegaM)=parameters[22]^2
  #rlb=parameters[23]
  rlb=0.25
  #Initialization
  X1<-matrix(NA,nrow=3,ncol=T+1)
  X2<-matrix(NA,nrow=3,ncol=T+1)  #Here we assume X1 and X2 are one the SAME time horizon.
  X1[,1]=0  #X1 starts at 0.
  
  
  #Compute aJ and bJ
  JJ=seq(1,max(J))
  bn=cbind(rhoQ1^JJ,rhoQ2^JJ,JJ*rhoQ2^(JJ-1)) #rho^J in Appendix A.
  Bn=rbind(c(1,1,0),bn)       
  Bn<-apply(Bn,2,cumsum)      #Summation of rho^J
  Bn<-Bn[-nrow(Bn),]          #Summation up to (n-1)
  delta1=c(1,1,0)
  aJ=c()
  SigmaJ=c()
  for(i in 1:length(J))
  {
    aJ[i]<-delta0-0.5*t(Bn[J[i],])%*%((sigma%*%t(sigma))%*%Bn[J[i],])/1200
  }
  
  #Compute SigmaJ.
  cn=rbind(c(1,1,0),bn)
  cn=cn[-nrow(cn),]
  Sigma_step=c()
  SigmaJ=c()
  
  for(i in 1:length(J))
  {
    for(j in 1:Maturity[i])
    {
      Sigma_step[j]<-t(cn[j,])%*%omega%*%cn[j,]
    }
    SigmaJ[i]=sqrt(sum(Sigma_step))
  }
  
  bJ=cbind(rhoQ1^J,rhoQ2^J,J*rhoQ2^(J-1)) 
  V1<-array(NA,dim=c(3,3,T+1))  #V1 is a 3-D arry of covariance at each time
  V2<-array(NA,dim=c(3,3,T+1))  #So is V2.
  V1[,,1]=0
  diag(V1[,,1])=100
  V1[3,3,1]<-V1[3,3,1]/144
  loglikvec<-c()
  
  for(i in 1:T)
  {
    F_h=aJ+bJ%*%X1[,i]  #Prediction of forward rate based on Factors. 
    Z1_temp=(F_h-rlb)/SigmaJ
    Z2_temp=rlb+(F_h-rlb)*pnorm(Z1_temp)+SigmaJ*dnorm(Z1_temp)
    
    err=forwardrates[,i]-Z2_temp #Error in prediction.
    
    H=rep(pnorm(Z1_temp))*bJ
    
    
    S=H%*%V1[,,i]%*%t(H)+omegaM
    
    if (rcond(S) < 1e-8 | is.na(rcond(S)))   return(10^8)
    
    flag=det(S)
    if(!is.finite(flag) | (flag<=0)) return(10^9)
    
    InvS=solve(S)
    
    K=V1[,,i]%*%t(H)%*%InvS #Kalman Gain using guessed value of V1
    
    loglikvec[i] = length(J)*log(2*pi) + log(det(S)) + t(err)%*%InvS%*%err 
    loglikvec[i] = -1/2*loglikvec[i]
    
    X2[,i+1]=X1[,i]+K%*%err #Correct the guess with Kalman Gain.
    X1[,i+1]=muP+rhoP%*%X2[,i+1] #Update next value of state variable(X1).
    I3=matrix(0,ncol=3,nrow=3) 
    diag(I3)=1
    
    V2[,,i+1]=(I3-K%*%H)%*%V1[,,i] #Update the covariance matrix.
    V1[,,i+1]=rhoP%*%V2[,,i+1]%*%t(rhoP)+sigma%*%t(sigma)  #Predict the covariance for next step.
  }
  llf=sum(loglikvec)
  X2<-X2[,-1] #Remove NA column.
  SR=t(delta0+matrix(c(1,1,0),nrow=1)%*%X2) #1st column of X2 is NA. Third row is muliply by 0.
  Xf=X2 #Factors.
  if(index==1) result=SR
  else if(index==2)
  {F_shadow = aJ + bJ%*%Xf #Repeat aJ.
  z1_temp = (F_shadow - rlb)/SigmaJ
  y = rlb + (F_shadow - rlb)*pnorm(z1_temp) + SigmaJ*dnorm(z1_temp)
  result=y
  }
  else result=Xf
  return(result)
}


#with Wu and Xia's parameters
par(mfrow=c(1,1))
rateswx=EKF_Shadow_Rate(parameters,1)
plot(rateswx,type='l')
# factorswx=EKF_Shadow_Rate(parameters,0)
# par(mfrow=c(3,1))
# plot(factorswx[1,],type='l')
# plot(factorswx[2,],type='l')
# plot(factorswx[3,],type='l')
# 
# #with our optimized parameters (with random starting values)
# par(mfrow=c(1,1))
# optimpara=unlist(test3)[1:22]
# rates=EKF_Shadow_Rate(optimpara,1)
# plot(rates,type='l')
# factors=EKF_Shadow_Rate(optimpara,0)
# par(mfrow=c(3,1))
# plot(factors[1,],type='l')
# plot(factors[2,],type='l')
# plot(factors[3,],type='l')

#Forwardrates Fitting for new Maturities.

# Maturity=seq(1,120)
# forwardrates<-MonthlyForwardRate(Fed,Maturity)
# forwardrates<-t(forwardrates)
# forwardrates_fitted<-EKF_Shadow_Rate(parameters_2,2)
# selected_month<-c(3,6,12,24,60,84,120)
# fr_dot<-rowMeans(forwardrates[,264:276])[selected_month]
# plot(seq(0,10,length.out=120),rowMeans(forwardrates_fitted[,264:276]),type='l',lwd=3,main="Fitted Average Forward Rate Curve in 2012",xlab="Maturity/month",ylab="Forward Rate/%")
# points(selected_month/12,fr_dot,pch=8,cex=2)
