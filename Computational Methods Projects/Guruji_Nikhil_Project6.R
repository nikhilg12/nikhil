#Nikhil Guruji (Cohort 1) Project 6

#Stock paths generator using Euler method
stock_path_euler=function(m1=100000,n1=100,r,sig,T,s0){
  
  # simulate S: m1 times
  # number of values for each S: n1
  st=matrix(nrow=m1,ncol=n1) # st will be a matrix which will store m1 simulations of S
  z=matrix(rnorm(m1*n1),nrow=m1,ncol=n1) #generate a matrix of random normals
  
  d=T/n1 #dividing t into n1 samples (d is the step size)
  st[,1]=s0 #start with s0
  
  progbar=txtProgressBar(min=2,max=n1,style=3)
  
  for(i in 2:n1){
    #considering the max(Vt,0) for sqrt terms as well as the second term of Vt (full truncation method)
    st[,i]=st[,(i-1)] + r*st[,(i-1)]*d + sig*sqrt(d)*st[,(i-1)]*z[,(i-1)]
    setTxtProgressBar(progbar,i)
  }
  st
}


#1

vol=seq(0.12,0.48,by=0.04)
strike=100




lookback=function(sig){
  payoff_call=vector(length=length(vol))
  payoff_put=vector(length=length(vol))
  stocks=stock_path_euler(m1=100000,n1=200,r=0.03,sig=sig,T=1,s0=98)
  smax=apply(stocks,1,max)
  smin=apply(stocks,1,min)
  
  payoff_call=mean(pmax((smax-strike),0))
  payoff_put=mean(pmax((strike-smin),0))
  c(payoff_call,payoff_put)
}

putcall=matrix(nrow=2,ncol=10) #2 rows because call and put, 10 columns because 10 values of volatility
putcall=sapply(vol,lookback) #runs 10 times
rownames(putcall)=c("Call","Put")
prices=exp(-0.03*1)*putcall
rownames(prices)=c("Call","Put")
plot(vol,prices[1,],type='o',xlab="Volatility",ylab="Call option Price",main="Q1 (Lookback Call)",col="red")
plot(vol,prices[2,],type='o',xlab="Volatility",ylab="Put option Price",main="Q1 (Lookback Put)",col="darkgreen")
#legend("bottomright",bty='n',legend=c("Fixed Strike Lookback Call","Fixed Strike Lookback Put"),lty=c(1,1),
 #      col=c("red","darkgreen"))


#2

default_option=function(lam1=0.2,lam2=0.4,T=5){
  
  
V0=20000
L0=22000
mu=-0.1
sig=0.2
gam=-0.4
#lam1=0.2
#T=5
r0=0.02
del=0.25
#lam2=0.4
alp=0.7
eps=0.95
bet=(eps-alp)/T

n=T*12
R=r0+del*lam2
r=R/12
pmt=L0*r/(1-(1/((1+r)^n)))
a=pmt/r
b=pmt/(r*(1+r)^n)
c=1+r



#first, simulating the value process with jumps for n number of paths


paths=1000
timesteps=n
delta=T/timesteps


Z_2 = matrix(rnorm(paths*n), paths,n)
W_2 = sqrt(delta) * Z_2
P <- matrix(rpois(paths*n, lam1 * delta), paths, n)

vtpaths = matrix(0,paths,n)
vtpaths[,1] <- V0

for(j in 1:(n-1)){
  vtpaths[,(j+1)] = vtpaths[,j] + mu * vtpaths[,j] * delta + sig * vtpaths[,j] * W_2[,j] + gam * vtpaths[,j] * P[,j]
}

#close(pb)
#vtpaths=vtpaths[-1,] #omitting the first row because it was only for the sake of declaration
#length(vtpaths)
tseq=seq(0,T,by=delta)
#plot(tseq,vtpaths[7,],type='l')
#second, simulating the loan process
Lt=a-b*c^(12*tseq)
qt=alp+bet*tseq

benchmark=qt*Lt
#Stopping time Q
Q_time <- c()
for (k in 1:paths){
  for (l in 1:n)
  {
    if (vtpaths[k,l] <= benchmark[l])
    {
      Q_time[k] = l
      break
    }
    else
    {
      Q_time[k] = 100
    }
  }
}

#Stopping time S
S_time <- c()
P_S <- matrix(rpois(paths*n, lam2 * delta), paths, n)
for (k in 1:paths){
  for (l in 1:n)
  {
    if (P_S[k,l] > 0)
    {
      S_time[k] = l
      break
    }
    else
    {
      S_time[k] = 100
    }
  }
}

#Payoff of the borrower
payoffQ <- c()
payoffS <- c()
payoff <- c()
for (z in 1:paths){
  if (Q_time[z] < S_time[z])
  {
    if (Q_time[z] == 100)
      payoff[z] <- 0
    else
    {
      payoffQ[z] <- max(Lt[Q_time[z]] - eps * vtpaths[z,Q_time[z]] , 0)
      payoff[z] <- exp(-r0 * delta * Q_time[z]) * payoffQ[z]
    }
    
  }
  else
  {
    if (S_time[z] == 100)
    {
      payoff[z] <- 0
    }
    
    else
    {
      payoffS[z] <- abs(Lt[S_time[z]] - eps * vtpaths[z,S_time[z]])
      payoff[z] <- exp(-r0 * delta * S_time[z]) * payoffS[z]
      
    }
  }
}

#Default Probability
tau <- pmin(Q_time,S_time)
prob_default <- (paths - length(which(tau == 100)))/paths

#Expected Stopping Time
expected_ex_time <- (sum(tau) - length(which(tau == 100))*100) / (paths - length(which(tau == 100))) / 12

final_answer=c(mean(payoff),prob_default,expected_ex_time)
names(final_answer)=c("Value of the option","Probability of Default","Expected Exercise Time given tau<T")
final_answer

}





#There will be warnings because Infinity is also generated
#There might be "replacement length zero" errors. Please run the code again. This occurs rarely if simultaneous 
#jumps are too low

default_option()


par(mfrow=c(2,1))

#a
lam2a1=0.4
lam1a1=seq(0.05,0.4,by=0.05)
Ta=seq(3,8,by=1)
basevalue=default_option()

color8=c("red","blue","chartreuse1","black","burlywood3","cyan2","gold1","deeppink")
#drawing a blank plot to add lines
plot(Ta,rep(basevalue[1],length(Ta)),type='n',ylab="Value of option",xlab="Time")



#Note: In the rare case of more than two jumps right at the beginning (time=0), the code may throw replacement length zero
# error. Please run the code again if that occurs. 
# this error can be resolved by increasing the number of timesteps and thus decreasing the step size. But, speed matters too!

valuea1=list(length=length(Ta))

for(l in 1:length(lam1a1)){
  #l=1
  valuea1[[l]]=sapply(Ta,default_option,lam1=lam1a1[l],lam2=lam2a1)
  lines(Ta,valuea1[[l]][1,],col=color8[l],lwd=2)
}
legend("bottomright",bty='n',legend=paste0("lambda1= ",lam1a1),lwd=rep(2,length(lam1a1)),
       col=color8)
#we get warnings because Infinity is generated for times when option is not exercised


lam2a2=seq(0,0.8,by=0.1)
lam1a2=0.2
Ta=seq(3,8,by=1)


#drawing a blank plot to add lines
plot(Ta,rep(basevalue[1],length(Ta)),type='n',ylab="Value of option",xlab="Time",ylim=c(3000,5500))

color9=c("red","blue","chartreuse1","black","burlywood3","cyan2","gold1","deeppink","violet")


#Note: In the rare case of more than two jumps right at the beginning (time=0), the code may throw replacement length zero
# error. Please run the code again if that occurs. 
# this error can be resolved by increasing the number of timesteps and thus decreasing the step size. But, speed matters too!

valuea2=list(length=length(Ta))

for(l in 1:length(lam2a2)){
  #l=1
  valuea2[[l]]=sapply(Ta,default_option,lam1=lam1a2,lam2=lam2a2[l])
  lines(Ta,valuea2[[l]][1,],col=color9[l],lwd=2)
}
legend("bottomright",bty='n',legend=paste0("lambda2= ",lam2a2),lwd=rep(2,length(lam2a2)),
       col=color9)
#we get warnings because Infinity is generated for times when option is not exercised

#b
lam2b1=0.4
lam1b1=seq(0.05,0.4,by=0.05)

#drawing a blank plot to add lines
plot(Ta,rep(basevalue[2],length(Ta)),type='n',ylab="Probability of Default",xlab="Time")


for(l in 1:length(lam1b1)){
  #l=1
  lines(Ta,(valuea1[[l]][2,]),col=color8[l],lwd=2)
}
legend("bottomright",bty='n',legend=paste0("lambda1= ",lam1b1),lwd=rep(2,length(lam1b1)),
       col=color8)


lam2b2=seq(0,0.8,by=0.1)
lam1b2=0.2


#drawing a blank plot to add lines
plot(Ta,rep(basevalue[2],length(Ta)),type='n',ylab="Probability of Default",xlab="Time",ylim=c(0.1,1))


for(l in 1:length(lam2b2)){
  #l=1
  lines(Ta,valuea2[[l]][2,],col=color9[l],lwd=2)
}
legend("bottomright",bty='n',legend=paste0("lambda2= ",lam2b2),lwd=rep(2,length(lam2b2)),
       col=color9)
#c
lam2c1=0.4
lam1c1=seq(0.05,0.4,by=0.05)

#drawing a blank plot to add lines
plot(Ta,rep(basevalue[3],length(Ta)),type='n',ylab="Expected Exercise Time",xlab="Time")


for(l in 1:length(lam1c1)){
  #l=1
  lines(Ta,valuea1[[l]][3,],col=color8[l],lwd=2)
}
legend("topright",bty='n',legend=paste0("lambda1= ",lam1c1),lwd=rep(2,length(lam1c1)),
       col=color8)


lam2c2=seq(0,0.8,by=0.1)
lam1c2=0.2


#drawing a blank plot to add lines
plot(Ta,rep(basevalue[3],length(Ta)),type='n',ylab="Expected Exercise Time",xlab="Time")


for(l in 1:length(lam2c2)){
  #l=1
  lines(Ta,valuea2[[l]][3,],col=color9[l],lwd=2)
}
legend("topright",bty='n',legend=paste0("lambda2= ",lam2c2),lwd=rep(2,length(lam2c2)),
       col=color9)

