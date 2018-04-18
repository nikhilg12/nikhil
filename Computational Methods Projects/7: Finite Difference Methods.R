#Nikhil Guruji Cohort 1 Project 7


euput_finite_diff=function(s0_start=10,r=0.04,sig=0.2,strike=10,T_mat=0.5,method,dtmul=1){
    
  dt=0.002
  # dtmul is a multiplier of dt for dxt
  dxt=sig*sqrt(dtmul*dt)

  s0min=0.0001 #small value but not zero so that log is not infinity
  s0max=s0_start+20
  xlower=seq(log(s0_start),log(s0min),by=-dxt)
  xupper=seq(log(s0max),(log(s0_start)-dxt),by=-dxt)

  x0=c(xupper,xlower)

  j_index=seq(length(xupper),(-length(xlower)+1),by=-1)
  i_index=seq(0,T_mat,by=dt)
  s0=exp(x0)
  s0_grid=matrix(0,nrow=length(j_index),ncol=length(i_index))
  #s0_grid=s0

  s0_grid[,(1:ncol(s0_grid))]=s0#*u^(j_index)

  payoffeuput=pmax((strike-s0_grid[,ncol(s0_grid)]),0)


  if(method=="E"){
    
    #Explicit Finite Difference Method
    
    
    pu=dt*(((sig^2)/(2*dxt^2)) + ((r- 0.5*sig^2)/(2*dxt)))
    pd=dt*(((sig^2)/(2*dxt^2)) - ((r- 0.5*sig^2)/(2*dxt)))
    pm=1- dt*(sig^2)/(dxt^2) - r*dt
    
    #creating matrix A
    A=matrix(0,nrow=length(x0),ncol=length(x0))
    diag(A)=c(pu,rep(pm,length(x0)-2),pd)
    A[which(A==pm)+1]=pu
    A[which(A==pm)-1]=pd
    A[nrow(A),((ncol(A)-2):ncol(A))]=c(pu,pm,pd)
    A[1,(1:3)]=c(pu,pm,pd)
    A[2,1]=pu
    A[(nrow(A)-1),ncol(A)]=pd
    
    #creating matrix B
    B=matrix(rep(0,length(x0)),nrow=length(x0))
    B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
    price=matrix(nrow=nrow(s0_grid),ncol=ncol(s0_grid))
    price[,ncol(price)]=payoffeuput
    for (i in (ncol(s0_grid)-1):1){
      #i=8
      price[,i]=A%*%price[,(i+1)] + B
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else if(method=="I"){
    
    #Implicit Finite Difference Method
    
    
    pu=-dt*(((sig^2)/(2*dxt^2)) + ((r- 0.5*sig^2)/(2*dxt))) #just the negative of pu for explicit method
    pd=-dt*(((sig^2)/(2*dxt^2)) - ((r- 0.5*sig^2)/(2*dxt))) #just the negative of pd for explicit method
    pm=1 + dt*(sig^2)/(dxt^2) + r*dt
    
    #creating matrix A
    A=matrix(0,nrow=length(x0),ncol=length(x0))
    diag(A)=c(pu,rep(pm,length(x0)-2),pd)
    A[which(A==pm)+1]=pu
    A[which(A==pm)-1]=pd
    A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,-1,1) #different from the A matrix for call option
    A[1,(1:3)]=c(1,-1,0)
    A[2,1]=pu
    A[(nrow(A)-1),ncol(A)]=pd
    
    #creating matrix B
    B=matrix(rep(0,length(x0)),nrow=length(x0))
    B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
    
    price=matrix(nrow=nrow(s0_grid),ncol=ncol(s0_grid))
    price[,ncol(price)]=payoffeuput
    A_inv=solve(A)
    for (i in (ncol(s0_grid)-1):1){
      B[(2:(nrow(B)-1)),1]=price[(2:(nrow(price)-1)),(i+1)]
      price[,i]=A_inv%*%B
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else if(method=="C"){ 
    
    #Crank Nicolson method
    
    pu=-0.5*dt*(((sig^2)/(2*dxt^2)) + ((r- 0.5*sig^2)/(2*dxt))) #just the 0.5*pu for implicit method
    pd=-0.5*dt*(((sig^2)/(2*dxt^2)) - ((r- 0.5*sig^2)/(2*dxt))) #just the 0.5*pd for implicit method
    pm=1 + dt*(sig^2)/(2*dxt^2) + 0.5*r*dt
    
    #creating matrix A (will be the same as that of implicit method)
    A=matrix(0,nrow=length(x0),ncol=length(x0))
    diag(A)=c(pu,rep(pm,length(x0)-2),pd)
    A[which(A==pm)+1]=pu
    A[which(A==pm)-1]=pd
    A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,-1,1) #different from the A matrix for call option
    A[1,(1:3)]=c(1,-1,0)
    A[2,1]=pu
    A[(nrow(A)-1),ncol(A)]=pd
    
    #creating matrix B
    B=matrix(rep(0,length(x0)),nrow=length(x0))
    
    A_intermediate=A #creating an intermediate matrix to create B
    
    #replacing the values
    A_intermediate[which(A_intermediate==pu)]=-pu
    A_intermediate[which(A_intermediate==pd)]=-pd
    A_intermediate[which(A_intermediate==pm)]=-(pm-2)
    
    price=matrix(nrow=nrow(s0_grid),ncol=ncol(s0_grid))
    price[,ncol(price)]=payoffeuput
    A_inv=solve(A)
    for (i in (ncol(s0_grid)-1):1){
      
      #calculating B matrix
      B=A_intermediate%*%price[,(i+1)]
      
      #replacing the extreme values of B
      B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
      B[1,1]=0
      
      price[,i]=A_inv%*%B
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else{
    print("Invalid method")
  }
}

black_scholes_eu_put_pricing=function(s0,t=0.5,k=10,r=0.04,s=0.2){
  
  #using the Black Scholes formula
  d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
  d2= d1 - s*sqrt(t)
  
  p_bs=-s0*pnorm(-d1) + k*exp(-r*t)*pnorm(-d2)
  p_bs
}

Pa1=euput_finite_diff(method="E",dtmul=1)
Pa3=euput_finite_diff(method="E",dtmul=3)
Pa4=euput_finite_diff(method="E",dtmul=4)

Pb1=euput_finite_diff(method="I",dtmul=1)
Pb3=euput_finite_diff(method="I",dtmul=3)
Pb4=euput_finite_diff(method="I",dtmul=4)

Pc1=euput_finite_diff(method="C",dtmul=1)
Pc3=euput_finite_diff(method="C",dtmul=3)
Pc4=euput_finite_diff(method="C",dtmul=4)

putvalues=matrix(c(Pa1,Pa3,Pa4,Pb1,Pb3,Pb4,Pc1,Pc3,Pc4),nrow=3)
rownames(putvalues)=c("mutiplier=1","multiplier=3","multiplier=4")
colnames(putvalues)=c("Explicit Finite Difference Method","Implicit Finite Difference Method","Crank Nicolson Method")
#View(putvalues)
test_s0=c(4:16)

putvaluesa1=sapply(test_s0,euput_finite_diff,method="E",dtmul=1)
putvaluesa3=sapply(test_s0,euput_finite_diff,method="E",dtmul=3)
putvaluesa4=sapply(test_s0,euput_finite_diff,method="E",dtmul=4)

putvaluesb1=sapply(test_s0,euput_finite_diff,method="I",dtmul=1)
putvaluesb3=sapply(test_s0,euput_finite_diff,method="I",dtmul=3)
putvaluesb4=sapply(test_s0,euput_finite_diff,method="I",dtmul=4)

putvaluesc1=sapply(test_s0,euput_finite_diff,method="C",dtmul=1)
putvaluesc3=sapply(test_s0,euput_finite_diff,method="C",dtmul=3)
putvaluesc4=sapply(test_s0,euput_finite_diff,method="C",dtmul=4)

bs_values=sapply(test_s0,black_scholes_eu_put_pricing)


errors=function(putvalues1,putvalues3,putvalues4){
  errors_=cbind((putvalues1-bs_values),(putvalues3-bs_values),
                             (putvalues4-bs_values))
  rownames(errors_)=test_s0
  colnames(errors_)=c("Explicit Method","Implicit Method","Crank Nicolson Method")
  errors_
}


errors_1=errors(putvaluesa1,putvaluesb1,putvaluesc1)
errors_3=errors(putvaluesa3,putvaluesb3,putvaluesc3)
errors_4=errors(putvaluesa4,putvaluesb4,putvaluesc4)


#2

am_option_finite_diff=function(s0_start=10,r=0.04,sig=0.2,strike=10,T_mat=0.5,method,ds=0.5,otype){
  
  dt=0.002
  
  if(ds!=1.5){
    s0=seq((s0_start+20),0,by=-ds)
  }else{
    s0=rev(seq(s0_start%%3,(s0_start+20),by=ds)) #by doing this, s0_start will always be included in this sample
  }
  
  j_index=seq(1,(length(s0)))
  i_index=seq(0,T_mat,by=dt)
  s0_grid=matrix(0,nrow=length(j_index),ncol=length(i_index))
  
  N=length(j_index)
  M=length(i_index)
  
  s0_grid[,(1:ncol(s0_grid))]=s0

  if(otype=="P"){
    payoff=pmax((strike-s0_grid),0)
  } else {payoff=pmax((s0_grid-strike),0)}
  
  if(method=="E"){
    
    #Explicit Finite Difference Method
    
    
    pu=0.5*dt*(r*j_index[-N] + sig^2*j_index[-N]^2)
    pd=0.5*dt*(-r*j_index[-N] + sig^2*j_index[-N]^2)
    pm=1- dt*(j_index[-N]^2*sig^2 + r)
    
    #creating matrix A
    A=matrix(0,nrow=N+1,ncol=N+1)
    diag(A)=c(pu[N-1],rev(pm),pd[1]) #setting the diagonal elements
    
    A[which(A %in% pm)+1]=c(rev(pu[-N+1]),pu[1]) #setting elements below the diagonal
    A[which(A %in% pu[-N+1])-2]=c(pm[N-1],rev(pd[-1])) #setting elements above the diagonal
    A[nrow(A),((ncol(A)-2):ncol(A))]=A[(nrow(A)-1),((ncol(A)-2):ncol(A))]=c(pu[1],pm[1],pd[1]) #last two rows
    A[1,(1:3)]=A[2,(1:3)]=c(pu[N-1],pm[N-1],pd[N-1]) #first two rows
    
    #creating matrix B
    B=matrix(rep(0,N+1),nrow=N+1)
    if(otype=="P"){
      B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
    }else{B[1,1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]}
    price=matrix(nrow=(N+1),ncol=M)
    price[,ncol(price)]=c(payoff[,ncol(payoff)],payoff[nrow(payoff),ncol(payoff)])
    for (i in (ncol(s0_grid)-1):1){
      price[,i]=A%*%price[,(i+1)] + B
      price[,i]=pmax(price[,i],c(payoff[,i],payoff[nrow(payoff),ncol(payoff)]))
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else if(method=="I"){
    
    #Implicit Finite Difference Method
    
    pu=-0.5*dt*(r*j_index[-N] + sig^2*j_index[-N]^2) #just the negative of pu for explicit method
    pd=-0.5*dt*(-r*j_index[-N] + sig^2*j_index[-N]^2) #just the negative of pd for explicit method
    pm=1 + dt*(j_index[-N]^2*sig^2 + r)
    
    #creating matrix A
    A=matrix(0,nrow=N+1,ncol=N+1)
    diag(A)=c(pu[N-1],rev(pm),pd[1]) #setting the diagonal elements
    
    A[which(A %in% pm)+1]=c(rev(pu[-N+1]),pu[1]) #setting elements below the diagonal
    A[which(A %in% pu[-N+1])-2]=c(pm[N-1],rev(pd[-1])) #setting elements above the diagonal
    A[(nrow(A)-1),((ncol(A)-2):ncol(A))]=c(pu[1],pm[1],pd[1]) #second last row
    A[2,(1:3)]=c(pu[N-1],pm[N-1],pd[N-1]) #second row
    
    A[1,(1:3)]=c(1,-1,0) #first row
    if(otype=="C"){
      A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,1,-1)
    }else {A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,-1,1)} #last row different for call and put option
    
    #creating matrix B
    B=matrix(rep(0,N+1),nrow=N+1)
    if(otype=="P"){
      B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
    }else{B[1,1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]}
    
    price=matrix(nrow=(N+1),ncol=M)
    price[,ncol(price)]=c(payoff[,ncol(payoff)],payoff[nrow(payoff),ncol(payoff)])
    
    A_inv=solve(A)
    
    for (i in (ncol(s0_grid)-1):1){
      B[(2:(nrow(B)-1)),1]=price[(2:(nrow(price)-1)),(i+1)]
      price[,i]=A_inv%*%B
      price[,i]=pmax(price[,i],c(payoff[,i],payoff[nrow(payoff),ncol(payoff)]))
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else if(method=="C"){ 
    
    #Crank Nicolson method
    
    pu=-0.5*0.5*dt*(r*j_index[-N] + sig^2*j_index[-N]^2) #just the 0.5*pu for implicit method
    pd=-0.5*0.5*dt*(-r*j_index[-N] + sig^2*j_index[-N]^2) #just the 0.5*pd for implicit method
    pm=1 + 0.5*dt*(j_index[-N]^2*sig^2 + r)
    
    #creating matrix A (will be the same as that of implicit method)
    A=matrix(0,nrow=N+1,ncol=N+1)
    diag(A)=c(pu[N-1],rev(pm),pd[1]) #setting the diagonal elements
    
    A[which(A %in% pm)+1]=c(rev(pu[-N+1]),pu[1]) #setting elements below the diagonal
    A[which(A %in% pu[-N+1])-2]=c(pm[N-1],rev(pd[-1])) #setting elements above the diagonal
    A[(nrow(A)-1),((ncol(A)-2):ncol(A))]=c(pu[1],pm[1],pd[1]) #second last row
    A[2,(1:3)]=c(pu[N-1],pm[N-1],pd[N-1]) #second row
    
    A[1,(1:3)]=c(1,-1,0) #first row
    if(otype=="C"){
      A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,1,-1)
    }else {A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,-1,1)} #last row different for call and put option
    
    #creating matrix B
    B=matrix(rep(0,N+1),nrow=N+1)

    
    #creating an intermediate matrix to create B
    pu_intermediate=-pu
    pd_intermediate=-pd
    pm_intermediate=-(pm-2)
    A_intermediate=matrix(0,nrow=N+1,ncol=N+1)
    diag(A_intermediate)=c(pu_intermediate[N-1],rev(pm_intermediate),pd_intermediate[1])
    
    A_intermediate[which(A_intermediate %in% pm_intermediate)+1]=c(rev(pu_intermediate[-N+1]),pu_intermediate[1])
    A_intermediate[which(A_intermediate %in% pu_intermediate[-N+1])-2]=c(pm_intermediate[N-1],rev(pd_intermediate[-1]))
    A_intermediate[(nrow(A_intermediate)-1),((ncol(A_intermediate)-2):ncol(A_intermediate))]=c(pu_intermediate[1],pm_intermediate[1],pd_intermediate[1])
    A_intermediate[2,(1:3)]=c(pu_intermediate[N-1],pm_intermediate[N-1],pd_intermediate[N-1]) #second row
    
    A_intermediate[1,(1:3)]=c(1,-1,0) #first row
    if(otype=="C"){
      A_intermediate[nrow(A_intermediate),((ncol(A_intermediate)-2):ncol(A_intermediate))]=c(0,1,-1)
    }else {A_intermediate[nrow(A_intermediate),((ncol(A_intermediate)-2):ncol(A_intermediate))]=c(0,-1,1)} 

    
    price=matrix(nrow=(N+1),ncol=M)
    price[,ncol(price)]=c(payoff[,ncol(payoff)],payoff[nrow(payoff),ncol(payoff)])
    
    A_inv=solve(A)
    for (i in (ncol(s0_grid)-1):1){
      
      #calculating B matrix
      B=A_intermediate%*%price[,(i+1)]
      
      #replacing the extreme values of B
      if(otype=="P"){
        B[nrow(B),1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]
      }else{B[1,1]=s0_grid[(nrow(s0_grid)-1),1] - s0_grid[nrow(s0_grid),1]}
      
      price[,i]=A_inv%*%B
    }
    price[which(s0_grid[,1]==as.character(s0_start))]
  } else{
    print("Invalid method")
  }
}

call_exp_05=sapply(test_s0,am_option_finite_diff,method="E",otype="C")
call_imp_05=sapply(test_s0,am_option_finite_diff,method="I",otype="C")
call_crank_05=sapply(test_s0,am_option_finite_diff,method="C",otype="C")

put_exp_05=sapply(test_s0,am_option_finite_diff,method="E",otype="P")
put_imp_05=sapply(test_s0,am_option_finite_diff,method="I",otype="P")
put_crank_05=sapply(test_s0,am_option_finite_diff,method="C",otype="P")


call_exp_1=sapply(test_s0,am_option_finite_diff,method="E",otype="C",ds=1)
call_imp_1=sapply(test_s0,am_option_finite_diff,method="I",otype="C",ds=1)
call_crank_1=sapply(test_s0,am_option_finite_diff,method="C",otype="C",ds=1)

put_exp_1=sapply(test_s0,am_option_finite_diff,method="E",otype="P",ds=1)
put_imp_1=sapply(test_s0,am_option_finite_diff,method="I",otype="P",ds=1)
put_crank_1=sapply(test_s0,am_option_finite_diff,method="C",otype="P",ds=1)


call_exp_15=sapply(test_s0,am_option_finite_diff,method="E",otype="C",ds=1.5)
call_imp_15=sapply(test_s0,am_option_finite_diff,method="I",otype="C",ds=1.5)
call_crank_15=sapply(test_s0,am_option_finite_diff,method="C",otype="C",ds=1.5)

put_exp_15=sapply(test_s0,am_option_finite_diff,method="E",otype="P",ds=1.5)
put_imp_15=sapply(test_s0,am_option_finite_diff,method="I",otype="P",ds=1.5)
put_crank_15=sapply(test_s0,am_option_finite_diff,method="C",otype="P",ds=1.5)

callput_05=cbind(call_exp_05,call_imp_05,call_crank_05,put_exp_05,put_imp_05,put_crank_05)
colnames(callput_05)=c("Call (Explicit Method)","Call (Implicit Method)","Call (Crank Nicolson Method)",
                       "Put (Explicit Method)","Put (Implicit Method)","Put (Crank Nicolson Method)")
rownames(callput_05)=test_s0
#View(callput_05)

callput_1=cbind(call_exp_1,call_imp_1,call_crank_1,put_exp_1,put_imp_1,put_crank_1)
colnames(callput_1)=c("Call (Explicit Method)","Call (Implicit Method)","Call (Crank Nicolson Method)",
                       "Put (Explicit Method)","Put (Implicit Method)","Put (Crank Nicolson Method)")
rownames(callput_1)=test_s0
#View(callput_1)

callput_15=cbind(call_exp_15,call_imp_15,call_crank_15,put_exp_15,put_imp_15,put_crank_15)
colnames(callput_15)=c("Call (Explicit Method)","Call (Implicit Method)","Call (Crank Nicolson Method)",
                      "Put (Explicit Method)","Put (Implicit Method)","Put (Crank Nicolson Method)")
rownames(callput_15)=test_s0
#View(callput_15)

#graph 1
plot(test_s0,call_exp_05,type='l',main="Call for ds=0.5",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,call_imp_05,col="chartreuse3")
lines(test_s0,call_crank_05,col="blue")
legend("topleft",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

#graph 2
plot(test_s0,put_exp_05,type='l',main="Put for ds=0.5",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,put_imp_05,col="chartreuse3")
lines(test_s0,put_crank_05,col="blue")
legend("topright",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

#graph 3
plot(test_s0,call_exp_1,type='l',main="Call for ds=1",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,call_imp_1,col="chartreuse3")
lines(test_s0,call_crank_1,col="blue")
legend("topleft",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

#graph 4
plot(test_s0,put_exp_1,type='l',main="Put for ds=1",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,put_imp_1,col="chartreuse3")
lines(test_s0,put_crank_1,col="blue")
legend("topright",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

#graph 5
plot(test_s0,call_exp_15,type='l',main="Call for ds=1.5",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,call_imp_15,col="chartreuse3")
lines(test_s0,call_crank_15,col="blue")
legend("topleft",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

#graph 6
plot(test_s0,put_exp_15,type='l',main="Put for ds=1.5",ylab="Price",xlab="Stock price",col="red")
lines(test_s0,put_imp_15,col="chartreuse3")
lines(test_s0,put_crank_15,col="blue")
legend("topright",legend=c("Explicit Method","Implicit Method","Crank Nicolson Method"),
       lty=c(1,1,1),col=c("red","chartreuse3","blue"))

