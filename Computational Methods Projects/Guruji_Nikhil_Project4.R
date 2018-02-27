#Nikhil Guruji Cohort 1 MFE 2018


#1
t1=0.5
r1=0.05
sig1=0.24
s01=32
k1=30
n1=c(10,20,40,80,100,200,500)
delta1=t1/n1

n_period_eu_call_binom=function(u,d,p,r,s0,k,t,n,delta){
  sum_n_c_k=0
  for(i in 0:n){
    sum_n_c_k=sum_n_c_k + choose(n,i) * p^(i) * (1-p)^(n-i) * pmax(((s0 * u^(i) * d^(n-i)) - k),0)
  } #choose(a,b) gives the value of a choose b
  
  return(exp(-r*n*delta)*sum_n_c_k)
}

#a

c0a=c()
for(ni in 1:length(n1)){
  N=n1[ni]
  del=delta1[ni]
  ctemp=0.5*(exp(-r1*del) + exp((r1 + sig1^2)*del))
  d=ctemp - sqrt(ctemp^2 - 1)
  u=1/d
  p=(exp(r1*del)-d)/(u-d)
  c0a[ni]=n_period_eu_call_binom(u=u,d=d,p=p,r=r1,s0=s01,k=k1,t=t1,n=N,delta=del)
}

plot(n1,c0a,type='l')

#b
c0b=c()
for(ni in 1:length(n1)){
  N=n1[ni]
  del=delta1[ni]
  d=exp(r1*del)*(1-sqrt(exp(del*(sig1^2)) - 1))
  u=exp(r1*del)*(1+sqrt(exp(del*(sig1^2)) - 1))
  p=0.5
  c0b[ni]=n_period_eu_call_binom(u=u,d=d,p=p,r=r1,s0=s01,k=k1,t=t1,n=N,delta=del)
}
plot(n1,c0b,type='l')

#c

c0c=c()
for(ni in 1:length(n1)){
  N=n1[ni]
  del=delta1[ni]
  d=exp((r1-0.5*(sig1^2))*del - sig1*sqrt(del))
  u=exp((r1-0.5*(sig1^2))*del + sig1*sqrt(del))
  p=0.5
  c0c[ni]=n_period_eu_call_binom(u=u,d=d,p=p,r=r1,s0=s01,k=k1,t=t1,n=N,delta=del)
}
plot(n1,c0c,type='l')

#d
c0d=c()
for(ni in 1:length(n1)){
  N=n1[ni]
  del=delta1[ni]
  d=exp(-sig1*sqrt(del))
  u=exp(sig1*sqrt(del))
  p=0.5 + 0.5*((r1-0.5*(sig1^2))*sqrt(del)/sig1)
  c0d[ni]=n_period_eu_call_binom(u=u,d=d,p=p,r=r1,s0=s01,k=k1,t=t1,n=N,delta=del)
}
plot(n1,c0d,type='l')


#combining all plots

plot(n1,c0a,type='l',col="red",ylim=c(3.7,3.78),
     ylab="Call price by different methods",xlab="Number of periods in the binomial model")
lines(n1,c0b,type='l',col="blue")
lines(n1,c0c,type='l',col="black")
lines(n1,c0d,type='l',col="green")
legend("topright",legend=c("Method (a)","Method (b)","Method (c)","Method (d)"),
       col=c("red","blue","black","green"),lty=c(1,1,1,1))

#2

#reading the data
library(httr)
library(quantmod)
library(xts)
url="https://github.com/nikhilg12/nikhil/raw/master/GOOG.csv"
GET(url, write_disk(tf <- tempfile(fileext = ".csv")))

historical_goog=read.csv(tf)

#assumptions: no stock splits occured during the last 5 years
#sorting by date
historical_goog$Date=as.Date(historical_goog$Date)
historical_goog=xts(historical_goog[,2:7],order.by = historical_goog[,1])

#calculating % daily returns
historical_goog$dailyret=100*dailyReturn(historical_goog$Close,type="log")
sig_goog_daily=sd(historical_goog$dailyret)/100 #converting volatility of daily returns to decimals

sig_goog=sig_goog_daily*sqrt(252) #annualizing volatility
current_price_goog=historical_goog[nrow(historical_goog),5] #current price as of February 6 
sig_goog #calculated annual volatility

strike=1.1*current_price_goog
strike=round(strike,-1) #rounding to -1th decimal place(i.e. nearest multiple of 10)
r2=0.02
n2=100
time_till_jan_18=(365-6-12)/365 #time in years till the expiry on Jan 18 2019
delta2=time_till_jan_18/n2

d2=exp((r2-0.5*(sig_goog^2))*delta2 - sig_goog*sqrt(delta2))
u2=exp((r2-0.5*(sig_goog^2))*delta2 + sig_goog*sqrt(delta2))
p2=0.5 + 0.5*((r2-0.5*(sig_goog^2))*sqrt(delta2)/sig_goog)
c0_goog=n_period_eu_call_binom(u=u2,d=d2,p=p2,r=r2,s0=current_price_goog
                               ,k=strike,t=time_till_jan_18,n=n2,delta=delta2)

c0_estimated=as.numeric(c0_goog)
c0_estimated
actual_c0_goog=67.51

#guessing a range for sigma
correct_imp_vol=c()
c0_goog_guess=c()
impvol_guess=seq(0.15,0.3,by=0.001) #15% to 30%
tolerance=0.1 #tolerance in the value of call. Implied vol is considered correct if price is between 67.41 and 67.61
for(i in 1:length(impvol_guess)){
  d2b=exp((r2-0.5*(impvol_guess[i]^2))*delta2 - impvol_guess[i]*sqrt(delta2))
  u2b=exp((r2-0.5*(impvol_guess[i]^2))*delta2 + impvol_guess[i]*sqrt(delta2))
  p2b=0.5 + 0.5*((r2-0.5*(impvol_guess[i]^2))*sqrt(delta2)/impvol_guess[i])
  c0_goog_guess[i]=n_period_eu_call_binom(u=u2b,d=d2b,p=p2b,r=r2,
                                          s0=current_price_goog,k=strike,t=time_till_jan_18,n=n2,delta=delta2)
  c0_goog_guess[i]=as.numeric(c0_goog_guess[i])
  if(abs(c0_goog_guess[i]-actual_c0_goog)<=tolerance){
    correct_imp_vol[i]=impvol_guess[i]
  }
}

correct_imp_vol=mean(correct_imp_vol[which(!is.na(correct_imp_vol))]) #whichever guess is stored
correct_imp_vol

#3

#setting initial values
s03=49
k3=50
r3=0.03
sig3=0.2
t3=0.3846
N3=100
delta3=t3/N3

#calculating u and d and risk neutral probabilities
d3=exp((r3-0.5*(sig3^2))*delta3 - sig3*sqrt(delta3))
u3=exp((r3-0.5*(sig3^2))*delta3 + sig3*sqrt(delta3))
p3=0.5 + 0.5*((r3-0.5*(sig3^2))*sqrt(delta3)/sig3)
c03=n_period_eu_call_binom(u=u3,d=d3,p=p3,r=r3,s0=s03,k=k3,t=t3,n=N3,delta=delta3)

# (i)
s0i=seq(20,80,by=2)
c03i=n_period_eu_call_binom(u=u3,d=d3,p=p3,r=r3,s0=s0i,k=k3,t=t3,n=N3,delta=delta3)

delta_est_i=c()
for(i in 1:(length(s0i))){
  delta_est_i[i]=(c03i[i+1]-c03i[i])/(s0i[i+1]-s0i[i]) #last value will be NA because there is no s0 beyond 80
} #estimating delta by the finite difference method

plot(s0i,delta_est_i,type='l') 


# (ii)
t_expir_ii=seq(0,t3,by=0.01)
delta3ii=t_expir_ii/N3
#s0ii=seq(20,80,by=((80-20)/(length(t_expir_ii)-1)))
epsilon=1

for(j in 1:length(delta3ii)){
  d3ii=exp(-sig3*sqrt(delta3ii))
  u3ii=exp(sig3*sqrt(delta3ii))
  p3ii=0.5 + 0.5*((r3-0.5*(sig3^2))*sqrt(delta3ii)/sig3)

  delta_est_ii=c()
  for(i in 1:(length(delta3ii)-1)){
    delta_est_ii[i]=(n_period_eu_call_binom(u=u3ii[i],d=d3ii[i],p=p3ii[i],r=r3,s0=(s03+epsilon),
                                            k=k3,t=t3,n=N3,delta=delta3ii[i]) 
                     - n_period_eu_call_binom(u=u3ii[i],d=d3ii[i],p=p3ii[i],r=r3,s0=(s03-epsilon),
                                              k=k3,t=t3,n=N3,delta=delta3ii[i]))/(2*epsilon)
  } #estimating delta with respect to time to maturity by central difference method
}  
plot(t_expir_ii[1:(length(t_expir_ii)-1)],delta_est_ii,type='l') 


# (iii)

s0iii=s0i
t_expir_iii=seq(0,t3,by=t3/(length(s0iii)-1)) #defining different points in time to estimate theta
delta3iii=t_expir_iii/N3

theta_mean=c()
for(j in 1:length(s0iii)){ #estimating the average theta for each value of stock price from 20 to 80
  d3iii=exp(- sig3*sqrt(delta3iii))
  u3iii=exp(sig3*sqrt(delta3iii))
  p3iii=0.5 + 0.5*((r3-0.5*(sig3^2))*sqrt(delta3iii)/sig3)

  c03iii=c()
  for (i in 1:length(delta3iii)){
    c03iii[i]=n_period_eu_call_binom(u=u3iii[i],d=d3iii[i],p=p3iii[i],r=r3,s0=s0iii[j]
                                    ,k=k3,t=t3,n=N3,delta=delta3iii[i])
  }
  theta_est_iii=c()
  for(i in 1:length(t_expir_iii)){
    theta_est_iii[i]=(c03iii[i]-c03iii[i+1])/(t_expir_iii[i+1]-t_expir_iii[i])
  }
  theta_mean[j]=mean(theta_est_iii[which(!is.na(theta_est_iii))])
}

plot(s0iii,theta_mean,type='l')

# (iv)
s0iv=seq(20,80,by=2)
gamma_est_iv=c()
for(i in 1:(length(s0iv)-2)){
  gamma_est_iv[i]=(delta_est_i[i+1] - delta_est_i[i])/(s0iv[i+1]-s0iv[i])
} #estimating gamma from delta available from part (i)
plot(s0iv[1:(length(s0iv)-2)],gamma_est_iv,type='l')

# (v)
s0v=s0iv
sig3v=seq(0.01,sig3,by=sig3/(length(s0v)-1)) #defining a volatility vector to estimate vega

#using the above vector to calculate u,d, and risk neutral probabilities
d3v=exp((r3-0.5*(sig3v^2))*delta3 - sig3v*sqrt(delta3))
u3v=exp((r3-0.5*(sig3v^2))*delta3 + sig3v*sqrt(delta3))
p3v=0.5 + 0.5*((r3-0.5*(sig3v^2))*sqrt(delta3)/sig3v)
vega_est=c()
c03v=c()
vega_means=c()
for(j in 1:length(s0v)){ #estimating the average vega for every value of stock price from 20 to 80
  c03v=c()
  for(i in 1:length(sig3v)){
    c03v[i]=n_period_eu_call_binom(u=u3v[i],d=d3v[i],p=p3v[i],r=r3,s0=s0v[j]
                              ,k=k3,t=t3,n=N3,delta=delta3)
  }
  for(i in 1:(length(sig3v)-1)){
  vega_est[i]=(c03v[i+1]-c03v[i])/(sig3v[i+1]-sig3v[i])  
  }
  vega_means[j]=mean(vega_est)
}
plot(s0v,vega_means,type='l')

# (vi)

s0vi=s0iv
r3vi=seq(0.001,r3,by=r3/(length(s0vi)-1)) #defining an interest rate vector to estimate rho

#using the above vector to calculate u,d, and risk neutral probabilities
d3vi=exp((r3vi-0.5*(sig3^2))*delta3 - sig3*sqrt(delta3))
u3vi=exp((r3vi-0.5*(sig3^2))*delta3 + sig3*sqrt(delta3))
p3vi=0.5 + 0.5*((r3vi-0.5*(sig3^2))*sqrt(delta3)/sig3)
rho_est=c()
c03vi=c()
rho_means=c()
for(j in 1:length(s0vi)){ #estimating the average rho for every stock price from 20 to 80
  c03v=c()
  for(i in 1:length(r3vi)){
    c03vi[i]=n_period_eu_call_binom(u=u3vi[i],d=d3vi[i],p=p3vi[i],r=r3vi[i],s0=s0v[j]
                                   ,k=k3,t=t3,n=N3,delta=delta3)
  }
  for(i in 1:(length(r3vi)-1)){
    rho_est[i]=(c03vi[i+1]-c03vi[i])/(r3vi[i+1]-r3vi[i])  
  }
  rho_means[j]=mean(rho_est)
}
plot(s0vi,rho_means,type='l')


#combining all the plots
par(mfrow=c(3,2))
plot(s0i,delta_est_i,type='l',col="red",xlab="Stock Price",ylab="Delta",main="(i) Delta")
plot(t_expir_ii[1:(length(t_expir_ii)-1)],delta_est_ii,type='l',col="darkred",xlab="Time to Maturity",ylab="Delta",main="(ii) Delta")
plot(s0iii,theta_mean,type='l',col="blue1",xlab="Stock Price",ylab="Theta",main="(iii) Theta")
plot(s0iv[1:(length(s0iv)-2)],gamma_est_iv,type='l',col="darkgoldenrod1",xlab="Stock Price",
     ylab="Gamma",main="(iv) Gamma")
plot(s0v,vega_means,type='l',col="darkorchid",xlab="Stock Price",ylab="Vega",main="(v) Vega")
plot(s0vi,rho_means,type='l',col="forestgreen",xlab="Stock Price",ylab="Rho",main="(vi) Rho")

par(mfrow=c(1,1)) #resetting the plots

#4

binomial_tree_generator=function(sig,t,n,r,s0,otype="C",extype="E",k){
  delta=t/n
  
  #calculating u,d, and risk neutral probability
  u=exp(sig*sqrt(delta))
  d=1/u
  p=(exp(r*(t/n)) - d)/(u-d)
  
  #defining a tree for the stock price
  stock_tree=matrix(nrow=(n+1),ncol=(n+1))
  stock_tree[1,1]=s0
  for(i in 2:(n+1)){
    #all the nodes with same row number get multiplied by u from previous stage
    stock_tree[(1:(i-1)),i]=stock_tree[(1:(i-1)),(i-1)]*u 
    stock_tree[i,i]=stock_tree[(i-1),(i-1)]*d #the one node remaining gets multiplied by d from previous stage
  }
  
  payoff_tree=matrix(nrow=(n+1),ncol=(n+1))
  if(otype=="C"){
    payoff_tree[,(n+1)]=pmax((stock_tree[,(n+1)] - k),0) #if option type is call
  } else if(otype=="P"){
    payoff_tree[,(n+1)]=pmax((k-stock_tree[,n+1]),0) #if option type is put
  }
  if(extype=="E"){#if option is European
    for(i in 1:n){
      payoff_tree[(1:(n+1-i)),(n+1-i)]=exp(-r*t/n)*(p*(payoff_tree[(1:(n+1-i)),(n+1-i+1)]) 
                                    + (1-p)*(payoff_tree[(2:(n+1-i+1)),(n+1-i+1)]))
    } 
  } else if(extype=="A" && otype=="C"){ #if option is American Call
      for(i in 1:n){
        payoff_tree[(1:(n+1-i)),(n+1-i)]=pmax(exp(-r*t/n)*(p*(payoff_tree[(1:(n+1-i)),(n+1-i+1)]) 
                                                  + (1-p)*(payoff_tree[(2:(n+1-i+1)),(n+1-i+1)])),
                                         pmax(stock_tree[(1:(n+1-i)),(n+1-i)]-k,0))
      }
  } else if(extype=="A" && otype=="P"){ #if option is American put
      for(i in 1:n){  
        payoff_tree[(1:(n+1-i)),(n+1-i)]=pmax(exp(-r*t/n)*(p*(payoff_tree[(1:(n+1-i)),(n+1-i+1)]) 
                                                       + (1-p)*(payoff_tree[(2:(n+1-i+1)),(n+1-i+1)])),
                                          pmax(k-stock_tree[(1:(n+1-i)),(n+1-i)],0))
      }
  }
  payoff_tree[1,1] #the option price
}

s04=seq(80,120,by=4)
euput4=c()
amput4=c()
for(i in 1:length(s04)){
  euput4[i]=binomial_tree_generator(t=1,sig=0.3,r=0.05,k=100,s0=s04[i],otype="P",extype="E",n=100)
  amput4[i]=binomial_tree_generator(t=1,sig=0.3,r=0.05,k=100,s0=s04[i],otype="P",extype="A",n=100)
}

plot(s04,euput4,type='l',col="blue",ylab="Put Price",xlab="Stock Price")
lines(s04,amput4,type='l',col="red")
legend("topright",legend=c("European Put","American Put")
          ,col=c("blue","red"),lty=c(1,1))

#5

trinomial_tree_generator=function(u,pu,pd,pm,n,s0,otype,k,r,t,extype){
  d=1/u #tree will be recombining because d=1/u
  stock_tree=matrix(nrow=(2*n+1),ncol=(n+1))
  stock_tree[1,1]=s0
  
  for(i in 1:n){
    stock_tree[(1:(2*i-1)),(i+1)]=stock_tree[(1:(2*i-1)),i]*u #goes up by u
    stock_tree[(i+1),(i+1)]=stock_tree[i,i] #nodes on middle line stays the same throughout
    stock_tree[2*i,(i+1)]=stock_tree[(2*i-1),i] #middle nodes following bottommost nodes in earlier periods 
    stock_tree[(2*i+1),(i+1)]=stock_tree[(2*i-1),i]*d #goes down by d
  }
  
  payoff_tree=matrix(nrow=(2*n+1),ncol=(n+1))
  if(otype=="C"){ #if option type is a call
    payoff_tree[,(n+1)]=pmax((stock_tree[,(n+1)] - k),0)
  } else if(otype=="P"){ #if option type is a put
    payoff_tree[,(n+1)]=pmax((k-stock_tree[,n+1]),0)
  }
  if(extype=="E"){ #if option is European
    for(i in 1:n){
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                    + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                    + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)]))
    }
  } else if(extype=="A" && otype=="C"){ #if option is American Call
    for(i in 1:n){
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=pmax(exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                         + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                        + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)])),
                                            pmax(stock_tree[(1:(2*n+1-2*i)),(n+1-i)]-k,0))
    }
  } else if(extype=="A" && otype=="P"){ #if option is American Put
    for(i in 1:n){  
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=pmax(exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                             + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                             + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)])),
      pmax(stock_tree[(1:(2*n+1-2*i)),(n+1-i)]-k,0))
    }
  }
  payoff_tree[1,1] #the option price
}

#a
r5=0.05
sig5=0.24
t5=0.5
k5=30
s05=32

n5=c(10,15,20,40,70,80,100,200,500)
delta5=t5/n5

d5=exp(-sig5*sqrt(3*delta5))
u5=1/d5 #tree will be recombining because u and d are reciprocals

pd5a=(r5*delta5*(1-u5) + (r5*delta5)**2 + (sig5**2)*delta5)/((u5-d5)*(1-d5))
pu5a=(r5*delta5*(1-d5) + (r5*delta5)**2 + (sig5**2)*delta5)/((u5-d5)*(u5-1)) 
pm5a=1-pu5a-pd5a

eu_call_5a=c()
for(i in 1:length(n5)){
  eu_call_5a[i]=trinomial_tree_generator(u=u5[i],pu=pu5a[i],pd=pd5a[i],pm=pm5a[i]
                                         ,n=n5[i],s0=s05,otype="C",k=k5,r=r5,extype="E",t=t5)
}

plot(n5,eu_call_5a,type='l')

#b

trinomial_tree_generator_log_stock=function(dxu,pu,pd,pm,n,s0,otype,k,r,t,extype){
  dxd=-dxu #tree will be recombining because dxd=-dxu
  stock_tree=matrix(nrow=(2*n+1),ncol=(n+1))
  stock_tree[1,1]=log(s0) #starting from the log of stock value and not the stock value itself
  
  for(i in 1:n){
    stock_tree[(1:(2*i-1)),(i+1)]=stock_tree[(1:(2*i-1)),i] + dxu #goes up by dxu
    stock_tree[(i+1),(i+1)]=stock_tree[i,i] #nodes on middle line stays the same throughout
    stock_tree[2*i,(i+1)]=stock_tree[(2*i-1),i] #middle nodes following bottommost nodes in earlier periods 
    stock_tree[(2*i+1),(i+1)]=stock_tree[(2*i-1),i] + dxd #goes down by dxd
  }
  
  payoff_tree=matrix(nrow=(2*n+1),ncol=(n+1))
  if(otype=="C"){ #call
    payoff_tree[,(n+1)]=pmax((exp(stock_tree[,(n+1)]) - k),0) #taking the exponential to convert log stock to stock
  } else if(otype=="P"){ #put
    payoff_tree[,(n+1)]=pmax((k-exp(stock_tree[,n+1])),0)
  }
  if(extype=="E"){ #European
    for(i in 1:n){
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                        + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                        + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)]))
    }
  } else if(extype=="A" && otype=="C"){ #American Call
    for(i in 1:n){
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=pmax(exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                             + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                             + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)])),
                                                pmax(exp(stock_tree[(1:(2*n+1-2*i)),(n+1-i)])-k,0))
    }#everytime we evaluate the payoff at any node, we take the exponential of the logprice to obtain the stock price
  } else if(extype=="A" && otype=="P"){ #American Put
    for(i in 1:n){  
      payoff_tree[(1:(2*n+1-2*i)),(n+1-i)]=pmax(exp(-r*t/n)*(pu*(payoff_tree[(1:(2*n+1-2*i)),(n+1-i+1)]) 
                                                             + pm*(payoff_tree[(2:(2*n+1-2*i+1)),(n+1-i+1)])
                                                             + pd*(payoff_tree[(3:(2*n+1-2*i+2)),(n+1-i+1)])),
                                                pmax(exp(stock_tree[(1:(2*n+1-2*i)),(n+1-i)])-k,0))
    }
  }
  payoff_tree[1,1] #the option price
}

dxu=sig5*sqrt(3*delta5)
dxd=-dxu
pd5b=0.5*(((sig5**2)*delta5 + ((r5-0.5*sig5**2)**2)*delta5**2)/(dxu**2) - 
            (r5-0.5*sig5**2)*delta5/dxu) 
pu5b=0.5*(((sig5**2)*delta5 + ((r5-0.5*sig5**2)**2)*delta5**2)/(dxu**2) + 
            (r5-0.5*sig5**2)*delta5/dxu)
pm5b=1-pu5b-pd5b

eu_call_5b=c()
for(i in 1:length(n5)){
  eu_call_5b[i]=trinomial_tree_generator_log_stock(dxu=dxu[i],pu=pu5b[i],pd=pd5b[i],pm=pm5b[i]
                                         ,n=n5[i],s0=s05,otype="C",k=k5,r=r5,extype="E",t=t5)
}

plot(n5,eu_call_5b,type='l')

#combining plots

plot(n5,eu_call_5a,type='l',col="red",xlab="Number of trinomial periods",
     ylab="Price of European Call option")
lines(n5,eu_call_5b,type='l',col="blue")
legend("bottomright",legend=c("Using Stock Price Process","Using Log Stock Price Process"),
       col=c("red","blue"),lty=c(1,1))
#6

q6=function(s0,k,t,r,sig,n,b1,b2){
  #Halton function from project 3


  #function to generate one halton number
  one_H=function(k,m){
    # m^r <= k < m^(r+1)
    # r*log(m) <= log(k) < (r+1)*log(m)
    # r <= log(k)/log(m) < r+1
    # r = floor(log(k)/log(m))
  
    r=floor(log(k)/log(m))
    coeffs=c()
    remainders=c()
    remainders[1]=k
    for(i in 1:(r+1)){
      #coeffs vector will be traversed in a reverse order(coeffs[1]=ar, coeffs[r]=a1, coeffs[r+1]=a0)
      remainders[i+1]=remainders[i]%%(m^(r-(i-1))) # remainder a_(i-1)=ai%%m^i
      coeffs[i]=floor(remainders[i]/(m^(r-(i-1)))) #quotient = a0,a1,...,ai,...ar
    }
    coeffs=coeffs[length(coeffs):1] #reversing the order of the vector so we get coeffs[1]=a0 and so on...
    h=0
    for(i in 1:length(coeffs)){
      h=h+coeffs[i]/m^i
    }
    h
  }

  #function to generate n Halton numbers
  H=function(n,m){
    hn=c()
    for(i in 1:n){
      hn[i]=one_H(i,m)
    }
    hn
  }

  # box muller method from project 1 with halton sequences instead of uniforms

  box_muller_with_halton=function(n=n,base1=b1,base2=b2){
    h1=H(n,base1)
    h2=H(n,base2)
    z=rep(NA,n)
    for(i in 0:((n/2)-1)){
      j=2*i+1
      z[j]=sqrt(-2*log(h1[i+1]))*cos(2*pi*h2[i+1])
      z[j+1]=sqrt(-2*log(h1[i+1]))*sin(2*pi*h2[i+1])
    }
    z
  }

  # a=box_muller_with_halton(100000,1234,2,7)
  # hist(a) #verifying that a is normal

  # pricing call using Monte Carlo from project 2
  
  z=box_muller_with_halton(n,b1,b2)
  
  st=s0*exp((r-0.5*sig^2)*t + sig*sqrt(t)*z)
  payoff=pmax(st-k,0)
  c_mc=exp(-r*t)*(payoff)
  mean(c_mc) #value of call option through monte carlo
} 

#verifying solution with project 2 Q4 european call
q6(s0=88,k=100,r=0.04,sig=0.2,t=5,b1=2,b2=5,n=10000)

