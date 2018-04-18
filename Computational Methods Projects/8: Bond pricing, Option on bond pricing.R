#Nikhil Guruji Cohort 1 Project 8

#1
vasicek_montecarlo_zcb=function(par,T_bond=0.5,r0=0.05,kappa=0.82,r_bar=0.05,sig=0.1,paths=10000){
  
  timesteps=ceiling(T_bond*365) #ensuring that the number of timesteps= no of days in that period
  delta=T_bond/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(delta)*z[,(j-1)]
  }
  R=delta*rowSums(r_t)
  price_zcb=par*mean(exp(-R))
  
  price_zcb
}

#a 
price_1a=vasicek_montecarlo_zcb(T_bond=0.5,par=1000)


#b 

vasicek_coupon_bond=function(couptimes,coupvalues,r_start,Npaths){
  prices=mapply(vasicek_montecarlo_zcb,T_bond=couptimes,par=coupvalues,r0=r_start,paths=Npaths) 
  #mapply used to apply one function to two argument vectors
  
  price_cb=sum(prices)
  price_cb
}

couptimes_1b=seq(0.5,4,by=0.5)
coupvalues_1b=c(rep(30,7),1030)

price_1b=vasicek_coupon_bond(couptimes = couptimes_1b,coupvalues = coupvalues_1b,r_start=0.05,Npaths=5000)

#c 
vasicek_explicit_price=function(T_,t_,r_t,r0=0.05,kappa=0.82,r_bar=0.05,sig=0.1,par){
  BtT=(1/kappa)*(1- exp(-kappa*(T_-t_)))
  AtT=exp((r_bar - 0.5*(sig^2/kappa^2))*(BtT - (T_ - t_) ) - BtT^2*sig^2/(4*kappa))
  price_explicit=par*AtT*exp(-BtT*r_t)
  price_explicit
}
par=1000
strike=980

vasicek_montecarlo_option_explicit=function(strike,par=0,cv=0,ct=0,btype="Z",
                                      T_option=0.25,T_bond=0.5,r0=0.05,kappa=0.82,r_bar=0.05,sig=0.1,paths=10000){
  
  timesteps=ceiling(T_option*365)
  delta=T_option/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(delta)*z[,(j-1)]
  }
  
  
  r_t_option=r_t[,ncol(r_t)]
  #call the function to calculate explicit price for given parameters
  
  if(btype=="C"){
    price_exp=matrix(nrow=length(r_t_option),ncol=length(cv))
    for(i in 1:length(cv)){
      price_exp[,i]=vasicek_explicit_price(par=cv[i],T_=ct[i],t_=T_option,r_t=r_t_option) 
      #pass the original values of time and coupon to the function (not actual values - 0.25)
    }
    price_explicit=rowSums(price_exp)
  }else{
    price_explicit=vasicek_explicit_price(T_=T_bond,t_=T_option,par=par,r_t=r_t_option)
  }
  payoff=pmax((price_explicit-strike),0)
  R=rowSums(r_t)*delta
  disc=exp(-R)*payoff
  price_option=mean(disc)
  price_option
}
price_1c=vasicek_montecarlo_option_explicit(par=par,strike=strike)

#d 
coupvalues_025=coupvalues_1b
couptimes_025=couptimes_1b-0.25

vasicek_montecarlo_option_montecarlo=function(strike,par,cv,ct,
                                              T_option=0.25,T_bond=0.5,r0=0.05,kappa=0.82,r_bar=0.05,sig=0.1,paths=100){
  timesteps=ceiling(T_option*365)
  delta=T_option/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(delta)*z[,(j-1)]
  }
  rT=r_t[,ncol(r_t)]
  price_T=vector(length=length(rT))
  for(i in 1:length(rT)){
    price_T[i]=vasicek_coupon_bond(couptimes=ct,coupvalues = cv,r_start = rT[i],Npaths=100)
  }
  payoff=pmax((price_T-strike),0)
  R=rowSums(r_t)*delta
  disc=exp(-R)*payoff
  price_option=mean(disc)
  price_option
}

price_1d=vasicek_montecarlo_option_montecarlo(strike=strike,par=par,cv=coupvalues_025,ct=couptimes_025)


#e
vasicek_call_price_explicit=function(t_.=0.25,couptimes=couptimes_1b,coupvalues=coupvalues_1b,par.=par,r0=0.05,
                                     kappa=0.82,r_bar=0.05,sig=0.1,strike=980){
  r_star_func=function(r_star){
    prices=mapply(vasicek_explicit_price,t_=t_.,T_=couptimes,r_t=r_star,par=par.)/par.
    prices=prices*coupvalues
    sum(prices)-strike
  }
  r_star_root=uniroot(r_star_func,interval=c(0.1,0.15))$root
  
  Ki=mapply(vasicek_explicit_price,T_=couptimes,t_=t_.,r_t=r_star_root,par=par.)/par.
  PtT=vasicek_explicit_price(T_=t_.,t_=t_.,r_t=r0,par=par.)/par.
  PtTi=mapply(vasicek_explicit_price,t_=t_.,T_=couptimes,r_t=r0,par=par.)/par.
  sigmap <- sig * (1-exp(-kappa * (couptimes - t_.)))/kappa * sqrt((1 - exp(-2 * kappa * (t_.- 0)))/(2*kappa))
  d1 <- (1/sigmap)*log(PtTi/(Ki * PtT)) + 0.5*sigmap
  d2=d1-sigmap
  call_exp <- PtTi * pnorm(d1) - Ki * PtT * pnorm(d2)
  sum(call_exp * coupvalues)
}


price_1e=vasicek_call_price_explicit()
# another method
# vasicek_option_explicit=function(t_=0,T_=0.25,S_=0.5,strike=980,cv,ct,
#                                  r0=0.05,kappa=0.82,r_bar=0.05,sig=0.1,par=1000){
#   
#   coupons_ytm=function(r_star){
#     dcf=vector(length=length(cv)+1)
#     for(i in 1:length(cv)){
#       dcf[i]=cv*exp(-r_star*i)
#     }
#     dcf[length(dcf)]=-strike
#     sum(dcf)
#   }
#   
#   r_star=uniroot(vasicek_explicit_price,t_=0.25,T_=0.5,par=1000,interval=c(0.1,0.15))$root
#   
#   e_power_minus_r_star=unique(abs(polyroot(c(-strike,cv))))
#   r_star=-log(e_power_minus_r_star)
#   PtT=mapply(vasicek_explicit_price,par=cv,t_=T_,T_=ct,r_t=r_star)
#   
#   PtT=vasicek_explicit_price(T_=T_,t_=t_,r_t=r0,par=par)/par
#   PtS=vasicek_explicit_price(T_=S_,t_=t_,r_t=r0,par=par)
#   
#   sigp=sqrt((1-exp(-2*kappa*(T_-t_)))/(2*kappa))*(1-exp(-kappa*(S_-T_)))*sig/kappa
#   
#   d1=log(PtS/(strike*PtT))/sigp + 0.5*sigp
#   d2=d1-sigp
#   
#   callprice=PtS*pnorm(d1) - strike*PtT*pnorm(d2)
#   callprice
# }

vasicek_option_explicit()
############################################################
############################################################

#2
#a
cir_montecarlo_zcb=function(par,T_bond=1,r0=0.05,kappa=0.92,r_bar=0.055,sig=0.12,paths=10000){
  
  timesteps=ceiling(T_bond*365) #ensuring that the number of timesteps= no of days in that period
  delta=T_bond/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(r_t[,(j-1)])*sqrt(delta)*z[,(j-1)]
  }
  R=delta*rowSums(r_t)
  price_zcb=par*mean(exp(-R))
  
  price_zcb
}
cir_montecarlo_zcb(par=1000)
cir_montecarlo_call=function(strike,par,T_option=0.5,T_bond=1,r0=0.05,kappa=0.92,r_bar=0.055,sig=0.12,paths=1000){
  
  timesteps=ceiling(T_option*365) #evaluating rates from 0 to option expiry
  delta=T_option/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(r_t[,(j-1)])*sqrt(delta)*z[,(j-1)]
  }
  rT=r_t[,ncol(r_t)]
  price_T=vector(length=length(rT))
  for(i in 1:length(rT)){
    price_T[i]=cir_montecarlo_zcb(par=par,r0=rT[i],T_bond=(T_bond-T_option),paths=1000)
  }
  payoff=pmax((price_T-strike),0)
  R=rowSums(r_t)*delta
  disc=exp(-R)*payoff
  price_option=mean(disc)
  price_option
}

price_2a=cir_montecarlo_call(strike=980,par=1000)


#b
cir_explicit_zcb=function(T_,t_,r_t,kappa=0.92,r_bar=0.055,sig=0.12,par){
  h1=sqrt(kappa^2 + 2*sig^2)
  h2=(kappa+h1)/2
  h3=2*kappa*r_bar/sig^2
  
  BtT=(exp(h1*(T_ - t_))-1)/(h2*(exp(h1*(T_ - t_))-1)+h1)
  
  AtT=((h1*exp(h2*(T_ - t_)))/(h2*(exp(h1*(T_ - t_))-1)+h1))^h3
  
  price_explicit=par*AtT*exp(-BtT*r_t)
  c(mean(price_explicit),AtT,BtT)
}

#checking
cir_explicit_zcb(r_t=0.05,t_=0,T_=1,par=1000)

#Implicit Finite Difference Method
cir_implicit=function(r_bar=0.055,kappa=0.92,sig=0.12,strike=980,par=1000,r0=0.05,dr=0.0005,dt=1/365){
  
  r_t=seq(0.1,0,by=-dr)
  dt=1/365
  N=length(r_t)
  M=ceiling(0.5/dt)
  pu=-0.5*dt*(kappa*(r_bar-r_t[-c(1,N)])/dr + sig^2*r_t[-c(1,N)]/dr^2) #just the negative of pu for explicit method
  pd=0.5*dt*(kappa*(r_bar-r_t[-c(1,N)])/dr - sig^2*r_t[-c(1,N)]/dr^2) #just the negative of pd for explicit method
  pm=1 + dt*(r_t[-c(1,N)]*sig^2/dr^2 + r_t[-c(1,N)])
  
  #creating matrix A
  
  pupmpd=cbind(pu,pm,pd)
  
  A=matrix(0,nrow=N,ncol=N)
  
  A[(2:(nrow(A)-1)),1:3]=pupmpd
  
  library(data.table)
  for(i in 2:(nrow(A)-1)){
    A[i,]=shift(A[i,],(i-2),fill=0) #shift the values of A towards the right by i-2 and fill the holes created by 0
  }
  
  # brute force to generate A matrix
  # for(i in 2:(nrow(A)-1)){
  #   for(j in 1:ncol(A)){
  #     if(i==j){
  #       A[i,((j-1):(j+1))]=pupmpd[i-1,]
  #     }
  #   }
  # }
  
  A[1,(1:3)]=c(1,-1,0) #first row
  
  A[nrow(A),((ncol(A)-2):ncol(A))]=c(0,1,-1)#last row
  
  #creating matrix B
  B=matrix(rep(0,N),nrow=N)
  
  bond_prices=sapply(r_t,cir_explicit_zcb,T_=1,t_=0.5,par=1000)[1,]
  payoff=pmax((bond_prices-strike),0)
  price=matrix(nrow=N,ncol=M)
  price[,ncol(price)]=payoff
  B[nrow(B),1]=bond_prices[nrow(B)-1]-bond_prices[nrow(B)-2]
  
  A_inv=solve(A)
  
  for (i in (ncol(price)-1):1){
    B[(2:(nrow(B)-1)),1]=price[(2:(nrow(price)-1)),(i+1)]
    price[,i]=A_inv%*%B
  }
  price[which(r_t==as.character(r0))]
}

price_2b=cir_implicit()

#c


cir_explicit_call=function(strike,T_,t_=0,S_,r0=0.05,kappa=0.92,r_bar=0.055,sig=0.12,par=1000){
  
  #assuming t_=0
  PtT=cir_explicit_zcb(r_t=r0,t_=t_,T_=T_,par=1000)[1]/par
  PtS=cir_explicit_zcb(r_t=r0,t_=t_,T_=S_,par=1000)[1]
  
  ATS=cir_explicit_zcb(r_t=0.05,t_=T_,T_=S_,par=1000)[2] #note: price will not be accurate because r_t should be simulated
  BTS=cir_explicit_zcb(r_t=0.05,t_=T_,T_=S_,par=1000)[3] #but ATS and BTS will be accurate because they don't depend on r_t
  
  theta=sqrt(kappa^2 + 2*sig^2)
  phi=2*theta/(sig^2*(exp(theta*(T_ - t_))-1))
  psi=(kappa+theta)/sig^2
  r_star=log(par*ATS/strike)/BTS
  
  call=PtS*pchisq((2*r_star*(phi+psi+BTS)),df=(4*kappa*r_bar/sig^2),ncp=((2*phi^2*r0*exp(theta*(T_-t_)))/(phi+psi+BTS)))-
    strike*PtT*pchisq((2*r_star*(phi+psi)),df=(4*kappa*r_bar/sig^2),ncp=((2*phi^2*r0*exp(theta*(T_-t_)))/(phi+psi)))
  call
}

price_2c=cir_explicit_call(strike=980,T_=0.5,S_=1)
####################################################################################
####################################################################################


#3


g2_montecarlo_zcb=function(par=1000,T_bond=1,x0,y0,r0=0.03,
                           rho=0.7,a=0.1,b=0.3,sig=0.03,eta=0.08,paths=10000){
  phi0=r0
  phit=phi0 #assuming phit is constant
  
  
  timesteps=ceiling(T_bond*365) #ensuring that the number of timesteps= no of days in that period
  delta=T_bond/timesteps
  
  z1=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate two matrices of random normals
  z2=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps)
  
  xt=yt=matrix(nrow=paths,ncol=timesteps)
  
  #generating xt and yt
  xt[,1]=x0
  yt[,1]=y0
  
  for(j in 2:timesteps){
    xt[,j]=xt[,(j-1)] - a*xt[,(j-1)]*delta + sig*sqrt(delta)*z1[,(j-1)]
    yt[,j]=yt[,(j-1)] - b*yt[,(j-1)]*delta + eta*sqrt(delta)*(rho*z1[,(j-1)] + sqrt(1-rho^2)*z2[,(j-1)])
  }
  
  r_t=xt+yt+phit
  
  R=delta*rowSums(r_t)
  price_zcb=par*mean(exp(-R))
  
  price_zcb
}

temp=g2_montecarlo_zcb(x0=0,y0=0)

g2_montecarlo_put=function(strike,par,T_option=0.5,T_bond=1,x0=0,y0=0,r0=0.03,
                           rho=0.7,a=0.1,b=0.3,sig=0.03,eta=0.08,paths=1000){
  
  phit=phi0=r0
  timesteps=ceiling(T_option*365) #ensuring that the number of timesteps= no of days in that period
  delta=T_option/timesteps
  
  z1=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate two matrices of random normals
  z2=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps)
  
  xt=yt=matrix(nrow=paths,ncol=timesteps)
  
  xt[,1]=x0
  yt[,1]=y0
  
  for(j in 2:timesteps){
    xt[,j]=xt[,(j-1)] - a*xt[,(j-1)]*delta + sig*sqrt(delta)*z1[,(j-1)]
    yt[,j]=yt[,(j-1)] - b*yt[,(j-1)]*delta + eta*sqrt(delta)*(rho*z1[,(j-1)] + sqrt(1-rho^2)*z2[,(j-1)])
  }
  
  r_t=xt+yt+phit
  
  xT=xt[,ncol(xt)]
  yT=yt[,ncol(yt)]
  
  #values of x0 and y0 change but r0 stays same
  price_T=mapply(g2_montecarlo_zcb,x0=xT,y0=yT,T_bond=(T_bond-T_option),paths=100) 
  
  payoff=pmax((strike-price_T),0)
  R=rowSums(r_t)*delta
  disc=exp(-R)*payoff
  price_option=mean(disc)
  price_option
}

price_3_montecarlo=g2_montecarlo_put(strike=950,par=1000)


#explicit

g2_explicit_zcb=function(x0,y0,r0=0.03,rho=0.7,a=0.1,b=0.3,sig=0.03,eta=0.08,T_=1,t_=0.5,par=1000,paths=10000){
  
  phit=phi0=r0
  timesteps=ceiling(T_*365) #ensuring that the number of timesteps= no of days in that period
  delta=T_/timesteps
  z1=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  z2=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps)
  xt=yt=matrix(nrow=paths,ncol=timesteps)
  
  xt[,1]=x0
  yt[,1]=y0
  
  for(j in 2:timesteps){
    xt[,j]=xt[,(j-1)] - a*xt[,(j-1)]*delta + sig*sqrt(delta)*z1[,(j-1)]
    yt[,j]=yt[,(j-1)] - b*yt[,(j-1)]*delta + eta*sqrt(delta)*(rho*z1[,(j-1)] + sqrt(1-rho^2)*z2[,(j-1)])
  }
  term1v=(sig^2/a^2)*(T_ - t_ + (2/a)*exp(-a*(T_ - t_)) - (1/(2*a))*exp(-2*a*(T_ - t_)) - 3/(2*a))
  term2v=(eta^2/b^2)*(T_ - t_ + (2/b)*exp(-b*(T_ - t_)) - (1/(2*b))*exp(-2*b*(T_ - t_)) - 3/(2*b))
  term3v=(2*rho*sig*eta/(a*b))*(T_ - t_ + (exp(-a*(T_ - t_))-1)/a + (exp(-b*(T_ - t_))-1)/b - (exp(-(a+b)*(T_ - t_))-1)/(a+b))
  
  vtT=term1v + term2v + term3v 
  t_location=floor((t_/T_)*timesteps)+1
  PtT=par*exp(-phi0*(T_-t_) - (1-exp(-a*(T_ - t_)))*xt[,t_location]/a - (1-exp(-b*(T_ - t_)))*yt[,t_location]/b + 0.5*vtT)
  mean(PtT)
}

g2_explicit_put=function(strike,phi0=0.03,r0=0.03,rho=0.7,a=0.1,b=0.3,sig=0.03,eta=0.08,T_=0.5,t_=0,S_=1,par=1000){
  
  capital_sigma_sq=(1-exp(-a*(S_-T_)))^2*(1-exp(-2*a*(T_-t_)))*sig^2/(2*a^3) + 
    (1-exp(-b*(S_-T_)))^2*(1-exp(-2*b*(T_-t_)))*eta^2/(2*b^3) +
    (1-exp(-a*(S_-T_)))*(1-exp(-b*(S_-T_)))*(1-exp(-(a+b)*(T_-t_)))*2*rho*sig*eta/(a*b*(a+b))
  capital_sigma=sqrt(capital_sigma_sq)
  
  #ensuring the PtT is for par value 1$ and PtS is for par value given 
  PtT=g2_explicit_zcb(t_=0,T_=0.5,par=par,x0=0,y0=0)/par
  PtS=g2_explicit_zcb(t_=0,T_=1,par=par,x0=0,y0=0)
  
  d1=log(strike*PtT/PtS)/capital_sigma - 0.5*capital_sigma
  d2=d1 + capital_sigma
  put=-PtS*pnorm(d1) + PtT*strike*pnorm(d2)
  put
}
price_3_exp=g2_explicit_put(strike=950)
