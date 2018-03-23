#Nikhil Guruji Project 9


#function which gives the expected interest rate at a given time (in years)
cir_rate=function(T_,r0=0.078,kappa=0.6,r_bar=0.08,sig=0.12,paths=10000){
  timesteps=ceiling(T_*100) #ensuring that the number of timesteps= no of days in that period
  delta=T_/timesteps
  r_t=matrix(nrow=paths,ncol=timesteps) # r_t will be a matrix which will store simulations of r
  z=matrix(rnorm(timesteps*paths),nrow=paths,ncol=timesteps) #generate a matrix of random normals
  
  r_t[,1]=r0
  
  for(j in 2:timesteps){
    r_t[,j]=r_t[,(j-1)] + kappa*delta*(r_bar-r_t[,(j-1)]) + sig*sqrt(abs(r_t[,(j-1)]))*sqrt(delta)*z[,(j-1)]
  }
  
  colMeans(r_t)
}

#generating the rates for the input of 10 year risk free rate in the global environment
r_cir=cir_rate(T_=30)

#CIR explicit function from Project 8
cir_explicit_zcb=function(T_,t_,r_t,kappa=0.6,r_bar=0.08,sig=0.12,par){
  h1=sqrt(kappa^2 + 2*sig^2)
  h2=(kappa+h1)/2
  h3=2*kappa*r_bar/sig^2
  
  BtT=(exp(h1*(T_ - t_))-1)/(h2*(exp(h1*(T_ - t_))-1)+h1)
  
  AtT=((h1*exp(h2*(T_ - t_)))/(h2*(exp(h1*(T_ - t_))-1)+h1))^h3
  
  price_explicit=par*AtT*exp(-BtT*r_t)
  mean(price_explicit)
}

#creating a function which calculates the cashflow at time 
ct=function(pv_last,pv0=100000,r,N=360,t_,month_no,r0=0.078,kappa=0.6,r_bar=0.08,sig=0.12,modeltype="Numerix"){

  if(modeltype=="Numerix"){
    r_start=r_cir[ceiling((t_-1)*100/12)+1]
    #r_start=-log(cir_explicit_zcb(par=1,T_=(t_/12),t_=0,r_t=r0,kappa=kappa,r_bar=r_bar,sig=sig))
    r_10=-0.1*log(cir_explicit_zcb(par=1,T_=((t_/12)+10),t_=(t_/12),r_t=r_start,kappa=kappa,r_bar=r_bar,sig=sig))
  
    RIt=0.28+0.14*atan(-8.57 + 430*((r*12)-r_10))
    BUt=0.3+0.7*(pv_last/pv0)
    SGt=min(1,t_/30)
    SYt=c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98)
  
    CPR=RIt*BUt*SGt*SYt[month_no]
  }else if(modeltype=="PSA"){
    CPR=ifelse(t_<=30,0.002*t_,0.06)
  }else{CPR=0}
  
  MPt=(pv_last*r)/(1-(1+r)^(-N+(t_-1)))
  PPt=(pv_last - pv_last*r*(1/(1-(1+r)^(-N+t_-1)) - 1))*(1-(1-CPR)^(1/12))
  MPt+PPt
}

mbs_price=function(modeltype="Numerix",N=30,wac=0.08,pv0=100000,price_of="MBS",
                   r0=0.078,kappa=0.6,r_bar=0.08,sig=0.12){
  times=seq(1,(12*N))
  
  #starting from January
  month_no=rep(c(1:12),N)
  pv=vector(length=(length(times)+1))
  pv[1]=pv0
  
  #note: 1st element of pv vector corresponds to t=0 and 1st element of pmt and ip vectors corresponds to t=1
  
  ip=vector(length=length(times))
  pmt=vector(length=length(times))
  for(i in 1:length(pmt)){
    pmt[i]=ct(pv_last = pv[i],pv0=pv[1],r=wac/12,
              month_no=month_no[i],t_=times[i],r0=r0,kappa=kappa,r_bar=r_bar,sig=sig,modeltype=modeltype)
    ip[i]=pv[i]*wac/12
    pv[i+1]=pv[i]-(pmt[i]-ip[i])
  }
  
  if(price_of=="MBS"){
    cash_flow=pmt
  }else if(price_of=="PO"){
    cash_flow=pmt-ip
  }else if(price_of=="IO"){
    cash_flow=ip
  }else{cash_flow=0}
  
  #100 paths for the final price
  dcf=matrix(nrow=100,ncol=12*N)
  for(i in 1:nrow(dcf)){
    dcf[i,]=mapply(cir_explicit_zcb,par=cash_flow,T_=(times/12),t_=0,r_t=r0,kappa=kappa,r_bar=r_bar,sig=sig)
  }
  mean(rowSums(dcf))
}

#1

#a
price_numerix=mbs_price(modeltype = "Numerix")

#b
ks=seq(0.3,0.9,by=0.1)
prices_ks=mapply(mbs_price,kappa=ks)
plot(ks,prices_ks,type='o',xlab="Kappa values",ylab="MBS Prices",main="Price for different values of Kappa (Numerix)")

#fancy package to show latex symbols on plots
if(!require(latex2exp)){
  install.packages("latex2exp")
}
library(latex2exp)

#c
r_bars=seq(0.03,0.09,by=0.01)
prices_rs=mapply(mbs_price,r_bar=r_bars)
plot(r_bars,prices_rs,type='o',xlab=TeX("\\bar{r}"),ylab="MBS Prices",
     main=TeX("Price for different values of \\bar{r} (Numerix)"))

#2
#a
price_psa=mbs_price(modeltype = "PSA")

#b
prices_ks_psa=mapply(mbs_price,kappa=ks,modeltype="PSA")
plot(ks,prices_ks_psa,type='o',xlab="Kappa values",ylab="MBS Prices",main="Price for different values of Kappa (PSA)")

#just checking
price_intrinsic=mbs_price(modeltype="N")


#3
oas_dur_conv_calculator=function(modeltype="Numerix",mkt_price=110000,return_qty,
                        N=30,wac=0.08,pv0=100000,r0=0.078,kappa=0.6,r_bar=0.08,sig=0.12){
  times=seq(1,(12*N))
  
  #starting from January
  month_no=rep(c(1:12),N)
  pv=vector(length=(length(times)+1))
  pv[1]=pv0
  
  #note: 1st element of pv vector corresponds to t=0 and 1st element of pmt and ip vectors corresponds to t=1
  
  
  ip=vector(length=length(times))
  pmt=vector(length=length(times))
  for(i in 1:length(pmt)){
    pmt[i]=ct(pv_last = pv[i],pv0=pv[1],r=wac/12,
              month_no=month_no[i],t_=times[i],r0=r0,kappa=kappa,r_bar=r_bar,sig=sig,modeltype=modeltype)
    ip[i]=pv[i]*wac/12
    pv[i+1]=pv[i]-(pmt[i]-ip[i])
  }
  
  oas_func=function(x){
    dcf=mapply(cir_explicit_zcb,par=pmt,T_=(times/12),t_=0,r_t=r0,kappa=kappa,r_bar=r_bar,sig=sig)
    dcf_adjusted=dcf*exp(-x*(times/12))
    sum(dcf_adjusted)-mkt_price
  }
  oas=uniroot(oas_func,interval=c(-1,1),tol=1e-100)$root #should be negative if market price is more than what is observed in Q1
  
  #to calculate duration and convexity
  y=0.0005 #value of perturbation in OAS 
  p0=mkt_price
  p_minus=oas_func(x=(oas-y)) + mkt_price
  p_plus=oas_func(x=(oas+y)) + mkt_price
  
  if(return_qty=="OAS"){
    return(oas)
  }else if(return_qty=="DUR"){
    dur=(p_minus - p_plus)/(2*y*p0)
    return(dur)
  }else if(return_qty=="CONV"){
    conv=(p_plus + p_minus -2*p0)/(2*p0*y^2)
    return(conv)
  }else{return("ERROR")}
  
}
oas_3=oas_dur_conv_calculator(return_qty = "OAS")

#4
dur_4=oas_dur_conv_calculator(return_qty = "DUR")
conv_4=oas_dur_conv_calculator(return_qty = "CONV")


#5

price_po_numerix=mbs_price(price_of = "PO")
price_io_numerix=mbs_price(price_of = "IO")

prices_po_num=mapply(mbs_price,r_bar=r_bars,price_of="PO")

prices_io_num=mapply(mbs_price,r_bar=r_bars,price_of="IO")

plot(r_bars,prices_io_num,type='o',ylab="Prices",xlab=TeX("\\bar{r}"),
     main=TeX("Price for different values of \\bar{r} (Numerix)"),
       lty=5,
     col="blue",ylim=c(20000,120000))
lines(r_bars,prices_po_num,type='o',col="red",lty=5)
lines(r_bars,prices_po_num+prices_io_num,type='o',lty=1,col="violet")
legend(x=0.07,y=90000,legend=c("Price of PO Tranche","Price of IO Tranche","Price of the MBS"),
       col=c("red","blue","violet"),lty=c(5,5,1))
# prices_po_psa=mapply(mbs_price,r_bar=r_bars,price_of="PO",modeltype="PSA")
# plot(r_bars,prices_po_psa,type='o',xlab=TeX("\\bar{r}"),ylab="PO Prices",
#      main=TeX("Price for different values of \\bar{r} (PSA)"))
# 
# prices_io_psa=mapply(mbs_price,r_bar=r_bars,price_of="IO",modeltype="PSA")
# plot(r_bars,prices_io_psa,type='o',xlab=TeX("\\bar{r}"),ylab="IO Prices",
#      main=TeX("Price for different values of \\bar{r} (PSA)"))

