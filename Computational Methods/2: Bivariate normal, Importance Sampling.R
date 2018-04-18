#Nikhil Guruji Cohort 1 MFE 2018


#uniform number generator from project 1

uniform_generator_01=function(seed=1234,n=10000){
  a=7^5
  b=0
  m=2^31 -1
  
  x=rep(NA,n+100) # generating extra 100 values so that we can remove any trends by removing first 100 values
  u=rep(NA,n+100)
  x[1]=seed
  for (i in 1:(n+100)) {
    x[i+1]=(a*x[i]+b)%%m
    u[i]=x[i]/m
  }
  
  uf=u[101:length(u)]#selecting 101th to 10100th entries to have a truly random sample
  uf
}

# box muller method from project 1

box_muller=function(n,seed){
  u=uniform_generator_01(seed=seed,n=n)
  z=rep(NA,n)
  for(i in 0:((n/2)-1)){
    j=2*i+1
    z[j]=sqrt(-2*log(u[j]))*cos(2*pi*u[j+1])
    z[j+1]=sqrt(-2*log(u[j]))*sin(2*pi*u[j+1])
  }
  z
}
bivar_normal=function(seed,n=1000,a=-0.7){

z=box_muller(2*n,seed)
z1=z[1:n]
z2=z[(n+1):(2*n)]
x=z1 # x= mu1 + sigma1*z1, mu1=0, sigma1=1

varcov=matrix(c(1,a,a,1),nrow=2)
q=a #from lecture notes
p=sqrt(1-q**2) #from lecture notes

y=p*z2+q*z1 # y = mu2 + p*z2 + q*z1
rho_numerator=(1/(n-1))*sum((x-mean(x))*(y-mean(y)))
rho_denominator=(1/(n-1))*sqrt(sum((x-mean(x))**2)*sum((y-mean(y))**2))
rho=rho_numerator/rho_denominator
r=list(rho,x,y)
r #returning a list of correlation, X, and Y
}

output1=bivar_normal(seed=29187089603456) #random seed input
correl=unlist(output1)[1]
correl

#2
output2=bivar_normal(seed=12345,a=0.6)
X=unlist(output2)[2:1001] # first value is the correlation 2nd to 1001st is the vector
Y=unlist(output2)[1002:2001]

E2=X**3 + sin(Y) + (X**2) * (Y)
E2[which(E2<0)]=0 #taking the max of 0 and the above expression
mean(E2) #expectation

#3

#a

z=box_muller(n=10000,seed=1123)
w5=sqrt(5)*z # Wt= sqrt(t)*Z
Ea1=w5**2 +sin(w5)
mean(Ea1) #estimate
(sd(Ea1)**2)/length(Ea1) #variance

E32=function(t){
 wt=sqrt(t)*z
 p=exp(0.5*t)*cos(wt)
 E2m=mean(p)
 E2v=(sd(p)**2)/length(p)
 c(E2m,E2v)
}
Ea2=E32(0.5)[1]
Ea3=E32(3.2)[1]
Ea4=E32(6.5)[1]

Ea2v=E32(0.5)[2]
Ea3v=E32(3.2)[2]
Ea4v=E32(6.5)[2]
#b

#It is observed that the three integrals are close to zero but the variance is higher as t is 

#c
z1=rnorm(10000) # or we could use the box_muller function with different seeds for z1 and z2
z2=rnorm(10000)
w5_vr1=sqrt(5)*z1
w5_vr2=sqrt(5)*(-z2)

Ec1_vr1=(w5_vr1**2 +sin(w5_vr1))
Ec1_vr2=(w5_vr2**2 +sin(w5_vr2))

Ec1_vr=0.5*(Ec1_vr1+Ec1_vr2)
mean(Ec1_vr)
(sd(Ec1_vr)**2)/length(Ec1_vr)

E32_vr=function(t){
  t=6.5
  wt1=sqrt(t)*z1
  wt2=sqrt(t)*(-z2)
  E21 = exp(0.5*t)*cos(wt1)
  E22 = exp(0.5*t)*cos(wt2)
  E3c=0.5*(E21+E22)
  E3cm=mean(E3c)
  E3cv=(sd(E3c)**2)/length(E3c)
  c(E3cm,E3cv)
}

Ec2_vrm=E32_vr(0.5)[1]
Ec3_vrm=E32_vr(3.2)[1]
Ec4_vrm=E32_vr(6.5)[1]

Ec2_vrv=E32_vr(0.5)[2]
Ec3_vrv=E32_vr(3.2)[2]
Ec4_vrv=E32_vr(6.5)[2]

#4

#a

r=0.04
t=5
s=0.2
s0=88
k=100

z=box_muller(1000,212147)

st=s0*exp((r-0.5*s^2)*t + s*sqrt(t)*z)
payoff=st-k
payoff[which(payoff<0)]=0 # taking only positive values of st-k
c_mc=exp(-r*t)*(payoff)
mean(c_mc) #value of call option through monte carlo
(sd(c_mc)**2)/length(c_mc) #variance of call opton through monte carlo

#b

#using the Black Scholes formula
d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
d2= d1 - s*sqrt(t)

c_bs=s0*pnorm(d1) - k*exp(-r*t)*pnorm(d2)
c_bs

st_vr1=s0*exp((r-0.5*(s^2))*t + s*sqrt(t)*z1)
payoff_vr1=st_vr1-k
payoff_vr1[which(payoff_vr1<0)]=0 # taking only positive values of st-k
c_vr1=exp(-r*t)*payoff_vr1

st_vr2=s0*exp((r-0.5*(s^2))*t + s*sqrt(t)*(-z2))
payoff_vr2=st_vr2-k
payoff_vr2[which(payoff_vr2<0)]=0 # taking only positive values of st-k
c_vr2=exp(-r*t)*payoff_vr2

c_vr=0.5*(c_vr1 + c_vr2)
mean(c_vr)
sd(c_vr)**2/length(c_vr) #variance of the price of call option through antithetic variance reduction

#5

#a
r=0.04 
s=0.18 # I will come back here and change this value to 0.35 to plot the graph and rechange it back to 0.18
s0=88
st=as.list(seq(from=1,to=9,by=1))# time points are defined first (1 to 9 because later I will append the S0 value)
z=box_muller(1000,22222) # 1000 random normal values for every time point
s_calc=function(t){
  s0*exp((r-0.5*s^2)*t + s*sqrt(t)*(z)) #1000 stock values at every point in time corresponding to 1000 normal values
}

stock_prices_1000=lapply(st,s_calc) # applying the function which calculates st given t
esn=c(s0,sapply(stock_prices_1000, mean)) # taking mean of the 1000 stock values for each point of time
plot(esn,xlab="Time steps (x1) ",type='o',ylab="Stock price",main="Plot for Q5(a)")

#simulated in steps of 0.01 for plotting together with the part(b) for part (c)
st1000=as.list(seq(from=0.01,to=10,by=0.01))
stock_prices_1000=lapply(st1000,s_calc)
esn1000=c(s0,sapply(stock_prices_1000, mean))

#b and c
six_esn=list() 
for_plots_max=s0 #this variable will be used for plot limits
for_plots_min=s0 #this variable will be used for plot limits so that all the points are captured
for(j in 1:6){
  st=c(s0) #every stock price path starts from s0
  delta=0.01
  point_in_time=seq(from=0.01,to=10,by=delta)  #time values
  z=rnorm(length(point_in_time))
  for(i in 1:length(point_in_time)){
    w=sqrt(delta)*(sum(z[1:i])) # Wt = (W1-W0) + (W2-W0) + ... + (Wt - W_(t-1))
    st[i+1]=s0*exp((r-0.5*s^2)*point_in_time[i] + s*w)
  }
  for_plots_min=c(for_plots_min,min(st)) #storing the minimum stock price for each path
  for_plots_max=c(for_plots_max,max(st)) #storing the maximum stock price for each path
  six_esn[[j]]=st
}

plotter=function(){
  for(i in 1:6){
    if(i==1){
      plot(six_esn[[1]],type='l',xlab="Time steps (x0.01)", col="blue", main="Plot for Q5 (c)",
           ylab="Stock price",xlim=c(0,1000),ylim=c(ceiling(min(for_plots_min)),ceiling(max(for_plots_max))))
    } # if it's the first plot, then use plot command to plot
    else{ # else use lines command to add the plots to existing plot with random colours
      lines(six_esn[[i]],col=rgb(runif(5),runif(5),runif(5)) )
    }
  }
  lines(esn1000,col=rgb(runif(5),runif(5),runif(5)))
}

plotter()

#d

# repeat the part (b) and (c) for s=0.35

#6

#a

f6a=function(x){
  sqrt(1-x**2)
} # the function to be integrated
trap6a=function(n){
  ans=0
  for(i in 1:(n-1)){
    ans=ans+f6a(i/n)
  }
  ans=ans+0.5*(f6a(0)+f6a(1))
  4*ans/n
} # applying trapezoidal rule

Ia=trap6a(1000) #taking large value of n
var(Ia)
actual=4*unlist(integrate(f6a,0,1))$value #actual value from integration (just for reference)

#b

u=uniform_generator_01(n=1000,seed=100*rnorm(1))
f=f6a(u)
Ib=4*mean(f)
var(4*f)

#c
a=0.74
t_x=function(x){
  if(x<=1 && x>=0){
    (1-(a*x**2))/(1-(a/3)) # from lecture notes
  }
  else{0}
} 

var_c= (1-(a/3))*((1/a)-((1-a)/a^(3/2))*atanh(a^(0.5))) - (pi/4)^2 #variance of the estimate from lecture notes
var_c

gf=f6a(u) ##should have drawn from the t function distribution and not uniform
h=t_x(u)

Ic=4*mean(gf/h) #estimate of pi through importance sampling
Ic_var=(sd(gf/h)**2)/length(gf) #variance of the estimate
Ic_var
var(gf/h) # variance using the var function
