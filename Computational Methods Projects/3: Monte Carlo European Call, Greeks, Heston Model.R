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

#1
q1=function(seed1=1111111,seed2=2222222,m1=1000,n1=1000,t){

  # simulate x and y: m1 times
  # number of values for each X and Y: n1
  xs1=matrix(nrow=m1,ncol=n1) # xs will be a matrix which will store m1 simulations of x
  ys1=matrix(nrow=m1,ncol=n1) # ys is a matrix which will store n1 simulations of x


  for(xyi in 1:m1){ #outer for loop to simulate m1 X's and Y's 
    #more number of simulations will be slow because there will be n1 simulations more to do with every 1 X and Y more
    d=t/n1 #dividing t into n1 samples (d is the step size)
    z1=box_muller(n1,xyi*seed1) #simulating n1 random normals for Wt with random seed everytime
    z2=box_muller(n1,(xyi*seed2 + xyi)) #simulating n1 random normals for Zt with random seed everytime
    ti=seq(d,t,by=d) #vector of points in time
  
    #simulating Wt for every step in time
    wt=c()
    for(i in 1:n1){
      wt[i]=sqrt(d)*sum(z1[1:i])
    }

    #simulating Zt for every step in time
    zt=c()
    for(i in 1:n1){
      zt[i]=sqrt(d)*sum(z2[1:i])
    }

    #simulating X
    x=c()
    x[1]=1
    for(i in 2:n1){
      x[i]=x[i-1] + ((1/5)-(x[i-1]/2))*d + (2/3)*(wt[i]-wt[i-1])
    }

    #simulating Y
    y=c()
    y[1]=3/4
    for(i in 2:n1){
      y[i]=y[i-1] + ((2/(ti[i-1]+1))*y[i-1] + (1+ti[i-1]**3)/3)*d + ((1+ti[i-1]**3)/3)*(zt[i]-zt[i-1])
    }
    xs1[xyi,]=x
    ys1[xyi,]=y
  }
  list(xs1,ys1)
}

ans1_t2=q1(m1=1000,n1=100,t=2)
Y2=unlist(ans1_t2[[2]])

X2=unlist(ans1_t2[[1]])
E1=X2[,ncol(X2)] #selecting the Xs at t=2 (the last column)
#I discovered a bug in R where cuberoot of a negative number in a vector returns an NaN
# So, I take the cuberoot of the absolute value, and then add the negative or positive sign
E1=((abs(E1))**(1/3))*sign(E1)

ans1_t3=q1(t=3,m1=1000,n1=100)
Y3=unlist(ans1_t3[[2]])

#final answers
prob_y2_greater_than_5=length(which(Y2[,ncol(Y2)]>5))/nrow(Y2) # rows which are above 5/ total rows
E1_1=mean(E1) 
E1_2=mean(Y3[,ncol(Y3)]) # last column corresponds to t=3
E1_3=mean(X2[which(X2[,ncol(X2)]>1),ncol(X2)]*Y2[which(X2[,ncol(X2)]>1),ncol(Y2)]) # last columns correspond to t=2


#2
q2=function(seed1=1111111,seed2=2222222,m2=1000,n2=1000,t){
    # simulate x and y: m2 times
  # number of values for each X and Y: n2
  xs2=matrix(nrow=m2,ncol=n2) # xs will be a matrix which will store m2 simulations of x
  ys2=matrix(nrow=m2,ncol=n2) # ys is a matrix which will store m2 simulations of x


  for(xyi in 1:m2){ #outer for loop to simulate m2 X's and Y's 
    #more number of simulations will be slow because there will be n2 simulations more to do with every 1 X and Y more
  
    d=t/n2 #dividing t into n2 samples (step size of d)
    z1=box_muller(n2,xyi*seed1) #simulating n2 random normals for Wt with random seed everytime
    z2=box_muller(n2,(xyi*seed2 + xyi)) #simulating n2 random normals for Zt with random seed everytime
    ti=seq(d,t,by=d) #vector of points in time from upto t
    

    #simulating Wt for every step in time
    wt=c()
    for(i in 1:n2){
      wt[i]=sqrt(d)*sum(z1[1:i])
    }
  
    #simulating Zt for every step in time
    zt=c()
    for(i in 1:n2){
      zt[i]=sqrt(d)*sum(z2[1:i])
    }
  
    #simulating X
    x=c()
    x[1]=1
    for(i in 2:n2){
      x[i]=x[i-1] + ((1/4)*x[i-1])*d + (1/3)*x[i-1]*(wt[i]-wt[i-1]) - (3/4)*x[i-1]*(zt[i]-zt[i-1])
    }
  
    #simulating Y
    y=c()
    y[1]=1
    for(i in 2:n2){
      y[i]=y[i-1] + (-0.08 + (1/18) + (9/32))*y[i-1]*d + (1/3)*y[i-1]*(wt[i]-wt[i-1]) + (3/4)*y[i-1]*(zt[i]-zt[i-1])
    }
    xs2[xyi,]=x
    ys2[xyi,]=y
  }
  list(xs2,ys2)  
}

ans2=q2(t=3,m2=1000,n2=100)
xs2=unlist(ans2[[1]])
ys2=unlist(ans2[[2]])

E21=1+xs2[,ncol(xs2)]
E21=((abs(E21))**(1/3))*sign(E21) #same as Q1 E1_1

E22=1+ys2[,ncol(ys2)]
E22=((abs(E22))**(1/3))*sign(E22)

#final answers

E2_1=mean(E21)
E2_2=mean(E22)


#3
#test values for five parameters : r,sigma, T, S0, K
s0=15
t=0.5
k=20
s=0.25
r=0.04

#a
monte_carlo_eu_call_pricing=function(seed=1234,s0,t=0.5,k=20,r=0.04,s=0.25){
  
  z1=box_muller(1000,seed)
  z2=box_muller(1000,seed=1234*seed)

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
}

C1=monte_carlo_eu_call_pricing(s0=s0,t=t,k=k,r=r,s=s)

#b

#numerical computation of N(.)
N=function(x){
  d1=0.0498673470
  d2=0.0211410061
  d3=0.0032776263
  d4=0.0000380036
  d5=0.0000488906
  d6=0.0000053830
  n=ifelse(x<0 # if true, then 1-N(-x)
           ,(1-(1-0.5*(1 + d1*(-x) + d2*(-x)**2 + d3*(-x)**3 + d4*(-x)**4 + d5*(-x)**5 + d6*(-x)**6)**(-16))) 
           ,(1-0.5*(1 + d1*x + d2*x**2 + d3*x**3 + d4*x**4 + d5*x**5 + d6*x**6)**(-16)))
  n
}

black_scholes_eu_call_pricing=function(seed=1234,s0,t=0.5,k=20,r=0.04,s=0.25){

  #using the Black Scholes formula
  d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
  d2= d1 - s*sqrt(t)

  c_bs=s0*N(d1) - k*exp(-r*t)*N(d2)
  c_bs
}

C2=black_scholes_eu_call_pricing(s0=s0,r=r,s=s,k=k,t=t) #verifying from previous assignment

#c 

#can be done using numerical methods, but there will be different answers for different step sizes
# so, estimation is done using formulae


small_n=function(t){ #defining function for n(.) (N(.) was defined earlier)
  (1/sqrt(2*pi))*exp(-0.5*t**2)
}

delta_eu_call=function(s0,t=0.5,k=20,r=0.04,s=0.25){
  d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
  N(d1)
}

gamma_eu_call=function(s0,t=0.5,k=20,r=0.04,s=0.25){
    d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
    small_n(d1)/(s0*s*sqrt(t))
}

theta_eu_call=function(s0,t=0.5,k=20,r=0.04,s=0.25){
    d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
    d2= d1 - s*sqrt(t)
    0-s0*s*small_n(d1)/(2*sqrt(t)) - r*k*exp(-r*t)*N(d2)
}

vega_eu_call=function(s0,t=0.5,k=20,r=0.04,s=0.25){
  d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
  s0*sqrt(t)*small_n(d1)
}

rho_eu_call=function(s0,t=0.5,k=20,r=0.04,s=0.25){
  d1= (log(s0/k)+(r+0.5*s^2)*t)/(s*sqrt(t))
  d2= d1 - s*sqrt(t)
  k*t*exp(-r*t)*N(d2)
}

s0=seq(15,25)
D=delta_eu_call(s0)
Th=theta_eu_call(s0)
V=vega_eu_call(s0)
G=gamma_eu_call(s0)
R=rho_eu_call(s0)


plot(s0,D,type='l',ylab="Delta")
plot(s0,G,type='l',ylab="Gamma")
plot(s0,Th,type='l',ylab="Theta")
plot(s0,V,type='l',ylab="Vega")
plot(s0,R,type='l',ylab="Rho")

#4

rho=-0.6
r=0.03
s0=48
v0=0.05
s=0.42 #not the volatility of stock but the volatility of variance process
alpha=5.8
beta=0.0625
seed1=111111

q4_full_trunc=function(seed1=111111,m1=1000,n1=1000,t){
  
  # simulate S and V: m1 times
  # number of values for each S and V: n1
  stock_paths=matrix(nrow=m1,ncol=n1) # stock_paths will be a matrix which will store n1 simulations of S
  
  
  for(si in 1:m1){ #outer for loop to simulate m1 S's and V's 
    #more number of simulations will be slow because there will be n1 simulations more to do with every 1 S and V more
    d=t/n1 #dividing t into n1 samples (d is the step size)
    z=box_muller(n1*2,si*seed1) #simulating 2*n1 random normals with random seed everytime
    z1=z[1:n1] # dividing z into two parts of length n1 each
    z2=z[(n1+1):(2*n1)]
    
    #simulating S and V
    st=c()
    st[1]=s0
    vt=c()
    vt[1]=v0
    
    for(i in 2:n1){
      #considering the max(Vt,0) for sqrt terms as well as the second term of Vt (full truncation method)
      st[i]=st[i-1] + r*st[i-1]*d + sqrt(max(vt[i-1],0)*d)*st[i-1]*z1[i]
      vt[i]=vt[i-1] + alpha*(beta-max(vt[i-1],0))*d + s*sqrt(d*max(vt[i-1],0))*(rho*z1[i] + sqrt(1-rho**2)*z2[i])
    }
    
    stock_paths[si,]=st
  }
  stock_paths
}

q4_part_trunc=function(seed1=1111111,m1=1000,n1=1000,t){
  
  # simulate S and V: m1 times
  # number of values for each S and V: n1
  stock_paths=matrix(nrow=m1,ncol=n1) # stock_paths will be a matrix which will store n1 simulations of S
  
  
  for(si in 1:m1){ #outer for loop to simulate m1 S's and V's 
    #more number of simulations will be slow because there will be n1 simulations more to do with every 1 S and V more
    d=t/n1 #dividing t into n1 samples (d is the step size)
    z=box_muller(n1*2,si*seed1) #simulating 2*n1 random normals with random seed everytime
    z1=z[1:n1] # dividing z into two parts of length n1 each
    z2=z[(n1+1):(2*n1)]
    
    #simulating S and V
    st=c()
    st[1]=s0
    vt=c()
    vt[1]=v0
    
    for(i in 2:n1){
      #considering the max(Vt,0) only for sqrt terms (partial truncation method)
      st[i]=st[i-1] + r*st[i-1]*d + sqrt(max(vt[i-1],0)*d)*st[i-1]*z1[i]
      vt[i]=vt[i-1] + alpha*(beta-vt[i-1])*d + s*sqrt(d*max(vt[i-1],0))*(rho*z1[i] + sqrt(1-rho**2)*z2[i])
    }
    
    stock_paths[si,]=st
  }
  stock_paths
}

q4_ref=function(seed1=111111,m1=1000,n1=1000,t){
  
  # simulate S and V: m1 times
  # number of values for each S and V: n1
  stock_paths=matrix(nrow=m1,ncol=n1) # stock_paths will be a matrix which will store n1 simulations of S
  
  
  for(si in 1:m1){ #outer for loop to simulate m1 S's and V's 
    #more number of simulations will be slow because there will be n1 simulations more to do with every 1 S and V more
    d=t/n1 #dividing t into n1 samples (d is the step size)
    z=box_muller(n1*2,si*seed1) #simulating 2*n1 random normals with random seed everytime
    z1=z[1:n1] # dividing z into two parts of length n1 each
    z2=z[(n1+1):(2*n1)]
    
    #simulating S and V
    st=c()
    st[1]=s0
    vt=c()
    vt[1]=v0
    
    for(i in 2:n1){
      #considering the absolute values of Vt (reflection method)
      st[i]=st[i-1] + r*st[i-1]*d + sqrt(abs(vt[i-1])*d)*st[i-1]*z1[i]
      vt[i]=abs(vt[i-1]) + alpha*(beta-abs(vt[i-1]))*d + s*sqrt(d*abs(vt[i-1]))*(rho*z1[i] + sqrt(1-rho**2)*z2[i])
    }
    
    stock_paths[si,]=st
  }
  stock_paths
}

#assuming strike = 50$ and time to maturity = 3 periods

k=50
t=3

#taking number of paths as the square (in terms of order) of number of simulations for each path
Sft=q4_full_trunc(t=t,m1=10000,n1=100) 
Spt=q4_part_trunc(t=t,m1=10000,n1=100)
Sref=q4_ref(t=t,m1=10000,n1=100)

#final stock price values (just for information)
ST_ft=mean(Sft[,ncol(Sft)])
ST_pt=mean(Spt[,ncol(Spt)])
ST_ref=mean(Sref[,ncol(Sref)])

#plotting and comparing
par(mfrow=c(3,1))
plot(Sft[241,],type='l',main="Full truncation")
plot(Spt[241,],type='l',main="Partial truncation")
plot(Sref[241,],type='l',main="Reflection")

#final answers
payoff_ft=Sft[,ncol(Sft)]-k #the last values correspond to S_T where T=3
payoff_ft[which(payoff_ft<0)]=0
cft=exp(-r*t)*mean(payoff_ft)

payoff_pt=Spt[,ncol(Spt)]-k
payoff_pt[which(payoff_pt<0)]=0
cpt=exp(-r*t)*mean(payoff_pt)

payoff_ref=Sref[,ncol(Sref)]-k
payoff_ref[which(payoff_ref<0)]=0
cref=exp(-r*t)*mean(payoff_ref)


#5

#a
seed=1111111
#generating 200 uniform numbers and separating them into two vectors
u=uniform_generator_01(seed=seed,n=200)
u1=u[1:100]
u2=u[101:200]
u2x2=cbind(u1,u2)

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

#b
p2=H(100,2) #generating 100 Halton numbers with base 2
p7=H(100,7)

pb2x2=cbind(p2,p7)

#c
p2=H(100,2)
p4=H(100,4)

pc2x2=cbind(p2,p4)

#d

par(mfrow=c(3,1))
plot(u1,u2,xlab="Random 100 Uniforms",ylab="Random 100 Uniforms")
plot(p2,p7,xlab="Random 100 Halton numbers (base = 2)",ylab="Random 100 Halton numbers (base = 7)")
plot(p2,p4,xlab="Random 100 Halton numbers (base = 2)",ylab="Random 100 Halton numbers (base = 4)")

#e

f=function(x,y){
  cos_cube_root=cos(2*pi*y)
  cos_cube_root=((abs(cos_cube_root))**(1/3))*sign(cos_cube_root)
  exp(-x*y)*(sin(6*pi*x)+cos_cube_root)
}

h2=H(10000,2)
h4=H(10000,4)
h5=H(10000,5)
h7=H(10000,7)


base_2_4=f(h2,h4)
I1=mean(base_2_4) #using average of 10000 samples

base_2_7=f(h2,h7)
I2=mean(base_2_7)

base_5_7=f(h5,h7)
I3=mean(base_5_7)

#It is observed that I1 is different from I2 and I3 which are very close. This is because 4 is not a prime base
#So, I1 is most likely an incorrect estimate of the double integral
