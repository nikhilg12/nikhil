# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 20:57:52 2018

@author: Nikhil
"""

from math import * #for functions like sqrt, log, exp, etc
from numpy import * #for functions like arange, asarray, etc
from matplotlib import * #for plots
from scipy import * #for normal density functions
from scipy.stats import *
import matplotlib.pyplot as plt

#From project 1

def uniform_generator_01(seed,n):
    a=7**5
    b=0
    m=2**31 -1

    x=[0]*n #declare an array of size n
    u=[0]*n
    x[0]=seed
    for i in arange(len(x)-1): #arange(len(x)) will give i=0,1,2,...,len(x)-1
      x[i+1]=(a*x[i]+b)%m
      u[i]=x[i]/m
    return(u)
   

def box_muller(n,seed):
  u=uniform_generator_01(seed=seed,n=n)
  z=[None]*n
  i=arange(0,n,2)
  for j in i:
    z[j]=sqrt(-2*log(u[j]))*cos(2*pi*u[j+1])
    z[j+1]=sqrt(-2*log(u[j]))*sin(2*pi*u[j+1])
  return(z)
  
  
#1
  
def bivar_normal(seed,n=1000,a=-0.7):

    z=box_muller(2*n,seed)
    z1=z[0:n] #takes the values from 0 to n-1 (not nth value)
    z2=z[n:2*n]
    x=z1 # x= mu1 + sigma1*z1, mu1=0, sigma1=1
    
    varcov=[[1,a],[a,1]]
    q=a #from lecture notes
    p=sqrt(1-q**2) #from lecture notes
    
    pz2=[p*l for l in z2] #multiplying each element in z2 by p
    qz1=[q*k for k in z1] #multiplying each element in z1 by q
    qz1=asarray(qz1) #converting into numpy array so + gives element wise addition and not concatenation
    pz2=asarray(pz2)
    y=pz2+qz1 # y = mu2 + p*z2 + q*z1
    rho_numerator=(1/(n-1))*sum((x-mean(x))*(y-mean(y)))
    rho_denominator=(1/(n-1))*sqrt(sum((x-mean(x))**2)*sum((y-mean(y))**2))
    rho=rho_numerator/rho_denominator
    r=[rho,x,y]
    return(r) #returning a list of correlation, X, and Y
    
out1=bivar_normal(seed=1234)
out1[0] #the first element is the correlation

#2


out2=bivar_normal(seed=123451409,a=0.6)
X=asarray(out2[1]) # first value is the correlation 2nd to 1001st is the vector
Y=asarray(out2[2])

E2=X**3 + sin(Y) + (X**2) * (Y)
E2=maximum(E2,0) #taking the max of 0 and the above expression
mean(E2) #expectation

#3
#a
z=asarray(box_muller(n=10000,seed=1123))
w5=sqrt(5)*z # Wt= sqrt(t)*Z
Ea1=w5**2 +sin(w5)
mean(Ea1) #estimate
(std(Ea1)**2)/len(Ea1) #variance

def E32(t):
     t=asarray(t)
     wt=sqrt(t)*z
     p=exp(0.5*t)*cos(wt)
     E2m=mean(p)
     E2v=(std(p)**2)/len(p)
     return([E2m,E2v])

Ea2=E32(0.5)[0]
Ea3=E32(3.2)[0]
Ea4=E32(6.5)[0]

Ea2v=E32(0.5)[1]
Ea3v=E32(3.2)[1]
Ea4v=E32(6.5)[1]

#c
z1=asarray(random.standard_normal(10000)) # or we could use the box_muller function with different seeds for z1 and z2
z2=asarray(random.standard_normal(10000))
w5_vr1=sqrt(5)*z1
w5_vr2=sqrt(5)*(-z2)

Ec1_vr1=(w5_vr1**2 +sin(w5_vr1))
Ec1_vr2=(w5_vr2**2 +sin(w5_vr2))

Ec1_vr=0.5*(Ec1_vr1+Ec1_vr2)
mean(Ec1_vr)
(std(Ec1_vr)**2)/len(Ec1_vr)

def E32_vr(t):
      wt1=sqrt(t)*z1
      wt2=sqrt(t)*(-z2)
      E21 = exp(0.5*t)*cos(wt1)
      E22 = exp(0.5*t)*cos(wt2)
      E3c=0.5*(E21+E22)
      E3cm=mean(E3c)
      E3cv=(std(E3c)**2)/len(E3c)
      return([E3cm,E3cv])

Ec2_vrm=E32_vr(0.5)[0]
Ec3_vrm=E32_vr(3.2)[0]
Ec4_vrm=E32_vr(6.5)[0]

Ec2_vrv=E32_vr(0.5)[1]
Ec3_vrv=E32_vr(3.2)[1]
Ec4_vrv=E32_vr(6.5)[1]

#4
#4

#a

r=0.04
t=5
s=0.2
s0=88
k=100

z=asarray(box_muller(1000,212147))

st=s0*exp((r-0.5*s**2)*t + s*sqrt(t)*z)
payoff=maximum(st-k,0) # taking only positive values of st-k
c_mc=exp(-r*t)*(payoff)
mean(c_mc) #value of call option through monte carlo
(std(c_mc)**2)/len(c_mc) #variance of call opton through monte carlo

#b

#using the Black Scholes formula
d1= (log(s0/k)+(r+0.5*s**2)*t)/(s*sqrt(t))
d2= d1 - s*sqrt(t)

c_bs=s0*norm.cdf(d1) - k*exp(-r*t)*norm.cdf(d2) #norm.cdf is a function from scipy.stats
c_bs

#using antithetic variance reduction

st_vr1=s0*exp((r-0.5*(s**2))*t + s*sqrt(t)*z1)
payoff_vr1=maximum(st_vr1-k,0) # taking only positive values of st-k
c_vr1=exp(-r*t)*payoff_vr1

st_vr2=s0*exp((r-0.5*(s**2))*t + s*sqrt(t)*(-z1))
payoff_vr2=maximum(st_vr2-k,0) # taking only positive values of st-k
c_vr2=exp(-r*t)*payoff_vr2

c_vr=0.5*(c_vr1 + c_vr2)
mean(c_vr)
std(c_vr)**2/len(c_vr) #variance of the price of call option through antithetic variance reduction

#5

#a
s0=88
r=0.04
s=0.18
n_values=range(1,11) #11 won't come
z=asarray(box_muller(1000,21247))

Esn=map(lambda t: mean(s0*exp((r-0.5*s**2)*t + s*sqrt(t)*(z))),n_values)
#This will apply the function mean(s0...) to each element of n_values
Esn=list(Esn) #convert the map to list where items can be accessed
Esn_new=[s0] + Esn

plt.plot(Esn_new,'go--')
plt.xlabel('n')
plt.ylabel('Expected value of Sn')

#b
paths=range(1,7)
n_values=arange(1,11,0.01)
stockpaths=[0]*6
for i in paths:
    z=asarray(random.standard_normal(1000)) #generate a vector of normals
    w=sqrt(0.01)*cumsum(z) #take cumsum because Wt = (W1-W0) + (W2-W1) + ... + (Wt - W_(t-1))
    #now, apply the function on each element of n_values and w together. Note that len(n_values)=len(w)
    stockpaths[i-1]=list(map(lambda t,brown: s0*exp((r-0.5*s**2)*t + s*brown),n_values,w))
    
#run the following lines together
    
plt.plot(n_values,stockpaths[0],'r')
plt.plot(n_values,stockpaths[1],'g')
plt.plot(n_values,stockpaths[2],'b')
plt.plot(n_values,stockpaths[3],'c')
plt.plot(n_values,stockpaths[4],'m')
plt.plot(n_values,stockpaths[5],'y')
plt.xlabel('time')
plt.ylabel('Stock values')