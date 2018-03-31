# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 10:40:15 2018

@author: Nikhil
"""


from math import *
from numpy import *
from matplotlib import *
import matplotlib.pyplot as plt

#1a

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

uni=uniform_generator_01(1234,10000)

#mean and std come from line 10
mean(uni)
std(uni)

#plots come from line 12
plt.hist(uni)
plt.plot(uni,'ro')

#1b
#this function comes from line 10
un=random.uniform(size=10000)

#1c

mean(un)
std(un)
plt.hist(un)

# it is observed that the mean and the standard deviation of output from part (a) and part (b) are very close
# so, the LGM method does a good job at generating random uniformly distributed numbers

#2a

p=[0.30,0.35,0.20,0.15]
xn=[-1,0,1,2]

def gdd_generator(p,xn,n,u):
  x=[77]*n
  for i in arange(n):
    for j in arange(len(p)):
      if x[i]==77:
        if u[i]<sum(p[0:(j+1)]): #taking the sum from 0:0 upto 0:4 because in the end we need sum as 1
          x[i]=xn[j]
  #x=asarray(x)
  #x=x-1
  return(x)

b=gdd_generator(p=p,xn=xn,n=10000,u=uni)

#2b

plt.hist(b)
mean(b)
std(b)

#3a

p=[0.64,0.36]
u=uniform_generator_01(seed=1571,n=44000)
bino=[None]*1000
for i in arange(1000):
  j=44*i
  b=gdd_generator(p=p,xn=[1,0],n=44,u=u[j:(44+j)]) #we input 44 random uniform numbers each time
  bino[i]=sum(b)


#3b

plt.hist(bino)
above_40=[i for i in bino if i>=40] #select all values from bino which are >= 40
p_above_40=len(above_40)/len(bino)
p_above_40
#it can be seen that none of the random variables are more than 40. So, P(x>=40)=0
# theoretically, P(X>=40)= 44C40(0.64)^40*(0.36)^4 + 44C41(0.64)^41*(0.36)^3 + 44C42(0.64)^42*(0.36)^2
# + 44C43(0.64)^43*(0.36)^1 + 44C44(0.64)^44*(0.36)^0

#it is verified that the probability is extremely small

#4a
lambd=1.5
u=[i for i in u if i>0]
e_=-lambd*log(u)

#4b
above_1=[i for i in e_ if i>=1]
above_4=[i for i in e_ if i>=4]
pe_above_1=len(above_1)/len(e_)
pe_above_4=len(above_4)/len(e_)
pe_above_1
pe_above_4

#4c

mean(e_)
std(e_)
plt.hist(e_,color='r')
plt.axvline(e_.mean(),linestyle='dashed')

#5a

u=uniform_generator_01(seed=2222,n=5000)

#5b
def box_muller(n):
  u=uniform_generator_01(seed=2222,n=n)
  z=[None]*n
  i=arange(0,n,2)
  for j in i:
    z[j]=sqrt(-2*log(u[j]))*cos(2*pi*u[j+1])
    z[j+1]=sqrt(-2*log(u[j]))*sin(2*pi*u[j+1])
  return(z)
zbm=box_muller(5000)

#5c
plt.hist(zbm)
plt.axvline(mean(zbm),color="r")
mean(zbm) #should be close to 0
std(zbm) #should be close to 1

#5d
def polar_mars(n):
  u=uniform_generator_01(seed=2222,n=n) #setting a random seed
  z=[None]*n
  v=[None]*n
  i=arange(0,n,2)
  for j in i:
    v[j]=2*u[j] - 1
    v[j+1]=2*u[j+1] - 1
    w=(v[j])**2 + (v[j+1])**2
    if w<=1:
      z[j]=v[j]*sqrt(-2*log(w)/w)
      z[j+1]=v[j+1]*sqrt(-2*log(w)/w)
  z=[k for k in z if k!=None] #removing NA values for which W>1
  return(z)
zpm=polar_mars(int(4*5000/pi)) # taking (4/pi) times the values needed = (4/pi)*5000 = 6367
len(zpm) #should be less than the original size because we have removed missing values

plt.hist(zpm,color='g')
plt.axvline(mean(zpm),color="r")
mean(zpm) #should be close to 0
std(zpm) #should be close to 1

#5f

from time import time
t0 = time()
box_muller(20000)
t1 = time()
polar_mars(int(4*20000/pi))
t2 = time()
bmtime=t1-t0
pmtime=t2-t1
print("function boxmuller takes",bmtime," seconds and function polarmars takes",pmtime," seconds")

#it is observed that the Polar Marsaglia Method is almost as fast as the Box Muller Method when the same
# number of uniform numbers are used. However, to generat the same number of normally distributed numbers,
# it is observed that the Box-Muller test is much faster.
