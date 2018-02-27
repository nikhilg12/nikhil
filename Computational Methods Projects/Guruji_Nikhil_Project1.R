#Nikhil Guruji - MFE 2018 Cohort 1


#1a
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

u=uniform_generator_01()
mean(u)
sd(u)
hist(u)

#1b

un=runif(10000)

#1c

mean(un)
sd(un)
hist(un)

# it is observed that the mean and the standard deviation of output from part (a) and part (b) are very close
# so, the LGM method does a good job at generating random uniformly distributed numbers

#2a

p=c(0.30,0.35,0.20,0.15)
xn=c(-1,0,1,2)

gdd_generator=function(p,xn,n,u){
  x=rep(NA,n)
  for(i in 1:n){
    for(j in 1:length(p)){
      if(is.na(x[i])){
        if(u[i]<sum(p[1:j])){
          x[i]=xn[j]
        }
      }
    }
  }
  x
}

b=gdd_generator(p,xn,10000,u)
#2b

hist(b,labels = TRUE)
abline(v=mean(b),col="red")
mean(b)
sd(b)

#3a

p=c(0.64,0.36)
u=uniform_generator_01(seed=1571,n=44000)
bino=rep(NA,1000)
for (i in 0:999){
  j=44*i+1
  b=gdd_generator(p=p,xn=c(1,0),n=44,u=u[j:(43+j)]) #we input 44 random uniform numbers each time
  bino[i+1]=sum(b)
}

#3b

hist(bino)
p_above_40=length(bino[bino>=40])/length(bino)
p_above_40
#it can be seen that none of the random variables are more than 40. So, P(x>=40)=0
# theoretically, P(X>=40)= 44C40(0.64)^40*(0.36)^4 + 44C41(0.64)^41*(0.36)^3 + 44C42(0.64)^42*(0.36)^2
# + 44C43(0.64)^43*(0.36)^1 + 44C44(0.64)^44*(0.36)^0

p_actual_above_40= ncol(combn(44,40))*(0.64)^40*(0.36)^4 + ncol(combn(44,41))*(0.64)^41*(0.36)^3 + ncol(combn(44,42))*(0.64)^42*(0.36)^2 + ncol(combn(44,43))*(0.64)^43*(0.36)^1 + ncol(combn(44,44))*(0.64)^44*(0.36)^0
p_actual_above_40

#it is verified that the probability is extremely small

#4a
u=uniform_generator_01()
lambda=1.5
e=-lambda*log(u)

#4b
pe_above_1=length(e[e>=1])/length(e)
pe_above_4=length(e[e>=4])/length(e)
pe_above_1
pe_above_4

#4c

mean(e)
sd(e)
hist(e)
abline(v=mean(e),col="red")

#5a

u=uniform_generator_01(seed=2222,n=5000)

#5b
box_muller=function(n){
  u=uniform_generator_01(seed=2222,n=n)
  z=rep(NA,n)
  for(i in 0:((n/2)-1)){
    j=2*i+1
    z[j]=sqrt(-2*log(u[j]))*cos(2*pi*u[j+1])
    z[j+1]=sqrt(-2*log(u[j]))*sin(2*pi*u[j+1])
  }
  z
}
zbm=box_muller(5000)

#5c
hist(zbm)
abline(v=mean(zbm),col="red")
mean(zbm) #should be close to 0
sd(zbm) #should be close to 1

#5d
polar_mars=function(n){
  u=uniform_generator_01(seed=2222,n=n) #setting a random seed
  z=rep(NA,n)
  v=rep(NA,n)
  for(i in 0:((n/2)-1)){
    j=2*i+1
    v[j]=2*u[j] - 1
    v[j+1]=2*u[j+1] - 1
    w=(v[j])^2 + (v[j+1])^2
  
    if(w<=1){
      z[j]=v[j]*sqrt(-2*log(w)/w)
      z[j+1]=v[j+1]*sqrt(-2*log(w)/w)
    }
    else{}
  }
  z=z[!is.na(z)]#removing NA values for which W>1
  z
}
zpm=polar_mars(ceiling(4*5000/pi)) # taking (4/pi) times the values needed = (4/pi)*5000 = 6367
length(zpm) #should be less than the original size because we have removed missing values

hist(zpm)
abline(v=mean(zpm),col="red")
mean(zpm) #should be close to 0
sd(zpm) #should be close to 1

#5f
library(microbenchmark)
microbenchmark(box_muller(5000),polar_mars(5000)) #to generate 5000 uniform numbers
microbenchmark(box_muller(5000),polar_mars(ceiling(4*5000/pi))) #to generate 5000 normal numbers (approx)

microbenchmark(box_muller(20000),polar_mars(20000)) #to generate 20000 uniform numbers
microbenchmark(box_muller(20000),polar_mars(26000)) #to generate 20000 normal numbers (approx)

#it is observed that the Polar Marsaglia Method is almost as fast as the Box Muller Method when the same
# number of uniform numbers are used. However, to generat the same number of normally distributed numbers,
# it is observed that the Box-Muller test is much faster.
