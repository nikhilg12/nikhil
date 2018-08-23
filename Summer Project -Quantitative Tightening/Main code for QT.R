#Quantitative Strategy
#Z_construction 950
library(readxl)
library(data.table)
library(lubridate)
library(zoo)
library(readr)
library(R.matlab)
library(matlib)
library(MASS)
library(optimx)
library(matrixcalc)
library(plyr)
library(dplyr)
library(tidyr)
library(bit64)

#####
#Load Parameters.

exact_parameters<-readMat("parameters_20100811.mat")

#Loading U.S.Treasury Data
UST<-fread("US Treasury.csv")
UST[,TMATDT:=mdy(TMATDT)][,MCALDT:=mdy(MCALDT)]


#Hanmilton & Wu's Z_public
z_p=readMat("z_public.mat")
z_2011=readMat("z_201101_public.mat")

Z_hw=rbind(z_p[[1]],z_2011[[2]])

Z_2_org=z_p[[2]]

#return the nearest Friday date
nextweekday <- function(date, wday) {
  date <- as.Date(date)
  diff <- wday - wday(date)
  return(date + diff)
}

UST[,FRIDT:=nextweekday(MCALDT,6)] #date of the nearest friday for each date
UST[,Maturity_d:=(TMATDT-FRIDT)][,Maturity_w:=as.integer(round(Maturity_d/7))]

temp=UST[,c("MCALDT","TMTOTOUT","Maturity_w")]
temp=temp[!is.na(TMTOTOUT),]
Z<-temp[,.(tot=sum(TMTOTOUT)),by=c("Maturity_w","MCALDT")] #adding the market value for similar maturities for every date
setcolorder(Z,c("MCALDT","Maturity_w","tot"))
setorder(Z,"MCALDT")


Z_2<-Z %>% spread(Maturity_w,tot) #reshape the data
Z_2[,c("1580","1585","1589"):=NULL] #Delete unused columns.
setnames(Z_2,"MCALDT","Date")

#Using Hanmilton&Wu's Treasury.
Date_90_03<-fread("Date_90_03.csv")

################
SOMA_old<-as.data.table(z_p[[1]]-z_p[[2]]) #SOMA from Hanmilton & Wu by Z-(Z-SOMA)
tempcnames=colnames(SOMA_old)
SOMA_old$Date<-Date_90_03$Date
SOMA_old[,Date:=dmy(paste("01",Date,sep="-"))]
setcolorder(SOMA_old,c("Date",tempcnames))


colnames(SOMA_old)<-c("Date",as.character(seq(0,1576)))
# Input Hanmilton & Wu's Z_public
Z_public_HW=z_p[[2]] ## Z_public from Hamilton & Wu data
Empty<-as.data.table(matrix(NA,1,ncol(SOMA_old)))  #Adding missing columns to Z by merge.
colnames(Empty)<-c("Date",as.character(seq(0,1576)))
Empty[,1]<-as.Date(c("2222-08-09"))
Z_2_merge<-merge(Empty,Z_2,by=colnames(Z_2),all.x=TRUE,all.y=TRUE)
Z_2_merge<-Z_2_merge[Date!="2222-08-09",]
Z_2_merge<-Z_2_merge[Date>="2011-02-01",]
Z_2_merge=as.matrix(Z_2_merge[,c(as.character(seq(0,1576)))]) #reordering the columns


Z_new=data.table(rbind(as.matrix(Z_hw),Z_2_merge))
Z_new[is.na(Z_new)]=0
# Calculate proportional percentage using Hanmilton & Wu Z_public
Z_new_Percent = as.data.table(matrix(numeric(),nrow(Z_new),ncol(Z_new)))


################

row_entry2 = matrix(0,nrow(Z_new_Percent),ncol(Z_new_Percent))
row_entry2[,1:2] = as.matrix(Z_new[,1:2])/rowSums(Z_new[,1:2]) # < 15 days
row_entry2[,3:13] = as.matrix(Z_new[,3:13])/rowSums(Z_new[,3:13]) # 15 to 90 days
row_entry2[,14:52] = as.matrix(Z_new[,14:52])/rowSums(Z_new[,14:52]) # 90 days to 1 years
row_entry2[,53:260] = as.matrix(Z_new[,53:260])/rowSums(Z_new[,53:260]) # 1 to 5 years
row_entry2[,261:520] = as.matrix(Z_new[,261:520])/rowSums(Z_new[,261:520]) # 5 to 10 years
row_entry2[,521:1577] = as.matrix(Z_new[,521:1577])/rowSums(Z_new[,521:1577]) # > 10 years
Z_new_Percent = as.data.table(row_entry2)

################
#Input SOMA six bins
SOMA<-data.table(read_xls("public debt data and Fed holdings_original.xls",sheet="SOMA by maturity - TIPS"))
setnames(SOMA,"X__1","Date")
SOMA[,Date:=ymd(Date)]
SOMA = SOMA[Date <"2018-07-01",] # Select dates from 1990 to 2009 to match Hamilton & Wu's data for comparison purpose

# Calculate expanded SOMA using proportional percentage above and our inputed SOMA six bins
SOMA_expanded = as.data.table(matrix(numeric(),nrow(Z_new_Percent),ncol(Z_new_Percent)))

################

row_entry3 = matrix(0,nrow(SOMA_expanded),ncol(SOMA_expanded))
row_entry3[,1:2] = as.matrix(Z_new_Percent[,1:2])*unlist(SOMA[,2]) # < 15 days
row_entry3[,3:13] = as.matrix(Z_new_Percent[,3:13])*unlist(SOMA[,3]) # 15 to 90 days
row_entry3[,14:52] = as.matrix(Z_new_Percent[,14:52])*unlist(SOMA[,4]) # 90 days to 1 years
row_entry3[,53:260] = as.matrix(Z_new_Percent[,53:260])*unlist(SOMA[,5]) # 1 to 5 years
row_entry3[,261:520] = as.matrix(Z_new_Percent[,261:520])*unlist(SOMA[,6]) # 5 to 10 years
row_entry3[,521:1577] = as.matrix(Z_new_Percent[,521:1577])*unlist(SOMA[,7]) # > 10 years
SOMA_expanded= as.data.table(row_entry3)
tempcnames=colnames(SOMA_expanded)
SOMA_expanded$Date<-Date_90_03$Date[1:342]
SOMA_expanded[,Date:=dmy(paste("01",Date,sep="-"))]
setcolorder(SOMA_expanded,c("Date",tempcnames))

colnames(SOMA_expanded)<-c("Date",as.character(seq(0,1576)))
SOMA_merge<-SOMA_expanded #

################

#Take same sample period as in paper(line 32 in code).
SOMA_new<-SOMA_merge[1:342,]

#Take the difference of Treasury and Fed.
Z_public<-Z_new-SOMA_new[,-c("Date")] 

Z_public[Z_public<0,]<-0


Z_public[,"Total":=rowSums(Z_public,na.rm=T)]
z_public<-copy(Z_public)

z_public[is.na(z_public)]<-0
z_public<-z_public/z_public$Total
z_public<-as.matrix(z_public)
z_public<-z_public[,-ncol(z_public)]


#Parameters.

rho_Q=exact_parameters[[1]]
Sigma=exact_parameters[[2]]
c=exact_parameters[[3]]
rho=exact_parameters[[4]]
b1=exact_parameters[[5]]
Lambda=exact_parameters[[6]]
a1=exact_parameters[[7]]
c_Q=exact_parameters[[8]]
lambda=exact_parameters[[9]]
N=ncol(z_public)

b_<-matrix(0,3,N)
b_[,1]<--b1
for(n in 2:N)
{
  b_[,n]<-b_[,(n-1)]%*%rho_Q-t(b1)
}  

a_=vector(length = N)
a_[1]=-a1
for(i in 1:(N-1)){
  a_[i+1]=a_[i]+t(b_[,i])%*%c_Q + 0.5*t(b_[,i])%*%Sigma%*%t(Sigma)%*%b_[,i] - a1
}

q<-matrix(0,3,nrow(z_public))
for(t in 1:nrow(z_public))
{
  temp=0
  for(n in 2:N)
  {
    temp<-temp+z_public[t,n]*b_[,n-1]
  }
  q[,t]<-100*Sigma%*%t(Sigma)%*%temp
}

#Factors.
yield_DGS6MO<-fread("DGS6MO.csv")
yield_DGS2<-fread("DGS2.csv")
yield_DGS10<-fread("DGS10.csv")

yield_d<-merge(merge(yield_DGS6MO,yield_DGS2,by="DATE"),yield_DGS10,by="DATE")
yield_d[,Year:=year(DATE)][,Month:=month(DATE)][,Day:=day(DATE)]
setorder(yield_d,DATE)
yield_d[,MaxDay:=max(Day),by=c("Year","Month")]  #Select the last line of each month.
yield_m<-yield_d[which(Day==MaxDay),]

yield_m[,DATE:=ymd(DATE)]
yield_m<-yield_m[DATE<'2018-07-31' & DATE>='1990-01-01',][,DGS6MO:=as.numeric(as.character(DGS6MO))][,DGS2:=as.numeric(as.character(DGS2))][,DGS10:=as.numeric(as.character(DGS10))]
yield_m[,F1:=(DGS6MO+DGS2+DGS10)/3][,F2:=(DGS10-DGS6MO)][,F3:=DGS10-DGS2*2+DGS6MO] #Factor1,2 &3.
yield_m[,F1_lag:=shift(F1)][,F2_lag:=shift(F2)][,F3_lag:=shift(F3)]
yield_m[,F1_dm:=F1_lag-mean(F1_lag,na.rm=TRUE)][,F2_dm:=F2_lag-mean(F2_lag,na.rm=TRUE)][,F3_dm:=F3_lag-mean(F3_lag,na.rm=TRUE)] #De-mean factors used.


q_m<-q[,1:nrow(yield_m)]
for(i in 1:3)
{q_m[i,]<-shift(q_m[i,])}


phi1<-summary(lm((yield_m$F1-mean(yield_m$F1,na.rm=TRUE))~yield_m$F1_dm+yield_m$F2_dm+yield_m$F3_dm+q_m[1,]+q_m[2,]+q_m[3,]))$coefficients[5:7]
phi2<-summary(lm((yield_m$F2-mean(yield_m$F2,na.rm=TRUE))~yield_m$F1_dm+yield_m$F2_dm+yield_m$F3_dm+q_m[1,]+q_m[2,]+q_m[3,]))$coefficients[5:7]
phi3<-summary(lm((yield_m$F3-mean(yield_m$F3,na.rm=TRUE))~yield_m$F1_dm+yield_m$F2_dm+yield_m$F3_dm+q_m[1,]+q_m[2,]+q_m[3,]))$coefficients[5:7]

phi<-matrix(c(phi1,phi2,phi3),byrow=TRUE,nrow=3)
delta<-c(0.026101,0.022712,-0.00780) #delta for QE
(phi)%*%(delta)
phi
#####

##Implementing new QT strategies.

#last row of SOMA consists of data for June 2018. We will forecast upto December 2019
SOMA_forecast_QT=matrix(0,nrow=18,ncol=1577)
SOMA_forecast_QT[1,]=as.matrix(SOMA_new[nrow(SOMA_new),-"Date"])

Z_forecast=matrix(0,nrow=18,ncol=1577)
Z_forecast[1:12,]=as.matrix(Z_new[(nrow(Z_new)-11):nrow(Z_new),])*1.05 #Assuming Treasury holding grows by 5% every year.
for(i in 13 : nrow(Z_forecast))
{
  Z_forecast[i,]=Z_forecast[(i-12),]*1.05 #Because 12 rows are defined in the previous line.
}


#With QT: reinvest with cap, assuming ALL maturities are available for purchase at the beginning of the month.
for(i in 2:18)
{
  Cap=ifelse(i<=4,24,30) #in billions. 
  Proceeds<-max(sum(SOMA_forecast_QT[i-1,1:4])-Cap*1000,0)  # in millions.
  Wts<-SOMA_forecast_QT[i-1,]/sum(SOMA_forecast_QT[i-1,])
  SOMA_forecast_QT[i,]=unlist(shift(SOMA_forecast_QT[i-1,],4,fill=0,type="lead"))
  SOMA_forecast_QT[i,]=SOMA_forecast_QT[i,]+Proceeds*Wts
}

Z_public_forecast_QT<-data.table(Z_forecast-SOMA_forecast_QT)
Z_public_forecast_QT[Z_public_forecast_QT<0,]<-0
Z_public_forecast_QT[,"Total":=rowSums(Z_public_forecast_QT,na.rm=T)]
Z_public_QT<-data.table(rbind(as.matrix(Z_public),as.matrix(Z_public_forecast_QT)))
colnames(Z_public_QT)<-c(as.character(seq(0,1576)),"Total")

#Without QT: reinvesting all the proceeds

# Removing QT effect from October 2017 to July 2018
# adding back 6 billion from Oct-Dec, 12 billion from Jan-Mar, 18 billion from Apr-Jun
# that means subtracting the above amounts from Z_public proportionately
SOMA_non_QT=copy(as.matrix(SOMA_new[,-"Date"]))

for(i in 334:342){ #Oct 2017 is 334th row, June 2018 is 342nd row
  cap=ifelse(i<=336,6000,ifelse(i<=339 & i>=337,12000,ifelse(i<=342 & i>=340,18000,0)))
  wts=as.numeric(SOMA_non_QT[i,]/sum(as.numeric(SOMA_new[i,])))
  addback=cap*wts
  SOMA_non_QT[i,]=as.numeric(SOMA_non_QT[i,])+addback
}

#new Z_public based on removing the effect of QT
Z_public_non_QT=copy(Z_public)
Z_public_non_QT$Total=NULL
Z_public_non_QT[334:342,]=Z_new[334:342,]-SOMA_non_QT[334:342,]
Z_public_non_QT[,"Total":=rowSums(Z_public_non_QT,na.rm=T)]

#last row of SOMA consists of data for June 2018. We will forecast upto December 2019
SOMA_forecast=matrix(0,nrow=18,ncol=1577)
SOMA_forecast[1,]=as.matrix(SOMA_non_QT[nrow(SOMA_non_QT),])

for(i in 2:18)
{
  Proceeds<-sum(SOMA_forecast[i-1,1:4]) # in millions.
  Wts<-SOMA_forecast[i-1,]/sum(SOMA_forecast[i-1,])
  SOMA_forecast[i,]=unlist(shift(SOMA_forecast[i-1,],4,fill=0,type="lead"))
  SOMA_forecast[i,]=SOMA_forecast[i,]+Proceeds*Wts
}

Z_public_forecast<-data.table(Z_forecast-SOMA_forecast)
Z_public_forecast[Z_public_forecast<0,]<-0
Z_public_forecast[,"Total":=rowSums(Z_public_forecast,na.rm=T)]
Z_public_normal<-data.table(rbind(as.matrix(Z_public_non_QT),as.matrix(Z_public_forecast)))
colnames(Z_public_normal)<-c(as.character(seq(0,1576)),"Total")

#################################
#########Calculating new q_t#####
#################################

z_publicA=as.data.table(Z_public_QT)
z_publicA<-z_publicA/z_publicA$Total
z_publicA<-as.matrix(z_publicA)
z_publicA<-z_publicA[,-ncol(z_publicA)]

q_A<-matrix(0,3,nrow(z_publicA))
for(t in 1:nrow(z_publicA))
{
  temp=0
  for(n in 2:N)
  {
    temp<-temp+z_publicA[t,n]*b_[,n-1]
  }
  q_A[,t]<-100*Sigma%*%t(Sigma)%*%temp
}

#Non-QT framework
z_publicB=as.data.table(Z_public_normal)
z_publicB<-z_publicB/z_publicB$Total
z_publicB<-as.matrix(z_publicB)
z_publicB<-z_publicB[,-ncol(z_publicB)]

q_B<-matrix(0,3,nrow(z_publicB))
for(t in 1:nrow(z_publicB))
{
  temp=0
  for(n in 2:N)
  {
    temp<-temp+z_publicB[t,n]*b_[,n-1]
  }
  q_B[,t]<-100*Sigma%*%t(Sigma)%*%temp
}


delta=rowSums((q_A-q_B)[,334:ncol(q_A)],na.rm=T)

phi%*%delta

#In line 218 to 230, 266 and 269, change 18 to 30 or 42 to implement QT for 2019, 2020, and 2021. 
