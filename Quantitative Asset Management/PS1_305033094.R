#Nikhil Guruji, Cohort 1
#Installing required libraries

if(!require(xts)){
  install.packages("xts")
}
if(!require(data.table)){
  install.packages("data.table")
}
if(!require(zoo)){
  install.packages("zoo")
}
if(!require(lubridate)){
  install.packages("lubridate")
}
if(!require(moments)){
  install.packages("moments")
}

library(data.table)
library(zoo)
library(xts) #for selecting dates
library(lubridate) #for ymd
library(moments) #for skewness and kurtosis

PS1_Q1=function(dataps1){

  dataps1$date=ymd(dataps1$date) #setting the date format

  #converting returns from factor to numeric (will generate NAs! Because some of the returns are "blank")
  dataps1$DLRET=as.numeric(as.character(dataps1$DLRET))
  dataps1$RET=as.numeric(as.character(dataps1$RET))
  
  ##calculating the total returns
  
  #setting all available RET values as TOTRET
  dataps1$TOTRET[which(!is.na(dataps1$RET))]=dataps1$RET[which(!is.na(dataps1$RET))]
  
  #setting all available DLRET values as TOTRET
  dataps1$TOTRET[which(!is.na(dataps1$DLRET))]=dataps1$DLRET[which(!is.na(dataps1$DLRET))]
  
  #takes rt=(1+rd)*(1+rh) - 1 for complete cases (when both are available)
  dataps1$TOTRET[which(complete.cases(dataps1$RET,dataps1$DLRET))]=
    (1+dataps1$DLRET[complete.cases(dataps1$RET,dataps1$DLRET)])*
    (1+dataps1$RET[complete.cases(dataps1$RET,dataps1$DLRET)]) - 1
  
  #calculate market equity
  dataps1$MKTEQ = abs(dataps1$PRC)*abs(dataps1$SHROUT)
  dataps1$MKTEQLAG=shift(dataps1$MKTEQ,fill=0) #calculate the lagged market equity
  
  #removing the rows where Total return or Lagged Market Equity is missing (doing this here and not earlier 
  # ensures values don't spillover to other stocks)
  dataps1=dataps1[which(complete.cases(dataps1$TOTRET,dataps1$MKTEQLAG)),] 
  dataps1=dataps1[SHRCD%in%c(10,11)] #selecting only stocks with sharecodes 10,11
  dataps1=dataps1[EXCHCD%in%c(1,2,3)] # selecting only stocks with exchange codes 1,2,3
  
  #function which calculates value weights given two vectors: market capitalization and returns
  vwrets=function(rets,mktcaps){
    vwts=mktcaps/sum(mktcaps)
    vw=sum(vwts*rets)
    vw
  }
  
  ews=dataps1[,lapply(.SD,mean),by="date",.SDcols="TOTRET"] #apply mean function to column TOTRET and group by date
  ews=ews[order(date)]
  
  mktlag=dataps1[,lapply(.SD,sum),by="date",.SDcols="MKTEQLAG"] #apply sum function to MKTEQLAG and group by date
  mktlag=mktlag[order(date)]

  vws=dataps1[,lapply(.SD,vwrets,MKTEQLAG),by="date",.SDcols="TOTRET"]
  vws=vws[order(date)]
  
  #merging the data together
  finaldata=merge(vws,ews,by="date")
  finaldata[,Stock_lag_MV:=mktlag$MKTEQLAG]
  
  #adding more columns
  finaldata[,Year:=year(date)]
  finaldata[,Month:=month(date)]
  finaldata[,YM:=paste(year(date),month(date))]
  
  #setting column names
  setnames(finaldata,c("date","Stock_Vw_Ret","Stock_Ew_Ret","Stock_lag_MV","Year","Month","YM"))
  
  finaldata
}

#1

#reading the CRSP monthly data
CRSP_Stocks=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/dataps1.csv")
Monthly_CRSP_Stocks=PS1_Q1(CRSP_Stocks)

#reading the Fama french 3 factor monthly data
FF_mkt=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/mktdata.csv")

#2
PS1_Q2=function(ffdata,q1data){

  #Converting returns to decimal
  ffdata$RF=ffdata$RF/100
  ffdata$`Mkt-RF`=ffdata$`Mkt-RF`/100
  
  ffdata$V1=paste(ffdata$V1,"01") #adding day to the date so that ymd workds
  ffdata$V1=ymd(ffdata$V1) #converting to date format
  ffdata$YM=paste(year(ffdata$V1),month(ffdata$V1)) #same YM column
  
  #outerjoin
  outdata=merge(q1data,ffdata[,c(2,5,6)],by="YM",sort=TRUE) #selecting Market-Rf, RF and YM columns from FFdata
  outdata=outdata[order(Year,Month)] #resorting by year and then month
  
  meanff=mean(outdata$`Mkt-RF`)
  sdff=sd(outdata$`Mkt-RF`)
  meanest=mean(outdata$Stock_Vw_Ret-outdata$RF)
  sdest=sd(outdata$Stock_Vw_Ret-outdata$RF)
  
  #annualizing mean and standard deviations
  annmeanff=(1+meanff)^12 - 1
  annsdff=sqrt((sdff^2 + (1+meanff)^2)^12 - (1+meanff)^24)
  annmeanest=(1+meanest)^12 - 1
  annsdest=sqrt((sdest^2 + (1+meanest)^2)^12 - (1+meanest)^24)
  
  out=matrix(c(annmeanff, #annualized mean of excess returns
               
               annsdff, #annualized sd of excess returns
               
               annmeanff/annsdff, #annualized Sharpe Ratio
               
               skewness(outdata$`Mkt-RF`), #skewness
               
               kurtosis(outdata$`Mkt-RF`)-3, #excess kurtosis
               
               annmeanest, #annualized mean of vw excess returns estimated
               
               annsdest, #annualized sd of vw excess returns estimated
               
               annmeanest/annsdest,
               
               skewness(outdata$Stock_Vw_Ret-outdata$RF),
               
               kurtosis(outdata$Stock_Vw_Ret-outdata$RF)-3 #excess kurtosis
               ),nrow=5,ncol=2)
  
  rownames(out)=c("Annualized Mean","Annualized SD","Annualized Sharpe Ratio","Skewness","Excess Kurtosis")
  colnames(out)=c("Market (FF) (Actual)","Value-Weighted Returns (Estimated)")
  out
}

Q2_out=PS1_Q2(FF_mkt,Monthly_CRSP_Stocks)

#3
PS1_Q3=function(ffdata,q1data){
  ffdata$V1=paste(ffdata$V1,"01") #adding day to the date so that ymd workds
  ffdata$V1=ymd(ffdata$V1) #converting to date format
  ffdata = xts(ffdata[,2:ncol(ffdata)],order.by = ffdata$V1)
  ffdata = ffdata["1926-07-01::2017-12-01"] #selecting desired dates
  ffdata = as.data.table(ffdata) #convert back to data table from xts (first column name has changed)

  ffdata$YM=paste(year(ffdata$index),month(ffdata$index)) #same YM column
  ffdata$MKT=(ffdata$`Mkt-RF`+ffdata$RF)/100 #converting market returns to decimal
  
  #outerjoin
  outdata=merge(q1data,ffdata[,c(6,7)],by="YM",sort=TRUE) #selecting market and YM columns from FFdata
  outdata=outdata[order(Year,Month)] #resorting by year and then month
  c(cor(outdata$Stock_Vw_Ret,outdata$MKT),max(abs(outdata$Stock_Vw_Ret-outdata$MKT)))
}

Q3_out=PS1_Q3(FF_mkt,Monthly_CRSP_Stocks)
names(Q3_out)=c("Correlation","Max Absolute Difference")
