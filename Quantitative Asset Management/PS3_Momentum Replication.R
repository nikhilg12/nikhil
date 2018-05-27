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
if(!require(tidyr)){
  install.packages("tidyr")
}
if(!require(moments)){
  install.packages("moments")
}


library(data.table)
library(zoo)
library(xts) #for selecting dates
library(lubridate) #for ymd
library(tidyr)
PS3_Q1<-function(dataps1){
  #set date format
  dataps1<-dataps1[, date := ymd(date)]
  
  #filter data by exchange codes and share codes
  dataps1<-dataps1[SHRCD %in% c(10,11) & EXCHCD %in% c(1,2,3)]
  
  #Setting the non-numeric returns as missing
  dataps1$RET[dataps1$RET %in% c("C","B")]<-NA
  dataps1$DLRET[dataps1$DLRET %in% c("A","S","P","T")]<-NA
  
  #calculating net returns
  dataps1<-dataps1[,TOTRET := ifelse(!is.na(as.numeric(dataps1$RET)) & !is.na(as.numeric(dataps1$DLRET)),(1+as.numeric(dataps1$RET)*(1+as.numeric(dataps1$DLRET))-1),
                                     ifelse(is.na(as.numeric(dataps1$RET)),as.numeric(dataps1$DLRET),as.numeric(dataps1$RET)))]
  
  #calculating market equity and lagging it
  dataps1<-dataps1[,MKTEQ := abs(as.numeric(dataps1$PRC))*as.numeric(dataps1$SHROUT)]
  dataps1<-dataps1[,MKTEQLAG := shift(MKTEQ),by=c("PERMNO")]
  
  #extracting the year and month
  dataps1[,c('Month','Year') := .(month(date),year(date))]
  dataps1<-dataps1[, prc13lag:=shift(PRC,13),by=.(PERMNO)]
  
  #creating lags to remove the ones which have less than 7 months of data
  dataps1<-dataps1[,c("Lag 2","Lag 3","Lag 4","Lag 5","Lag 6","Lag 7","Lag 8","Lag 9","Lag 10","Lag 11","Lag 12"):=
                     .(shift(TOTRET,2,fill = 0),shift(TOTRET,3,fill = 0),shift(TOTRET,4,fill = 0),shift(TOTRET,5,fill = 0),shift(TOTRET,6,fill = 0),
                       shift(TOTRET,7,fill = 0),shift(TOTRET,8,fill = 0),shift(TOTRET,9,fill = 0),shift(TOTRET,10,fill = 0),
                       shift(TOTRET,11,fill = 0),shift(TOTRET,12,fill = 0)),by=.(PERMNO)]
  
  #checking the number of returns available over an 11 month period
  dataps1<-dataps1[,Min8:=.(sum(abs(`Lag 2`)>0,abs(`Lag 3`)>0,abs(`Lag 4`)>0,abs(`Lag 5`)>0,abs(`Lag 6`)>0,abs(`Lag 7`)>0,abs(`Lag 8`)>0,
                                abs(`Lag 9`)>0,abs(`Lag 10`)>0,abs(`Lag 11`)>0,abs(`Lag 12`)>0)),by=.(Year,Month,PERMNO)]
  
  #calculating the ranking return
  dataps1<-dataps1[,rankret:=log(1+`Lag 2`)+log(1+`Lag 3`)+log(1+`Lag 4`)+log(1+`Lag 5`)+log(1+`Lag 6`)+log(1+`Lag 7`)+
                     log(1+`Lag 8`)+log(1+`Lag 9`)+log(1+`Lag 10`)+log(1+`Lag 11`)+log(1+`Lag 12`)]
  
  #selecting stock values >=8
  dataps13<-dataps1[Min8>=8]
  
  
  dataps1<-dataps1[!is.na(as.numeric(prc13lag))]
  dataps1<-dataps1[!abs(`Lag 2`)==0]
  out1<-data.table(date=dataps1$date,Year=dataps1$Year,Month=dataps1$Month,
                   PERMNO=dataps1$PERMNO,EXCHCD=dataps1$EXCHCD,
                   lag_Mkt_Cap=dataps1$MKTEQLAG,Ret=dataps1$TOTRET,
                   Ranking_Ret=dataps1$rankret)
  
  #selecting desired years
  out1<-out1[!Year %in% c(1925,1926)]
  out1<-out1[!is.na(`lag_Mkt_Cap`)]
  out1<-out1[order(Year,Month,PERMNO)]
  
  return(out1)
}
CRSP_Stocks=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/dataps1.csv")

CRSP_Stocks_Momentum=PS3_Q1(CRSP_Stocks)

PS3_Q2=function(neudata){  
  GetPortNums <- function(x, numports=10) {
    as.integer(cut(x,
                   unique(quantile(x, probs=0:numports/numports, na.rm = TRUE)),
                   include.lowest=TRUE))
  }
  
  #Calculate the stock rankings based on quantiles (equally distributed)
  neudata[,DM_decile:=lapply(.(Ranking_Ret),GetPortNums),by=(date)]
  
  GetPortNumsKRF <- function(x,y) { #I add -infinity and +infinity to 8 quantiles to make it 10
    breakpoints <- c(-50,quantile(x, probs=seq(0,1,0.1)[2:10],na.rm=TRUE),50)
    
    as.integer(cut(y,
                   breakpoints,
                   include.lowest=TRUE))
  }
  
  neudata[,KRF_decile:=mapply(GetPortNumsKRF,.(Ranking_Ret[which(EXCHCD==1)]),.(Ranking_Ret)),by=(date)]
  
  #remove NAs
  out2=neudata[complete.cases(neudata),-c("EXCHCD","Ranking_Ret")]
  out2
}

CRSP_Stocks_Momentum_decile=PS3_Q2(CRSP_Stocks_Momentum)

PS3_Q3=function(ffdata,newdata){
  newdata=newdata[!is.na(Ret)]
  #function which calculates value weights given two vectors: market capitalization and returns
  vwrets=function(rets,mktcaps){
    vwts=mktcaps/sum(mktcaps,na.rm=TRUE)
    vw=sum(vwts*rets,na.rm=TRUE)
    vw
  }
  
  reqdataDM=newdata[,c(1,5,6,7)] #selecting required columns (for DM_Deciles only)
  
  groupedDM=reqdataDM[,.(decileretsDM=lapply(.SD,vwrets,lag_Mkt_Cap)),.SDcols="Ret",by=c("date","DM_decile")]
  groupedDM=groupedDM[order(DM_decile,date)]
  
  
  reqdataKRF=newdata[,c(1,5,6,8)] #selecting required columns (for KRF_Deciles only)
  
  #take the value weighted returns
  groupedKRF=reqdataKRF[,.(decileretsKRF=lapply(.SD,vwrets,lag_Mkt_Cap)),.SDcols="Ret",by=c("date","KRF_decile")]
  groupedKRF=groupedKRF[order(KRF_decile,date)]
  
  out3=cbind(groupedDM,unlist(groupedKRF$decileretsKRF))
  out3$decileretsDM=unlist(out3$decileretsDM)
  riskfree=ffdata[V1>=192701 & V1<=201712,"RF"]/100
  
  out3=cbind(out3,riskfree)
  colnames(out3)=c("date","decile","DM_Ret","KRF_Ret","Rf")
  
  out3=out3[order(date,decile)]
  
  out3
}

ff=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW3/10_Portfolios_Prior_12_2.CSV")
FF_mkt=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/mktdata.CSV")
ffdata=FF_mkt

CRSP_Stocks_Momentum_returns=PS3_Q3(FF_mkt,CRSP_Stocks_Momentum_decile)

PS3_Q4=function(stockmomdata){
  #converting to excess
  stockmomdata=stockmomdata[complete.cases(stockmomdata)]
  
  stockmomdata[,3]=stockmomdata[,3] - stockmomdata[,5]
  stockmomdata[,4]=stockmomdata[,4] - stockmomdata[,5]
  
  #no need of RF column now that we have excess returns
  stockmomdata[,5]=NULL
  
  #separating DM and KRF datas
  stockmomdatadm=stockmomdata[,1:3]
  stockmomdatakrf=stockmomdata[,c(1,2,4)]
  
  #reshaping the data so that for each date, I see all the decile returns
  stockmomdatadm=stockmomdatadm %>% spread(decile,DM_Ret)
  stockmomdatakrf=stockmomdatakrf %>% spread(decile,KRF_Ret)
  
  #selecting returns only till March 2013
  matchingmomdata=stockmomdatadm[date<"2013-04-01"]
  
  #calculating winners minus losers portfolio returns
  matchingmomdata=cbind(matchingmomdata,matchingmomdata[,11]-matchingmomdata[,2])
  
  mdm=colMeans((matchingmomdata[,c(2:12)])*1200)
  sdm=apply(matchingmomdata[,c(2:12)],2,sd)*sqrt(12)*100
  srdm=mdm/sdm
  skdm=skewness(log(1+matchingmomdata[,c(2:12)]))
  
  out4=cbind(mdm,sdm,srdm,skdm)
  rownames(out4)=c(as.character(seq(1:10)),"WML")
  out4
}

q4=PS3_Q4(CRSP_Stocks_Momentum_returns)

#Daniel Moskowitz data
dir.create("dmport")
untar("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW3/DM_data_2017_03.tar.gz",exdir="dmport")
dm <- fread("dmport/m_m_pt_tot.txt")



PS3_Q5=function(stockmomdata,dmd,ffd){
  dmd <- dmd[,1:3] %>% spread(V2,V3)
  dmd$"WML"=dmd[,10]-dmd[,2]
  ffd$"WML"=ffd[,10]-ffd[,2]
  
  #converting to excess
  
  stockmomdata=stockmomdata[complete.cases(stockmomdata)]
  
  stockmomdata[,3]=stockmomdata[,3] - stockmomdata[,5]
  stockmomdata[,4]=stockmomdata[,4] - stockmomdata[,5]
  
  stockmomdata[,5]=NULL
  
  #separating DM and KRF data
  stockmomdatadm=stockmomdata[,1:3]
  stockmomdatakrf=stockmomdata[,c(1,2,4)]
  
  #reshaping
  stockmomdatadm=stockmomdatadm %>% spread(decile,DM_Ret)
  stockmomdatakrf=stockmomdatakrf %>% spread(decile,KRF_Ret)
  
  #WML returns
  stockmomdatadm$WML=stockmomdatadm[,11]-stockmomdatadm[,2]
  stockmomdatakrf$WML=stockmomdatakrf[,11]-stockmomdatakrf[,2]
  
  #taking correlations in the same timeframe
  dmcormat=diag(cor(stockmomdatadm[1:1080,c(2:12)],dmd[,-1]))
  ffcormat=diag(cor(stockmomdatakrf[,c(2:12)],ffd[1:1092,-1]))
  cbind(dmcormat,ffcormat)
}

q5=PS3_Q5(CRSP_Stocks_Momentum_returns,dm,ff)


PS3_Q6=function(stockmomdata){
  #converting to excess
  stockmomdata=stockmomdata[complete.cases(stockmomdata)]
  stockmomdata[,3]=stockmomdata[,3] - stockmomdata[,5]
  stockmomdatadm=stockmomdata[,1:3]
  stockmomdatadm=stockmomdatadm %>% spread(decile,DM_Ret)
  
  recent=stockmomdatadm[date>"2007-12-31",]
  WML=recent[,11]-recent[,2] #10-1
  WML=unlist(WML)
  dates=recent$date
  plot(dates,cumprod(1+WML)-1,type='l',main="Cumulative Momentum Returns for 2011-2017",xlab="Months",ylab="")
  skewness(WML)
}

q6=PS3_Q6(CRSP_Stocks_Momentum_returns)


