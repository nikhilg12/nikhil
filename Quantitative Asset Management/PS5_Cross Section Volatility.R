#Nikhil Guruji Cohort 1

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
if(!require(DataAnalytics)){
  install.packages("DataAnalytics")
}
if(!require(xtable)){
  install.packages("xtable")
}


library(data.table)
library(zoo)
library(xts) #for selecting dates
library(lubridate) #for ymd
library(tidyr)
library(DataAnalytics)
library(xtable)

dailydata <- fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW5/stocksdaily.csv")
monthlydata <- fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/dataps1.csv")


####cleaning

cleaner <- function(crspdata){
  #date format
  crspdata[,date:=ymd(date)]
  
  #filter data by exchange codes and share codes
  #crspdata<-crspdata[SHRCD %in% c(10,11) & EXCHCD %in% c(1,2,3)]
  crspdata<-crspdata[EXCHCD %in% c(1,2,3)]
  
  #Setting the non-numeric returns as missing
  crspdata$RET[crspdata$RET %in% c("C","B")]<-NA
  crspdata$DLRET[crspdata$DLRET %in% c("A","S","P","T")]<-NA
  
  crspdata$DLRET=as.numeric(crspdata$DLRET)
  crspdata$RET=as.numeric(crspdata$RET)
  
  #calculating net returns
  crspdata[is.na(DLRET) & !is.na(RET),
           TOTRET := RET][!is.na(DLRET) & is.na(RET),
                          TOTRET := DLRET][!is.na(DLRET) & !is.na(RET),
                                           TOTRET := (1 + RET)*(1 + DLRET) - 1]
  crspdata[,Year:=year(date)][,Month:=month(date)]
  
  return(crspdata)
}

dailydata <- cleaner(dailydata)
monthlydata <- cleaner(monthlydata)
monthlydata[,MKTEQ := abs(as.numeric(monthlydata$PRC))*as.numeric(monthlydata$SHROUT)/1000] #converting to millions
monthlydata[,lag_Mkt_Cap:=shift(MKTEQ,fill=NA),by=PERMNO]
monthlydata <- monthlydata[date<="2015-12-31" & date>="2012-12-01"]

cleanday<-copy(dailydata)
cleanmonth<-copy(monthlydata)
monthlydata<-copy(cleanmonth)
dailydata<-copy(cleanday)

#count the number of returns available
obscount <- dailydata[,.(Obs=length(which(!is.na(TOTRET)))),by=c("Year","Month","PERMNO")]

#calculate monthly volatility
monthlyvols <- dailydata[,lapply(.SD,sd,na.rm=TRUE),.SDcols="TOTRET",by=c("Year","Month","PERMNO")]
setnames(monthlyvols,c("Year","Month","PERMNO","Volatility"))

#final merge
maindata <- merge(monthlyvols,monthlydata,by=c("Year","Month","PERMNO"))

#filter data to include observations with at least 17 values of returns available
maindata <- merge(maindata,obscount,by=c("Year","Month","PERMNO"))

maindata=maindata[which(Obs>17),]

#create quintiles

GetPortNums <- function(x, numports=5) {
  as.integer(cut(x,
                 quantile(x, probs=0:numports/numports, na.rm = TRUE),
                 include.lowest=TRUE))
}

#Calculate the quintiles based on volatility

#lag the volatility
maindata[,lag_vol:=shift(Volatility,fill=NA),by=PERMNO]

#select the required subset of values
maindata <- maindata[date<="2015-12-31" & date>="2013-01-01"]
maindata<-maindata[!is.na(lag_vol)]

#calculating the quintiles
maindata[,Qtl_Tot_vol:=lapply(.(lag_vol),GetPortNums),by=(date)]
checkpt1=copy(maindata)
maindata=copy(checkpt1)

vwrets=function(rets,mktcaps){
  vwts=mktcaps/sum(mktcaps,na.rm=TRUE)
  vw=sum(vwts*rets,na.rm=TRUE)
  vw
}

#VW
Tot_Vol_ret <- maindata[,.(Qtl_Tot_Vol_Rets=lapply(.SD,vwrets,lag_Mkt_Cap)),.SDcols="TOTRET",by=c("date","Qtl_Tot_vol")]
#EW
#Tot_Vol_ret <- maindata[,.(Qtl_Tot_Vol_Rets=lapply(.SD,mean)),.SDcols="TOTRET",by=c("date","Qtl_Tot_vol")]

#Size
forsize=copy(maindata)
forsize=forsize[!is.na(MKTEQ)]
avg_size_tot_vol <- forsize[,.(size=lapply(log(.SD),mean)),.SDcols="MKTEQ",by="Qtl_Tot_vol"]
setorder(avg_size_tot_vol,Qtl_Tot_vol)
# % market share

perc_mkt_share <- forsize[,.(perc=lapply(.SD,sum)),.SDcols="MKTEQ",by="Qtl_Tot_vol"]
perc_mkt_share[,percent:=unlist(perc)/sum(unlist(perc))]
setorder(perc_mkt_share,Qtl_Tot_vol)

#mean and SD

tot_ret <- Tot_Vol_ret %>% spread(Qtl_Tot_vol,Qtl_Tot_Vol_Rets)
tot_ret_mat <- matrix(unlist(tot_ret[,-1]),ncol=5)
100*colMeans(tot_ret_mat)
apply(tot_ret_mat,2,sd)

#Fama French Monthly Data
ffdata <-fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW4/F-F_Research_Data_Factors.CSV")
ffdata<-ffdata[V1>=201301 & V1<=201512]

capmalpha=c()
ffalpha=c()
capmalphat=c()
ffalphat=c()
for(i in 1:5){
  capmalpha[i]=lm(tot_ret_mat[,i]~ffdata$`Mkt-RF`)$coefficients[1]*100 #converting to %
  ffalpha[i]=lm(tot_ret_mat[,i]~ffdata$`Mkt-RF`+ffdata$SMB+ffdata$HML)$coefficients[1]*100 #converting to %
  ct=lmSumm(lm(tot_ret_mat[,i]~ffdata$`Mkt-RF`),HAC=TRUE) #fifth element is the t-stat of intercept
  capmalphat[i]=ct$coef.table[5]
  ft=lmSumm(lm(tot_ret_mat[,i]~ffdata$`Mkt-RF`+ffdata$SMB+ffdata$HML),HAC=TRUE)
  ffalphat[i]=ft$coef.table[9]
}

#################
##idiosyncratic volatilities
############

#Fama French daily data
ffdaily<-fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW5/F-F_Research_Data_Factors_daily.CSV")
ffdaily<-ffdaily[V1>=20121201 & V1<=20151201]

#Take the required subset from daily data
dailysample<-dailydata[date<="2015-12-31" & date>="2012-12-01"]
dailysample[,Day:=day(date)]
ffdaily[,date:=ymd(V1)]
ffdaily[,Year:=year(date)][,Month:=month(date)][,Day:=day(date)]

#Merging by year,month and day so it is easier to regress
ff_daily_merge<-merge(dailysample,ffdaily[,-c("date","V1")],by=c("Year","Month","Day"))

#merging again to remove rows where returns < 17 are available for a month
ff_daily_merge_17 <- merge(ff_daily_merge,obscount,by=c("Year","Month","PERMNO"))
ff_daily_merge_17<-ff_daily_merge_17[Obs>17]
setorder(ff_daily_merge_17,PERMNO)

#remove the missing returns
small_sample<-ff_daily_merge_17[!is.na(TOTRET)]

#calculate the idiosyncratic volatilities
idiovols<-small_sample[,.(Idio_vol=summary(lm(TOTRET~`Mkt-RF`+SMB+HML))$sigma),
                       by=c("Year","Month","PERMNO")]

#Alternate method
#idiovols<-small_sample[,.(Idio_vol=sd(lm(TOTRET~`Mkt-RF`+SMB+HML)$residuals)),
#                    by=c("Year","Month","PERMNO")]

#merge these monthly idiosyncratic volatilities with the main monthly data
maindata_idio <- merge(idiovols,monthlydata,by=c("Year","Month","PERMNO"))

#calculating the quintiles based on lagged idiosyncratic volatilites
maindata_idio[,lag_idio_vol:=shift(Idio_vol,fill=NA),by=PERMNO]
maindata_idio<-maindata_idio[date<="2015-12-31" & date>="2013-01-01"]
maindata_idio<-maindata_idio[!is.na(lag_idio_vol)]
maindata_idio[,Qtl_Tot_vol:=lapply(.(lag_idio_vol),GetPortNums),by=(date)]
checkpt1_idio=copy(maindata_idio)
maindata_idio=copy(checkpt1_idio)


#VW
Tot_Vol_ret_idio <- maindata_idio[,.(Qtl_Tot_Vol_Rets=lapply(.SD,vwrets,lag_Mkt_Cap)),.SDcols="TOTRET",by=c("date","Qtl_Tot_vol")]

#Size
forsize_idio=copy(maindata_idio)
forsize_idio=forsize_idio[!is.na(MKTEQ)]
avg_size_tot_vol_idio <- forsize_idio[,.(size=lapply(log(.SD),mean)),.SDcols="MKTEQ",by="Qtl_Tot_vol"]
setorder(avg_size_tot_vol_idio,Qtl_Tot_vol)

# % market share

perc_mkt_share_idio <- forsize_idio[,.(perc=lapply(.SD,sum)),.SDcols="MKTEQ",by="Qtl_Tot_vol"]
perc_mkt_share_idio[,percent:=unlist(perc)/sum(unlist(perc))]
setorder(perc_mkt_share_idio,Qtl_Tot_vol)

#mean and SD

tot_ret_idio <- Tot_Vol_ret_idio %>% spread(Qtl_Tot_vol,Qtl_Tot_Vol_Rets)
tot_ret_mat_idio <- matrix(unlist(tot_ret_idio[,-1]),ncol=5)
100*colMeans(tot_ret_mat_idio)
apply(tot_ret_mat_idio,2,sd)

capmalpha_idio=c()
ffalpha_idio=c()
capmalphat_idio=c()
ffalphat_idio=c()
for(i in 1:5){
  capmalpha_idio[i]=lm(tot_ret_mat_idio[,i]~ffdata$`Mkt-RF`)$coefficients[1]*100 #converting to %
  ffalpha_idio[i]=lm(tot_ret_mat_idio[,i]~ffdata$`Mkt-RF`+ffdata$SMB+ffdata$HML)$coefficients[1]*100 #converting to %
  ct=lmSumm(lm(tot_ret_mat_idio[,i]~ffdata$`Mkt-RF`),HAC=TRUE) #fifth element is the t-stat of intercept
  capmalphat_idio[i]=ct$coef.table[5]
  ft=lmSumm(lm(tot_ret_mat_idio[,i]~ffdata$`Mkt-RF`+ffdata$SMB+ffdata$HML),HAC=TRUE)
  ffalphat_idio[i]=ft$coef.table[9]
}

outputA=data.table(Rank=c(1:5),Mean=100*colMeans(tot_ret_mat),Std.Dev=100*apply(tot_ret_mat,2,sd),
                   "%Mkt Share"=perc_mkt_share$percent*100,Size=avg_size_tot_vol$size,"CAPM Alpha"=capmalpha,
                   "FF-3 Alpha"=ffalpha,"t CAPM Alpha"=capmalphat,"t FF-3 Alpha"=ffalphat)

outputB=data.table(Rank=c(1:5),Mean=100*colMeans(tot_ret_mat_idio),Std.Dev=100*apply(tot_ret_mat_idio,2,sd),
                   "%Mkt Share"=perc_mkt_share_idio$percent*100,Size=avg_size_tot_vol_idio$size,
                   "CAPM Alpha"=capmalpha_idio,"FF-3 Alpha"=ffalpha_idio,
                   "t CAPM Alpha"=capmalphat_idio,"t FF-3 Alpha"=ffalphat_idio)

outA=xtable(outputA)
outB=xtable(outputB)

print.xtable(outA,digits=c(0,rep(2,8)))
print.xtable(outB,digits=c(0,rep(2,8)))