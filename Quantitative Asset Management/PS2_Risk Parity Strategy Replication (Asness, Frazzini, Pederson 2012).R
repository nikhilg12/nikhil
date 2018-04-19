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

PS2_Q1=function(bonddata){

  bonddata$MCALDT=mdy(bonddata$MCALDT) #setting the date format
  
  #first return of every date is -99
  bonddata[which(bonddata$TMRETNUA==-99),"TMRETNUA"]=NA
  
  #next I calculate the lagged market equity
  bonddata=bonddata[,TMTOTOUTLAG:=lapply(.SD,shift),.SDcols="TMTOTOUT",by="KYCRSPID"]  
  #take complete cases
  
  bonddata[which(is.na(bonddata$TMTOTOUTLAG)),"TMTOTOUTLAG"]=0
  bonddata[which(is.na(bonddata$TMRETNUA)),"TMRETNUA"]=0

  #function which calculates value weights given two vectors: market capitalization and returns
  vwrets=function(rets,mktcaps){
    vwts=mktcaps/sum(mktcaps)
    vw=sum(vwts*rets)
    vw
  }
  
  ews=bonddata[,lapply(.SD,mean),by="MCALDT",.SDcols="TMRETNUA"] #apply mean function to column TMRETNUA and group by date
  ews=ews[order(MCALDT)]
  
  mktlag=bonddata[,lapply(.SD,sum),by="MCALDT",.SDcols="TMTOTOUTLAG"] #apply sum function to TMTOTOUTLAG and group by date
  mktlag=mktlag[order(MCALDT)]
  
  
  vws=bonddata[,lapply(.SD,vwrets,TMTOTOUTLAG),by="MCALDT",.SDcols="TMRETNUA"]
  vws=vws[order(MCALDT)]
  
  #merging the data together
  finaldata=merge(vws,ews,by="MCALDT")
  finaldata[,Bond_lag_MV:=mktlag$TMTOTOUTLAG]
  
  #adding more columns
  finaldata[,Year:=year(MCALDT)]
  finaldata[,Month:=month(MCALDT)]
  finaldata[,YM:=paste(year(MCALDT),month(MCALDT))]
  
  #setting column names
  setnames(finaldata,c("date","Bond_Vw_Ret","Bond_Ew_Ret","Bond_lag_MV","Year","Month","YM"))
  finaldata 
}

#Copying function PS1_Q1 from PS1

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


#reading the CRSP monthly data
CRSP_Stocks=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW1/dataps1.csv")
CRSP_Bonds=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW2/bonddata.csv")
Monthly_CRSP_Stocks=PS1_Q1(CRSP_Stocks)
Monthly_CRSP_Bonds=PS2_Q1(CRSP_Bonds)
Monthly_CRSP_Riskless=fread("D:/MSFE/MFE/Q3/Quantitative Asset Management/HW2/riskfree.csv")


PS2_Q2=function(stockdata,q1outbond,risklessdata){
  
  risklessdata$caldt=ymd(risklessdata$caldt)
  risklessdata$YM=paste(year(risklessdata$caldt),month(risklessdata$caldt))
  
  #removing NAs, (if any)
  risklessdata=risklessdata[complete.cases(risklessdata$t90ret,risklessdata$t30ret),]
  
  #merging all three datasets
  combineddata=merge(stockdata,q1outbond,by="YM")
  combineddata=merge(combineddata,risklessdata,by="YM")
  
  outdata=data.table(Date=combineddata$date.x,
                     Year=combineddata$Year.x,
                     Month=combineddata$Month.x,
                     Stock_lag_MV=combineddata$Stock_lag_MV/1000,  #converting to millions
                     Stock_Excess_Vw_Ret=combineddata$Stock_Vw_Ret-combineddata$t30ret,
                     Bond_lag_MV=combineddata$Bond_lag_MV,
                     Bond_Excess_Vw_Ret=combineddata$Bond_Vw_Ret-combineddata$t30ret)
  outdata
}

Monthly_CRSP_Universe=PS2_Q2(Monthly_CRSP_Stocks,Monthly_CRSP_Bonds,Monthly_CRSP_Riskless)

PS2_Q3=function(universedata){
  #excess value weighted return
  universedata$Excess_Vw_Ret = (universedata$Stock_lag_MV*universedata$Stock_Excess_Vw_Ret + 
                                  universedata$Bond_lag_MV*universedata$Bond_Excess_Vw_Ret)/(
                                    universedata$Stock_lag_MV+universedata$Bond_lag_MV)
  #excess 60/40 return
  universedata$Excess_60_40_Ret = 0.6*universedata$Stock_Excess_Vw_Ret + 0.4*universedata$Bond_Excess_Vw_Ret
  
  setorder(universedata,Year,Month)
  
  #apply sd function at every 36 months starting from 1 to end
  universedata$sigma_inverse_stocks=shift(rollapply(universedata$Stock_Excess_Vw_Ret,36,sd,fill=NA,align="right"))^(-1)
  universedata$sigma_inverse_bonds=shift(rollapply(universedata$Bond_Excess_Vw_Ret,36,sd,fill=NA,align="right"))^(-1)
  universedata$sigma_inverse_vw=shift(rollapply(universedata$Excess_Vw_Ret,36,sd,fill=NA,align="right"))^(-1)
  
  #calculating k_unlevered
  universedata$k_unlevered=1/(universedata$sigma_inverse_stocks+universedata$sigma_inverse_bonds)
  
  #calculating unlevered excess returns as a weighted return
  universedata$unlevered_returns=(universedata$sigma_inverse_stocks*universedata$Stock_Excess_Vw_Ret + 
      universedata$sigma_inverse_bonds*universedata$Bond_Excess_Vw_Ret)*universedata$k_unlevered
  
  #calculating k_levered
  universedata$k_lev=sd(universedata$Excess_Vw_Ret)/sd(universedata$sigma_inverse_stocks*universedata$Stock_Excess_Vw_Ret + 
                                                       universedata$sigma_inverse_bonds*universedata$Bond_Excess_Vw_Ret,na.rm=TRUE)
  #calculating levered excess returns as a weighted return
  universedata$lev_returns=(universedata$sigma_inverse_stocks*universedata$Stock_Excess_Vw_Ret + 
                                    universedata$sigma_inverse_bonds*universedata$Bond_Excess_Vw_Ret)*universedata$k_lev
  

  #output required in the question
  outputQ3=universedata[,-c(4,6,12)]
  
  colnames(outputQ3)=c("Date",
                       "Year",
                       "Month",
                       "Stock_Excess_Vw_Ret",
                       "Bond_Excess_Vw_Ret",
                       "Excess_Vw_Ret",
                       "Excess_60_40_Ret",
                        "Stock_inverse_sigma_hat",
                        "Bond_inverse_sigma_hat",
                        "Unlevered_k",
                        "Excess_Unlevered_RP_Ret",
                        "Levered_k",
                        "Excess_Levered_RP_Ret")
  outputQ3
}

Port_Rets=PS2_Q3(Monthly_CRSP_Universe)

PS2_Q4=function(universedata){
  #converting to xts so it is easy to select the desired period
  universexts=xts(universedata,order.by = universedata$Date)["1929-01::2010-06"]
  
  output=matrix(c(mean(as.numeric(universexts$Stock_Excess_Vw_Ret))*1200,
                  mean(as.numeric(universexts$Bond_Excess_Vw_Ret))*1200,
                  mean(as.numeric(universexts$Excess_Vw_Ret))*1200,
                  mean(as.numeric(universexts$Excess_60_40_Ret))*1200,
                  mean(as.numeric(universexts$Excess_Unlevered_RP_Ret,na.rm=TRUE))*1200,
                  mean(as.numeric(universexts$Excess_Levered_RP_Ret,na.rm=TRUE))*1200,
                  t.test(as.numeric(universexts$Stock_Excess_Vw_Ret))$statistic,
                  t.test(as.numeric(universexts$Bond_Excess_Vw_Ret))$statistic,
                  t.test(as.numeric(universexts$Excess_Vw_Ret))$statistic,
                  t.test(as.numeric(universexts$Excess_60_40_Ret))$statistic,
                  t.test(as.numeric(universexts$Excess_Unlevered_RP_Ret))$statistic,
                  t.test(as.numeric(universexts$Excess_Levered_RP_Ret))$statistic,
                  sd(as.numeric(universexts$Stock_Excess_Vw_Ret))*sqrt(12)*100,
                  sd(as.numeric(universexts$Bond_Excess_Vw_Ret))*sqrt(12)*100,
                  sd(as.numeric(universexts$Excess_Vw_Ret))*sqrt(12)*100,
                  sd(as.numeric(universexts$Excess_60_40_Ret))*sqrt(12)*100,
                  sd(as.numeric(universexts$Excess_Unlevered_RP_Ret,na.rm=TRUE))*sqrt(12)*100,
                  sd(as.numeric(universexts$Excess_Levered_RP_Ret,na.rm=TRUE))*sqrt(12)*100,
                  mean(as.numeric(universexts$Stock_Excess_Vw_Ret))*12/
                    (sd(as.numeric(universexts$Stock_Excess_Vw_Ret))*sqrt(12)),
                  mean(as.numeric(universexts$Bond_Excess_Vw_Ret))*12/
                    (sd(as.numeric(universexts$Bond_Excess_Vw_Ret))*sqrt(12)),
                  mean(as.numeric(universexts$Excess_Vw_Ret))*12/
                    (sd(as.numeric(universexts$Excess_Vw_Ret))*sqrt(12)),
                  mean(as.numeric(universexts$Excess_60_40_Ret))*12/
                    (sd(as.numeric(universexts$Excess_60_40_Ret))*sqrt(12)),
                  mean(as.numeric(universexts$Excess_Unlevered_RP_Ret),na.rm = TRUE)*12/
                    (sd(as.numeric(universexts$Excess_Unlevered_RP_Ret),na.rm = TRUE)*sqrt(12)),
                  mean(as.numeric(universexts$Excess_Levered_RP_Ret),na.rm = TRUE)*12/
                    (sd(as.numeric(universexts$Excess_Levered_RP_Ret),na.rm = TRUE)*sqrt(12)),
                  skewness(as.numeric(universexts$Stock_Excess_Vw_Ret)),
                  skewness(as.numeric(universexts$Bond_Excess_Vw_Ret)),
                  skewness(as.numeric(universexts$Excess_Vw_Ret)),
                  skewness(as.numeric(universexts$Excess_60_40_Ret)),
                  skewness(as.numeric(universexts$Excess_Unlevered_RP_Ret,na.rm=TRUE)),
                  skewness(as.numeric(universexts$Excess_Levered_RP_Ret,na.rm=TRUE)),
                  kurtosis(as.numeric(universexts$Stock_Excess_Vw_Ret))-3,
                  kurtosis(as.numeric(universexts$Bond_Excess_Vw_Ret))-3,
                  kurtosis(as.numeric(universexts$Excess_Vw_Ret))-3,
                  kurtosis(as.numeric(universexts$Excess_60_40_Ret))-3,
                  kurtosis(as.numeric(universexts$Excess_Unlevered_RP_Ret,na.rm=TRUE))-3,
                  kurtosis(as.numeric(universexts$Excess_Levered_RP_Ret,na.rm=TRUE))-3),nrow=6)
  rownames(output)=c("CRSP Stocks","CRSP Bonds","Value-Weighted Portfolio","60/40 Portfolio",
                     "Unlevered RP","Levered RP")
  colnames(output)=c("Annualized Mean (%)","t-statistic","Annualized SD(%)",
                     "Annualized Sharpe Ratio","Skewness","Excess Kurtosis")
  output
}

Q4_out=PS2_Q4(Port_Rets)
