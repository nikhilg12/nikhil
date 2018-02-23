library(quantmod)
library(httr)
library(plyr)
library(dplyr)

url1="https://github.com/nikhilg12/nikhil/raw/master/corpfinhw6/S%26P%20Stocks%20Financial%20Data.csv"
GET(url1, write_disk(tf1 <- tempfile(fileext = ".csv")))
stocksdata=read.csv(tf1)

tickers=as.vector(as.character(unique(stocksdata$tic)))

getpriceswithNA=function(ticker){
  tryCatch((getQuote(ticker,src="yahoo",what=standardQuote())$Last), error=function(err) NA)
} #trycatch function will return NA if the price is not available

latestprices=sapply(tickers,getpriceswithNA)
names(latestprices)=tickers

stocksdata_table=as.data.table(stocksdata)
latesteps=stocksdata_table[,.(epslatest=last(epsfx)),by=tic] #group by tic and then take the last value of eps for each tic
epratio=latesteps$epslatest/latestprices

getlaggednetincome_3yearavg=function(netincome){
  netincome=netincome[!is.na(netincome)]
  l=length(netincome)
  if(l>3){
    mean(diff(netincome)[(l-3):(l-1)]/netincome[(l-2):l]) #taking mean of last three values of the net income
  } else{
    return(as.double(NA))
  }
}

perc_lagged_ni_3yravg=stocksdata_table[,.(AvgLaggedEarnings=getlaggednetincome_3yearavg(ni)),by=tic]

length(epratio)
nrow(perc_lagged_ni_3yravg)
plotdata=cbind(perc_lagged_ni_3yravg,epratio)
plot(pull(plotdata[,2])*100,pull(plotdata[,3]),xlab="3-year average of the lagged change in percent earnings (%)"
     ,ylab="E/P",main="All stocks")

ss=subset(plotdata,plotdata$AvgLaggedEarnings>-1 & plotdata$AvgLaggedEarnings<1)
nrow(ss)
plot(pull(ss[,2])*100,pull(ss[,3]),xlab="3-year average of the lagged change in percent earnings(%)"
     ,ylab="E/P",main="Stocks with x axis only from -100% to +100%")
