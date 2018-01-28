
library(readxl)
library(httr)
packageVersion("readxl")
url="https://github.com/nikhilg12/nikhil/raw/master/sixstocks_corpfin.xlsx"
GET(url, write_disk(tf <- tempfile(fileext = ".xlsx")))
stocks=read_excel(tf)


View(stocks)
msft=stocks[which(stocks[,3]=="MICROSOFT CORP"),]
msft[31,4]=NA #making the dividend paid in 2010 as NA
mun=stocks[which(stocks[,3]=="MANCHESTER UNITED PLC"),]
apple=stocks[which(stocks[,3]=="APPLE INC"),]
intel=stocks[which(stocks[,3]=="INTEL CORP"),]
jnj=stocks[which(stocks[,3]=="JOHNSON & JOHNSON"),]
nike=stocks[which(stocks[,3]=="NIKE INC"),]

library("dplyr")

div_plotter=function(x){
x=cbind(pull(x,3),pull(x,2),pull(x,4),pull(x,6),pull(x,7))
decl=x[which(!is.na(x[,3])),3] #all dates when dividends were announced

pos=matrix(nrow=length(decl),ncol=41)

for(i in 1:length(decl)){
  p=which(x[,2]==decl[i]) #finding the 0th day on which dividend was announced
  if(p<20){
    p=21 #lower limit of p (might induce distortions)
  }
  if(p>(nrow(x)-20)){
    p=nrow(x)-20 #upper limit of p (might induce distortions)
  }
  pos[i,]=as.numeric(x[((p-20):(p+20)),4]) - as.numeric(x[((p-20):(p+20)),5])  #taking the excess stock return
}

div=100*colMeans(pos) #taking the mean for each relative day and converting to %
barplot(div,main=paste("Average returns 20 days before dividend declaration to 20 days after for:",x[1,1]))
}
div_plotter(intel)
