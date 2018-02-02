data=read.csv("D:/MSFE/MFE/Q2/Corporate Finance/HW3/monthlystockdata.csv")
data$Mktcap=data$PRC*data$SHROUT #price x shares outstanding

data$year=substr(data$date,7,10) #extracting year from the date
data$year=as.numeric(data$year) #converting year to numeric type
min(data$year) #verifying that the minimum year is 1925 and maximum is 2016
max(data$year)

permnos=unique(data[,1]) # extracting all the unique PERMNOS (which uniquely identify each firm)
data$RET=as.numeric(as.character(data$RET)) #converting the returns to numeric type (already in decimals)

list_firms=list()
for(j in 1:2000){
  one_firm=data[which(data[,1]==permnos[j]),]
  years=unique(one_firm$year) #finding the years for which returns are available
  one_firm=one_firm[which(!is.na(one_firm$RET)),] #removing missing values from months in which data is not available
  #(years are not removed)
  one_firm$RET=1+one_firm$RET
  ann_ret=c()
  for(i in 1:length(years)){
    ann_ret[i]=prod(one_firm[which(one_firm$year==years[i]),6])-1
  }
  ann_ret=100*ann_ret #to convert into %
  ret=cbind(years,ann_ret)
  colnames(ret)=c("Year","Annual Returns in %")
  rownames(ret)=c(ret[,1])
  list_firms[[j]]=ret
}

#any stock you want to see the barplot of enter the position in the permnos vector (from 1 to 2000)
# can run the for loop for all the 31000 something stocks. But, takes too long to run

pos_in_permnos=6 #returns from 1925
permnos[pos_in_permnos] #PERMNO of the stock whose returns I plot a few lines later 

pos_in_permnos=1789 #returns from 1972 to 2016
permnos[pos_in_permnos] #PERMNO of the stock whose returns I plot a few lines later 

tickername=unique(as.character(data[which(data[,1]==permnos[pos_in_permnos]),3])) #gets the ticker name of the stock

pl=barplot(list_firms[[pos_in_permnos]][,2],main=paste("% annual returns of ticker: ",tickername[length(tickername)]))
