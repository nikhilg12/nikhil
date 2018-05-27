
#2.2 To Calculate Historical 99% 1 day VaR on each day of 2008 and backtesting at the end of 2008

library(quantmod)

#getting the data for CitiGroup ("C") from Yahoo! Finance
data=getSymbols("C",src="yahoo",from="2006-01-01",to="2008-12-31")

#Stock Price
plot(C$C.Close)

#Calculate Daily Return
C$DR=dailyReturn(C$C.Close)
plot(C$DR)


#the 99% VaR through historical daily returns is given by:

var99=c()
dr0607=C$DR["2006::2007"] #take the daily return values in 2006-2007
dr2008=C$DR["2008"] #take the values in 2008
pos2008jan=which(C$DR==dr2008[1])
for(i in 1:length(dr2008)){
  returns=c(dr0607,dr2008[1:i])
  #returns=returns[-last(returns)]
  var99[i]=quantile(returns,prob=0.01)
}
plot(var99,type='l')

#merging for plotting VaR with returns
varandreturns=merge(dr2008,var99)
plot(varandreturns$DR,type='l')
lines(varandreturns$var99,col="red")

exceptions=length(which(varandreturns$DR<varandreturns$var99))
#Goes out of bound a lot more than 1% of the time


#2.3 To calculate VaR of portfolio with 1 million in Goldman Sachs, JP Morgan,Barclays,Deutsche Bank, BNP Paribas
# and 2 million in UBS, Citigroup, Morgan Stanley, Bank of America, and Credit Suisse

# also, calculate the marginal VaR and the component VaR for each bank

tickers=c("GS","UBS","JPM","C","BCS","MS","DB","BAC","BNPQY","CS")
data2=getSymbols(tickers
                 ,src="yahoo",from="2006-01-01",to="2008-12-31")

data2=merge(GS,UBS,JPM,C,BCS,MS,DB,BAC,BNPQY,CS)

#calculate the daily returns
data2$GS.DR=dailyReturn(data2$GS.Adjusted)
data2$UBS.DR=dailyReturn(data2$UBS.Adjusted)
data2$JPM.DR=dailyReturn(data2$JPM.Adjusted)
data2$C.DR=dailyReturn(data2$C.Adjusted)
data2$BCS.DR=dailyReturn(data2$BCS.Adjusted)
data2$MS.DR=dailyReturn(data2$MS.Adjusted)
data2$DB.DR=dailyReturn(data2$DB.Adjusted)
data2$BAC.DR=dailyReturn(data2$BAC.Adjusted)
data2$BNPQY.DR=dailyReturn(data2$BNPQY.Adjusted)
data2$CS.DR=dailyReturn(data2$CS.Adjusted)

#create a vector of portfolio weights (1 million in odd positions, 2 million in even)
pfwts=rep(c(1e6,2e6),5)
pfwts=pfwts/sum(pfwts)

#extract a matrix of the daily returns for all 10 companies for all the dates
dailyretspf=as.matrix(data2[,((ncol(data2)-9):ncol(data2))])

#calculate the portfolio daily returns
dailyreturns=dailyretspf%*%pfwts

#the VaR of the portfolio
var99pf=sum(rep(c(1e6,2e6),5))-(1+quantile(dailyreturns,prob=0.01))*sum(rep(c(1e6,2e6),5))
var99pf


dvar=function(weights,value){
  dailrets=dailyretspf%*%weights
  varnew=value-(1+quantile(dailrets,prob=0.01))*value
  (varnew-var99pf) # increase in VaR due to 1$ increase in a position
}

cvar=function(weights,value){
  dailrets=dailyretspf%*%weights
  varnew=value-(1+quantile(dailrets,prob=0.01))*value
  (varnew-var99pf)/0.01 # increase in VaR due to 1% increase in a position
}

#calculate DVaR
dvars=c()
for(i in 1:10){
  wts=c(rep(c(1e6,2e6),5))
  wts[i]=wts[i]+1 #1$ increase in one of the position
  val=sum(wts)
  wts=wts/val
  dvars[i]=dvar(wts,val)
}

#Calculate CVaR
cvars=c()
for(i in 1:10){
  wts=c(rep(c(1e6,2e6),5))
  wts[i]=wts[i]*1.01 #1% increase in one of the positions
  val=sum(wts)
  wts=wts/val
  cvars[i]=cvar(wts,val)
}

output=cbind(dvars,cvars)
rownames(output)=tickers

#sum of CVaRs
sum(output[,2])

#Portfolio VaR
var99pf
