#Nikhil Guruji - Project 5 - Cohort 

#stock path generator function

stock_path_generator=function(s0,r,sig,paths,timesteps,T){
  
  stock_paths=matrix(nrow=paths,ncol=timesteps)
  delta=T/timesteps
  for(i1 in 1:(nrow(stock_paths)/2)){
    i=2*i1-1
    z_pos=rnorm(timesteps)
    z_neg=-z_pos #antithetic 
    stock_paths[i,] = s0*exp((r-0.5*sig^2)*T + sig*sqrt(delta)*cumsum(z_pos)) #simulating ith path
    stock_paths[(i+1),] = s0*exp((r-0.5*sig^2)*T + sig*sqrt(delta)*cumsum(z_neg)) #simulating (i+1)th path using antithetic values
  }
  
  stock_paths=cbind(s0,stock_paths) #making the first row as the stock price
}

#laguerre definition
laguerre_upto_k=function(x,k){
  f1=exp(-x/2)
  f2=f1*(1-x)
  f3=f1*(1 - 2*x + (0.5*x^2))
  f4=f1*(1 - 3*x + (1.5*x^2) - (x^3)/6)
  
  if(k==2){
    cbind(f1,f2)
  } else if(k==3){
    cbind(f1,f2,f3)
  } else if(k==4){
    cbind(f1,f2,f3,f4)
  }
}

hermite_upto_k=function(x,k){
  f1=1
  f2=2*x
  f3=4*(x^2) - 2
  f4=8*(x^3) - 12*x
  
  if(k==2){
    cbind(f1,f2)
  } else if(k==3){
    cbind(f1,f2,f3)
  } else if(k==4){
    cbind(f1,f2,f3,f4)
  }
}

monomials_upto_k=function(x,k){
  f1=1
  f2=x
  f3=x^2
  f4=x^3
  
  if(k==2){
    cbind(f1,f2)
  } else if(k==3){
    cbind(f1,f2,f3)
  } else if(k==4){
    cbind(f1,f2,f3,f4)
  }
}

#laguerre_upto_k(10,4)


#1
#using lecture notes method
american_put_price_lsmc=function(s0,r,sig,paths,timesteps,T,func,strike,k){
  stock=stock_path_generator(s0,r,sig,paths,timesteps,T)
  
  #defining index matrix
  index=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  
  #defining Y matrix (payoff)
  Y=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues=pmax((strike-stock),0)
  Y[,ncol(Y)]=exercisevalues[,ncol(stock)]
  
  index[which(Y[,ncol(Y)]>0),ncol(index)]=1
  #View(cbind(Y[,ncol(Y)],index[,ncol(index)]))
  
  delta=T/timesteps
  pb=txtProgressBar(min=0,max=(ncol(Y)-2),style=3)
  for(j in 1:(ncol(Y)-2)){

    #start from the end
    j_reverse=ncol(Y)-j
    
    #Y = ind1*exp(-rd)*ev1 + ind2*exp(-2rd)*ev2 + . . .
    Y[,j_reverse]=rowSums(as.matrix(index[,((j_reverse+1):ncol(index))]*exp(-r*(((j_reverse+1):ncol(index))-j_reverse)
                                                          *delta)*exercisevalues[,((j_reverse+1):ncol(index))]))
    
    #find all positions in the current column for which Y=0 (we might try exercising here, hoping to get something>0)
    exercise_positions=which(Y[,j_reverse]==0)
    
    #convert all index values to zero (for now)
    index[exercise_positions,j_reverse]=1
    
    #make all Ys equal to corresponding exercise values
    Y[exercise_positions,j_reverse]=exercisevalues[exercise_positions,j_reverse]
    
    #even after exercising, some of the values in the current column will remain zero if option is out of the money
    exercise_positions=which(Y[,j_reverse]==0)
    
    #we won't exercise here and change those index values back to 0 for which even after exercising we get 0
    index[exercise_positions,j_reverse]=0
    
    #intuitive but slower for loop
    
    
    #now take all the values for which option is in the money
    nodes=which(Y[,j_reverse]>0)
    ys=Y[nodes,j_reverse]
    
    #and corresponding stock price
    xs=stock[nodes,j_reverse]
    
    #choose the appropriate L function
    if(func=="L"){
      L_xs=laguerre_upto_k(xs,k)
    }else if(func=="H"){
      L_xs=hermite_upto_k(xs,k)
    }else if(func=="M"){
      L_xs=monomials_upto_k(xs,k)
    }
    
    A=matrix(nrow=k,ncol=k)
    A=t(L_xs)%*%(L_xs) #sum of all combinations of L function (summation Li(x)Lj(x))
     
    b=matrix(nrow=k,ncol=1)
    b=t(t(as.matrix(ys))%*%L_xs)
    
    a=solve(A,b,tol = 1e-300)
    ecv=rep(0,length=nrow(Y))
    ecv[nodes]=L_xs%*%a #get the continuation values
     
    ev=rep(0,length=nrow(Y))
    ev[nodes]=exercisevalues[nodes,j_reverse]
    
    index[which(ev>ecv),j_reverse]=1 #wherever exercise value > continuation value, make that index 1
    
    index[which(index[,j_reverse]==1),(j_reverse+1):ncol(index)]=0 #making all index values zero after 1 is seen
    
    #for loop to do the same operation as line above (slower)
    #for(i in 1:nrow(index)){
     # index[i,(which(index[i,]==1)[-1])]=0
    #}
    
    #put the maximum of exercise value and the value that was there before (discounted payoff from future)
    Y[,j_reverse]=pmax(ev,Y[,j_reverse])
    #View(cbind(ev,ecv,index[,j_reverse]))
    
    setTxtProgressBar(pb,j)
  }
  close(pb)
  
  
  #multiply index by exercise matrix to get only the exercise values where it is optimal
  #and discount by appropriate multiple of delta (col 1 will have 0 as multiple and so on)
  payoff_each_path=rowSums(exp(-r*seq(0,(ncol(Y)-1))*delta)*(index*exercisevalues))
  
  #take payoff by monte carlo

  p=mean(payoff_each_path)
  p
}




#a
laguerrek2=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="L",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="L",strike=40,k=2,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="L",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="L",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="L",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="L",strike=40,k=2,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(laguerrek2)=c("S0=36","S0=40","S0=44")
rownames(laguerrek2)=c("T=0.5","T=1","T=2")

laguerrek3=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="L",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="L",strike=40,k=3,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="L",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="L",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="L",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="L",strike=40,k=3,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(laguerrek3)=c("S0=36","S0=40","S0=44")
rownames(laguerrek3)=c("T=0.5","T=1","T=2")

laguerrek4=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="L",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="L",strike=40,k=4,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="L",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="L",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="L",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="L",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="L",strike=40,k=4,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(laguerrek4)=c("S0=36","S0=40","S0=44")
rownames(laguerrek4)=c("T=0.5","T=1","T=2")

#b
hermitek2=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="H",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="H",strike=40,k=2,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="H",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="H",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="H",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="H",strike=40,k=2,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(hermitek2)=c("S0=36","S0=40","S0=44")
rownames(hermitek2)=c("T=0.5","T=1","T=2")

hermitek3=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="H",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="H",strike=40,k=3,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="H",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="H",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="H",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="H",strike=40,k=3,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(hermitek3)=c("S0=36","S0=40","S0=44")
rownames(hermitek3)=c("T=0.5","T=1","T=2")

hermitek4=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="H",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="H",strike=40,k=4,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="H",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="H",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="H",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="H",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="H",strike=40,k=4,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(hermitek4)=c("S0=36","S0=40","S0=44")
rownames(hermitek4)=c("T=0.5","T=1","T=2")


#c

monomialk2=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="M",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="M",strike=40,k=2,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="M",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="M",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=2,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="M",strike=40,k=2,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="M",strike=40,k=2,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(monomialk2)=c("S0=36","S0=40","S0=44")
rownames(monomialk2)=c("T=0.5","T=1","T=2")

monomialk3=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="M",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="M",strike=40,k=3,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="M",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="M",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=3,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="M",strike=40,k=3,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="M",strike=40,k=3,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(monomialk3)=c("S0=36","S0=40","S0=44")
rownames(monomialk3)=c("T=0.5","T=1","T=2")


monomialk4=matrix(c(
  american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=1,func="M",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=36,r=0.06,sig=0.2,T=2,func="M",strike=40,k=4,paths=100000,timesteps = 200)  
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=1,func="M",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=40,r=0.06,sig=0.2,T=2,func="M",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=0.5,func="M",strike=40,k=4,paths=100000,timesteps = 126)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=1,func="M",strike=40,k=4,paths=100000,timesteps = 200)
  ,american_put_price_lsmc(s0=44,r=0.06,sig=0.2,T=2,func="M",strike=40,k=4,paths=100000,timesteps = 200))
  ,nrow=3,ncol=3)

colnames(monomialk4)=c("S0=36","S0=40","S0=44")
rownames(monomialk4)=c("T=0.5","T=1","T=2")


#2
#a

r=0.06
T_put=1
timesteps_put=100
t_put=0.2
stock_2a=stock_path_generator(s0=65,r=0.06,sig=0.2,paths=100000,timesteps=timesteps_put,T=T_put)
europut=exp(-r*T_put)*pmax((stock_2a[,(t_put*timesteps_put)]-stock_2a[,ncol(stock_2a)]),0)
europut_value=mean(europut)


#b
fwd_start_american_put_price_lsmc=function(s0,r,sig,paths,timesteps,T,func,k,t){

  stock=stock_path_generator(s0,r,sig,paths,timesteps,T)
  
  #defining index matrix
  index=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  
  strike=stock[,(t*timesteps)]
  #defining Y matrix (payoff)
  Y=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues[,(t*timesteps):ncol(exercisevalues)]=pmax((strike-stock[,(t*timesteps):ncol(exercisevalues)]),0)
  Y[,ncol(Y)]=exercisevalues[,ncol(stock)]
  
  index[which(Y[,ncol(Y)]>0),ncol(index)]=1
  #View(cbind(Y[,ncol(Y)],index[,ncol(index)]))
  
  delta=T/timesteps
  for(j_reverse in (ncol(Y)-1):(t*timesteps + 1)){ #exercising starts right after t=0.2
    #starting from the end
    
    #Y = ind1*exp(-rd)*ev1 + ind2*exp(-2rd)*ev2 + . . .
    Y[,j_reverse]=rowSums(as.matrix(index[,((j_reverse+1):ncol(index))]*exp(-r*(((j_reverse+1):ncol(index))-j_reverse)
                                                                            *delta)*exercisevalues[,((j_reverse+1):ncol(index))]))
    
    #find all positions in the current column for which Y=0 (we might try exercising here, hoping to get something>0)
    exercise_positions=which(Y[,j_reverse]==0)
    
    #convert all index values to zero (for now)
    index[exercise_positions,j_reverse]=1
    
    #make all Ys equal to corresponding exercise values
    Y[exercise_positions,j_reverse]=exercisevalues[exercise_positions,j_reverse]
    
    #even after exercising, some of the values in the current column will remain zero if option is out of the money
    exercise_positions=which(Y[,j_reverse]==0)
    
    #we won't exercise here and change those index values back to 0 for which even after exercising we get 0
    index[exercise_positions,j_reverse]=0
    
    
    #now take all the values for which option is in the money
    nodes=which(Y[,j_reverse]>0)
    ys=Y[nodes,j_reverse]
    
    #and corresponding stock price
    xs=stock[nodes,j_reverse]
    
    #choose the appropriate L function
    if(func=="L"){
      L_xs=laguerre_upto_k(xs,k)
    }else if(func=="H"){
      L_xs=hermite_upto_k(xs,k)
    }else if(func=="M"){
      L_xs=monomials_upto_k(xs,k)
    }
    
    A=matrix(nrow=k,ncol=k)
    A=t(L_xs)%*%(L_xs) #sum of all combinations of L function (summation Li(x)Lj(x))
    
    b=matrix(nrow=k,ncol=1)
    b=t(t(as.matrix(ys))%*%L_xs)
    
    a=solve(A,b,tol = 1e-30)
    ecv=rep(0,length=nrow(Y))
    ecv[nodes]=L_xs%*%a #get the continuation values
    
    ev=rep(0,length=nrow(Y))
    ev[nodes]=exercisevalues[nodes,j_reverse]
    
    index[which(ev>ecv),j_reverse]=1 #wherever exercise value > continuation value, make that index 1
    
    index[which(index[,j_reverse]==1),(j_reverse+1):ncol(index)]=0 #making all index values zero after 1 is seen
    
    #for loop to do the same operation as line above (slower)
    #for(i in 1:nrow(index)){
    # index[i,(which(index[i,]==1)[-1])]=0
    #}
    
    #put the maximum of exercise value and the value that was there before (discounted payoff from future)
    Y[,j_reverse]=pmax(ev,Y[,j_reverse])
    #View(cbind(ev,ecv,index[,j_reverse]))
    
  }

  #multiply index by exercise matrix to get only the exercise values where it is optimal
  #and discount by appropriate multiple of delta (col 1 will have 0 as multiple and so on)
  payoff_each_path=rowSums(exp(-r*seq(0,(ncol(Y)-1))*delta)*(index*exercisevalues))
  
  #take payoff by monte carlo

  p=mean(payoff_each_path)
  p
}
americanput=fwd_start_american_put_price_lsmc(s0=65,sig=0.2,r=0.06,k=3,func="M",t=0.2,timesteps=100000,paths=100,T=1)






#################################################################################################################
#################################################################################################################
#################################################################################################################

#here is a slower function that does the same thing for Q1. However, I observed that the 
# results are more accurate using this function (but only that it is slower)

slow_american_put_price_lsmc=function(s0,r,sig,paths,timesteps,T,func,strike,k){
  stock=stock_path_generator(s0,r,sig,paths,timesteps,T)
  
  #defining index matrix
  index=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  
  #defining Y matrix (payoff)
  Y=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues=matrix(0,nrow=nrow(stock),ncol=ncol(stock))
  exercisevalues=pmax((strike-stock),0)
  Y[,ncol(Y)]=exercisevalues[,ncol(stock)]
  
  index[which(Y[,ncol(Y)]>0),ncol(index)]=1
  #View(cbind(Y[,ncol(Y)],index[,ncol(index)]))
  
  delta=T/timesteps
  pb=txtProgressBar(min=0,max=(ncol(Y)-2),style=3)
  for(j in 1:(ncol(Y)-2)){
    #j=1
    j_reverse=ncol(Y)-j
    for(i in 1:nrow(Y)){
      if(length(which(index[i,]==1))>0){ #not considering out of the money paths
        exerc=which(index[i,]==1)[1] #taking the first value for which index=1
        Y[i,j_reverse]=index[i,exerc]*exp(-r*(exerc-j_reverse)*delta)*exercisevalues[i,exerc]
      }else {
        Y[i,j_reverse]=exercisevalues[i,j_reverse]
        if(Y[i,j_reverse]>0){
          index[i,j_reverse]=1
        }
      }
    }
    nodes=which(Y[,j_reverse]>0)
    ys=Y[nodes,j_reverse]
    xs=stock[nodes,j_reverse]
    if(func=="L"){
      L_xs=laguerre_upto_k(xs,k)
    }else if(func=="H"){
      L_xs=hermite_upto_k(xs,k)
    }else if(func=="M"){
      L_xs=monomials_upto_k(xs,k)
    }
    
    A=matrix(nrow=k,ncol=k)
    A=t(L_xs)%*%(L_xs) #sum of all combinations of L function
    
    b=matrix(nrow=k,ncol=1)
    b=t(t(as.matrix(ys))%*%L_xs)
    
    a=solve(A,b,tol = 1e-30)
    ecv=rep(0,length=nrow(Y))
    ecv[nodes]=L_xs%*%a
    
    ev=rep(0,length=nrow(Y))
    ev[nodes]=exercisevalues[nodes,j_reverse]
    
    index[which(ev>ecv),j_reverse]=1
    for(i in 1:nrow(index)){
      index[i,(which(index[i,]==1)[-1])]=0
    }
    Y[,j_reverse]=pmax(ev,ecv)
    #View(cbind(ev,ecv,index[,j_reverse]))
    
    setTxtProgressBar(pb,j)
  }
  close(pb)
  
  payoff_each_path=vector(length = nrow(Y))
  
  for(i in 1:length(payoff_each_path)){
    ind=which(index[i,]==1)
    if(length(ind)>0){
      payoff_each_path[i]=exp(-r*ind*delta)*exercisevalues[i,ind]
    }
  }
  p=mean(payoff_each_path)
  p
}



#################################################################################################################
#################################################################################################################
#################################################################################################################
