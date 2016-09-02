install.packages("tmvtnorm")
install.packages("matrixStats")
install.packages("PerformanceAnalytics")
install.packages("ggplot")
library("tmvtnorm")
library("MASS")
library("matrixStats")
library('nloptr')
library("PerformanceAnalytics")
library("ggplot2")
library("plyr")

data<-read.table("dataset20142016.txt")#,header=F)
data<-apply(data,2,as.numeric) #508*11
DATA<-matrix(0,nrow=508,ncol=21)
DATA[,1]<-seq(1,508,by=1)
DATA[,2:11]<-data[,2:11]
DATA[2:508,12:21]<-t(sapply(2:508,function(x){data[x,-1]/data[(x-1),-1]-rep(1,10)}))
DATA.train<-DATA[1:409,]
dim(DATA.train)#409*21
DATA.test<-DATA[409:508,]
dim(DATA.test)#100*21
#ftse 100 iondex as risk free rate 
ftse.price.data<-read.table("bench100date.txt")
ftse.price<-as.numeric(ftse.price.data[,2])
risk.free.rate<-ftse.price/ftse.price[1]-rep(1,100)

#risk.free.rate.seires<-data.frame(risk.free.rate,row.names=ftse.price.data[,1])
risk.free.rate.seires<-data.frame(risk.free.rate,row.names=seq(as.Date('2011-01-01'),as.Date('2011-04-10'),by = 1))
#------------------------------------------------------------------------------------------------------#
#-------------------------------------parameter estimation: EM algortihm--------------------------------#
#-----particle filter SIR-----#
########prepare function 1: SIR###################
SIR<-function(Rho,Tau,Sigma,N,T1,ass){
  #result sample
  res1<-numeric()#E[(x_n)^2]#-----#sufficient stats 1, 3
  res2<-numeric()#E[(x_n)*(x_n_1)]#-----#sufficient statistic 2
  res2[1]<-0
  res3<-numeric()#E[(y_n-x_n)^2]#-----#sufficient statistic 4
  #matrix for x and y
  mx<-matrix(0,nrow=(T1+1),ncol=N)
  #matrix for w and W
  my<-matrix(0,nrow=(T1+1),ncol=N)
  mw<-matrix(0,nrow=(T1+1),ncol=N)
  mW<-matrix(0,nrow=(T1+1),ncol=N)
  #initial values of x
  mx[1,]<-rnorm(N,mean=0,sd=1)
  #values of y synthesised in the beginning
  my[1,]<-rep(DATA.train[2,ass+11],N)
  for(i in 1:N){
    mw[1,i]<-dnorm(my[1,i],mean=mx[1,i],sd=Sigma)
  }
  mW[1,]<-mw[1,]/sum(mw[1,])
  ####resample
  mx[1,]<-sample(mx[1,], N, mW[1,],replace=T)
  res1[1]<-sum((mx[1,])^2)/N
  res3[1]<-sum((my[1,]-mx[1,])^2)/N
  for(j in 1:T1){
    #generate new x
    mx[j+1,]<-Rho*mx[j,]+Tau*rnorm(N)
    #use y values synthesised at the beginning
    my[j+1,]<-rep(DATA.train[j+2,ass+11],N)
    
    for(k in 1:N){
      mw[j+1,k]<- dnorm(my[j+1,k],mean=mx[j+1,k],sd=Sigma)
    }
    #normalize weights to sum to 1
    mW[j+1,]<-mw[j+1,]/sum(mw[j+1,])
    mx[j+1,]<-sample(mx[j+1,], N, mW[j+1,],replace=T)
    res1[j+1]<-sum((mx[j+1,])^2)/N
    res2[j+1]<-sum(mx[j+1,]*mx[j,])/N
    res3[j+1]<-sum((my[j+1,]-mx[j+1,])^2)/N
  }
  list(res1,res2,res3,mw)
}

##
em<-function(N,T,int,asset){
  loglik<-numeric()
  marginal.log.lik<-numeric()
  rho.res<-numeric()
  tau2.res<-numeric()
  sigma2.res<-numeric()
  rho.res[1]<-1
  tau2.res[1]<-0.0005
  sigma2.res[1]<-0.0005
  for (i in 1:int){
    sir.res<-SIR(rho.res[i],sqrt(tau2.res[i]),sqrt(sigma2.res[i]),N,T,asset)
    s1<-sum(sir.res[[1]])
    s2<-sum(sir.res[[2]])
    s3<-sum(sir.res[[1]][-(T+1)])
    s4<-sum(sir.res[[3]])
    loglik[i]<- -(T+1)/2*log(4*(pi^2)*sigma2.res[i]*tau2.res[i])-1/(2*tau2.res[i])*s1+rho.res[i]/tau2.res[i]*s2-rho.res[i]^2*s3-1/(2*sigma2.res[i])*s4
    rho.res[i+1]<-s2/s3
    tau2.res[i+1]<-(s1-2*s2*rho.res[i+1]+s3*(rho.res[i+1])^2)/(T+1)
    sigma2.res[i+1]<-s4/(T+1)
  }
  list(rho.res,tau2.res,sigma2.res)
}



#-------------------------------------------numerical estimation of parameters-------------------------------------#
set.seed(1998)
#for asset 1 #
estimate.model2.a1<-em(N=500,T=407,int=500,asset=1)
#asset 2
estimate.model2.a2<-em(N=500,T=407,int=500,asset=2)
#asset 3
estimate.model2.a3<-em(N=500,T=407,int=500,asset=3)
#asset 4
estimate.model2.a4<-em(N=500,T=407,int=500,asset=4)
#asset 5
estimate.model2.a5<-em(N=500,T=407,int=500,asset=5)
#asset 6
estimate.model2.a6<-em(N=500,T=407,int=500,asset=6)
#asset 7
estimate.model2.a7<-em(N=500,T=407,int=500,asset=7)
#asset 8
estimate.model2.a8<-em(N=500,T=407,int=500,asset=8)
#asset 9
estimate.model2.a9<-em(N=500,T=407,int=500,asset=9)
#asset 10
estimate.model2.a10<-em(N=500,T=407,int=500,asset=10)



#1000 particles
estimate.model2.a1.1000<-em(N=1000,T=407,int=500,asset=1)
estimate.model2.a2.1000<-em(N=1000,T=407,int=500,asset=2)
estimate.model2.a3.1000<-em(N=1000,T=407,int=500,asset=3)
estimate.model2.a4.1000<-em(N=1000,T=407,int=500,asset=4)
estimate.model2.a5.1000<-em(N=1000,T=407,int=500,asset=5)
estimate.model2.a6.1000<-em(N=1000,T=407,int=500,asset=6)
estimate.model2.a7.1000<-em(N=1000,T=407,int=500,asset=7)
estimate.model2.a8.1000<-em(N=1000,T=407,int=500,asset=8)
estimate.model2.a9.1000<-em(N=1000,T=407,int=500,asset=9)
estimate.model2.a10.1000<-em(N=1000,T=407,int=500,asset=10)







#parameter matrix#each column representing three parameters for each asset
pm.model2<-matrix(0,nrow=3,ncol=10)
for(j.model in 1:3){
  pm.model2[j.model,1]<-median(estimate.model2.a1.1000[[j.model]])
  pm.model2[j.model,2]<-median(estimate.model2.a2.1000[[j.model]])
  pm.model2[j.model,3]<-median(estimate.model2.a3.1000[[j.model]])
  pm.model2[j.model,4]<-median(estimate.model2.a4.1000[[j.model]])
  pm.model2[j.model,5]<-median(estimate.model2.a5.1000[[j.model]])
  pm.model2[j.model,6]<-median(estimate.model2.a6.1000[[j.model]])
  pm.model2[j.model,7]<-median(estimate.model2.a7.1000[[j.model]])
  pm.model2[j.model,8]<-median(estimate.model2.a8.1000[[j.model]])
  pm.model2[j.model,9]<-median(estimate.model2.a9.1000[[j.model]])
  pm.model2[j.model,10]<-median(estimate.model2.a10.1000[[j.model]])
}










#--------------------------------------------------generate future stock price-----------------------------------------------#
#function for generating one sample of 10 asset over the following 100 time periods
#Tfu is the future time periods up to which we want to predict the return for
#res1.return is a matrix, where each col represent the predciting return for each aseet. It is a time series for one specific asset
pred.return.model2<-function(Tfu){
  res1.return<-matrix(0,nrow=Tfu,ncol=10)
  res1.price<-matrix(0,nrow=Tfu,ncol=10)
  for (as in 1:10){
    x.real<-numeric()
    #x0<-rnorm(0,1)
    x0<-DATA.train[409,as+11]-sqrt(pm.model2[3,as])*rnorm(1)
    x.real[1]<-x0
    res1.return[1,as]<-DATA.train[409,as+11]
    res1.price[1,as]<-DATA.train[409,as+1]
    for(itfu in 2: (Tfu)){
      x.real[itfu]<-pm.model2[1,as]*x.real[itfu-1]+sqrt(pm.model2[2,as])*rnorm(1)
      res1.return[itfu,as]<-x.real[itfu]+sqrt(pm.model2[3,as])*rnorm(1)
      res1.price[itfu,as]<-(1+res1.return[itfu,as])*res1.price[(itfu-1),as]
    }
  }
  return(list("predicted return"=res1.return,"predicted price"=res1.price))
}

set.seed(998)
price.list.model2<-list()
for(jj in 1:1000){
  price.list.model2[[jj]]<-pred.return.model2(100)[[2]]
}

#----choose one matrix in the result list as predicted price----#
#now it is a function, input x indicates which sample it is. Totally 1000 samples.
Predict.price.model2<-function(x){
  price.list.model2[[x]]
}

#----mean and variance prepared for return----#
Mean.price.model2<-Reduce("+", price.list.model2) / length(price.list.model2)
price.var.list.model2<-list()
for(jjj in 1:1000){
  price.var.list.model2[[jjj]]<-(price.list.model2[[jjj]]-Mean.price.model2)^2
}
Var.price.model2<-Reduce("+", price.var.list.model2) / (length(price.var.list.model2)-1)

#-----------------------------perforamnce measure--------------------------------------------#
#######################################################################
rb<-diff(ftse.price)/ftse.price[1:99]
RB<-data.frame(rb,row.names=seq(as.Date('2016-03-24'),as.Date('2016-06-30'),by = 1))

performance<-function(Y){
  SR<-c();ADD<-c();MDD<-c();NDD<-c();MT<-c();LS<-c();RE<-c();UI<-c()
  for(i in 1:1000){
    y<-Y[i,]
    ret<-diff(y)/y[1:99]
    RET<-data.frame(ret,row.names=seq(as.Date('2016-03-24'),as.Date('2016-06-30'),by = 1))
    SR[i]<-SharpeRatio(RET,RB,FUN="StdDev")
    ADD[i]<-AverageDrawdown(RET)
    MDD[i]<-maxDrawdown(RET)
    FD<-findDrawdowns(RET)
    NDD[i]<-length(FD$from)
    MT[i]<-max(FD$to-FD$from)
    #final return wrt initial wealth 100000
    RE[i]<-(y[100]-100000)/100000
    UI[i]<-UlcerIndex(RET)
    #MAR[i]<-MartinRatio(RET,RB,FUN="StdDev")
  }
  return(list("SR"=mean(SR),"MDD"=mean(MDD),"ADD"=mean(ADD),"NDD"=mean(NDD),"MT"=mean(MT),"R"=mean(RE),"V"=sd(RE),"UI"<-mean(UI)))
}
