install.packages("coda")
install.packages("stochvol")

library("tmvtnorm")
library("MASS")
library("matrixStats")
library('nloptr')
library("PerformanceAnalytics")
library("ggplot2")
library("plyr")
library("coda")
library("stochvol")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

data<-read.table("dataset20142016.txt",header=F)
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
ftse.price.data<-read.table("bench100date.txt")#,header=F)
ftse.price<-as.numeric(ftse.price.data[,2])
#risk.free.rate<-ftse.price/ftse.price[1]-rep(1,100)


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


#------------asset 1----------------#
ret.asset1 <- logret(DATA.train[-1,2], demean = TRUE)
res.asset1 <- svsample(ret.asset1, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset1, showlatent = FALSE)
#------------asset 2----------------#
ret.asset2 <- logret(DATA.train[-1,3], demean = TRUE)
res.asset2 <- svsample(ret.asset2, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset2, showlatent = FALSE)
#------------asset 3----------------#
ret.asset3 <- logret(DATA.train[-1,4], demean = TRUE)
res.asset3 <- svsample(ret.asset3, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset3, showlatent = FALSE)
#------------asset 4----------------#
ret.asset4 <- logret(DATA.train[-1,5], demean = TRUE)
res.asset4 <- svsample(ret.asset4, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset4, showlatent = FALSE)
#------------asset 5----------------#
ret.asset5 <- logret(DATA.train[-1,6], demean = TRUE)
res.asset5 <- svsample(ret.asset5, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset5, showlatent = FALSE)
#------------asset 6----------------#
ret.asset6 <- logret(DATA.train[-1,7], demean = TRUE)
res.asset6 <- svsample(ret.asset6, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset6, showlatent = FALSE)
#------------asset 7----------------#
ret.asset7 <- logret(DATA.train[-1,8], demean = TRUE)
res.asset7 <- svsample(ret.asset7, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset7, showlatent = FALSE)
#------------asset 8----------------#
ret.asset8 <- logret(DATA.train[-1,9], demean = TRUE)
res.asset8 <- svsample(ret.asset8, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset8, showlatent = FALSE)
#------------asset 9----------------#
ret.asset9 <- logret(DATA.train[-1,10], demean = TRUE)
res.asset9 <- svsample(ret.asset9, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset9, showlatent = FALSE)
#------------asset 10----------------#
ret.asset10 <- logret(DATA.train[-1,11], demean = TRUE)
res.asset10 <- svsample(ret.asset10, priormu = c(0, 100), priorphi = c(5, 1.5),  priorsigma = 1)
summary(res.asset10, showlatent = FALSE)


pm.model3<-matrix(0,nrow=3,ncol=10)
pm.model3[,1]<-c(-8.375,0.957,0.262)
pm.model3[,2]<-c(-8.891,0.876,0.239)
pm.model3[,3]<-c(-8.794,0.908,0.275)
pm.model3[,4]<-c(-8.950,0.813,0.300)
pm.model3[,5]<-c(-8.874,0.607,0.469)
pm.model3[,6]<-c(-9.041,0.885,0.579)
pm.model3[,7]<-c(-8.480,0.845,0.515)
pm.model3[,8]<-c(-8.615,0.637,0.492)
pm.model3[,9]<-c(-8.698,0.710,0.456)
pm.model3[,10]<-c(-8.953,0.586,0.225)

pred.return.model3<-function(Tfu){
  res1.return<-matrix(0,nrow=Tfu,ncol=10)
  res1.price<-matrix(0,nrow=Tfu,ncol=10)
  res1.price[1,]<-DATA.test[1,2:11]
  for (as in 1:10){
    res1.return[,as]<-svsim(Tfu, mu = pm.model3[1,as], phi = pm.model3[2,as], sigma = pm.model3[3,as])[[1]]
    for (period in 2:Tfu){
      res1.price[period,as]<-(1+res1.return[period,as])*res1.price[period-1,as]
    }
  }
  return(list("predicted return"=res1.return,"predicted price"=res1.price))
}

##----simulate x samples----representing as a list-MC#
set.seed(998)
price.list.model3<-list()
for(jj in 1:1000){
  price.list.model3[[jj]]<-pred.return.model3(100)[[2]]
}



Predict.price.model3<-function(x){
  price.list.model3[[x]]
}
#

#----mean and variance prepared for return----#
Mean.price.model3<-Reduce("+", price.list.model3) / length(price.list.model3)
price.var.list.model3<-list()
for(jjj in 1:1000){
  price.var.list.model3[[jjj]]<-(price.list.model3[[jjj]]-Mean.price.model3)^2
}
Var.price.model3<-Reduce("+", price.var.list.model3) / (length(price.var.list.model3)-1)


#-------------------------------------R1 S1----------------------------------------------#
ub<-function(x){
  rep(100000*0.5,10)/abs(Predict.price.model3(x)[1,])
}
lb<-function(x){
  -rep(100000*0.5,10)/abs(Predict.price.model3(x)[1,])
}


#----------------correct----------------------------------------------------------#
opt.initial.R1.S1<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*(x[i]^2)*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S1<-function(gam,w.prev,L.length,intervalind,samind){
  #end of the interval
  time<-L.length*intervalind
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_eq<-function(x){
    constr<-0
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    for(i in 1:10){
      if(x[i]>w.prev[i]){
        constr<-constr+(1+0.003)*(x[i]-w.prev[i])*pri[i]
      }
      else{constr<-constr+(1-0.003)*(x[i]-w.prev[i])*pri[i]
      }
    }
    g<-c()
    for(i in 1:10){
      if(x[i]>w.prev[i]){
        g[i]<-(1+0.003)*pri[i]
      }
      else{g[i]<-(1-0.003)*pri[i]}
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-runif(10,lb(samind),ub(samind))
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}


#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R1.S1<-function(ntime,gammm,intervallength,sam){
  #RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S1(gammm,100000,intervallength,sam)
  for(i in 2:ntime){
    res[i,]<- opt.consecutive.R1.S1(gammm,res[i-1,],intervallength,i,sam)
  }
  res
}
#----------------final answer for the optimised weight matrix
#just input 
AW.model3.R1.S1<-function(ga,totalinte,intelength,whichsam){
  #res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R1.S1(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      #res[i2,]<-res.inter[[1]][i1,]
      res[i2,]<-res.inter[i1,]
    }
  }
  res
}


#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R1.S1<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R1.S1(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[n,]-AW.sam[n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[x,]%*%DATA.test[x,2:11]})
    res.matrix[i,]<-res.matrix[i,]-tc
    print(i)
  }
  res.matrix
}

#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=10-----------------#
model3.L20.R1.S1.gam1<-AV.model3.R1.S1(L=20,paragamma=10,nosample=1000)
#--------------gamma=5-----------------#
model3.L20.R1.S1.gam2<-AV.model3.R1.S1(L=20,paragamma=5,nosample=1000)
#--------------gamma=1-----------------#
model3.L20.R1.S1.gam3<-AV.model3.R1.S1(L=20,paragamma=1,nosample=1000)
#--------------gamma=0.1-----------------#
model3.L20.R1.S1.gam4<-AV.model3.R1.S1(L=20,paragamma=0.1,nosample=1000)
#--------------gamma=0.05-----------------#
model3.L20.R1.S1.gam5<-AV.model3.R1.S1(L=20,paragamma=0.05,nosample=1000)
#--------------gamma=0.01-----------------#
model3.L20.R1.S1.gam6<-AV.model3.R1.S1(L=20,paragamma=0.01,nosample=1000)
#--------------gamma=0.005-----------------#
model3.L20.R1.S1.gam7<-AV.model3.R1.S1(L=20,paragamma=0.005,nosample=1000)
#--------------gamma=0.001-----------------#
model3.L20.R1.S1.gam8<-AV.model3.R1.S1(L=20,paragamma=0.001,nosample=1000)
#--------------gamma=0.0005-----------------#
model3.L20.R1.S1.gam9<-AV.model3.R1.S1(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model3.L20.R1.S1.gam10<-AV.model3.R1.S1(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model3.L20.R1.S1.gam11<-AV.model3.R1.S1(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model3.L20.R1.S1.gam12<-AV.model3.R1.S1(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model3.L20.R1.S1.gam13<-AV.model3.R1.S1(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model3.L20.R1.S1.gam14<-AV.model3.R1.S1(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model3.L20.R1.S1.gam15<-AV.model3.R1.S1(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model3.L20.R1.S1.gam16<-AV.model3.R1.S1(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model3.L20.R1.S1.gam17<-AV.model3.R1.S1(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model3.L20.R1.S1.gam18<-AV.model3.R1.S1(L=20,paragamma=0.00001,nosample=1000)
#----------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#
per.model3.L20.R1.S1.gam1<-performance(model3.L20.R1.S1.gam1)
per.model3.L20.R1.S1.gam2<-performance(model3.L20.R1.S1.gam2)
per.model3.L20.R1.S1.gam3<-performance(model3.L20.R1.S1.gam3)
per.model3.L20.R1.S1.gam4<-performance(model3.L20.R1.S1.gam4)
per.model3.L20.R1.S1.gam5<-performance(model3.L20.R1.S1.gam5)
per.model3.L20.R1.S1.gam6<-performance(model3.L20.R1.S1.gam6)
per.model3.L20.R1.S1.gam7<-performance(model3.L20.R1.S1.gam7)
per.model3.L20.R1.S1.gam8<-performance(model3.L20.R1.S1.gam8)
per.model3.L20.R1.S1.gam9<-performance(model3.L20.R1.S1.gam9)
per.model3.L20.R1.S1.gam10<-performance(model3.L20.R1.S1.gam10)
per.model3.L20.R1.S1.gam11<-performance(model3.L20.R1.S1.gam11)
per.model3.L20.R1.S1.gam12<-performance(model3.L20.R1.S1.gam12)
per.model3.L20.R1.S1.gam13<-performance(model3.L20.R1.S1.gam13)
per.model3.L20.R1.S1.gam14<-performance(model3.L20.R1.S1.gam14)
per.model3.L20.R1.S1.gam15<-performance(model3.L20.R1.S1.gam15)
per.model3.L20.R1.S1.gam16<-performance(model3.L20.R1.S1.gam16)
per.model3.L20.R1.S1.gam17<-performance(model3.L20.R1.S1.gam17)
per.model3.L20.R1.S1.gam18<-performance(model3.L20.R1.S1.gam18)

#--------------------------------------r1s2-------------------------
#----------------for the first time interval------------------------#
opt.initial.R1.S2<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*(x[i]^2)*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S2<-function(gam,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  #objective function
  eval_f<-function(x){
    dif.vec<-Predict.price.model3(samind)[time,]-Predict.price.model3(samind)[L.length*(intervalind-1),]
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]
    }
    g<-0
    for(i in 1:10){
      g<-g-dif.vec[i]*Mean.price.model3[time,i]+2*gam*dif.vec[i]*Var.price.model3[time,i]*(x*dif.vec[i]+1)
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    dif.vec<-Predict.price.model3(samind)[time,]-Predict.price.model3(samind)[L.length*(intervalind-1),]
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    constr<-0
    for(i in 1:10){
      constr<-constr+x*(dif.vec[i])*pri[i]+0.003*x*abs(dif.vec[i])*pri[i]
    }
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+0.003*abs(dif.vec[i])*pri[i]
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-0
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=0,
              ub=Inf,
              eval_g_ineq=eval_g_ineq,
              opts=opts)
  res$solution
}

#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R1.S2<-function(ntime,gammm,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S2(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  for(i in 2:ntime){
    dif.vec<-Predict.price.model3(sam)[i*intervallength,]-Predict.price.model3(sam)[(i-1)*intervallength,]
    inter<-opt.consecutive.R1.S2(gammm,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}

#----------------final answer for the optimised weight matrix
#just input 
AW.model3.R1.S2<-function(ga,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R1.S2(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}

#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R1.S2<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R1.S2(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    for(j in 2:(100/L)){
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*(DATA.test[j*L,2:11]-DATA.test[(j-1)*L,2:11])%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(DATA.test[j*L,2:11]-DATA.test[(j-1)*L,2:11])%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}


model3.L20.R1.S2.gam9<-AV.model3.R1.S2(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model3.L20.R1.S2.gam10<-AV.model3.R1.S2(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model3.L20.R1.S2.gam11<-AV.model3.R1.S2(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model3.L20.R1.S2.gam12<-AV.model3.R1.S2(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model3.L20.R1.S2.gam13<-AV.model3.R1.S2(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model3.L20.R1.S2.gam14<-AV.model3.R1.S2(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model3.L20.R1.S2.gam15<-AV.model3.R1.S2(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model3.L20.R1.S2.gam16<-AV.model3.R1.S2(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model3.L20.R1.S2.gam17<-AV.model3.R1.S2(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model3.L20.R1.S2.gam18<-AV.model3.R1.S2(L=20,paragamma=0.00001,nosample=1000)


per.model3.L20.R1.S2.gam9<-performance(model3.L20.R1.S2.gam9[[1]])
per.model3.L20.R1.S2.gam10<-performance(model3.L20.R1.S2.gam10[[1]])
per.model3.L20.R1.S2.gam11<-performance(model3.L20.R1.S2.gam11[[1]])
per.model3.L20.R1.S2.gam12<-performance(model3.L20.R1.S2.gam12[[1]])
per.model3.L20.R1.S2.gam13<-performance(model3.L20.R1.S2.gam13[[1]])
per.model3.L20.R1.S2.gam14<-performance(model3.L20.R1.S2.gam14[[1]])
per.model3.L20.R1.S2.gam15<-performance(model3.L20.R1.S2.gam15[[1]])
per.model3.L20.R1.S2.gam16<-performance(model3.L20.R1.S2.gam16[[1]])
per.model3.L20.R1.S2.gam17<-performance(model3.L20.R1.S2.gam17[[1]])
per.model3.L20.R1.S2.gam18<-performance(model3.L20.R1.S2.gam18[[1]])

#--------------------------------------------r1 s3--------------------------
opt.initial.R1.S3<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S3<-function(gam,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  price.mat<-Predict.price.model3(samind)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  dif.vec<-colMeans(dif.mat[(L.length*(intervalind-1)+1):time,])
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]
    }
    g<-0
    for(i in 1:10){
      g<-g-dif.vec[i]*Mean.price.model3[time,i]+2*gam*dif.vec[i]*Var.price.model3[time,i]*(x*dif.vec[i]+1)
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    constr<-0
    for(i in 1:10){
      #constr<-constr+x*(dif.vec[i])*pri[i]+0.003*x*abs(dif.vec[i])*pri[i]-(1+rfr)*cash
      constr<-constr+x*(dif.vec[i])*pri[i]+0.003*x*abs(dif.vec[i])*pri[i]
    }
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+0.003*abs(dif.vec[i])*pri[i]
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-0
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=0,
              ub=Inf,
              eval_g_ineq=eval_g_ineq,
              opts=opts)
  res$solution
}

opt.consecutive.R1.S3(gam=3,w.prev=tt1,L.length=20,intervalind=2,samind=1000)

#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R1.S3<-function(ntime,gammm,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S3(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  price.mat<-Predict.price.model3(sam)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  for(i in 2:ntime){
    dif.vec<-colMeans(dif.mat[(intervallength*(i-1)+1):(i*intervallength),])
    inter<-opt.consecutive.R1.S3(gammm,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}
#-----------------------------------------------------------------------------------------------------------------#
AW.model3.R1.S3<-function(ga,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R1.S3(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}

#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R1.S3<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  #
  price.mat<-DATA.test[,2:11]
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  #
  for (i in 1:nosample){
    AW.sam<-AW.model3.R1.S3(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    for(j in 2:(100/L)){
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*colMeans(dif.mat[((j-1)*L):(j*L),])%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(colMeans(dif.mat[((j-1)*L):(j*L),]))%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}

#------------------------------------------------L=20---------------------------------------------------#
#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=0.0005-----------------#
model3.L20.R1.S3.gam9<-AV.model3.R1.S3(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model3.L20.R1.S3.gam10<-AV.model3.R1.S3(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model3.L20.R1.S3.gam11<-AV.model3.R1.S3(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model3.L20.R1.S3.gam12<-AV.model3.R1.S3(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model3.L20.R1.S3.gam13<-AV.model3.R1.S3(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model3.L20.R1.S3.gam14<-AV.model3.R1.S3(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model3.L20.R1.S3.gam15<-AV.model3.R1.S3(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model3.L20.R1.S3.gam16<-AV.model3.R1.S3(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model3.L20.R1.S3.gam17<-AV.model3.R1.S3(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model3.L20.R1.S3.gam18<-AV.model3.R1.S3(L=20,paragamma=0.00001,nosample=1000)
#-------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#
per.model3.L20.R1.S3.gam9<-performance(model3.L20.R1.S3.gam9[[1]])
per.model3.L20.R1.S3.gam10<-performance(model3.L20.R1.S3.gam10[[1]])
per.model3.L20.R1.S3.gam11<-performance(model3.L20.R1.S3.gam11[[1]])
per.model3.L20.R1.S3.gam12<-performance(model3.L20.R1.S3.gam12[[1]])
per.model3.L20.R1.S3.gam13<-performance(model3.L20.R1.S3.gam13[[1]])
per.model3.L20.R1.S3.gam14<-performance(model3.L20.R1.S3.gam14[[1]])
per.model3.L20.R1.S3.gam15<-performance(model3.L20.R1.S3.gam15[[1]])
per.model3.L20.R1.S3.gam16<-performance(model3.L20.R1.S3.gam16[[1]])
per.model3.L20.R1.S3.gam17<-performance(model3.L20.R1.S3.gam17[[1]])
per.model3.L20.R1.S3.gam18<-performance(model3.L20.R1.S3.gam18[[1]])

#----------------------------------R1 S4
opt.initial.R1.S4<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*(x[i]^2)*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}


eps<-0.0005      
#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S4<-function(gam,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  price.mat<-Predict.price.model3(samind)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  dif.vec.indicator<-colMeans(dif.mat[(L.length*(intervalind-1)+1):time,])
  dif.vec<-c()
  for(j in 1:10){
    if(abs(dif.vec.indicator[j])>eps){dif.vec[j]<-dif.vec.indicator[j]}
    else{dif.vec[j]<-0}
  }
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]
    }
    g<-0
    for(i in 1:10){
      g<-g-dif.vec[i]*Mean.price.model3[time,i]+2*gam*dif.vec[i]*Var.price.model3[time,i]*(x*dif.vec[i]+1)
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    constr<-0
    for(i in 1:10){
      #constr<-constr+x*(dif.vec[i])*pri[i]+0.003*x*abs(dif.vec[i])*pri[i]-(1+rfr)*cash
      constr<-constr+x*(dif.vec[i])*pri[i]+0.003*x*abs(dif.vec[i])*pri[i]
    }
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+0.003*abs(dif.vec[i])*pri[i]
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  #
  x0<-0.1
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG","xtol_rel" = 1.0e-7,"maxeval" = 1000, "local_opts" = local_opts )
  res<-nloptr(x0=x0,eval_f=eval_f,lb=0, ub=Inf,eval_g_ineq=eval_g_ineq,opts=opts)
  res$solution
}

#------------------------------
RES.model3.R1.S4<-function(ntime,gammm,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S4(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  price.mat<-Predict.price.model3(sam)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  for(i in 2:ntime){
    dif.vec.indicator<-colMeans(dif.mat[(intervallength*(i-1)+1):(i*intervallength),])
    dif.vec<-c()
    for(j in 1:10){
      if(abs(dif.vec.indicator[j])>eps){dif.vec[j]<-dif.vec.indicator[j]}
      else{dif.vec[j]<-0}
    }
    inter<-opt.consecutive.R1.S4(gammm,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}

#----------------final answer for the optimised weight matrix
AW.model3.R1.S4<-function(ga,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R1.S4(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}
############################
AV.model3.R1.S4<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  #
  price.mat<-DATA.test[,2:11]
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  #
  for (i in 1:nosample){
    AW.sam<-AW.model3.R1.S4(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    for(j in 2:(100/L)){
      dif.vec.indicator<-colMeans(dif.mat[((j-1)*L):(j*L),])
      dif.vec<-c()
      for(j1 in 1:10){
        if(abs(dif.vec.indicator[j1])>eps){dif.vec[j1]<-dif.vec.indicator[j1]}
        else{dif.vec[j1]<-0}
      }
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*dif.vec%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(dif.vec)%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}
#------------------------------------------------L=20---------------------------------------------------#
#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=0.0005-----------------#
model3.L20.R1.S4.gam9<-AV.model3.R1.S4(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model3.L20.R1.S4.gam10<-AV.model3.R1.S4(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model3.L20.R1.S4.gam11<-AV.model3.R1.S4(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model3.L20.R1.S4.gam12<-AV.model3.R1.S4(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model3.L20.R1.S4.gam13<-AV.model3.R1.S4(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model3.L20.R1.S4.gam14<-AV.model3.R1.S4(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model3.L20.R1.S4.gam15<-AV.model3.R1.S4(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model3.L20.R1.S4.gam16<-AV.model3.R1.S4(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model3.L20.R1.S4.gam17<-AV.model3.R1.S4(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model3.L20.R1.S4.gam18<-AV.model3.R1.S4(L=20,paragamma=0.00001,nosample=1000)
#-------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#

per.model3.L20.R1.S4.gam9<-performance(model3.L20.R1.S4.gam9[[1]])
per.model3.L20.R1.S4.gam10<-performance(model3.L20.R1.S4.gam10[[1]])
per.model3.L20.R1.S4.gam11<-performance(model3.L20.R1.S4.gam11[[1]])
per.model3.L20.R1.S4.gam12<-performance(model3.L20.R1.S4.gam12[[1]])
per.model3.L20.R1.S4.gam13<-performance(model3.L20.R1.S4.gam13[[1]])
per.model3.L20.R1.S4.gam14<-performance(model3.L20.R1.S4.gam14[[1]])
per.model3.L20.R1.S4.gam15<-performance(model3.L20.R1.S4.gam15[[1]])
per.model3.L20.R1.S4.gam16<-performance(model3.L20.R1.S4.gam16[[1]])
per.model3.L20.R1.S4.gam17<-performance(model3.L20.R1.S4.gam17[[1]])
per.model3.L20.R1.S4.gam18<-performance(model3.L20.R1.S4.gam18[[1]])

#----------------------------------R2 S1
#----------------for the first time interval------------------------#
opt.initial.R2.S1<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R2.S1<-function(gam,theta,w.prev,L.length,intervalind,samind){
  #end of the interval
  time<-L.length*intervalind
  time_1<-L.length*(intervalind-1)+1
  #objective function 
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]+theta*(x[i]-w.prev[i])^2*Mean.price.model3[time_1,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]+2*theta*(x[i]-w.prev[i])*Mean.price.model3[time_1,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_eq<-function(x){
    constr<-0
    pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
    for(i in 1:10){
      constr<-constr+x[i]*pri[i]+0.003*abs(x[i]-w.prev[i])*pri[i]-w.prev[i]*pri[i]
    }
    g<-c()
    for(i in 1:10){
      if(x[i]-w.prev[i]>0){
        g[i]<-(1+0.003)*pri[i]
      }
      else{g[i]<-(1-0.003)*pri[i]}
    }
    return(list("constraints"=constr,"jacobian"=g))
  }

  x0<-runif(10,lb(samind),ub(samind))
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R2.S1<-function(ntime,gammm,thet,intervallength,sam){
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S1(gammm,100000,intervallength,sam)
  for(i in 2:ntime){
    res[i,]<- opt.consecutive.R2.S1(gammm,thet,res[i-1,],intervallength,i,sam)
  }
  res
}
#-----------------------------------------------------------------------------------------------------------------#
AW.model3.R2.S1<-function(ga,the,totalinte,intelength,whichsam){
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R2.S1(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[i1,]
    }
  }
  res
}

#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R2.S1<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R2.S1(paragamma,paramthet,100/L,L,i)
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[x,]%*%DATA.test[x,2:11]})
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[n,]-AW.sam[n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-res.matrix[i,]-tc
    print(i)
  }
  
  res.matrix
}



model3.L20.R2.S1.theta1<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=1,nosample=1000)
model3.L20.R2.S1.theta2<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
model3.L20.R2.S1.theta3<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
model3.L20.R2.S1.theta4<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
model3.L20.R2.S1.theta5<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
model3.L20.R2.S1.theta6<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=0,nosample=1000)
model3.L20.R2.S1.theta7<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model3.L20.R2.S1.theta8<-AV.model3.R2.S1(L=20,paragamma=0.1,paramthet=10,nosample=1000)


#---------------L=20--------------------#
per.model3.L20.R2.S1.theta1<-performance(model3.L20.R2.S1.theta1)
per.model3.L20.R2.S1.theta2<-performance(model3.L20.R2.S1.theta2)
per.model3.L20.R2.S1.theta3<-performance(model3.L20.R2.S1.theta3)
per.model3.L20.R2.S1.theta4<-performance(model3.L20.R2.S1.theta4)
per.model3.L20.R2.S1.theta5<-performance(model3.L20.R2.S1.theta5)
per.model3.L20.R2.S1.theta6<-performance(model3.L20.R2.S1.theta6)
per.model3.L20.R2.S1.theta7<-performance(model3.L20.R2.S1.theta7)
per.model3.L20.R2.S1.theta8<-performance(model3.L20.R2.S1.theta8)

#----------------------------R2S2
opt.initial.R2.S2<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }

  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG","xtol_rel" = 1.0e-7,"maxeval" = 1000, "local_opts" = local_opts )
  res<-nloptr(x0=x0,eval_f=eval_f,lb=lb(samind),ub=ub(samind), eval_g_eq=eval_g_eq,opts=opts)
  res$solution
}



#-------------for the consecutive time intervals-------------#
opt.consecutive.R2.S2<-function(gam,theta,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  time_1<-L.length*(intervalind-1)+1
  dif.vec<-Predict.price.model3(samind)[time,]-Predict.price.model3(samind)[L.length*(intervalind-1),]
  pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]+theta*(x^2)*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    g<-0
    for(i in 1:10){
      g<- g-dif.vec[i]*Mean.price.model3[time,i]+gam*2*(w.prev[i]+x*dif.vec[i])*dif.vec[i]*Var.price.model3[time,i]+2*theta*x*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    constr<-0
    for(i in 1:10){
      constr<-constr+x*dif.vec[i]*pri[i]+x*abs(dif.vec[i])*pri[i]*0.003
    }
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+abs(dif.vec[i])*pri[i]*0.003
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-0
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=0,
              ub=Inf,
              eval_g_ineq=eval_g_ineq,
              opts=opts)
  res$solution
}


#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R2.S2<-function(ntime,gammm,thet,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S2(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  for(i in 2:ntime){
    dif.vec<-Predict.price.model3(sam)[i*intervallength,]-Predict.price.model3(sam)[(i-1)*intervallength,]
    inter<-opt.consecutive.R2.S2(gammm,thet,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}

#----------------final answer for the optimised weight matrix
#just input 
AW.model3.R2.S2<-function(ga,the,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R2.S2(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}
#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R2.S2<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R2.S2(paragamma,paramthet,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    for(j in 2:(100/L)){
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*(DATA.test[j*L,2:11]-DATA.test[(j-1)*L,2:11])%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(DATA.test[j*L,2:11]-DATA.test[(j-1)*L,2:11])%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}

#--------------gamma=20, theta=60-----------------#
model3.L20.R2.S2.theta1<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=1,nosample=1000)
model3.L20.R2.S2.theta2<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
model3.L20.R2.S2.theta3<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
model3.L20.R2.S2.theta4<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
model3.L20.R2.S2.theta5<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
model3.L20.R2.S2.theta6<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=0,nosample=1000)
model3.L20.R2.S2.theta7<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model3.L20.R2.S2.theta8<-AV.model3.R2.S2(L=20,paragamma=0.1,paramthet=10,nosample=1000)




#---------------L=20--------------------#
per.model3.L20.R2.S2.theta1<-performance(model3.L20.R2.S2.theta1[[1]])
per.model3.L20.R2.S2.theta2<-performance(model3.L20.R2.S2.theta2[[1]])
per.model3.L20.R2.S2.theta3<-performance(model3.L20.R2.S2.theta3[[1]])
per.model3.L20.R2.S2.theta4<-performance(model3.L20.R2.S2.theta4[[1]])
per.model3.L20.R2.S2.theta5<-performance(model3.L20.R2.S2.theta5[[1]])
per.model3.L20.R2.S2.theta6<-performance(model3.L20.R2.S2.theta6[[1]])
per.model3.L20.R2.S2.theta7<-performance(model3.L20.R2.S2.theta7[[1]])
per.model3.L20.R2.S2.theta8<-performance(model3.L20.R2.S2.theta8[[1]])



#-----------------------------------------R2S3
#----------------for the first time interval------------------------#
opt.initial.R2.S3<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }
  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}

#-------------for the consecutive time intervals-------------#
opt.consecutive.R2.S3<-function(gam,theta,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  time_1<-L.length*(intervalind-1)+1
  pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
  price.mat<-Predict.price.model3(samind)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  dif.vec<-colMeans(dif.mat[(L.length*(intervalind-1)+1):time,])
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]+theta*(x^2)*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    g<-0
    for(i in 1:10){
      g<- g-dif.vec[i]*Mean.price.model3[time,i]+gam*2*(w.prev[i]+x*dif.vec[i])*dif.vec[i]*Var.price.model3[time,i]+2*theta*x*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    constr<-0
    for(i in 1:10){
      constr<-constr+x*dif.vec[i]*pri[i]+x*abs(dif.vec[i])*pri[i]*0.003
    }
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+abs(dif.vec[i])*pri[i]*0.003
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-0
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=0,
              ub=Inf,
              eval_g_ineq=eval_g_ineq,
              opts=opts)
  res$solution
}


#-----------------------------------------------------------------------------------------------------------------#
RES.model3.R2.S3<-function(ntime,gammm,thet,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S3(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  price.mat<-Predict.price.model3(sam)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  for(i in 2:ntime){
    dif.vec<-colMeans(dif.mat[(intervallength*(i-1)+1):(i*intervallength),])
    inter<-opt.consecutive.R2.S3(gammm,thet,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}

AW.model3.R2.S3<-function(ga,the,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R2.S3(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}
#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
AV.model3.R2.S3<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  for (i in 1:nosample){
    AW.sam<-AW.model3.R2.S3(paragamma,paramthet,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    #
    price.mat<-DATA.test[,2:11]
    dif.mat<-matrix(0,nrow=100,ncol=10)
    dif.mat[-1,]<-diff(price.mat)
    #
    for(j in 2:(100/L)){
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*(colMeans(dif.mat[((j-1)*L):(j*L),]))%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(colMeans(dif.mat[((j-1)*L):(j*L),]))%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}

#--------------gamma=20, theta=60-----------------#
model3.L20.R2.S3.theta1<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=1,nosample=1000)
#--------------gamma=20, theta=40-----------------#
model3.L20.R2.S3.theta2<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
#--------------gamma=20, theta=30-----------------#
model3.L20.R2.S3.theta3<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
#--------------gamma=20, theta=20-----------------#
model3.L20.R2.S3.theta4<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
#--------------gamma=20, theta=10-----------------#
model3.L20.R2.S3.theta5<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
#--------------gamma=20, theta=5-----------------#
model3.L20.R2.S3.theta6<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=0,nosample=1000)
model3.L20.R2.S3.theta7<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model3.L20.R2.S3.theta8<-AV.model3.R2.S3(L=20,paragamma=0.1,paramthet=10,nosample=1000)




#---------------L=20--------------------#
per.model3.L20.R2.S3.theta1<-performance(model3.L20.R2.S3.theta1[[1]])
per.model3.L20.R2.S3.theta2<-performance(model3.L20.R2.S3.theta2[[1]])
per.model3.L20.R2.S3.theta3<-performance(model3.L20.R2.S3.theta3[[1]])
per.model3.L20.R2.S3.theta4<-performance(model3.L20.R2.S3.theta4[[1]])
per.model3.L20.R2.S3.theta5<-performance(model3.L20.R2.S3.theta5[[1]])
per.model3.L20.R2.S3.theta6<-performance(model3.L20.R2.S3.theta6[[1]])
per.model3.L20.R2.S3.theta7<-performance(model3.L20.R2.S3.theta7[[1]])
per.model3.L20.R2.S3.theta8<-performance(model3.L20.R2.S3.theta8[[1]])


#-------------------------R2S4
#reward function 2 #strategy 4
#----------------for the first time interval------------------------#
opt.initial.R2.S4<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model3[time,i]+gam*x[i]^2*Var.price.model3[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model3[time,i]+gam*2*x[i]*Var.price.model3[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model3(samind)[1,]-iw)
    grad<-Predict.price.model3(samind)[1,]
    return(list("constraints"=constr,"jacobian"=grad))
  }

  x0<-runif(10,lb(samind),ub(samind))
  
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb(samind),
              ub=ub(samind),
              eval_g_eq=eval_g_eq,
              opts=opts)
  res$solution
}


eps<-0.0005  
#-------------for the consecutive time intervals-------------#
opt.consecutive.R2.S4<-function(gam,theta,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  time_1<-L.length*(intervalind-1)+1
  pri<-Predict.price.model3(samind)[L.length*(intervalind-1)+1,]
  price.mat<-Predict.price.model3(samind)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  dif.vec.indicator<-colMeans(dif.mat[(L.length*(intervalind-1)+1):time,])
  dif.vec<-c()
  for(j in 1:10){
    if(abs(dif.vec.indicator[j])>eps){dif.vec[j]<-dif.vec.indicator[j]}
    else{dif.vec[j]<-0}
  }
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model3[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model3[time,i]+theta*(x^2)*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    g<-0
    for(i in 1:10){
      g<- g-dif.vec[i]*Mean.price.model3[time,i]+gam*2*(w.prev[i]+x*dif.vec[i])*dif.vec[i]*Var.price.model3[time,i]+2*theta*x*(dif.vec[i]^2)*Mean.price.model3[time_1,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    constr<-0
    for(i in 1:10){
      constr<-constr+x*dif.vec[i]*pri[i]+x*abs(dif.vec[i])*pri[i]*0.003
    }
    #constr<-constr-cash
    g<-0
    for(i in 1:10){
      g<-g+dif.vec[i]*pri[i]+abs(dif.vec[i])*pri[i]*0.003
    }
    return(list("constraints"=constr,"jacobian"=g))
  }
  x0<-0
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 1000,
                "local_opts" = local_opts )
  res<-nloptr(x0=x0,
              eval_f=eval_f,
              lb=0,
              ub=Inf,
              eval_g_ineq=eval_g_ineq,
              opts=opts)
  res$solution
}


RES.model3.R2.S4<-function(ntime,gammm,thet,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S4(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  price.mat<-Predict.price.model3(sam)
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  for(i in 2:ntime){
    dif.vec.indicator<-colMeans(dif.mat[(intervallength*(i-1)+1):(i*intervallength),])
    dif.vec<-c()
    for(j in 1:10){
      if(abs(dif.vec.indicator[j])>eps){dif.vec[j]<-dif.vec.indicator[j]}
      else{dif.vec[j]<-0}
    }
    inter<-opt.consecutive.R2.S4(gammm,thet,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model3(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}
#----------------final answer for the optimised weight matrix
#just input 
AW.model3.R2.S4<-function(ga,the,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model3.R2.S4(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}


#-------------------write a function, only have to input frequency L,interval length,gamma and total number of samples
#-------------------and the return is the sequence of average value of the portfolio, ie sequence of length 100
#-----------answer to the value of the portfolio
AV.model3.R2.S4<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  #
  price.mat<-DATA.test[,2:11]
  dif.mat<-matrix(0,nrow=100,ncol=10)
  dif.mat[-1,]<-diff(price.mat)
  #
  for (i in 1:nosample){
    AW.sam<-AW.model3.R2.S4(paragamma,paramthet,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[[1]][n,]-AW.sam[[1]][n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[[1]][x,]%*%DATA.test[x,2:11]})
    CA.vector<-rep(0,100/L)
    CA<-rep(0,100)
    for(j in 2:(100/L)){
      dif.vec.indicator<-colMeans(dif.mat[((j-1)*L):(j*L),])
      dif.vec<-c()
      for(j1 in 1:10){
        if(abs(dif.vec.indicator[j1])>eps){dif.vec[j1]<-dif.vec.indicator[j1]}
        else{dif.vec[j1]<-0}
      }
      CA.vector[j]<-(1+rfr)*CA.vector[j-1]-AW.sam[[2]][j-1]*dif.vec%*%DATA.test[(j-1)*L+1,2:11]-0.003*AW.sam[[2]][j-1]*abs(dif.vec)%*%DATA.test[(j-1)*L+1,2:11]
      
      for(k2 in ((j-1)*L+1):(j*L)){
        CA[k2]<-CA.vector[j]
      }
    }
    res.matrix[i,]<-res.matrix[i,]+CA-tc
    beta.matrix[i,]<-AW.sam[[2]]
    CA.matrix[i,]<-CA.vector
    print(i)
  }
  return(list(res.matrix,beta.matrix,CA.matrix))
}

#--------------gamma=20, theta=60-----------------#
model3.L20.R2.S4.theta1<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=1,nosample=1000)
#--------------gamma=20, theta=40-----------------#
model3.L20.R2.S4.theta2<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
#--------------gamma=20, theta=30-----------------#
model3.L20.R2.S4.theta3<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
#--------------gamma=20, theta=20-----------------#
model3.L20.R2.S4.theta4<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
#--------------gamma=20, theta=10-----------------#
model3.L20.R2.S4.theta5<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
#--------------gamma=20, theta=5-----------------#
model3.L20.R2.S4.theta6<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=0,nosample=1000)
model3.L20.R2.S4.theta7<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model3.L20.R2.S4.theta8<-AV.model3.R2.S4(L=20,paragamma=0.1,paramthet=10,nosample=1000)




#---------------L=20--------------------#
per.model3.L20.R2.S4.theta1<-performance(model3.L20.R2.S4.theta1[[1]])
per.model3.L20.R2.S4.theta2<-performance(model3.L20.R2.S4.theta2[[1]])
per.model3.L20.R2.S4.theta3<-performance(model3.L20.R2.S4.theta3[[1]])
per.model3.L20.R2.S4.theta4<-performance(model3.L20.R2.S4.theta4[[1]])
per.model3.L20.R2.S4.theta5<-performance(model3.L20.R2.S4.theta5[[1]])
per.model3.L20.R2.S4.theta6<-performance(model3.L20.R2.S4.theta6[[1]])
per.model3.L20.R2.S4.theta7<-performance(model3.L20.R2.S4.theta7[[1]])
per.model3.L20.R2.S4.theta8<-performance(model3.L20.R2.S4.theta8[[1]])

