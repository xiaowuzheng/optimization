#reward function 2 with transaction costs
#----------------for the first time interval------------------------#
opt.initial.R2.S1<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model2[time,i]+gam*x[i]^2*Var.price.model2[time,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model2[time,i]+gam*2*x[i]*Var.price.model2[time,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  
  #constriant functions
  eval_g_eq<-function(x){
    constr<-c(x%*%Predict.price.model2(samind)[1,]-iw)
    grad<-Predict.price.model2(samind)[1,]
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
      fun<-fun-x[i]*Mean.price.model2[time,i]+gam*x[i]^2*Var.price.model2[time,i]+theta*(x[i]-w.prev[i])^2*Mean.price.model2[time_1,i]
    }
    g<-c()
    for(i in 1:10){
      g[i]<- -Mean.price.model2[time,i]+gam*2*x[i]*Var.price.model2[time,i]+2*theta*(x[i]-w.prev[i])*Mean.price.model2[time_1,i]
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_eq<-function(x){
    constr<-0
    pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
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
RES.model2.R2.S1<-function(ntime,gammm,thet,intervallength,sam){
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S1(gammm,100000,intervallength,sam)
  for(i in 2:ntime){
    res[i,]<- opt.consecutive.R2.S1(gammm,thet,res[i-1,],intervallength,i,sam)
  }
  res
}
#-----------------------------------------------------------------------------------------------------------------#
AW.model2.R2.S1<-function(ga,the,totalinte,intelength,whichsam){
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R2.S1(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[i1,]
    }
  }
  res
}

#-----------------------------------------------------------------------------------------------------------------#
AV.model2.R2.S1<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  for (i in 1:nosample){
    AW.sam<-AW.model2.R2.S1(paragamma,paramthet,100/L,L,i)
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[x,]%*%DATA.test[x,2:11]})
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[n,]-AW.sam[n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-res.matrix[i,]-tc
    print(i)
  }
  
  res.matrix
}


model2.L20.R2.S1.theta1<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=60,nosample=1000)
model2.L20.R2.S1.theta2<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=40,nosample=1000)
model2.L20.R2.S1.theta3<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=30,nosample=1000)
model2.L20.R2.S1.theta4<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=20,nosample=1000)
model2.L20.R2.S1.theta5<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=10,nosample=1000)
model2.L20.R2.S1.theta6<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model2.L20.R2.S1.theta7<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=1,nosample=1000)
model2.L20.R2.S1.theta12<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
model2.L20.R2.S1.theta11<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
model2.L20.R2.S1.theta10<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
model2.L20.R2.S1.theta9<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
model2.L20.R2.S1.theta8<-AV.model2.R2.S1(L=20,paragamma=0.1,paramthet=0,nosample=1000)
#
#---------------L=20--------------------#
per.model2.L20.R2.S1.theta1<-performance(model2.L20.R2.S1.theta1)
per.model2.L20.R2.S1.theta2<-performance(model2.L20.R2.S1.theta2)
per.model2.L20.R2.S1.theta3<-performance(model2.L20.R2.S1.theta3)
per.model2.L20.R2.S1.theta4<-performance(model2.L20.R2.S1.theta4)
per.model2.L20.R2.S1.theta5<-performance(model2.L20.R2.S1.theta5)
per.model2.L20.R2.S1.theta6<-performance(model2.L20.R2.S1.theta6)
per.model2.L20.R2.S1.theta7<-performance(model2.L20.R2.S1.theta7)
per.model2.L20.R2.S1.theta8<-performance(model2.L20.R2.S1.theta8)
per.model2.L20.R2.S1.theta9<-performance(model2.L20.R2.S1.theta9)
per.model2.L20.R2.S1.theta10<-performance(model2.L20.R2.S1.theta10)
per.model2.L20.R2.S1.theta11<-performance(model2.L20.R2.S1.theta11)
per.model2.L20.R2.S1.theta12<-performance(model2.L20.R2.S1.theta12)

#--------------------------------ggplot for profit and loss----------------------------#
pandl.R2.S1<-c(colMeans(model2.L20.R2.S1.theta5)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta6)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta7)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta12)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta11)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta10)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta9)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta8)/100000-rep(1,100))
para.R2.S1<-c(rep("10",100),rep("5",100),rep("1",100),rep("0.8",100),rep("0.5",100),rep("0.3",100),rep("0.1",100),rep("0",100))
df.R2.S1<-data.frame(time=rep(1:100,8),return=(pandl.R2.S1),theta=para.R2.S1)
ggplot(data=df.R2.S1,aes(x=time,y=return))+geom_line(aes(colour=theta))+ scale_colour_hue()

#--------------------------------ggplot for profit and loss----------------------------#
pandl.R2.S1<-c(colMeans(model2.L20.R2.S1.theta1)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta2)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta3)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta4)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta5)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta6)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta7)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta8)/100000-rep(1,100),colMeans(model2.L20.R2.S1.theta9)/100000-rep(1,100))
para.R2.S1<-c(rep("60",100),rep("40",100),rep("30",100),rep("20",100),rep("10",100),rep("5",100),rep("1",100),rep("0",100),rep("0.1",100))
df.R2.S1<-data.frame(time=rep(1:100,9),return=(pandl.R2.S1),theta=para.R2.S1)
ggplot(data=df.R2.S1,aes(x=time,y=return))+geom_line(aes(colour=theta))+ scale_colour_hue()

#---------------------------------ggplot for density of final wealth----------------------------#
den.R2.S1<-c(model2.L20.R2.S1.theta1[,100],model2.L20.R2.S1.theta2[,100],model2.L20.R2.S1.theta3[,100],model2.L20.R2.S1.theta4[,100],model2.L20.R2.S1.theta5[,100],model2.L20.R2.S1.theta6[,100],model2.L20.R2.S1.theta7[,100],model2.L20.R2.S1.theta8[,100],model2.L20.R2.S1.theta9[,100])
para.sam.R2.S1<-c(rep("60",1000),rep("40",1000),rep("30",1000),rep("20",1000),rep("10",1000),rep("5",1000),rep("1",1000),rep("0",1000),rep("0.1",1000))
df.den.R2.S1<-data.frame(wealth=den.R2.S1,theta=para.sam.R2.S1)
mu<-ddply(df.den.R2.S1,"theta",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R2.S1,aes(x=wealth,color=theta))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(110000,114500))+geom_vline(data=mu,aes(xintercept=grp.mean,color=theta),linetype="dashed")
