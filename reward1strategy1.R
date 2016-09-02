#upper and lower bound
ub<-function(x){
  rep(100000*0.5,10)/abs(Predict.price.model2(x)[1,])
}
lb<-function(x){
  -rep(100000*0.5,10)/abs(Predict.price.model2(x)[1,])
}

#---------------initial interval-----------------------------------#
opt.initial.R1.S1<-function(gam,iw,time,samind){
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-x[i]*Mean.price.model2[time,i]+gam*(x[i]^2)*Var.price.model2[time,i]
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




#-------------for the consecutive time intervals--------------------#
opt.consecutive.R1.S1<-function(gam,w.prev,L.length,intervalind,samind){
  #end of the interval
  time<-L.length*intervalind
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
    constr<-0
    pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
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


#-------------allocation result for each interval------------------#
RES.model2.R1.S1<-function(ntime,gammm,intervallength,sam){
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S1(gammm,100000,intervallength,sam)
  for(i in 2:ntime){
    res[i,]<- opt.consecutive.R1.S1(gammm,res[i-1,],intervallength,i,sam)
  }
  res
}
#-------------allocation result for entire horizon------------------#
AW.model2.R1.S1<-function(ga,totalinte,intelength,whichsam){
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R1.S1(totalinte,ga,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[i1,]
    }
  }
  res
}

#-------------value of portfolio------------------#
AV.model2.R1.S1<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  for (i in 1:nosample){
    AW.sam<-AW.model2.R1.S1(paragamma,100/L,L,i)
    tc<-c(0,0.003*sapply(2:100,function(n){sum(abs(AW.sam[n,]-AW.sam[n-1,])%*%DATA.test[n,2:11])}))
    res.matrix[i,]<-sapply(1:100,function(x){AW.sam[x,]%*%DATA.test[x,2:11]})
    res.matrix[i,]<-res.matrix[i,]-tc
    print(i)
  }
  res.matrix
}

#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=10-----------------#
model2.L20.R1.S1.gam1<-AV.model2.R1.S1(L=20,paragamma=10,nosample=1000)
#--------------gamma=5-----------------#
model2.L20.R1.S1.gam2<-AV.model2.R1.S1(L=20,paragamma=5,nosample=1000)
#--------------gamma=1-----------------#
model2.L20.R1.S1.gam3<-AV.model2.R1.S1(L=20,paragamma=1,nosample=1000)
#--------------gamma=0.1-----------------#
model2.L20.R1.S1.gam4<-AV.model2.R1.S1(L=20,paragamma=0.1,nosample=1000)
#--------------gamma=0.05-----------------#
model2.L20.R1.S1.gam5<-AV.model2.R1.S1(L=20,paragamma=0.05,nosample=1000)
#--------------gamma=0.01-----------------#
model2.L20.R1.S1.gam6<-AV.model2.R1.S1(L=20,paragamma=0.01,nosample=1000)
#--------------gamma=0.005-----------------#
model2.L20.R1.S1.gam7<-AV.model2.R1.S1(L=20,paragamma=0.005,nosample=1000)
#--------------gamma=0.001-----------------#
model2.L20.R1.S1.gam8<-AV.model2.R1.S1(L=20,paragamma=0.001,nosample=1000)
#--------------gamma=0.0005-----------------#
model2.L20.R1.S1.gam9<-AV.model2.R1.S1(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model2.L20.R1.S1.gam10<-AV.model2.R1.S1(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model2.L20.R1.S1.gam11<-AV.model2.R1.S1(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model2.L20.R1.S1.gam12<-AV.model2.R1.S1(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model2.L20.R1.S1.gam13<-AV.model2.R1.S1(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model2.L20.R1.S1.gam14<-AV.model2.R1.S1(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model2.L20.R1.S1.gam15<-AV.model2.R1.S1(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model2.L20.R1.S1.gam16<-AV.model2.R1.S1(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model2.L20.R1.S1.gam17<-AV.model2.R1.S1(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model2.L20.R1.S1.gam18<-AV.model2.R1.S1(L=20,paragamma=0.00001,nosample=1000)
#----------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#
per.model2.L20.R1.S1.gam1<-performance(model2.L20.R1.S1.gam1)
per.model2.L20.R1.S1.gam2<-performance(model2.L20.R1.S1.gam2)
per.model2.L20.R1.S1.gam3<-performance(model2.L20.R1.S1.gam3)
per.model2.L20.R1.S1.gam4<-performance(model2.L20.R1.S1.gam4)
per.model2.L20.R1.S1.gam5<-performance(model2.L20.R1.S1.gam5)
per.model2.L20.R1.S1.gam6<-performance(model2.L20.R1.S1.gam6)
per.model2.L20.R1.S1.gam7<-performance(model2.L20.R1.S1.gam7)
per.model2.L20.R1.S1.gam8<-performance(model2.L20.R1.S1.gam8)
per.model2.L20.R1.S1.gam9<-performance(model2.L20.R1.S1.gam9)
per.model2.L20.R1.S1.gam10<-performance(model2.L20.R1.S1.gam10)
per.model2.L20.R1.S1.gam11<-performance(model2.L20.R1.S1.gam11)
per.model2.L20.R1.S1.gam12<-performance(model2.L20.R1.S1.gam12)
per.model2.L20.R1.S1.gam13<-performance(model2.L20.R1.S1.gam13)
per.model2.L20.R1.S1.gam14<-performance(model2.L20.R1.S1.gam14)
per.model2.L20.R1.S1.gam15<-performance(model2.L20.R1.S1.gam15)
per.model2.L20.R1.S1.gam16<-performance(model2.L20.R1.S1.gam16)
per.model2.L20.R1.S1.gam17<-performance(model2.L20.R1.S1.gam17)
per.model2.L20.R1.S1.gam18<-performance(model2.L20.R1.S1.gam18)

#--------------------------------ggplot for profit and loss----------------------------#
pandl.R1.S1<-c(colMeans(model2.L20.R1.S1.gam1)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam1)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam3)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam5)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam4)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam6)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam7)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam8)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam9)/100000-rep(1,100),colMeans(model2.L20.R1.S1.gam10)/100000-rep(1,100))
para.R1.S1<-c(rep("0.0005",100),rep("0.0004",100),rep("0.0003",100),rep("0.0002",100),rep("0.0001",100),rep("0.00008",100),rep("0.00007",100),rep("0.00005",100),rep("0.00003",100),rep("0.00001",100))
df.R1.S1<-data.frame(time=rep(1:100,10),return=(pandl.R1.S1),gamma=para.R1.S1)
ggplot(data=df.R1.S1,aes(x=time,y=return))+geom_line(aes(colour=gamma))+ scale_colour_hue()


#---------------------------------ggplot for density of final wealth----------------------------#
den.R1.S1<-c(model2.L20.R1.S1.gam1[,100],model2.L20.R1.S1.gam2[,100],model2.L20.R1.S1.gam3[,100],model2.L20.R1.S1.gam4[,100],model2.L20.R1.S1.gam5[,100],model2.L20.R1.S1.gam6[,100],model2.L20.R1.S1.gam7[,100],model2.L20.R1.S1.gam8[,100],model2.L20.R1.S1.gam9[,100],model2.L20.R1.S1.gam10[,100])
para.sam.R1.S1<-c(rep("0.0005",1000),rep("0.0004",1000),rep("0.0003",1000),rep("0.0002",1000),rep("0.0001",1000),rep("0.00008",1000),rep("0.00007",1000),rep("0.00005",1000),rep("0.00003",1000),rep("0.00001",1000))
df.den.R1.S1<-data.frame(wealth=den.R1.S1,gamma=para.sam.R1.S1)
mu<-ddply(df.den.R1.S1,"gamma",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R1.S1,aes(x=wealth,color=gamma))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(108000,130000))+geom_vline(data=mu,aes(xintercept=grp.mean,color=gamma),linetype="dashed")


