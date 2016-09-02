#reward function 2 with transaction costs
#----------------for the first time interval------------------------#
opt.initial.R2.S2<-function(gam,iw,time,samind){
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
  dif.vec<-Predict.price.model2(samind)[time,]-Predict.price.model2(samind)[L.length*(intervalind-1),]
  pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
  #objective function
  eval_f<-function(x){
    fun<-0
    for(i in 1:10){
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model2[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model2[time,i]+theta*(x^2)*(dif.vec[i]^2)*Mean.price.model2[time_1,i]
    }
    g<-0
    for(i in 1:10){
      g<- g-dif.vec[i]*Mean.price.model2[time,i]+gam*2*(w.prev[i]+x*dif.vec[i])*dif.vec[i]*Var.price.model2[time,i]+2*theta*x*(dif.vec[i]^2)*Mean.price.model2[time_1,i]
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
RES.model2.R2.S2<-function(ntime,gammm,thet,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S2(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  for(i in 2:ntime){
    dif.vec<-Predict.price.model2(sam)[i*intervallength,]-Predict.price.model2(sam)[(i-1)*intervallength,]
    inter<-opt.consecutive.R2.S2(gammm,thet,res[i-1,],intervallength,i,sam,ca[i-1])
    res[i,]<- res[i-1,]+inter*dif.vec
    betapara[i-1]<-inter
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model2(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model2(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}

#-----------------------------------------------------------------------------------------------------------------#
AW.model2.R2.S2<-function(ga,the,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R2.S2(totalinte,ga,the,intelength,whichsam)
  for(i1 in 1:totalinte){
    for(i2 in ((i1-1)*intelength+1):(i1*intelength)){
      res[i2,]<-res.inter[[1]][i1,]
    }
  }
  res.list[[1]]<-res#weight
  res.list[[2]]<-res.inter[[2]]#beta parameter
  res.list
}
#-----------------------------------------------------------------------------------------------------------------#
AV.model2.R2.S2<-function(L,paragamma,paramthet,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  for (i in 1:nosample){
    AW.sam<-AW.model2.R2.S2(paragamma,paramthet,100/L,L,i)
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
#----------------------------------------------------------------------------------#
#------------------------------------------------L=20---------------------------------------------------#
#########################################################################################################
model2.L20.R2.S2.theta1<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=1,nosample=1000)
model2.L20.R2.S2.theta2<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
model2.L20.R2.S2.theta3<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
model2.L20.R2.S2.theta4<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
model2.L20.R2.S2.theta5<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
model2.L20.R2.S2.theta6<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=0,nosample=1000)
model2.L20.R2.S2.theta7<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=5,nosample=1000)
model2.L20.R2.S2.theta8<-AV.model2.R2.S2(L=20,paragamma=0.1,paramthet=10,nosample=1000)
#---------------L=20--------------------#
per.model2.L20.R2.S2.theta1<-performance(model2.L20.R2.S2.theta1[[1]])
per.model2.L20.R2.S2.theta2<-performance(model2.L20.R2.S2.theta2[[1]])
per.model2.L20.R2.S2.theta3<-performance(model2.L20.R2.S2.theta3[[1]])
per.model2.L20.R2.S2.theta4<-performance(model2.L20.R2.S2.theta4[[1]])
per.model2.L20.R2.S2.theta5<-performance(model2.L20.R2.S2.theta5[[1]])
per.model2.L20.R2.S2.theta6<-performance(model2.L20.R2.S2.theta6[[1]])
per.model2.L20.R2.S2.theta7<-performance(model2.L20.R2.S2.theta7[[1]])
per.model2.L20.R2.S2.theta8<-performance(model2.L20.R2.S2.theta8[[1]])



pandl.R2.S2<-c(colMeans(model2.L20.R2.S2.theta1[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta2[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta3[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta4[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta5[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta6[[1]])/100000-rep(1,100))#,colMeans(model2.L20.R2.S2.theta7[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S2.theta8.1[[1]])/100000-rep(1,100))
para.R2.S2<-c(rep("1",100),rep("0.8",100),rep("0.5",100),rep("0.3",100),rep("0.1",100),rep("0",100))
df.R2.S2<-data.frame(time=rep(1:100,6),return=(pandl.R2.S2),theta=para.R2.S2)
ggplot(data=df.R2.S2,aes(x=time,y=return))+geom_line(aes(colour=theta))+ scale_colour_hue()


#---------------------------------ggplot for density of final wealth----------------------------#
den.R2.S2<-c(model2.L20.R2.S2.theta1.1[[1]][,100],model2.L20.R2.S2.theta2.1[[1]][,100],model2.L20.R2.S2.theta3.1[[1]][,100],model2.L20.R2.S2.theta4.1[[1]][,100],model2.L20.R2.S2.theta5.1[[1]][,100],model2.L20.R2.S2.theta6.1[[1]][,100],model2.L20.R2.S2.theta7.1[[1]][,100],model2.L20.R2.S2.theta8.1[[1]][,100])
para.sam.R2.S2<-c(rep("60",1000),rep("40",1000),rep("30",1000),rep("20",1000),rep("10",1000),rep("5",1000),rep("1",1000),rep("0",1000))
df.den.R2.S2<-data.frame(wealth=den.R2.S2,theta=para.sam.R2.S2)
mu<-ddply(df.den.R2.S2,"theta",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R2.S2,aes(x=wealth,color=theta))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(100000,122500))+geom_vline(data=mu,aes(xintercept=grp.mean,color=theta),linetype="dashed")

