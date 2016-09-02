#reward function 2 #strategy 4
#----------------for the first time interval------------------------#
opt.initial.R2.S4<-function(gam,iw,time,samind){
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

eps<-0.0005  
#-------------for the consecutive time intervals-------------#
opt.consecutive.R2.S4<-function(gam,theta,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  time_1<-L.length*(intervalind-1)+1
  pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
  price.mat<-Predict.price.model2(samind)
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
RES.model2.R2.S4<-function(ntime,gammm,thet,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R2.S4(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  price.mat<-Predict.price.model2(sam)
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
    tc[i]<-0.003*inter*abs(dif.vec)%*%Predict.price.model2(sam)[(i-1)*intervallength+1,]
    ca[i]<-ca[i-1]*(1+rfr)-tc[i]-betapara[i-1]*dif.vec%*%Predict.price.model2(sam)[(i-1)*intervallength+1,]
  }
  RES.list[[1]]<-res
  RES.list[[2]]<-betapara
  RES.list
}
#-----------------------------------------------------------------------------------------------------------------#
AW.model2.R2.S4<-function(ga,the,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R2.S4(totalinte,ga,the,intelength,whichsam)
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
AV.model2.R2.S4<-function(L,paragamma,paramthet,nosample){
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
    AW.sam<-AW.model2.R2.S4(paragamma,paramthet,100/L,L,i)
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





#########################################################################################################
#--------------gamma=20, theta=60-----------------#
model2.L20.R2.S4.theta1<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=1,nosample=1000)
#--------------gamma=20, theta=40-----------------#
model2.L20.R2.S4.theta2<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=0.8,nosample=1000)
#--------------gamma=20, theta=30-----------------#
model2.L20.R2.S4.theta3<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=0.5,nosample=1000)
#--------------gamma=20, theta=20-----------------#
model2.L20.R2.S4.theta4<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=0.3,nosample=1000)
#--------------gamma=20, theta=10-----------------#
model2.L20.R2.S4.theta5<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=0.1,nosample=1000)
#--------------gamma=20, theta=5-----------------#
model2.L20.R2.S4.theta6<-AV.model2.R2.S4(L=20,paragamma=0.1,paramthet=0,nosample=1000)

#---------------L=20--------------------#
per.model2.L20.R2.S4.theta1<-performance(model2.L20.R2.S4.theta1[[1]])
per.model2.L20.R2.S4.theta2<-performance(model2.L20.R2.S4.theta2[[1]])
per.model2.L20.R2.S4.theta3<-performance(model2.L20.R2.S4.theta3[[1]])
per.model2.L20.R2.S4.theta4<-performance(model2.L20.R2.S4.theta4[[1]])
per.model2.L20.R2.S4.theta5<-performance(model2.L20.R2.S4.theta5[[1]])
per.model2.L20.R2.S4.theta6<-performance(model2.L20.R2.S4.theta6[[1]])

#--------------------------------ggplot for profit and loss----------------------------#
pandl.R2.S4<-c(colMeans(model2.L20.R2.S4.theta1[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta2[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta3[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta4[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta5[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta6[[1]])/100000-rep(1,100))#,colMeans(model2.L20.R2.S4.theta7[[1]])/100000-rep(1,100),colMeans(model2.L20.R2.S4.theta8.1[[1]])/100000-rep(1,100))
para.R2.S4<-c(rep("1",100),rep("0.8",100),rep("0.5",100),rep("0.3",100),rep("0.1",100),rep("0",100))
df.R2.S4<-data.frame(time=rep(1:100,6),return=(pandl.R2.S4),theta=para.R2.S4)
ggplot(data=df.R2.S4,aes(x=time,y=return))+geom_line(aes(colour=theta))+ scale_colour_hue()