#reward function 1 #strategy 4

opt.initial.R1.S4<-function(gam,iw,time,samind){
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

opt.initial.R1.S4(3,100000,20,1)

eps<-0.0005      
#tt1<-opt.initial.R1.S4(3,100000,20,1)
#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S4<-function(gam,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
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
      fun<-fun-(w.prev[i]+x*dif.vec[i])*Mean.price.model2[time,i]+gam*(w.prev[i]+x*dif.vec[i])^2*Var.price.model2[time,i]
    }
    g<-0
    for(i in 1:10){
      g<-g-dif.vec[i]*Mean.price.model2[time,i]+2*gam*dif.vec[i]*Var.price.model2[time,i]*(x*dif.vec[i]+1)
    }
    return(list("objective"=fun,"gradient"=g))
  }
  #constriant functions
  eval_g_ineq<-function(x){
    pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
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
  #
  x0<-0.1
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG","xtol_rel" = 1.0e-7,"maxeval" = 1000, "local_opts" = local_opts )
  res<-nloptr(x0=x0,eval_f=eval_f,lb=0, ub=Inf,eval_g_ineq=eval_g_ineq,opts=opts)
  res$solution
}

#------------------------------
RES.model2.R1.S4<-function(ntime,gammm,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S4(gammm,100000,intervallength,sam)
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
    inter<-opt.consecutive.R1.S4(gammm,res[i-1,],intervallength,i,sam,ca[i-1])
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
AW.model2.R1.S4<-function(ga,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R1.S4(totalinte,ga,intelength,whichsam)
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
AV.model2.R1.S4<-function(L,paragamma,nosample){
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
    AW.sam<-AW.model2.R1.S4(paragamma,100/L,L,i)
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
#--------------gamma=10-----------------#
model2.L20.R1.S4.gam1<-AV.model2.R1.S4(L=20,paragamma=10,nosample=1000)
#--------------gamma=5-----------------#
model2.L20.R1.S4.gam2<-AV.model2.R1.S4(L=20,paragamma=5,nosample=1000)
#--------------gamma=1-----------------#
model2.L20.R1.S4.gam3<-AV.model2.R1.S4(L=20,paragamma=1,nosample=1000)
#--------------gamma=0.1-----------------#
model2.L20.R1.S4.gam4<-AV.model2.R1.S4(L=20,paragamma=0.1,nosample=1000)
#--------------gamma=0.05-----------------#
model2.L20.R1.S4.gam5<-AV.model2.R1.S4(L=20,paragamma=0.05,nosample=1000)
#--------------gamma=0.01-----------------#
model2.L20.R1.S4.gam6<-AV.model2.R1.S4(L=20,paragamma=0.01,nosample=1000)
#--------------gamma=0.005-----------------#
model2.L20.R1.S4.gam7<-AV.model2.R1.S4(L=20,paragamma=0.005,nosample=1000)
#--------------gamma=0.001-----------------#
model2.L20.R1.S4.gam8<-AV.model2.R1.S4(L=20,paragamma=0.001,nosample=1000)
#--------------gamma=0.0005-----------------#
model2.L20.R1.S4.gam9<-AV.model2.R1.S4(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model2.L20.R1.S4.gam10<-AV.model2.R1.S4(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model2.L20.R1.S4.gam11<-AV.model2.R1.S4(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model2.L20.R1.S4.gam12<-AV.model2.R1.S4(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model2.L20.R1.S4.gam13<-AV.model2.R1.S4(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model2.L20.R1.S4.gam14<-AV.model2.R1.S4(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model2.L20.R1.S4.gam15<-AV.model2.R1.S4(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model2.L20.R1.S4.gam16<-AV.model2.R1.S4(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model2.L20.R1.S4.gam17<-AV.model2.R1.S4(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model2.L20.R1.S4.gam18<-AV.model2.R1.S4(L=20,paragamma=0.00001,nosample=1000)
#-------------------------------------------------------------------------------------#
#---------------L=20 performance--------------------#
per.model2.L20.R1.S4.gam1<-performance(model2.L20.R1.S4.gam1[[1]])
per.model2.L20.R1.S4.gam2<-performance(model2.L20.R1.S4.gam2[[1]])
per.model2.L20.R1.S4.gam3<-performance(model2.L20.R1.S4.gam3[[1]])
per.model2.L20.R1.S4.gam4<-performance(model2.L20.R1.S4.gam4[[1]])
per.model2.L20.R1.S4.gam5<-performance(model2.L20.R1.S4.gam5[[1]])
per.model2.L20.R1.S4.gam6<-performance(model2.L20.R1.S4.gam6[[1]])
per.model2.L20.R1.S4.gam7<-performance(model2.L20.R1.S4.gam7[[1]])
per.model2.L20.R1.S4.gam8<-performance(model2.L20.R1.S4.gam8[[1]])
per.model2.L20.R1.S4.gam9<-performance(model2.L20.R1.S4.gam9[[1]])
per.model2.L20.R1.S4.gam10<-performance(model2.L20.R1.S4.gam10[[1]])
per.model2.L20.R1.S4.gam11<-performance(model2.L20.R1.S4.gam11[[1]])
per.model2.L20.R1.S4.gam12<-performance(model2.L20.R1.S4.gam12[[1]])
per.model2.L20.R1.S4.gam13<-performance(model2.L20.R1.S4.gam13[[1]])
per.model2.L20.R1.S4.gam14<-performance(model2.L20.R1.S4.gam14[[1]])
per.model2.L20.R1.S4.gam15<-performance(model2.L20.R1.S4.gam15[[1]])
per.model2.L20.R1.S4.gam16<-performance(model2.L20.R1.S4.gam16[[1]])
per.model2.L20.R1.S4.gam17<-performance(model2.L20.R1.S4.gam17[[1]])
per.model2.L20.R1.S4.gam18<-performance(model2.L20.R1.S4.gam18[[1]])

#--------------------------------ggplot for profit and loss----------------------------#
pandl.R1.S4<-c(colMeans(model2.L20.R1.S4.gam9[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam10[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam11[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam12[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam13[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam14[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam15[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam16[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam17[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam18[[1]])/100000-rep(1,100))
para.R1.S4<-c(rep("0.0005",100),rep("0.0004",100),rep("0.0003",100),rep("0.0002",100),rep("0.0001",100),rep("0.00008",100),rep("0.00007",100),rep("0.00005",100),rep("0.00003",100),rep("0.00001",100))
df.R1.S4<-data.frame(time=rep(1:100,10),return=(pandl.R1.S4),gamma=para.R1.S4)
ggplot(data=df.R1.S4,aes(x=time,y=return))+geom_line(aes(colour=gamma))+ scale_colour_hue()+coord_cartesian(ylim=c(-0.09,0.18))#+ ggtitle("Policy 4")

#---------------------------------ggplot for density of final wealth----------------------------#
den.R1.S4<-c(model2.L20.R1.S4.gam9[[1]][,100],model2.L20.R1.S4.gam10[[1]][,100],model2.L20.R1.S4.gam11[[1]][,100],model2.L20.R1.S4.gam12[[1]][,100],model2.L20.R1.S4.gam13[[1]][,100],model2.L20.R1.S4.gam14[[1]][,100],model2.L20.R1.S4.gam15[[1]][,100],model2.L20.R1.S4.gam16[[1]][,100],model2.L20.R1.S4.gam17[[1]][,100],model2.L20.R1.S4.gam18[[1]][,100])
para.sam.R1.S4<-c(rep("0.0005",1000),rep("0.0004",1000),rep("0.0003",1000),rep("0.0002",1000),rep("0.0001",1000),rep("0.00008",1000),rep("0.00007",1000),rep("0.00005",1000),rep("0.00003",1000),rep("0.00001",1000))
df.den.R1.S4<-data.frame(wealth=den.R1.S4,gamma=para.sam.R1.S4)
mu<-ddply(df.den.R1.S4,"gamma",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R1.S4,aes(x=wealth,color=gamma))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(110000,215000))+geom_vline(data=mu,aes(xintercept=grp.mean,color=gamma),linetype="dashed")



#--------------------------------ggplot for profit and loss----------------------------#
pandl.R1.S4<-c(colMeans(model2.L20.R1.S4.gam1[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam2[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam3[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam4[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam5[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam6[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam7[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam8[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam9[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S4.gam10[[1]])/100000-rep(1,100))
para.R1.S4<-c(rep("0.0005",100),rep("0.0004",100),rep("0.0003",100),rep("0.0002",100),rep("0.0001",100),rep("0.00008",100),rep("0.00007",100),rep("0.00005",100),rep("0.00003",100),rep("0.00001",100))
df.R1.S4<-data.frame(time=rep(1:100,10),return=(pandl.R1.S4),gamma=para.R1.S4)
ggplot(data=df.R1.S4,aes(x=time,y=return))+geom_line(aes(colour=gamma))+ scale_colour_hue()+coord_cartesian(ylim=c(-0.09,0.18))#+ ggtitle("Policy 4")
