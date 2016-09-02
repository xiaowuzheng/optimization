#reward 1 strategy 2
#----------------for the first time interval------------------------#
opt.initial.R1.S2<-function(gam,iw,time,samind){
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

rfr<-0.0002
#-------------for the consecutive time intervals-------------#
opt.consecutive.R1.S2<-function(gam,w.prev,L.length,intervalind,samind,cash){
  #end of the interval
  time<-L.length*intervalind
  #objective function
  eval_f<-function(x){
    dif.vec<-Predict.price.model2(samind)[time,]-Predict.price.model2(samind)[L.length*(intervalind-1),]
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
    dif.vec<-Predict.price.model2(samind)[time,]-Predict.price.model2(samind)[L.length*(intervalind-1),]
    pri<-Predict.price.model2(samind)[L.length*(intervalind-1)+1,]
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

#-----------------------------------------------------------------------------------------------------------------#
RES.model2.R1.S2<-function(ntime,gammm,intervallength,sam){
  RES.list<-list()
  res<-matrix(0,nrow=ntime,ncol=10)
  res[1,]<-opt.initial.R1.S2(gammm,100000,intervallength,sam)
  betapara<-c()#4
  tc<-c()#5
  tc[1]<-0
  ca<-c()#5
  ca[1]<-0#cash not used in investment at t=1 (initial wealth=initial investment)
  for(i in 2:ntime){
    dif.vec<-Predict.price.model2(sam)[i*intervallength,]-Predict.price.model2(sam)[(i-1)*intervallength,]
    inter<-opt.consecutive.R1.S2(gammm,res[i-1,],intervallength,i,sam,ca[i-1])
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
AW.model2.R1.S2<-function(ga,totalinte,intelength,whichsam){
  res.list<-list()
  res<-matrix(0,nrow=100,ncol=10)
  res.inter<-RES.model2.R1.S2(totalinte,ga,intelength,whichsam)
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
AV.model2.R1.S2<-function(L,paragamma,nosample){
  #each row is a series of value of portfolio over 100 trading days for one sample
  res.matrix<-matrix(0,nrow=nosample,ncol=100)
  beta.matrix<-matrix(0,nrow=nosample,ncol=100/L-1)
  CA.matrix<-matrix(0,nrow=nosample,ncol=100/L)
  for (i in 1:nosample){
    AW.sam<-AW.model2.R1.S2(paragamma,100/L,L,i)
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

#------------------------------------------------L=20---------------------------------------------------#
#--------------gamma=10-----------------#
model2.L20.R1.S2.gam1<-AV.model2.R1.S2(L=20,paragamma=10,nosample=1000)
#--------------gamma=5-----------------#
model2.L20.R1.S2.gam2<-AV.model2.R1.S2(L=20,paragamma=5,nosample=1000)
#--------------gamma=1-----------------#
model2.L20.R1.S2.gam3<-AV.model2.R1.S2(L=20,paragamma=1,nosample=1000)
#--------------gamma=0.1-----------------#
model2.L20.R1.S2.gam4<-AV.model2.R1.S2(L=20,paragamma=0.1,nosample=1000)
#--------------gamma=0.05-----------------#
model2.L20.R1.S2.gam5<-AV.model2.R1.S2(L=20,paragamma=0.05,nosample=1000)
#--------------gamma=0.01-----------------#
model2.L20.R1.S2.gam6<-AV.model2.R1.S2(L=20,paragamma=0.01,nosample=1000)
#--------------gamma=0.005-----------------#
model2.L20.R1.S2.gam7<-AV.model2.R1.S2(L=20,paragamma=0.005,nosample=1000)
#--------------gamma=0.001-----------------#
model2.L20.R1.S2.gam8<-AV.model2.R1.S2(L=20,paragamma=0.001,nosample=1000)
#--------------gamma=0.0005-----------------#
model2.L20.R1.S2.gam9<-AV.model2.R1.S2(L=20,paragamma=0.0005,nosample=1000)
#--------------gamma=0.0004-----------------#
model2.L20.R1.S2.gam10<-AV.model2.R1.S2(L=20,paragamma=0.0004,nosample=1000)
#--------------gamma=0.0003-----------------#
model2.L20.R1.S2.gam11<-AV.model2.R1.S2(L=20,paragamma=0.0003,nosample=1000)
#--------------gamma=0.0002-----------------#
model2.L20.R1.S2.gam12<-AV.model2.R1.S2(L=20,paragamma=0.0002,nosample=1000)
#--------------gamma=0.0001-----------------#
model2.L20.R1.S2.gam13<-AV.model2.R1.S2(L=20,paragamma=0.0001,nosample=1000)
#--------------gamma=0.00008-----------------#
model2.L20.R1.S2.gam14<-AV.model2.R1.S2(L=20,paragamma=0.00008,nosample=1000)
#--------------gamma=0.00007-----------------#
model2.L20.R1.S2.gam15<-AV.model2.R1.S2(L=20,paragamma=0.00007,nosample=1000)
#--------------gamma=0.00005-----------------#
model2.L20.R1.S2.gam16<-AV.model2.R1.S2(L=20,paragamma=0.00005,nosample=1000)
#--------------gamma=0.00003-----------------#
model2.L20.R1.S2.gam17<-AV.model2.R1.S2(L=20,paragamma=0.00003,nosample=1000)
#--------------gamma=0.00001-----------------#
model2.L20.R1.S2.gam18<-AV.model2.R1.S2(L=20,paragamma=0.00001,nosample=1000)

#---------------L=20 performance--------------------#
per.model2.L20.R1.S2.gam1<-performance(model2.L20.R1.S2.gam1[[1]])
per.model2.L20.R1.S2.gam2<-performance(model2.L20.R1.S2.gam2[[1]])
per.model2.L20.R1.S2.gam3<-performance(model2.L20.R1.S2.gam3[[1]])
per.model2.L20.R1.S2.gam4<-performance(model2.L20.R1.S2.gam4[[1]])
per.model2.L20.R1.S2.gam5<-performance(model2.L20.R1.S2.gam5[[1]])
per.model2.L20.R1.S2.gam6<-performance(model2.L20.R1.S2.gam6[[1]])
per.model2.L20.R1.S2.gam7<-performance(model2.L20.R1.S2.gam7[[1]])
per.model2.L20.R1.S2.gam8<-performance(model2.L20.R1.S2.gam8[[1]])
per.model2.L20.R1.S2.gam9<-performance(model2.L20.R1.S2.gam9[[1]])
per.model2.L20.R1.S2.gam10<-performance(model2.L20.R1.S2.gam10[[1]])
per.model2.L20.R1.S2.gam11<-performance(model2.L20.R1.S2.gam11[[1]])
per.model2.L20.R1.S2.gam12<-performance(model2.L20.R1.S2.gam12[[1]])
per.model2.L20.R1.S2.gam13<-performance(model2.L20.R1.S2.gam13[[1]])
per.model2.L20.R1.S2.gam14<-performance(model2.L20.R1.S2.gam14[[1]])
per.model2.L20.R1.S2.gam15<-performance(model2.L20.R1.S2.gam15[[1]])
per.model2.L20.R1.S2.gam16<-performance(model2.L20.R1.S2.gam16[[1]])
per.model2.L20.R1.S2.gam17<-performance(model2.L20.R1.S2.gam17[[1]])
per.model2.L20.R1.S2.gam18<-performance(model2.L20.R1.S2.gam18[[1]])


#--------------------------------ggplot for profit and loss----------------------------#
pandl.R1.S2<-c(colMeans(model2.L20.R1.S2.gam9[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam10[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam11[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam12[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam13[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam14[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam15[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam16[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam17[[1]])/100000-rep(1,100),colMeans(model2.L20.R1.S2.gam18[[1]])/100000-rep(1,100))
para.R1.S2<-c(rep("0.0005",100),rep("0.0004",100),rep("0.0003",100),rep("0.0002",100),rep("0.0001",100),rep("0.00008",100),rep("0.00007",100),rep("0.00005",100),rep("0.00003",100),rep("0.00001",100))
df.R1.S2<-data.frame(time=rep(1:100,10),return=(pandl.R1.S2),gamma=para.R1.S2)
ggplot(data=df.R1.S2,aes(x=time,y=return))+geom_line(aes(colour=gamma))+ scale_colour_hue()+coord_cartesian(ylim=c(-0.05,0.18))#+ ggtitle("Policy 2")

#---------------------------------ggplot for density of final wealth----------------------------#
den.R1.S2<-c(model2.L20.R1.S2.gam9[[1]][,100],model2.L20.R1.S2.gam10[[1]][,100],model2.L20.R1.S2.gam11[[1]][,100],model2.L20.R1.S2.gam12[[1]][,100],model2.L20.R1.S2.gam13[[1]][,100],model2.L20.R1.S2.gam14[[1]][,100],model2.L20.R1.S2.gam15[[1]][,100],model2.L20.R1.S2.gam16[[1]][,100],model2.L20.R1.S2.gam17[[1]][,100],model2.L20.R1.S2.gam18[[1]][,100])
para.sam.R1.S2<-c(rep("0.0005",1000),rep("0.0004",1000),rep("0.0003",1000),rep("0.0002",1000),rep("0.0001",1000),rep("0.00008",1000),rep("0.00007",1000),rep("0.00005",1000),rep("0.00003",1000),rep("0.00001",1000))
df.den.R1.S2<-data.frame(wealth=den.R1.S2,gamma=para.sam.R1.S2)
mu<-ddply(df.den.R1.S2,"gamma",summarise,grp.mean=mean(wealth))
ggplot(data=df.den.R1.S2,aes(x=wealth,color=gamma))+stat_density(geom='line',position = 'identity')+coord_cartesian(xlim=c(90000,120000))+geom_vline(data=mu,aes(xintercept=grp.mean,color=gamma),linetype="dashed")


beta1.int1<-which(model2.L20.R1.S2.gam18[[2]][,1]!=0)
model2.L20.R1.S2.gam18[[2]][c(139, 161, 525, 674, 700),1]
beta1.int2<-which(model2.L20.R1.S2.gam18[[2]][,2]!=0)
model2.L20.R1.S2.gam18[[2]][c(131,172,202,268,587,852,907,944,1000),2]
beta1.int3<-which(model2.L20.R1.S2.gam18[[2]][,3]!=0)
model2.L20.R1.S2.gam18[[2]][c(127,384,451,493,802,943),3]
beta1.int4<-which(model2.L20.R1.S2.gam18[[2]][,4]!=0)
model2.L20.R1.S2.gam18[[2]][c(26,31,66,70, 102, 136, 143, 381, 442, 509, 652, 655, 679, 774),4]
