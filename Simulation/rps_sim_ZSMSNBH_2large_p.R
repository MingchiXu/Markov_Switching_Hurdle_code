###################################################################
# RPS for Markov Switching models: ZS-MSNBH
# second large reporting rate: p0=0.8
###################################################################

args=commandArgs(trailingOnly = TRUE)

library(nimble)
library(readxl)
library(coda)
library(utils)



if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}


# the hyper parameter
T0<- as.numeric(args[1])
print(T0)
K<-4 
Tsim<- 84
# set seed T0
set.seed(T0)



# read in the data files 

cases<-read.csv('/lustre03/project/6003552/mingchi/cases_80.csv')
N_matrix<-read.csv("/lustre03/project/6003552/mingchi/spatial_matrix.csv", 
                   header=TRUE, stringsAsFactors=FALSE)
cases<-cases[,-1]
NM<- N_matrix


#Read in the covariates: HDI, pop  
area.level<-read.csv("/lustre03/project/6003552/mingchi/DistrictCovariates.csv")
load("/lustre03/project/6003552/mingchi/final_data_cleaned.Rdata")

# a spatial covariate: HDI (159 districts)
HDI<- area.level$HDI
# a temporal covariate: Temperature (284 weeks)
temp <- d$MeanMaxTemp[1:84]



count <- rep(NA,159)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}


######################################################
#---------- 1. rps for ZS-MSNBH model 
#####################################################


library(nimble)


# def functions
dhurdleNegBinMarg <- nimbleFunction(
  run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    if(x==0){
      logProb <- log(1-p)
    }else if(x>0){
      new.p<- r/(r+lambda)
      logProb <- log(p)+dnbinom(x, size=r ,prob= new.p, log=TRUE)-log(1-(new.p)^r)
    }else{
      logProb <- -Inf
    }
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

#now the rhurdlePoisMarg
rhurdleNegBinMarg <- nimbleFunction(
  run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
    returnType(double(0))
    if(n != 1) print("rmixunif only allows n = 1; using n = 1.")
    
    xt <- rbinom(n=1,size=1,prob=p)
    new.p<- r/(r+lambda)
    
    if(xt==1){
      xsim <- rnbinom(n=1, size=r, prob= new.p)
      while(xsim==0){
        xsim <- rnbinom(n=1, size=r, prob=new.p)
      }
    }else if(xt==0){
      xsim <- 0
    }
    return(xsim)
  })

assign('dhurdleNegBinMarg', dhurdleNegBinMarg, envir = .GlobalEnv)
assign('rhurdleNegBinMarg', rhurdleNegBinMarg, envir = .GlobalEnv)
dhurdleNegBinMargV <- Vectorize(dhurdleNegBinMarg)
rhurdleNegBinMargV <- Vectorize(rhurdleNegBinMarg)
assign('dhurdleNegBinMargV', dhurdleNegBinMargV, envir = .GlobalEnv)
assign('rhurdleNegBinMargV', rhurdleNegBinMargV, envir = .GlobalEnv)


# define X
cases<- as.matrix(cases)
X<- matrix(rep(0,159*Tsim),nrow=159)
X[cases>0] <- 1


# This is the key: we def a new variable out of nimble
NeiXsum<- matrix(rep(0, 159*T0), nrow = 159)
for (i in 1:159){
  for (t in 2:T0){
    for (j in which(NM[i,]==1)){
      NeiXsum[i,t]<-  NeiXsum[i,t]+X[j,t-1]
    }
  }
}



dengeeConsts <- list(N=159,
                     T=T0,
                     HDI=HDI,
                     temp=temp,
                     NeiXsum=NeiXsum)
dengeeData <- list(y=cases,X=X) 




parallel_nimble <- function(dengeeData, dengeeConsts,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin,N,T) {
  
  require(nimble)
  # first we define 2 functions for NB models
  
  
  dhurdleNegBinMarg <- nimbleFunction(
    run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                   log = integer(0, default = 0)) {
      
      returnType(double(0))
      if(x==0){
        logProb <- log(1-p)
      }else if(x>0){
        new.p<- r/(r+lambda)
        logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE)-log(1-(new.p)^r)
      }else{
        logProb <- -Inf
      }
      if(log) return(logProb)
      else return(exp(logProb)) 
    })
  
  #now the rhurdlePoisMarg
  rhurdleNegBinMarg <- nimbleFunction(
    run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
      returnType(double(0))
      if(n != 1) print("rmixunif only allows n = 1; using n = 1.")
      
      xt <- rbinom(n=1,size=1,prob=p)
      new.p<- r/(r+lambda)
      
      if(xt==1){
        xsim <- rnbinom(n=1, size=r, prob= new.p)
        while(xsim==0){
          xsim <- rnbinom(n=1, size=r, prob=new.p)
        }
      }else if(xt==0){
        xsim <- 0
      }
      return(xsim)
    })
  
  assign('dhurdleNegBinMarg', dhurdleNegBinMarg, envir = .GlobalEnv)
  assign('rhurdleNegBinMarg', rhurdleNegBinMarg, envir = .GlobalEnv)
  dhurdleNegBinMargV <- Vectorize(dhurdleNegBinMarg)
  rhurdleNegBinMargV <- Vectorize(rhurdleNegBinMarg)
  assign('dhurdleNegBinMargV', dhurdleNegBinMargV, envir = .GlobalEnv)
  assign('rhurdleNegBinMargV', rhurdleNegBinMargV, envir = .GlobalEnv)
  
  
  # we have changes in the model defining
  dengeeCode <- nimbleCode({
    #priors
    
    beta0 ~ dnorm(0,sd=100)
    beta1 ~ dnorm(0,sd=100)
    beta2 ~ dnorm(0, sd=100)
    
    alpha1 ~dnorm(0,sd=5)
    alpha2 ~dnorm(0,sd=5)
    alpha3 ~dnorm(0,sd=5)
    alpha4 ~dnorm(0,sd=5)
    alpha5 ~dnorm(0,sd=5)
    alpha6 ~dnorm(0,sd=5)
    
    r~dunif(0, 50)
    
    gamma1 ~ dnorm(0,sd=100)
    gamma2 ~ dnorm(0, sd=100)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        # logistics regs: p01, p11 and r
        logit(p01[i,t])<- alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
          alpha3*(temp[t]-mean(temp[1:T]))+
          gamma1*NeiXsum[i,t]
        
        logit(p11[i,t])<- alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ 
          alpha6*(temp[t]-mean(temp[1:T]))+
          gamma2*NeiXsum[i,t]
        
        
        # mu
        mup[i,t]<- exp(beta0+
                         beta1*(HDI[i]-mean(HDI[1:159]))+
                         beta2*(temp[t]-mean(temp[1:T])))
        # y
        y[i,t] ~ dhurdleNegBinMarg(p=p01[i,t]*(1-X[i,t-1])+p11[i,t]*X[i,t-1], 
                                   lambda= mup[i,t], r=r)
      }
    }
    
  })
  
  
  
  inits <- list(beta0 = rnorm(n=1,mean=0,sd=1),
                beta1 = rnorm(n=1,mean=0,sd=1),
                beta2 = rnorm(n=1,mean=0,sd=1) ,
                alpha1 = rnorm(n=1,mean=0,sd=1),
                alpha2 = rnorm(n=1,mean=0,sd=1) ,
                alpha3 = rnorm(n=1,mean=0,sd=1) ,
                alpha4 = rnorm(n=1,mean=0,sd=1),
                alpha5 = rnorm(n=1,mean=0,sd=1) ,
                alpha6 = rnorm(n=1,mean=0,sd=1) ,
                gamma1= rnorm(n=1,mean=0,sd=1) ,
                gamma2= rnorm(n=1,mean=0,sd=1), 
                r=runif(n=1, min=0, max=10))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  print(dengeeConf)
  dengeeConf$addMonitors(c("beta0","beta1","beta2",
                           "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6",
                           "gamma1","gamma2", "r"))
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),
                                   beta1 = rnorm(n=1,mean=0,sd=1),
                                   beta2 = rnorm(n=1,mean=0,sd=1) ,
                                   alpha1 = rnorm(n=1,mean=0,sd=1),
                                   alpha2 = rnorm(n=1,mean=0,sd=1) ,
                                   alpha3 = rnorm(n=1,mean=0,sd=1) ,
                                   alpha4 = rnorm(n=1,mean=0,sd=1),
                                   alpha5 = rnorm(n=1,mean=0,sd=1) ,
                                   alpha6 = rnorm(n=1,mean=0,sd=1) ,
                                   gamma1= rnorm(n=1,mean=0,sd=1) ,
                                   gamma2= rnorm(n=1,mean=0,sd=1), 
                                   r=runif(n=1, min=0, max=10))
  samples.ZSMSNBH <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                             ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNBH)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile=paste0("/lustre03/project/6003552/mingchi/ZSMSNBH_sim_2large_p_",as.character(T0),".txt"))
doParallel::registerDoParallel(cl)
clusterExport(cl, list( "parallel_nimble"),
              envir = globalenv())
fit_par <- foreach::foreach(i = 1:3,
                            .packages = c("nimble", "Rcpp")) %dopar% {
                              parallel_nimble(dengeeData, dengeeConsts,
                                              n_iterations = 80000,
                                              n_chains = 1,
                                              warm_up =30000,
                                              n_thin = 5,N=159,T=T0)
                            }
parallel::stopCluster(cl)


library(coda)
combine_par_mcmc <- function(fit_par) {
  n_par <- length(fit_par)
  n_chain <- 1
  out_lst <- list()
  if (n_chain == lengths(fit_par)[1]) {
    index <- 0
    for (i in 1:n_par) {
      for (j in 1:n_chain) {
        index <- index + 1
        trace_fit <- as.data.frame(do.call(rbind, fit_par[[i]][j]))
        out_lst[[index]] <- coda::mcmc(trace_fit)
      } }
    val <- coda::mcmc.list(out_lst)
  } else {
    val <- coda::mcmc.list(fit_par)
  }
  
  return(val)
}

ZSMSNBH.samples <- combine_par_mcmc(fit_par)


gelman.matrix<-gelman.diag(ZSMSNBH.samples[,c("beta0","beta1","beta2",
                                              "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6",
                                              "gamma1","gamma2", "r")],multivariate = FALSE)


write.csv(gelman.matrix$psrf,paste0("/lustre03/project/6003552/mingchi/rps_result/ZSMSNBH_gelman_sim_2large_p_",as.character(T0),".csv"))


####################
# calculate rps  ###
####################



samps.zsmsnbh <- data.frame(rbind(ZSMSNBH.samples[[1]],
                                  ZSMSNBH.samples[[2]],
                                  ZSMSNBH.samples[[3]]))


# 1. For ZS-MSNBH model 
M <- nrow(samps.zsmsnbh)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.zsmsnbh<- matrix(NA,nrow=159,ncol=K)


# rps
rps_zsmsnbh<- matrix(NA, nrow=159, ncol=K)
logs_zsmsnbh<- matrix(NA, nrow=159, ncol=K)


# spatial constants
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
num=as.numeric(as.carAdjacency(N_matrix)$num)
count=count
# X
X.h<- array(0,dim=c(159,K+1,M))# this is X[t-1], which is Xs of the previous 1 time 

NeiXsum_pred <- array(0,dim=c(159,K,M))


for(i in 1:159){
  X.h[i,1,]<- ifelse(rep(cases[i,T0],M)>0,1,0)# X of "T0"
}


# iterations: 
for (k in 1:K){
  print(k)
  #for t=T0+1,T0+2,...T0+K  
  for (i in 1:159){
    # the neighborhoods: 
    for (j in which(NM[i,]==1)){
      NeiXsum_pred[i,k,]<-  NeiXsum_pred[i,k,]+X.h[j,k,]
    }
    
    # parameters: 
    beta0<- as.numeric(unlist(samps.zsmsnbh[paste0("beta0")]))
    beta1<- as.numeric(unlist(samps.zsmsnbh[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zsmsnbh[paste0("beta2")]))
    alpha1<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha6")]))
    gamma1<- as.numeric(unlist(samps.zsmsnbh[paste0("gamma1")]))
    gamma2<- as.numeric(unlist(samps.zsmsnbh[paste0("gamma2")]))
    r<- as.numeric(unlist(samps.zsmsnbh[paste0("r")]))
    
    # mup
    mup.pred.zsmsnbh<- exp(beta0+
                             beta1*(HDI[i]-mean(HDI[1:159]))+
                             beta2*(temp[T0+k]-mean(temp[1:Tsim])))
    
    # p01
    p01.h<- expit(alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
                    alpha3*(temp[T0+k]-mean(temp[1:Tsim]))+gamma1*NeiXsum_pred[i,k,])
    # p11
    p11.h<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+
                    alpha6*(temp[T0+k]-mean(temp[1:Tsim]))+gamma2*NeiXsum_pred[i,k,])
    # r: 
    r_nb<- r
    
    # p: 
    p_x<- p01.h*(1-X.h[i,k,])+p11.h*X.h[i,k,]
    
    
    # X: 
    X.h[i,k+1,]<- rbinom(n=rep(1,M),
                         size=1,
                         p=p_x)
    
    # simulated case
    for(m in 1:M){
      if (X.h[i,k+1,m]==0){
        simulated.pred[i,k,m]<- 0
      }else{
        case.if<- rnbinom(1, prob=r_nb[m]/(r_nb[m]+mup.pred.zsmsnbh[m]), size=r_nb[m])
        while(case.if==0){
          case.if<- rnbinom(1, prob=r_nb[m]/(r_nb[m]+mup.pred.zsmsnbh[m]), size=r_nb[m])
        }
        simulated.pred[i,k,m]<-case.if
      }
    }
    # logscores
    logs_zsmsnbh[i,k] <- -log(mean(dhurdleNegBinMargV(x=rep(cases[i,T0+k],M),p = p_x,
                                                      lambda=mup.pred.zsmsnbh,r = r_nb,log = FALSE)))
    # rps 
    P_pred <- ecdf(simulated.pred[i,k,])(seq(0,10000,by=1))
    rps_zsmsnbh[i,k] <- sum((P_pred -ifelse(cases[i,T0+k]<=seq(0,10000,by=1),1,0))^2)
    
  }
}



# rps 
K1<-cbind(mean(rps_zsmsnbh[,1]))
K2<-cbind(mean(rps_zsmsnbh[,2]))
K3<-cbind(mean(rps_zsmsnbh[,3]))
K4<-cbind(mean(rps_zsmsnbh[,4]))

# logscores
L1<-cbind(mean(logs_zsmsnbh[,1]))
L2<-cbind(mean(logs_zsmsnbh[,2]))
L3<-cbind(mean(logs_zsmsnbh[,3]))
L4<-cbind(mean(logs_zsmsnbh[,4]))


dat.rps.h<-data.frame(rbind(K1,K2,K3,K4))
dat.logs.h<-data.frame(rbind(L1,L2,L3,L4))


row.names(dat.rps.h)<- c("K=1","K=2","K=3","K=4")
row.names(dat.logs.h)<- c("D=1","D=2","D=3","D=4")

write.csv(dat.rps.h,paste0("/lustre03/project/6003552/mingchi/rps_result/ZSMSNBH_rps_sim_2large_p_",as.character(T0),".csv"))
write.csv(dat.logs.h,paste0("/lustre03/project/6003552/mingchi/dss_result/ZSMSNBH_logs_sim_2large_p_",as.character(T0),".csv"))
