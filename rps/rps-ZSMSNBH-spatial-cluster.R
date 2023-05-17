#########################################
#  rps for ZS-MSNBH with spatial and r  #
#########################################


# rps cluster for ZS-MSNBH model

# cluster path: /lustre03/project/6003552/mingchi # change 
# git: store the code
args=commandArgs(trailingOnly = TRUE)
# packages
library(nimble)
library(readxl)
library(coda)
library(utils)

if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}

#sim_num<- args[1]
K<-4 
T0<-as.numeric(args[1])

print(T0)

# 1. load datasets
CHIKV2015 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2015.xlsx")
CHIKV2016 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2016.xlsx")
CHIKV2017 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2017.xlsx")
CHIKV2018 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2018.xlsx")
CHIKV2019 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2019.xlsx")
CHIKV2020 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2020.xlsx")
CHIKV2021 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2021.xlsx")
CHIKV2022 <- read_excel("/lustre03/project/6003552/mingchi/CHIKV2022.xlsx")
area.level<-read.csv("/lustre03/project/6003552/mingchi/DistrictCovariates.csv")

# then merge these data
datCHIKV<-merge.data.frame(CHIKV2015,CHIKV2016,by.x ="District",by.y="District",
                           suffixes = c("/15","/16") ,sort = F)
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2017,by.x ="District",by.y="District",sort = F)
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2018,by.x ="District",by.y="District",
                           suffixes = c("/17","/18") ,sort = F)
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2019,by.x ="District",by.y="District",sort = F)
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2020,by.x ="District",by.y="District",
                           suffixes = c("/19","/20") ,sort = F)
colnames(datCHIKV)[length(colnames(datCHIKV))]<-"53/20"
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2021,by.x ="District",by.y="District", 
                           suffixes=c("/20","/21"),sort = F)
datCHIKV<-merge.data.frame(datCHIKV,CHIKV2022,by.x ="District",by.y="District",
                           suffixes = c("/21","/22") ,sort = F)
datCHIKV<-datCHIKV[,-c(386:418)]
colnames(datCHIKV)[c(334:366)]<-c(paste0(seq(from=20,to=52,by=1),"/21"))
dataCHIKV<-datCHIKV[-c(14),]
cases<-matrix(seq(1,159*384,by=1),nrow=159)
cases<-as.data.frame(cases)
for (i in 1:159){
  cases[i,]<-as.numeric(dataCHIKV[i,][-1])
}
cases<-as.matrix(cases)
# change cases to 1:T0
#cases<-cases[,c(1:T0)]



# spatial matrix
N_matrix<-read.csv("/lustre03/project/6003552/mingchi/spatial_matrix.csv", 
                   header=TRUE, stringsAsFactors=FALSE)
#N_matrix<-read.csv("/Users/xumingchi/Desktop/paper/spatial_matrix.csv", 
#header=TRUE, stringsAsFactors=FALSE)

count <- rep(NA,159)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}




#sum cases in neighboring areas
#log_sum_nei_cases <- matrix(nrow=159,ncol=384)
#for(i in 1:159){
#  for(t in 2:384){
#    nei <- which(N_matrix[i,]==1)
#    log_sum_nei_cases[i,t] <- 0
#    for(j in nei){
#      log_sum_nei_cases[i,t] <- log_sum_nei_cases[i,t]+log(cases[j,t-1]+1)
#    }
#  }
#}


#hist(log_sum_nei_cases)
#mean_log_sum_nei_cases <- mean(log_sum_nei_cases,na.rm = TRUE)

#sum of prevalence in neighboring areas
#log_sum_nei_prev <- matrix(nrow=159,ncol=384)
#for(i in 1:159){
#  for(t in 2:384){
#    nei <- which(N_matrix[i,]==1)
#    log_sum_nei_prev[i,t] <- 0
#    for(j in nei){
#      log_sum_nei_prev[i,t] <- log_sum_nei_prev[i,t]+log((cases[j,t-1]/pop[j])+1)
#    }
#  }
#}
#hist(log_sum_nei_prev)
#mean_log_sum_nei_prev <- mean(log_sum_nei_prev,na.rm = TRUE)



#need to calculate log of prevalence
#prev <- matrix(nrow = 159,ncol=384)
#for(i in 1:159){
#  for(t in 1:384){
#    prev[i,t] <- cases[i,t]/pop[i]
#  }
#}

#hist(log(prev+1))


#---------------------------------------
#  ZS-MSNBH with spatial, modeling r
#---------------------------------------

# covariates
HDI <- area.level$HDI
#N_i 
N <- area.level$Npop
#population in thousands people
pop <- N/1000
# the green area
greenarea<- area.level$green_area


# time and sin/cos of time
time<-seq(1,384,by=1)
sin.time<-sin(time*2*pi/52)
cos.time<-cos(time*2*pi/52)


mcos.time <- mean(cos.time[1:384])
msin.time <- mean(sin.time[1:384])

#log(cases) and mean of log(cases)
lpsi <- log(cases+1)
mlpsi <- mean(log(cases+1))


# functions:

# def functions
dhurdleNegBinMarg <- nimbleFunction(
  run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    if(x==0){
      logProb <- log(1-p)
    }else if(x>0){
      new.p<- r/(r+lambda)
      logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE)
      - log(1-(new.p)^r)
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
      #xsim <- rpois(n=1,lambda=lambda)
      xsim <- rnbinom(n=1, size=r, prob= new.p)
      while(xsim==0){
        #xsim <- rpois(n=1,lambda=lambda)
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



# set seed T0
set.seed(T0)

# define X
X<- matrix(rep(0,159*384),nrow=159,ncol=384)
X[cases>0] <- 1


# This is the key: we def a new variable out of nimble

NM<-N_matrix

NeiXsum<- matrix(rep(0, 159*384), nrow = 159)
for (i in 1:159){
  for (t in 2:384){
    for (j in which(NM[i,]==1)){
      NeiXsum[i,t]<-  NeiXsum[i,t]+X[j,t-1]
    }
  }
}




#X.star<-X[which(is.na(X))]

dengeeConsts <- list(N=159,T=T0,
                     pop=pop,
                     HDI=HDI,
                     greenarea=greenarea,
                     NeiXsum=NeiXsum,
                     psi =cases,
                     mpsi=mean(cases),
                     mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time,
                     msin.time=msin.time,
                     mcos.time=mcos.time,
                     adj=as.numeric(as.carAdjacency(N_matrix)$adj),
                     num=as.numeric(as.carAdjacency(N_matrix)$num),
                     count=count)
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
        logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE)
        - log(1-(new.p)^r)
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
        #xsim <- rpois(n=1,lambda=lambda)
        xsim <- rnbinom(n=1, size=r, prob= new.p)
        while(xsim==0){
          #xsim <- 1 if goes too long
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
    beta0 ~ dnorm(0, sd=100)
    beta1 ~ dnorm(0, sd=100)
    beta2 ~ dnorm(0, sd=100)
    beta3 ~ dnorm(0, sd=100)
    beta4 ~ dnorm(0, sd=100)
    beta5 ~ dnorm(0, sd=100)
    beta6 ~ dnorm(0, sd=100)
    
    alpha1 ~ dnorm(0, sd=100)
    alpha2 ~ dnorm(0, sd=100)
    alpha3 ~ dnorm(0, sd=100)
    alpha4 ~ dnorm(0, sd=100)
    alpha5 ~ dnorm(0, sd=100)
    alpha6 ~ dnorm(0, sd=100)# for previous cases
    alpha7 ~ dnorm(0, sd=100)
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~dnorm(0, sd=100)
    alpha11 ~dnorm(0, sd=100)
    gamma1 ~ dnorm(0, sd=100)# green area
    gamma2 ~ dnorm(0, sd=100)
    gamma3 ~ dnorm(0, sd=100)
    delta1 ~ dnorm(0, sd=100)# spatial coef
    delta2 ~ dnorm(0, sd=100)
    rho~dnorm(0,sd=100)
    
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,prec_b)
    }
    prec_b ~ dgamma(.1,.1)
    sigma_b <- 1/sqrt(prec_b)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        
        # logistics regression: p01, p11 and r
        logit(p01[i,t])<- alpha1+alpha2*(HDI[i]-mean(HDI[1:N]))+
          alpha3*(pop[i]-mean(pop[1:N]))+gamma1*(greenarea[i]-mean(greenarea[1:N]))+
          delta1*NeiXsum[i,t]
        
        logit(p11[i,t])<- alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+alpha6*(pop[i]-mean(pop[1:N]))+
          alpha7*(log(psi[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:N]))+
          delta2*NeiXsum[i,t]
        
        log(r[i,t])<- alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
          alpha10*(pop[i]-mean(pop[1:N]))+
          alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N]))
        
        # mu
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]- msin.time )+
                                               beta2*(cos.time[t]- mcos.time ))
        # y
        y[i,t] ~ dhurdleNegBinMarg(p=p01[i,t]*(1-X[i,t-1])+p11[i,t]*X[i,t-1], 
                                   lambda= mup[i,t], r=r[i,t])
       
      }
    }
    
  })
  inits <- list(beta0 = rnorm(n=1,mean=0,sd=1),
                beta1 = rnorm(n=1,mean=0,sd=.01),
                beta2 = rnorm(n=1,mean=0,sd=.01) ,
                beta3 = rnorm(n=1,mean=0,sd=.01) ,
                beta4 = rnorm(n=1,mean=0,sd=.01) ,
                beta5 = rnorm(n=1,mean=0,sd=.01) ,
                beta6 = rnorm(n=1,mean=0,sd=.01) ,
                alpha1 = rnorm(n=1,mean=0,sd=.01),
                alpha2 = rnorm(n=1,mean=0,sd=.01) ,
                alpha3 = rnorm(n=1,mean=0,sd=.01) ,
                alpha4 = rnorm(n=1,mean=0,sd=.01),
                alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                gamma1= rnorm(n=1,mean=0,sd=.01) ,
                gamma2= rnorm(n=1,mean=0,sd=.01) ,
                gamma3= rnorm(n=1,mean=0,sd=.01) ,
                delta1= rnorm(n=1,mean=0,sd=.01) ,
                delta2= rnorm(n=1,mean=0,sd=.01) ,
                prec_b = 1,prec_b0 = 1,
                rho=rnorm(n=1,mean=0,sd=1),
                b=rnorm(n=159,mean=0,sd=.01),
                b0=rnorm(n=159,mean=0,sd=.01))# X size 384/T0
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  #dengeeConf$resetMonitors()
  
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                           "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                           "gamma3","delta1","delta2","sigma_b0",
                           "rho","b","sigma_b","b0"))# "X_final"
  
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),
                                   beta1 = rnorm(n=1,mean=0,sd=.01),
                                   beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta3 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta4 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta5 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta6 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha1 = rnorm(n=1,mean=0,sd=.01),
                                   alpha2 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha3 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha4 = rnorm(n=1,mean=0,sd=.01),
                                   alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                                   gamma1= rnorm(n=1,mean=0,sd=.01) ,
                                   gamma2= rnorm(n=1,mean=0,sd=.01) ,
                                   gamma3= rnorm(n=1,mean=0,sd=.01) ,
                                   delta1= rnorm(n=1,mean=0,sd=.01) ,
                                   delta2= rnorm(n=1,mean=0,sd=.01) ,
                                   prec_b = 1,prec_b0 = 1,
                                   rho=rnorm(n=1,mean=0,sd=1),
                                   b=rnorm(n=159,mean=0,sd=.01),
                                   b0=rnorm(n=159,mean=0,sd=.01))# X size 384/T0
  
  samples.ZSMSNBH <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                              ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNBH)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile=paste0("/lustre03/project/6003552/mingchi/ZSMSNBH_parallel_dynamic_",T0,".txt"))

#cl <- parallel::makeCluster(3,
#                            outfile=paste0("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_parallel_dynamic_",T0,".txt"))


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

ZSMSNBH.samples<- combine_par_mcmc(fit_par)


gelman.matrix<-gelman.diag(ZSMSNBH.samples[,c("beta0","beta1","beta2","beta3","beta4","beta5","beta6","sigma_b","sigma_b0","rho",
                                                     "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                                                     "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                                                     "gamma3","delta1","delta2")],multivariate = FALSE)


write.csv(gelman.matrix$psrf,paste0("/lustre03/project/6003552/mingchi/rps_result/gelman_ZSMSNBH_T0_",as.character(T0),".csv"))




samps.zsmsnbh <- data.frame(rbind(ZSMSNBH.samples[[1]],
                                  ZSMSNBH.samples[[2]],
                                  ZSMSNBH.samples[[3]]))


# 1. For ZS-MSNBH model 
M <- nrow(samps.zsmsnbh)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.zsmsnbh<- matrix(0,nrow=159,ncol=K)


# rps
#LS_zip<- matrix(0,nrow=159,ncol=K) 
rps_zsmsnbh<- matrix(0, nrow=159, ncol=K)



# spatial constants
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
num=as.numeric(as.carAdjacency(N_matrix)$num)
count=count
# X
X.h<- array(0,dim=c(159,K,M))# this is X[t-1], which is Xs of the previous 1 time 



for (i in 1:159){
  print(i)
  
 
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K){
    print(k)
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.zsmsnbh[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.zsmsnbh[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.zsmsnbh[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zsmsnbh[paste0("beta2")]))
    alpha1<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.zsmsnbh[paste0("alpha11")]))
    
     
    gamma1<- as.numeric(unlist(samps.zsmsnbh[paste0("gamma1")]))
    gamma2<- as.numeric(unlist(samps.zsmsnbh[paste0("gamma2")]))
    gamma3<- as.numeric(unlist(samps.zsmsnbh[paste0("gamma3")]))
    
    delta1<- as.numeric(unlist(samps.zsmsnbh[paste0("delta1")]))
    delta2<- as.numeric(unlist(samps.zsmsnbh[paste0("delta2")]))
    
    
    # mup
    if (k== 1){
      mup.pred.zsmsnbh<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred-msin.time )+beta2*(cos.time.pred- mcos.time ))
    }else {
      mup.pred.zsmsnbh<-  exp(b0)*simulated.pred[i,k-1,]+exp(b+beta1*(sin.time.pred-msin.time )+beta2*(cos.time.pred-mcos.time ))
    }
    # p01 and p11 should be changed into these, and r: 
    
    
    
    # p01
    p01.h<- expit(alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
                  alpha3*(pop[i]-mean(pop[1:159]))+gamma1*(greenarea[i]-mean(greenarea[1:159]))+
                  delta1*NeiXsum[i,k+T0])
    
    
    # p11: change the parameter: 
    if(k==1){
      p11.h<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(cases[i,T0]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*NeiXsum[i,k+T0])
    }else{
      p11.h<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*NeiXsum[i,k+T0])
    }
    
    # r: change the parameters: 
    
    if (k==1){
      r_nb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(cases[i,T0]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      r_nb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    # X: 
    if (k==1){
      X.h[i,1,]<- ifelse(rep(cases[i,T0],M)>0,1,0)# X of "T0"
    }else{
      X.h[i,k,]<-ifelse(simulated.pred[i,k-1,]>0,1,0)# X of "T0+1",...,"TO+K-1"
    }
    
    p_nb<- p01.h*(1-X.h[i,k,])+p11.h*X.h[i,k,]
    
    # simulated y
    simulated.pred[i,k,]<-rhurdleNegBinMargV(n=rep(1,M), 
                                           p =p_nb, 
                                           lambda= mup.pred.zsmsnbh,
                                           r=r_nb) 
    
    
    # rps
    P_pred <- ecdf(simulated.pred[i,k,])
    rps_zsmsnbh[i,k] <- sum((P_pred(seq(0,1000,by=1))-ifelse(cases[i,T0+k]<=seq(0,1000,by=1),1,0))^2)#T=mT
    
    
  }
}

 

K1<-cbind(mean(rps_zsmsnbh[,1]))
K2<-cbind(mean(rps_zsmsnbh[,2]))
K3<-cbind(mean(rps_zsmsnbh[,3]))
K4<-cbind(mean(rps_zsmsnbh[,4]))

dat.rps<-data.frame(rbind(K1,K2,K3,K4))


row.names(dat.rps)<- c("K=1","K=2","K=3","K=4")

write.csv(dat.rps,paste0("/lustre03/project/6003552/mingchi/rps_result/ZSMSNBH_T0_",as.character(T0),".csv"))
