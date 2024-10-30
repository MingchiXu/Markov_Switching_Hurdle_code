
#########################################################################
######  Fits the models for 5 mods to T0=231 vs. the actual cases  ######
#########################################################################


#-------------(1). For ZSMSNB model with fitting up to T0=231 -------------------
T0<-231
set.seed(T0)


# neighborhood matrix 
library(nimble)
library(ggplot2)
library(gridExtra)
library(boot)
library(latex2exp)
library(ggpubr)
library(spdep)
library(maptools)
library(leaflet)
library(MASS)


# compare to my order: 
library(readxl)
CHIKV2015 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2015.xlsx")
CHIKV2016 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2016.xlsx")
CHIKV2017 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2017.xlsx")
CHIKV2018 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2018.xlsx")
CHIKV2019 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2019.xlsx")
CHIKV2020 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2020.xlsx")
CHIKV2021 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2021.xlsx")
CHIKV2022 <- read_excel("/Users/xumingchi/Desktop/chikv-csv/CHIKV2022.xlsx")
# extra covariates data
area.level<-read.csv("/Users/xumingchi/Desktop/chikv-csv/DistrictCovariates.csv")


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


############

# we choose to delete "Paqueta", then 
dataCHIKV<-datCHIKV[-c(14),]
cases<-matrix(seq(1,159*384,by=1),nrow=159)
cases<-as.data.frame(cases)
for (i in 1:159){
  cases[i,]<-as.numeric(dataCHIKV[i,][-1])
}
cases<-as.matrix(cases)


library(nimble)

HDI <- area.level$HDI

#N_i 
N <- area.level$Npop

#population in thousands 
pop <- N/1000

# time and sin/cos of time
time<-seq(1,384,by=1)
sin.time<-sin(time*2*pi/52)
cos.time<-cos(time*2*pi/52)



N_matrix<-read.csv("/Users/xumingchi/Desktop/paper/spatial_matrix.csv", 
                   header=TRUE, stringsAsFactors=FALSE)

NM<-N_matrix

count <- rep(NA,159)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}



######################################
############# ZS-MSNB model ##########
######################################


library(nimble)

HDI <- area.level$HDI

#N_i 
N <- area.level$Npop

#population in thousands 
pop <- N/1000

# green area
greenarea<- area.level$green_area

# time and sin/cos of time
time<-seq(1,384,by=1)
sin.time<-sin(time*2*pi/52)
cos.time<-cos(time*2*pi/52)


#log(cases) and mean of log(cases)
lpsi <- log(cases+1)
mlpsi <- mean(log(cases+1))


# define X
X<- matrix(rep(NA,159*384),nrow=159,ncol=384)# change 0 to na
X[cases>0] <- 1

X.star<-X[which(is.na(X))]
length(X.star)


dengeeConsts <- list(N=159,T=T0,
                     pop=pop,HDI=HDI, greenarea=greenarea,
                     psi =cases,mpsi=mean(cases),mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time,
                     # add here for spatial correlation
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
  
  
  # we have changes in the model defining
  dengeeCode <- nimbleCode({
    #priors
    beta0 ~ dnorm(0,sd=100)
    beta1 ~dnorm(0,sd=100)
    beta2 ~ dnorm(0, sd=100)
    beta3 ~ dnorm(0, sd=100)
    beta4 ~ dnorm(0, sd=100)
    beta5 ~ dnorm(0, sd=100)
    beta6 ~ dnorm(0, sd=100)
    
    alpha1 ~dnorm(0,sd=100)
    alpha2 ~ dnorm(0, sd=100)
    alpha3 ~ dnorm(0, sd=100)
    alpha4 ~dnorm(0,sd=100)
    alpha5 ~ dnorm(0, sd=100)
    alpha6 ~ dnorm(0, sd=100)# for previous cases
    alpha7 ~ dnorm(0, sd=100)
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~ dnorm(0, sd=100)
    alpha11 ~ dnorm(0, sd=100)
    gamma1 ~ dnorm(0, sd=100)# green area
    gamma2 ~ dnorm(0, sd=100)
    gamma3 ~ dnorm(0, sd=100)
    delta1 ~ dnorm(0, sd=100)# spatial coef
    delta2 ~ dnorm(0, sd=100)
    # gen r for NB
    #for (i in 1:N){
    #  lr[i] ~ dnorm(lru[i],sd=lrsd[i])
    #  r[i]<- exp(lr[i]) 
    #}
    #r ~ dunif(0,100)# add this
    rho~dnorm(0,sd=100)
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)# add betas
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,sd=sigma_b)
    }
    sigma_b ~ dunif(0,10)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      X[i,1] ~ dbern(0.5)
      for (t in 2:T){
        # cal Neighbor sum for t=t-1
        Snei[num[i],i,(t-1)] <- X[adj[count[i]],(t-1)]
        for(j in 2:num[i]){
          Snei[num[i]-j+1,i,(t-1)] <- Snei[num[i]-j+2,i,(t-1)] + X[adj[count[i]+j-1],(t-1)]
        }
        
        # logistics regs: p01, p11 and r
        logit(p01[i,t])<- alpha1+alpha2*(HDI[i]-mean(HDI[1:N]))+
          alpha3*(pop[i]-mean(pop[1:N]))+gamma1*(greenarea[i]-mean(greenarea[1:N]))+
          delta1*Snei[1,i,(t-1)]
        
        logit(p11[i,t])<- alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+alpha6*(pop[i]-mean(pop[1:N]))+
          alpha7*(log(psi[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:N]))+
          delta2*Snei[1,i,(t-1)]
        
        X[i,t]~ dbern(X[i,t-1]*p11[i,t]+(1-X[i,t-1])*p01[i,t])# by markov chain
        
        log(r[i,t])<- alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
          alpha10*(pop[i]-mean(pop[1:N]))+
          alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N]))
      }
    }
    X_final[1:159] <- X[1:159,T]
    for(i in 1:159) {
      for(t in 2:T){
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
        
        
        p[i,t] <- r[i,t]/(r[i,t]+(X[i,t])*mup[i,t]) - 1e-10*(1-X[i,t])
        y[i,t] ~ dnegbin(prob=p[i,t],size=r[i,t]) 
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
                sigma_b = .01,prec_b0 = 1,
                rho=rnorm(n=1,mean=0,sd=1),
                b=rnorm(n=159,mean=0,sd=.01),
                b0=rnorm(n=159,mean=0,sd=.01),
                #r= runif(n=1,min=0,max=100),
                X=matrix(rbinom(n=159*384,size=1,prob=.5),nrow=159))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                           "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                           "gamma3","delta1","delta2","sigma_b0","X_final",
                           "rho","b","sigma_b","b0"))
  
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
                                   gamma1= rnorm(n=1,mean=0,sd=.01),
                                   gamma2= rnorm(n=1,mean=0,sd=.01),
                                   gamma3= rnorm(n=1,mean=0,sd=.01),
                                   delta1= rnorm(n=1,mean=0,sd=.01),
                                   delta2= rnorm(n=1,mean=0,sd=.01),
                                   sigma_b = .01,prec_b0 = 1,
                                   rho=rnorm(n=1,mean=0,sd=1),
                                   b=rnorm(n=159,mean=0,sd=.01),
                                   b0=rnorm(n=159,mean=0,sd=.01),
                                   #r= runif(n=1,min=0,max=100),
                                   X=matrix(rbinom(n=159*384,size=1,prob=.5),nrow=159))
  samples.ZSMSNB.d <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                              ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNB.d)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNB_spatial_231.txt")
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

ZSMSNB.samples.dynamic <- combine_par_mcmc(fit_par)

save(ZSMSNB.samples.dynamic, file="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNB_231.RData")
















#----------------------------------(2). For ZSMSNBH model with fitting up to T0=231 -------------------


######################################
############# ZS-MSNBH model #########
######################################


set.seed(T0)
# def functions
dhurdleNegBinMarg <- nimbleFunction(
  run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    if(x==0){
      logProb <- log(1-p)
    }else if(x>0){
      new.p<- r/(r+lambda)
      logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE) -pnbinom(q = 0,size = r,prob = new.p,lower.tail = FALSE,log.p = TRUE)
    }else{
      logProb <- -Inf
    }
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

#now the rhurdlePoisMarg
rhurdleNegBinMarg <- nimbleFunction(
  run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
    returnType(integer(0))
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



# define X
X<- matrix(rep(0,159*384),nrow=159,ncol=384)
X[cases>0] <- 1

# This is the key: we def a new variable out of nimble
NeiXsum<- matrix(rep(0, 159*384), nrow = 159)
for (i in 1:159){
  for (t in 2:384){
    for (j in which(NM[i,]==1)){
      NeiXsum[i,t]<-  NeiXsum[i,t]+X[j,t-1]
    }
  }
}



dengeeConsts <- list(N=159,T=T0,
                     pop=pop,HDI=HDI, greenarea=greenarea,NeiXsum=NeiXsum,
                     psi =cases,mpsi=mean(cases),mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time,
                     # add here for spatial correlation
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
        logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE) -pnbinom(q = 0,size = r,prob = new.p,lower.tail = FALSE,log.p = TRUE)
      }else{
        logProb <- -Inf
      }
      if(log) return(logProb)
      else return(exp(logProb)) 
    })
  
  #now the rhurdlePoisMarg
  rhurdleNegBinMarg <- nimbleFunction(
    run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
      returnType(integer(0))
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
  
  
  # we have changes in the model defining
  dengeeCode <- nimbleCode({
    #priors
    beta0 ~ dnorm(0,sd=100)
    beta1 ~dnorm(0,sd=100)
    beta2 ~ dnorm(0, sd=100)
    beta3 ~ dnorm(0, sd=100)
    beta4 ~ dnorm(0, sd=100)
    beta5 ~ dnorm(0, sd=100)
    beta6 ~ dnorm(0, sd=100)
    
    alpha1 ~dnorm(0,sd=100)
    alpha2 ~ dnorm(0, sd=100)
    alpha3 ~ dnorm(0, sd=100)
    alpha4 ~dnorm(0,sd=100)
    alpha5 ~ dnorm(0, sd=100)
    alpha6 ~ dnorm(0, sd=100)
    alpha7 ~ dnorm(0, sd=100)
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~ dnorm(0, sd=100)
    alpha11 ~ dnorm(0, sd=100)
    gamma1 ~ dnorm(0, sd=100)# green area
    gamma2 ~ dnorm(0, sd=100)
    gamma3 ~ dnorm(0, sd=100)
    delta1 ~ dnorm(0, sd=100)# spatial coef
    delta2 ~ dnorm(0, sd=100)
    rho~dnorm(0,sd=100)
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)# add betas
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,prec_b)
    }
    prec_b ~ dgamma(.1,.1)
    sigma_b <- 1/sqrt(prec_b)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        # logistics regs: p01, p11 and r
        logit(p01[i,t])<- alpha1+alpha2*(HDI[i]-mean(HDI[1:N]))+
          alpha3*(pop[i]-mean(pop[1:N]))+gamma1*(greenarea[i]-mean(greenarea[1:N]))+
          delta1*NeiXsum[i,t]
        
        logit(p11[i,t])<- alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+ alpha6*(pop[i]-mean(pop[1:N]))+
          alpha7*(log(psi[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:N]))+
          delta2*NeiXsum[i,t]
        
        
        log(r[i,t])<- alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
          alpha10*(pop[i]-mean(pop[1:N]))+
          alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N]))
        
        # mu
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
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
                b0=rnorm(n=159,mean=0,sd=.01))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                           "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                           "gamma3","delta1","delta2","sigma_b0",
                           "rho","b","sigma_b","b0"))
  
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
                                   b0=rnorm(n=159,mean=0,sd=.01))
  samples.ZSMSNBH <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                             ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNBH)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_parallel_231.txt")
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
#end_time <- Sys.time()
#end_time - start_time
#run_time <- end_time-start_time

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

save(ZSMSNBH.samples, file="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_231.RData")







#----------------------------------(3). For NB model with fitting up to T0=231 -------------------




# set seed T0
set.seed(T0)
# define X, 384/ T0

dengeeConsts <- list(N=159,T=T0,
                     pop=pop,
                     HDI=HDI,
                     greenarea=greenarea,
                     psi =cases,
                     mpsi=mean(cases),
                     mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time)
dengeeData <- list(y=cases) 


parallel_nimble <- function(dengeeData, dengeeConsts,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin,N,T) {
  
  require(nimble)
  
  
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
    rho ~ dnorm(0,sd=100)
    
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~dnorm(0, sd=100)
    alpha11 ~dnorm(0, sd=100)
    
    gamma3 ~ dnorm(0, sd=100)
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,sd=sigma_b)
    }
    
    sigma_b ~ dunif(0,10)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        # r
        r[i,t]<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
                       alpha10*(pop[i]-mean(pop[1:N]))+
                       alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N])))
        # for the mean
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
        
        p[i,t] <- r[i,t]/(r[i,t]+mup[i,t]) 
        y[i,t] ~ dnegbin(prob=p[i,t],size=r[i,t]) 
        
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
                rho=rnorm(n=1,mean=0,sd=1),
                alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                gamma3= rnorm(n=1,mean=0,sd=.01) ,
                sigma_b = .01,prec_b0 = 1,
                b=rnorm(n=159,mean=0,sd=.01),
                b0=rnorm(n=159,mean=0,sd=.01))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha8","alpha9","alpha10","alpha11",
                           "gamma3","sigma_b0","rho","b","sigma_b","b0"))
  
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),
                                   beta1 = rnorm(n=1,mean=0,sd=.01),
                                   beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta3 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta4 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta5 = rnorm(n=1,mean=0,sd=.01) ,
                                   beta6 = rnorm(n=1,mean=0,sd=.01) ,
                                   rho=rnorm(n=1,mean=0,sd=1),
                                   alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                                   alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                                   gamma3= rnorm(n=1,mean=0,sd=.01) ,
                                   sigma_b = .01,prec_b0 = 1,
                                   b=rnorm(n=159,mean=0,sd=.01),
                                   b0=rnorm(n=159,mean=0,sd=.01))
  
  samples.NB <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                        ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.NB)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/NB_parallel_231.txt")
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

NB.samples<- combine_par_mcmc(fit_par)

save(NB.samples, file="/Users/xumingchi/Desktop/hurdle-codes/NB_231.RData")












#----------------------------------(4). For NBH model with fitting up to T0=231 -------------------


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


# def functions
dhurdleNegBinMarg <- nimbleFunction(
  run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    if(x==0){
      logProb <- log(1-p)
    }else if(x>0){
      new.p<- r/(r+lambda)
      logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE) -pnbinom(q = 0,size = r,prob = new.p,lower.tail = FALSE,log.p = TRUE)
    }else{
      logProb <- -Inf
    }
    if(log) return(logProb)
    else return(exp(logProb)) 
  })

#now the rhurdlePoisMarg
rhurdleNegBinMarg <- nimbleFunction(
  run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
    returnType(integer(0))
    #if(n != 1) print("rmixunif only allows n = 1; using n = 1.")
    
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



set.seed(T0)

dengeeConsts <- list(N=159,T=T0,
                     pop=pop,
                     HDI=HDI,
                     greenarea=greenarea,
                     psi =cases,
                     mpsi=mean(cases),
                     mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time)
dengeeData <- list(y=cases) 


parallel_nimble <- function(dengeeData, dengeeConsts,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin,N,T) {
  
  require(nimble)
  
  # def functions
  dhurdleNegBinMarg <- nimbleFunction(
    run = function(x = double(0), p = double(0), lambda = double(0), r= double(0),
                   log = integer(0, default = 0)) {
      
      returnType(double(0))
      if(x==0){
        logProb <- log(1-p)
      }else if(x>0){
        new.p<- r/(r+lambda)
        logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE) -pnbinom(q = 0,size = r,prob = new.p,lower.tail = FALSE,log.p = TRUE)
      }else{
        logProb <- -Inf
      }
      if(log) return(logProb)
      else return(exp(logProb)) 
    })
  
  #now the rhurdlePoisMarg
  rhurdleNegBinMarg <- nimbleFunction(
    run = function(n = integer(0),  p = double(0), lambda = double(0), r=double(0)) {
      returnType(integer(0))
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
    rho ~ dnorm(0,sd=100)
    
    alpha4 ~ dnorm(0, sd=100)
    alpha5 ~ dnorm(0, sd=100)
    alpha6 ~dnorm(0, sd=100)
    alpha7 ~dnorm(0, sd=100)
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~dnorm(0, sd=100)
    alpha11 ~dnorm(0, sd=100)
    
    gamma2 ~ dnorm(0, sd=100)
    gamma3 ~ dnorm(0, sd=100)
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,sd=sigma_b)
    }
    
    sigma_b ~ dunif(0,10)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        # r
        r[i,t]<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
                       alpha10*(pop[i]-mean(pop[1:N]))+
                       alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N])))
        # p11
        p11[i,t]<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+ alpha6*(pop[i]-mean(pop[1:N]))+
                           alpha7*(log(psi[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:N])))
        
        
        # for the mean
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
        #p[i,t] <- r[i,t]/(r[i,t]+mup[i,t]) 
        y[i,t] ~ dhurdleNegBinMarg(p=p11[i,t], 
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
                rho=rnorm(n=1,mean=0,sd=1),
                alpha4 = rnorm(n=1,mean=0,sd=.01) ,
                alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                gamma2= rnorm(n=1,mean=0,sd=.01) ,
                gamma3= rnorm(n=1,mean=0,sd=.01) ,
                sigma_b = .01,prec_b0 = 1,
                b=rnorm(n=159,mean=0,sd=.01),
                b0=rnorm(n=159,mean=0,sd=.01))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  #dengeeConf$resetMonitors()
  
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha4","alpha5","alpha6","alpha7", 
                           "alpha8","alpha9","alpha10","alpha11", "gamma2",
                           "gamma3","sigma_b0","rho","b","sigma_b","b0"))
  
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function()  list(beta0 = rnorm(n=1,mean=0,sd=1),
                                    beta1 = rnorm(n=1,mean=0,sd=.01),
                                    beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta3 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta4 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta5 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta6 = rnorm(n=1,mean=0,sd=.01) ,
                                    rho=rnorm(n=1,mean=0,sd=1),
                                    alpha4 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                                    gamma2= rnorm(n=1,mean=0,sd=.01) ,
                                    gamma3= rnorm(n=1,mean=0,sd=.01) ,
                                    sigma_b = .01,prec_b0 = 1,
                                    b=rnorm(n=159,mean=0,sd=.01),
                                    b0=rnorm(n=159,mean=0,sd=.01))
  
  samples.NBH <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                         ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.NBH)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/NBH_parallel_231.txt")
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

NBH.samples<- combine_par_mcmc(fit_par)

save(NBH.samples, file="/Users/xumingchi/Desktop/hurdle-codes/NBH_231.RData")










#----------------------------------(5). For ZINB model with fitting up to T0=231 -------------------



# time and sin/cos of time
time<-seq(1,384,by=1)
sin.time<-sin(time*2*pi/52)
cos.time<-cos(time*2*pi/52)



# define X, 384/ T0
X<- matrix(rep(NA,159*384),nrow=159,ncol=384)
X[cases>0] <- 1


dengeeConsts <- list(N=159,T=T0,
                     pop=pop,
                     HDI=HDI,
                     greenarea=greenarea,
                     psi =cases,
                     mpsi=mean(cases),
                     mlpsi=mean(log(cases+1)),lpsi=log(cases+1),
                     sin.time=sin.time,
                     cos.time=cos.time)
dengeeData <- list(y=cases, X=X) 


parallel_nimble <- function(dengeeData, dengeeConsts,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin,N,T) {
  
  require(nimble)
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
    rho ~ dnorm(0,sd=100)
    
    alpha4 ~ dnorm(0, sd=100)
    alpha5 ~ dnorm(0, sd=100)
    alpha6 ~dnorm(0, sd=100)
    alpha7 ~dnorm(0, sd=100)
    alpha8 ~ dnorm(0, sd=100)
    alpha9 ~ dnorm(0, sd=100)
    alpha10 ~dnorm(0, sd=100)
    alpha11 ~dnorm(0, sd=100)
    
    gamma2 ~ dnorm(0, sd=100)
    gamma3 ~ dnorm(0, sd=100)
    
    for(i in 1:N){
      b0[i]~dnorm(beta0+beta3*(pop[i]-mean(pop[1:N]))+beta4*(HDI[i]-mean(HDI[1:N]))+
                    beta6*(greenarea[i]-mean(greenarea[1:N])),prec_b0)
      b[i]~dnorm(rho+(log(pop[i])-mean(log(pop[1:N])))*beta5,sd=sigma_b)
    }
    
    sigma_b ~ dunif(0,10)
    prec_b0 ~ dgamma(.1,.1)
    sigma_b0 <- 1/sqrt(prec_b0)
    
    #likelihood
    for (i in 1:N){
      for (t in 2:T){
        # r
        r[i,t]<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
                       alpha10*(pop[i]-mean(pop[1:N]))+
                       alpha11*(log(psi[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N])))
        # p11
        p11[i,t]<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+ alpha6*(pop[i]-mean(pop[1:N]))+
                           alpha7*(log(psi[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:N])))
        
        # X
        X[i,t]~ dbern(p11[i,t])
        # X_final
        X_final[1:159] <- X[1:159,T]
        # for the mean
        mup[i,t]<- exp(b0[i])*psi[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
        
        p[i,t] <- r[i,t]/(r[i,t]+(X[i,t])*mup[i,t]) - 1e-10*(1-X[i,t])
        
        # y
        y[i,t] ~ dnegbin(prob=p[i,t],size=r[i,t]) 
        
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
                rho=rnorm(n=1,mean=0,sd=1),
                alpha4 = rnorm(n=1,mean=0,sd=.01) ,
                alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                gamma2= rnorm(n=1,mean=0,sd=.01) ,
                gamma3= rnorm(n=1,mean=0,sd=.01) ,
                sigma_b = .01,prec_b0 = 1,
                b=rnorm(n=159,mean=0,sd=.01),
                b0=rnorm(n=159,mean=0,sd=.01),
                X=matrix(rbinom(n=159*384,size=1,prob=.5),nrow=159))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  dengeeConf$removeSampler(c("sigma_b"))
  dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  
  print(dengeeConf)
  #dengeeConf$resetMonitors()
  
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha4","alpha5","alpha6","alpha7", 
                           "alpha8","alpha9","alpha10","alpha11", "gamma2",
                           "gamma3","sigma_b0","rho","b","sigma_b","b0","X_final"))# "X_final"
  
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function()  list(beta0 = rnorm(n=1,mean=0,sd=1),
                                    beta1 = rnorm(n=1,mean=0,sd=.01),
                                    beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta3 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta4 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta5 = rnorm(n=1,mean=0,sd=.01) ,
                                    beta6 = rnorm(n=1,mean=0,sd=.01) ,
                                    rho=rnorm(n=1,mean=0,sd=1),
                                    alpha4 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha5 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha6 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha7 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha8 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha9 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha10 = rnorm(n=1,mean=0,sd=.01) ,
                                    alpha11 = rnorm(n=1,mean=0,sd=.01) ,
                                    gamma2= rnorm(n=1,mean=0,sd=.01) ,
                                    gamma3= rnorm(n=1,mean=0,sd=.01) ,
                                    sigma_b = .01,prec_b0 = 1,
                                    b=rnorm(n=159,mean=0,sd=.01),
                                    b0=rnorm(n=159,mean=0,sd=.01),
                                    X=matrix(rbinom(n=159*384,size=1,prob=.5),nrow=159))
  
  samples.ZINB <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                          ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZINB)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/ZINB_parallel_231.txt")
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

ZINB.samples<- combine_par_mcmc(fit_par)

save(ZINB.samples, file="/Users/xumingchi/Desktop/hurdle-codes/ZINB_231.RData")


