
library(nimble)
library(readxl)
library(coda)
library(utils)



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



set.seed(1993)



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


psi =cases
mpsi=mean(cases)
mlpsi=mean(log(cases+1))
lpsi=log(cases+1)

dengeeConsts <- list(N=159,T=384,
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
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/NBH_parallel.txt")
doParallel::registerDoParallel(cl)
clusterExport(cl, list( "parallel_nimble"),
              envir = globalenv())
fit_par <- foreach::foreach(i = 1:3,
                            .packages = c("nimble", "Rcpp")) %dopar% {
                              parallel_nimble(dengeeData, dengeeConsts,
                                              n_iterations = 80000,
                                              n_chains = 1,
                                              warm_up =30000,
                                              n_thin = 5,N=159,T=384)
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

save(NBH.samples, file="/Users/xumingchi/Desktop/hurdle-codes/NBH.RData")




#---------------------------------
# 2. WAIC for NBH model
#---------------------------------

load("/Users/xumingchi/Desktop/hurdle-codes/NBH.RData")


library(RMKdiscrete)
library(boot)

samps.NBH<- data.frame(rbind(NBH.samples[[1]],
                             NBH.samples[[2]],
                             NBH.samples[[3]]))


X<- matrix(rep(0,159*384),nrow=159,ncol=384)
X[cases>0] <- 1


lppd.h.nbh <- 0
pwaic.h.nbh <- 0
for(i in 1:159){
  print(i)
  for(t in 2:384){
    # state the parameters
    b0<- as.numeric(unlist(samps.NBH[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.NBH[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.NBH[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.NBH[paste0("beta2")]))
    
    alpha4<- as.numeric(unlist(samps.NBH[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.NBH[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.NBH[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.NBH[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.NBH[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.NBH[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.NBH[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.NBH[paste0("alpha11")]))
    
    gamma2<-  as.numeric(unlist(samps.NBH[paste0("gamma2")]))
    gamma3<-  as.numeric(unlist(samps.NBH[paste0("gamma3")]))
    # mup
    mup<- exp(b0)*cases[i,t-1]+exp(b+beta1*(sin.time[t]-mean(sin.time[]))+
                                     beta2*(cos.time[t]-mean(cos.time[])))
    
    # r 
    r.samp<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(cases[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    
    # p11
    p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ alpha6*(pop[i]-mean(pop[1:159]))+
                      alpha7*(log(cases[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159])))
    
    # waic
    lppd.h.nbh <- lppd.h.nbh + log(mean(dhurdleNegBinMarg(x = cases[i,t], p=p11, 
                                                          lambda= mup, r=r.samp )))
    pwaic.h.nbh<- pwaic.h.nbh+ var(log(dhurdleNegBinMarg(x = cases[i,t],p=p11,
                                                          lambda= mup, r=r.samp )))
    
  }
}

waic.nbh <- -2*(lppd.h.nbh-pwaic.h.nbh) 

# this is the end
