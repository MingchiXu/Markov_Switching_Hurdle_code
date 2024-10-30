#=============================
# Missing data in time series 
#=============================

rm(list=ls())

library(gridExtra)
library(boot)
library(latex2exp)
library(ggpubr)
library(spdep)
library(maptools)
library(leaflet)
library(MASS)
library(nimble)
library(readxl)

set.seed(1998)
# Load Data ---------------------------------------------------------------

CHIKV2015 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2015.xlsx")
CHIKV2016 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2016.xlsx")
CHIKV2017 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2017.xlsx")
CHIKV2018 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2018.xlsx")
CHIKV2019 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2019.xlsx")
CHIKV2020 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2020.xlsx")
CHIKV2021 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2021.xlsx")
CHIKV2022 <- read_excel("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/true_case/data/CHIKV2022.xlsx")

# extra covariates data
area.level<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/DistrictCovariates.csv")


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

# we choose to delete "Paqueta", then 
dataCHIKV<-datCHIKV[-c(14),]
cases<-matrix(seq(1,159*384,by=1),nrow=159)
cases<-as.data.frame(cases)
for (i in 1:159){
  cases[i,]<-as.numeric(dataCHIKV[i,][-1])
}

# convert to matrix
cases<-as.matrix(cases)

T0<- 80


# the covariates 

HDI <- area.level$HDI

#N_i 
N <- area.level$Npop

#population in thousands 
pop <- N/1000

# green area
greenarea<- area.level$green_area

# time and sin/cos of time
time<-seq(1,80,by=1)
sin.time<-sin(time*2*pi/52)
cos.time<-cos(time*2*pi/52)



#cases<-read.csv('C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases.csv')
N_matrix<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/spatial_matrix.csv", 
                   header=TRUE, stringsAsFactors=FALSE)

NM<- as.matrix(N_matrix)


count <- rep(NA,159)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}



# def functions
dhurdleNegBinMarg <- nimbleFunction(
  run = function(x = integer(0), p = double(0), lambda = double(0), r= double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    if(x==0){
      logProb <- log(1-p)
    }else if(x>0){
      new.p<- r/(r+lambda)
      logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE)-pnbinom(q = 0,size = r, prob = new.p, lower.tail = FALSE, log.p = TRUE)
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




# indicators customized functions 

indicatorY <- nimbleFunction(
  run = function(yvalue = integer(0)) {
    returnType(double(0))
    if (yvalue<0){
      print('The count must be positive.')
    }else{
      if(yvalue>0){
        ind=1
      }else{
        ind=0
      }
    }
    return(ind)
  })

assign('indicatorY', indicatorY, envir = .GlobalEnv)
indicatorYV <- Vectorize(indicatorY)


########################################################################################
#select_index<-sample(1:384, 3, replace=FALSE)
#cases[,select_index]<-NA

# Assign Missing values in the middle part of the time serial 
# we first cut the time series to make it shorter (to increase the fitting time) to be 200 points 
# we treat the 68th to the 70th time points as missing 

shortcases<- cases[,1:80] # 0, 1, 3 are missing cases (roughly at 30%, 50%, 75%)
#shortcases[5, 25]<- NA # we want to do inference on the missing values on these points 
#shortcases[5, 46]<- NA # we want to do inference on the missing values on these points 
#shortcases[5, 62]<- NA # we want to do inference on the missing values on these points 
shortcases[5, 57:59]<- NA
shortcases[5, 61:62]<- NA
shortcases[5, 71:72]<- NA


# the initial value 
y.init<- ifelse(is.na(shortcases), 0, shortcases)


dengeeConsts <- list(N=159,T=T0,
                     pop=pop,HDI=HDI, greenarea=greenarea,
                     mlpsi=mean(log(na.omit(shortcases)+1)),
                     sin.time=sin.time,
                     cos.time=cos.time,
                     # add here for spatial correlation
                     adj=as.numeric(as.carAdjacency(N_matrix)$adj),
                     num=as.numeric(as.carAdjacency(N_matrix)$num),
                     count=count)

dengeeData <- list(y=shortcases) 




parallel_nimble <- function(dengeeData, dengeeConsts,
                            n_iterations,
                            n_chains,
                            warm_up,
                            n_thin,N,T) {
  
  require(nimble)
  # first we define 2 functions for NB models
  dhurdleNegBinMarg <- nimbleFunction(
    run = function(x = integer(0), p = double(0), lambda = double(0), r= double(0),
                   log = integer(0, default = 0)) {
      
      returnType(double(0))
      if(x==0){
        logProb <- log(1-p)
      }else if(x>0){
        new.p<- r/(r+lambda)
        logProb <- log(p)+ dnbinom(x, size=r ,prob= new.p, log=TRUE)-pnbinom(q = 0,size = r, prob = new.p, lower.tail = FALSE, log.p = TRUE)
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
  # indicator for Y
  indicatorY <- nimbleFunction(
    run = function(yvalue = integer(0)) {
      returnType(double(0))
      if (yvalue<0){
        print('The count must be positive.')
      }else{
        if(yvalue>0){
          ind=1
        }else{
          ind=0
        }
      }
      return(ind)
    })
  
  assign('indicatorY', indicatorY, envir = .GlobalEnv)
  assign('dhurdleNegBinMarg', dhurdleNegBinMarg, envir = .GlobalEnv)
  assign('rhurdleNegBinMarg', rhurdleNegBinMarg, envir = .GlobalEnv)
  dhurdleNegBinMargV <- Vectorize(dhurdleNegBinMarg)
  rhurdleNegBinMargV <- Vectorize(rhurdleNegBinMarg)
  indicatorYV <- Vectorize(indicatorY)
  assign('dhurdleNegBinMargV', dhurdleNegBinMargV, envir = .GlobalEnv)
  assign('rhurdleNegBinMargV', rhurdleNegBinMargV, envir = .GlobalEnv)
  assign('indicatorYV', indicatorYV, envir = .GlobalEnv)
  
  
  # let nimble know we use the discrete distribution 
  registerDistributions(list(
    dhurdleNegBinMarg = list(
      BUGSdist = "dhurdleNegBinMarg(p, lambda, r)",
      types = c('value = integer(0)', 'p = double(0)', 'lambda = double(0)', 'r = double(0)'),
      discrete = TRUE
    )
  ))
  
  # we have changes in the model defining
  dengeeCode <- nimbleCode({
    #priors
    beta0 ~ dnorm(0,sd=10)
    beta1 ~dnorm(0,sd=10)
    beta2 ~ dnorm(0, sd=10)
    beta3 ~ dnorm(0, sd=10)
    beta4 ~ dnorm(0, sd=10)
    beta5 ~ dnorm(0, sd=10)
    beta6 ~ dnorm(0, sd=10)
    
    alpha1 ~dnorm(0,sd=5)
    alpha2 ~ dnorm(0, sd=5)
    alpha3 ~ dnorm(0, sd=5)
    alpha4 ~dnorm(0,sd=5)
    alpha5 ~ dnorm(0, sd=5)
    alpha6 ~ dnorm(0, sd=5)
    alpha7 ~ dnorm(0, sd=5)
    alpha8 ~ dnorm(0, sd=5)
    alpha9 ~ dnorm(0, sd=5)
    alpha10 ~ dnorm(0, sd=5)
    alpha11 ~ dnorm(0, sd=5)
    gamma1 ~ dnorm(0, sd=10)# green area
    gamma2 ~ dnorm(0, sd=10)
    gamma3 ~ dnorm(0, sd=10)
    delta1 ~ dnorm(0, sd=10)# spatial coef
    delta2 ~ dnorm(0, sd=10)
    
    #r~dunif(0, 50)
    rho~dnorm(0,sd=10)
    #X[5, 68]~ dbern(0.5)
    #X[5, 69]~ dbern(0.5)
    #X[5, 70]~ dbern(0.5)
  
    
    # the random effect 
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
    
    # MC for the X
    #X[5, 25] ~ dbern(0.5)
    #X[5, 46] ~ dbern(0.5)
    #X[5, 62] ~ dbern(0.5)
    
    
    #y[5, 25] ~ dpois(1)
    #y[5, 46] ~ dpois(1)
    #y[5, 62] ~ dpois(2)
    # the missing indicators
   
    
    #likelihood before y[5,68]
    for (i in 1:159){
      for (t in 2:T){
        # cal Neighbor sum for t=t-1
        Snei[num[i],i,(t-1)] <- indicatorY(y[adj[count[i]],(t-1)])
        for(j in 2:num[i]){
          Snei[num[i]-j+1,i,(t-1)] <- Snei[num[i]-j+2,i,(t-1)] + indicatorY( y[adj[count[i]+j-1],(t-1)])
        }
        # logistics regs: p01, p11 and r
        logit(p01[i,t])<- alpha1+alpha2*(HDI[i]-mean(HDI[1:N]))+
          alpha3*(pop[i]-mean(pop[1:N]))+gamma1*(greenarea[i]-mean(greenarea[1:N]))+
          #delta1*NeiXsum[i,t]
          delta1*Snei[1,i,(t-1)]
        
        logit(p11[i,t])<- alpha4+alpha5*(HDI[i]-mean(HDI[1:N]))+ alpha6*(pop[i]-mean(pop[1:N]))+
          alpha7*(log(y[i,t-1]+1)- mlpsi )+gamma2*(greenarea[i]-mean(greenarea[1:N]))+
          #delta2*NeiXsum[i,t]
          delta2*Snei[1,i,t-1]
        
        log(r[i,t])<- alpha8+ alpha9*(HDI[i]-mean(HDI[1:N]))+
          alpha10*(pop[i]-mean(pop[1:N]))+
          alpha11*(log(y[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:N]))
        
        # mu
        mup[i,t]<- exp(b0[i])*y[i,t-1]+exp(b[i]+beta1*(sin.time[t]-mean(sin.time[]))+
                                               beta2*(cos.time[t]-mean(cos.time[])))
        
        # y
        y[i,t] ~ dhurdleNegBinMarg(p=p01[i,t]*(1- indicatorY(y[i,t-1]))+p11[i,t]* indicatorY(y[i,t-1]), 
                                   lambda= mup[i,t], 
                                   r=r[i,t])
        
      }
    }
    
    #  }
    #}
    
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
                b0=rnorm(n=159,mean=0,sd=.01),
                y=y.init 
                #X= X.init
                )
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  #sample sd on log scale
  #dengeeConf$removeSampler(c("sigma_b"))
  #dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
  # add a sampler named 'y'
  #dengeeConf$removeSampler(c("y[5, 25]", "y[5, 46]", "y[5, 62]" ))
  
  
  #dengeeConf$addSampler(target = "y[5, 25]", type = "slice")
  #dengeeConf$addSampler(target = "y[5, 46]", type = "slice")
  #dengeeConf$addSampler(target = "y[5, 62]", type = "slice")
  
  print(dengeeConf)
  dengeeConf$addMonitors(c("beta0","beta1","beta2","beta3","beta4","beta5","beta6",
                           "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                           "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                           "gamma3","delta1","delta2","sigma_b0",
                           "rho","b","sigma_b","b0",
                           #"y[5, 25]","y[5, 46]","y[5, 62]"
                           "y[5, 57]", 'y[5, 58]', 'y[5, 59]', 
                           'y[5, 61]', 'y[5, 62]', 'y[5, 71]', 'y[5, 72]'
                           ))
  
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
                                   b0=rnorm(n=159,mean=0,sd=.01),
                                   y=y.init 
                                   #X= X.init
                                   )
  samples.ZSMSNBH <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                             ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNBH)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile=paste0("C:/Users/xumc7/OneDrive/Desktop/Simulation/rps/ZSMSNBH_Missing_middle.txt"))
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

ZSMSNBH.samples.miss <- combine_par_mcmc(fit_par)



save(ZSMSNBH.samples.miss, file="C:/Users/xumc7/OneDrive/Desktop/SiM_v1/Missing_middle/hurdle-new.RData")

load("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/Missing_middle/hurdle-new.RData")



samps.MSH<- data.frame(rbind(ZSMSNBH.samples.miss[[1]],
                                            ZSMSNBH.samples.miss[[2]],
                                            ZSMSNBH.samples.miss[[3]]))



# the samples 
y_5_57<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",57,".")]))
y_5_58<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",58,".")]))
y_5_59<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",59,".")]))

y_5_61<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",61,".")]))
y_5_62<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",62,".")]))

y_5_71<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",71,".")]))
y_5_72<-as.numeric(unlist(samps.MSH[paste0("y.",5,"..",72,".")]))

# posterior mean and 95% credible intervals 
mean.5.57<-mean(y_5_57)
upper.5.57<- quantile(y_5_57, 0.975)
lower.5.57<- quantile(y_5_57, 0.025)
mean.5.57# 0.0189

mean.5.58<-mean(y_5_58)
upper.5.58<- quantile(y_5_58, 0.975)
lower.5.58<- quantile(y_5_58, 0.025)
mean.5.58 # 0.01636667

mean.5.59<-mean(y_5_59)
upper.5.59<- quantile(y_5_59, 0.975)
lower.5.59<- quantile(y_5_59, 0.025)
mean.5.59 # 0.01636667


mean.5.61<-mean(y_5_61)
upper.5.61<- quantile(y_5_61, 0.975)
lower.5.61<- quantile(y_5_61, 0.025)
mean.5.61 # 0.01636667

mean.5.62<-mean(y_5_62)
upper.5.62<- quantile(y_5_62, 0.975)
lower.5.62<- quantile(y_5_62, 0.025)
mean.5.62 # 0.7403667

mean.5.71<-mean(y_5_71)
upper.5.71<- quantile(y_5_71, 0.975)
lower.5.71<- quantile(y_5_71, 0.025)
mean.5.71 # 0.7403667


mean.5.72<-mean(y_5_72)
upper.5.72<- quantile(y_5_72, 0.975)
lower.5.72<- quantile(y_5_72, 0.025)
mean.5.72 # 0.7403667




# make a plot to show the mean and 95% CI 
# build a data frame to summarize the data 
timepts<- seq(1, 80, by=1)
true.val<- cases[5,1:80]
# posterior mean
posterior.mean<- c(cases[5, 1:56],
                   mean.5.57,
                   mean.5.58,
                   mean.5.59,
                   cases[5, 60], 
                   mean.5.61,
                   mean.5.62, 
                   cases[5, 63:70],
                   mean.5.71,
                   mean.5.72,
                   cases[5, 73:80])
# upper 95% 
upper.95<- c(cases[5, 1:56],
             upper.5.57,
             upper.5.58,
             upper.5.59,
             cases[5, 60], 
             upper.5.61,
             upper.5.62, 
             cases[5, 63:70],
             upper.5.71,
             upper.5.72,
             cases[5, 73:80])
# lower 95% 
lower.95<- c(cases[5, 1:56],
             lower.5.57,
             lower.5.58,
             lower.5.59,
             cases[5, 60], 
             lower.5.61,
             lower.5.62, 
             cases[5, 63:70],
             lower.5.71,
             lower.5.72,
             cases[5, 73:80])
# Missing index 
Missing<- factor(c(rep('Observed', 56), rep('Not Observed', 3), 'Observed',
                   rep('Not Observed', 2), rep('Observed', 8), rep('Not Observed',2), 
                   rep('Observed', 8)))

# the data frame 
miss_5<- data.frame(timepts,
                    true.val,
                    posterior.mean,
                    upper.95,
                    lower.95, 
                    Missing)

colnames(miss_5)<- c('weeks',"true", 'post_mean', 'upper', 'lower', 'Missing')
# a selected smaller ones 
short_miss_5<- miss_5[40:80,]
short_miss_5$solid_dash<- c(rep(1, 16), 2,2, 3, 4, 4, 
                             rep(5, 1), 10, 10, 
                            rep(6, 6), 7,7, 8,8, 
                            rep(9, 7))
short_miss_5$extra<- c(rep(1, 17), rep(2, 3), 
                       7, 7, rep(5, 1), 
                       rep(3, 8), rep(6, 2), 
                       rep(4, 8))

short_miss_5$group_connect<- c(rep(1, 17), rep(2, 3), 
                            7, rep(5, 2), 
                            rep(3, 8), rep(6, 2), 
                            rep(4, 8))

short_miss_5$add_mean<- c(rep(0, 17),1, 1, 1, 0, 1, 1, rep(0, 8), 1, 1, rep(0, 8) )

# plots 
library(ggplot2)
library(cowplot)
Missing_plt<-ggplot(data=short_miss_5)+
  # add the points 
  geom_point(data=subset(short_miss_5, add_mean==1), aes(x=weeks, y=post_mean, fill='Posterior Mean'), 
             color='black', size=5, shape=4)+
  # add the observed 
  geom_point(aes(x=weeks, y=true, shape=Missing), size=4)+
  geom_line(data=subset(short_miss_5, group_connect==1), aes(x=weeks, y=true), linetype='solid')+
  geom_line(data=subset(short_miss_5, group_connect==3), aes(x=weeks, y=true), linetype='solid')+
  geom_line(data=subset(short_miss_5, group_connect==4), aes(x=weeks, y=true), linetype='solid')+
  
  # upper
  geom_line(data=subset(short_miss_5, group_connect==2), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, group_connect==5), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, group_connect==6), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, solid_dash==2), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, solid_dash==4), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, solid_dash==7), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, solid_dash==8), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, solid_dash==10), aes(x=weeks, y=upper), linetype='dashed')+
  geom_line(data=subset(short_miss_5, extra==7), aes(x=weeks, y=upper), linetype='dashed')+
  # lower 
  geom_line(data=subset(short_miss_5, group_connect==2), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, group_connect==5), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, group_connect==6), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, solid_dash==2), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, solid_dash==4), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, solid_dash==7), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, solid_dash==8), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, solid_dash==10), aes(x=weeks, y=lower), linetype='longdash')+
  geom_line(data=subset(short_miss_5, extra==7), aes(x=weeks, y=lower), linetype='dashed')+
  
  scale_shape_manual(values=c(1,16))+
  labs(title='Missing in Centro district', x='Weeks', y='Cases')+
  theme_cowplot()+
  theme(axis.text = element_text(size=25),
        axis.title = element_text(size=25),
        text=element_text(size=25))+
  theme(legend.position=c(0.08, 0.95), 
        legend.justification = c(0.01, 1))+
  theme(legend.text=element_text(size=25))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 5.2)))+
  theme(legend.margin = margin(-0.5,0,0,0, unit="cm"))

Missing_plt

ggsave(filename = 'C:/Users/xumc7/OneDrive/Desktop/Simulation/Figures/Figure_Missing.eps',
       Missing_plt, 
       width=20, height=11, dpi=600, , unit='in', device=cairo_ps, 
       fallback_resolution=600)






