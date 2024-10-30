#####################################################
#########  ZS-MSNBH with spatial correlation  #######
#####################################################
#rm
rm(list=ls())

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


#create map
#read in shape file for map
shape_rj=readShapePoly("/Users/xumingchi/Desktop/paper/map/shape-RJ2-2.shp") 
shape_rj$oo <- seq(1:160)

nb2mat(poly2nb(shape_rj), style = "B",zero.policy = TRUE)# B is the basic binary coding
# this is the default order:
shape_rj$NOME

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


##############

# this is the order used in my project
datCHIKV$District

datCHIKV$District[59]<- 'Freguesia (Ilha)'
datCHIKV$District[120]<- 'Parque Colúmbia'
datCHIKV$District[105]<- 'Turiaçú'
datCHIKV$District[75]<- 'Tomás Coelho'
datCHIKV$District[108]<- 'Osvaldo Cruz'
datCHIKV$District[83]<- 'Todos os Santos'
datCHIKV$District[149]<- 'Gericinó'
#datCHIKV$District[14]<- 'Lapa' # this place does not appears in my dataset


# match 
trans.index<- rep(0, 160)

for (i in 1:160){
  trans.index[i]<-match(shape_rj$NOME[i], datCHIKV$District)
}


#seq(1,160,by=1), sort(trans.index)

# check the different names
datCHIKV$District[setdiff(seq(1,160,by=1), sort(trans.index))]
shape_rj$NOME[c(which(is.na(trans.index)))]


mat1<-nb2mat(poly2nb(shape_rj), style = "B",zero.policy = TRUE)
#mat<- mat1[1:159,1:159]
my.matrix<- matrix(rep(0,160*160), nrow = 160)

#[-14] Paqueta

# to change order, we deliberately set the 160 elemtnt to be 14. Take care to delete this
trans.index.art<- c(trans.index[-160],14)


ser<- seq(1,160,by=1)
ser1<- ser

for (k in 1:160){
  ser1[k]<- ser[trans.index.art[k]]
}



for (i in 1:160 ){
  mat2<-mat1[i,]
  mat3<- mat2
  # reorder mat3
  for (j in 1:160){
    mat3[j]<-mat2[which(trans.index.art==j)]
  }
  my.matrix[trans.index.art[i],]<- mat3
}


# my matrix is the ranked matrix. 
# randomly check some elements:
my.matrix[58,57]
my.matrix[57,58]
my.matrix[58,59]
my.matrix[59,58]

# then we delete one place, which is 14
my.matrix.adj<- my.matrix[-14,]
my.matrix.adj<- my.matrix.adj[,-14]

# Input our matrix to the neighborhood matrix. 
N_matrix<- my.matrix.adj




#check symmetric
#sum(N_matrix[1:160,2:161] != t(N_matrix[1:160,2:161]))
#N_matrix <- N_matrix[1:160,2:161]
sum(N_matrix[1:159,1:159] != t(N_matrix[1:159,1:159]))
min(rowSums(N_matrix))
#symmetric and every location has at least one neighbor (except for "Cidade Universitária")

NM <- as.matrix(N_matrix)

NM[65,67]<-1
NM[67,65]<-1
min(rowSums(NM))
N_matrix<- NM

count <- rep(NA,159)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}




#sum cases in neighboring areas
log_sum_nei_cases <- matrix(nrow=159,ncol=384)
for(i in 1:159){
  for(t in 2:384){
    nei <- which(N_matrix[i,]==1)
    log_sum_nei_cases[i,t] <- 0
    for(j in nei){
      log_sum_nei_cases[i,t] <- log_sum_nei_cases[i,t]+log(cases[j,t-1]+1)
    }
  }
}


hist(log_sum_nei_cases)
mean_log_sum_nei_cases <- mean(log_sum_nei_cases,na.rm = TRUE)

#sum of prevalence in neighboring areas
log_sum_nei_prev <- matrix(nrow=159,ncol=384)
for(i in 1:159){
  for(t in 2:384){
    nei <- which(N_matrix[i,]==1)
    log_sum_nei_prev[i,t] <- 0
    for(j in nei){
      log_sum_nei_prev[i,t] <- log_sum_nei_prev[i,t]+log((cases[j,t-1]/pop[j])+1)
    }
  }
}
hist(log_sum_nei_prev)
mean_log_sum_nei_prev <- mean(log_sum_nei_prev,na.rm = TRUE)



#need to calculate log of prevalence
prev <- matrix(nrow = 159,ncol=384)
for(i in 1:159){
  for(t in 1:384){
    prev[i,t] <- cases[i,t]/pop[i]
  }
}

hist(log(prev+1))

#log(cases) and mean of log(cases)
lpsi <- log(cases+1)
mlpsi <- mean(log(cases+1))


######################################
############# ZS-MSNBH model #########
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



#N_matrix<-read.csv("/Users/xumingchi/Desktop/paper/spatial_matrix.csv", 
#                   header=TRUE, stringsAsFactors=FALSE) 

set.seed(1993)
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



dengeeConsts <- list(N=159,T=384,
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
                            outfile="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_parallel_spatial.txt")
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

save(ZSMSNBH.samples, file="/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_spatial.RData")


load("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_spatial.RData")# load the samples





# check convergence
library(coda)
effectiveSize(ZSMSNBH.samples[,c("beta0","beta1","beta2","beta3","beta4","beta5","beta6","sigma_b","sigma_b0","rho",
                                        "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7")])
min(effectiveSize(ZSMSNBH.samples[,c("beta0","beta1","beta2","beta3","beta4","beta5","beta6","sigma_b","sigma_b0","rho")]))

gelman.diag(ZSMSNBH.samples[,c("beta0","beta1","beta2","beta3","beta4","beta5","beta6","sigma_b","sigma_b0","rho",
                                      "alpha1","alpha2","alpha3","alpha4","alpha5","alpha7")],multivariate = FALSE)
#randomly check location specific parameters
gelman.diag(ZSMSNBH.samples[,c("b[1]","b[2]","b[3]","b[4]",
                                      "b[21]","b[39]","b[112]","b[55]",
                                      "b[23]","b[29]","b[15]","b[74]",
                                      "b[41]","b[123]","b[16]","b[59]")],multivariate = FALSE)

gelman.diag(ZSMSNBH.samples[,c("b0[1]","b0[2]","b0[3]","b0[4]",
                                      "b0[21]","b0[39]","b0[112]","b0[55]",
                                      "b0[23]","b0[29]","b0[15]","b0[74]",
                                      "b0[41]","b0[123]","b0[16]","b0[59]")],multivariate = FALSE)






# summary


summary(ZSMSNBH.samples[,c("beta0","beta1","beta2","beta3","beta4","beta5","beta6","sigma_b","sigma_b0","rho",
                                  "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6","alpha7",
                                  "alpha8","alpha9","alpha10","alpha11","gamma1","gamma2",
                                  "gamma3","delta1","delta2")])



# persistence
summary(ZSMSNBH.samples[,c("alpha4","alpha5","alpha6","gamma2", "delta2","alpha7")])
# reemergence
summary(ZSMSNBH.samples[,c("alpha1","alpha2","alpha3","gamma1", "delta1")])






#---------------------------------
# 2. WAIC for ZS-MSNBH model
#---------------------------------


library(RMKdiscrete)
library(boot)

samps.ZSMSNBH<- data.frame(rbind(ZSMSNBH.samples[[1]],
                                ZSMSNBH.samples[[2]],
                                ZSMSNBH.samples[[3]]))

X<- matrix(rep(0,159*384),nrow=159,ncol=384)
X[cases>0] <- 1



lppd.h <- 0
pwaic.h <- 0
for(i in 1:2){
  print(i)
  for(t in 2:384){
    # state the parameters
    b0<- as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.ZSMSNBH[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))
    
    alpha1<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.ZSMSNBH[paste0("alpha11")]))
    
    gamma1<-  as.numeric(unlist(samps.ZSMSNBH[paste0("gamma1")]))
    gamma2<-  as.numeric(unlist(samps.ZSMSNBH[paste0("gamma2")]))
    gamma3<-  as.numeric(unlist(samps.ZSMSNBH[paste0("gamma3")]))
    delta1<-  as.numeric(unlist(samps.ZSMSNBH[paste0("delta1")]))
    delta2<-  as.numeric(unlist(samps.ZSMSNBH[paste0("delta2")]))
    # mu
    mup<- exp(b0)*cases[i,t-1]+exp(b+beta1*(sin.time[t]-mean(sin.time[]))+
                                     beta2*(cos.time[t]-mean(cos.time[])))
    
    # r 
    r.samp<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(cases[i,t-1]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    
    # p01 and p11
    p01<- inv.logit(alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
      alpha3*(pop[i]-mean(pop[1:159]))+gamma1*(greenarea[i]-mean(greenarea[1:159]))+
      delta1*NeiXsum[i,t])
    
    p11<- inv.logit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
      alpha7*(log(cases[i,t-1]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
      delta2*NeiXsum[i,t])
    
    # WAIC
    lppd.h <- lppd.h + log(mean(dhurdleNegBinMarg(x = cases[i,t],
                                        p=p01*(1-X[i,t-1])+p11*X[i,t-1],
                                        lambda= mup, r=r.samp )))
    pwaic.h<- pwaic.h+ var(log(dhurdleNegBinMarg(x = cases[i,t],
                                                 p=p01*(1-X[i,t-1])+p11*X[i,t-1],
                                                 lambda= mup, r=r.samp )))
    
  }
}

waic.zsmsnbh <- -2*(lppd.h-pwaic.h)




# posterior summary 

summary(ZSMSNBH.samples[,c("beta0","beta3","beta4","beta6","beta1","beta2","rho","beta5")])



# persistence

# odds for alpha5, alpha6, gamma2, delta2, 
round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha5")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha5")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha5")]))), 0.975),3)


round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha6")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha6")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha6")]))), 0.975),3)

round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma2")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma2")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma2")]))), 0.975),3)


round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta2")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta2")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta2")]))), 0.975),3)


round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha7")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha7")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha7")]))), 0.975),3)

# for alpha4 prob: 
round(mean(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha4")]+
                                     samps.ZSMSNBH[paste0("alpha7")]*(log(2)-0.277)
)))),3)

round(quantile(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha4")]+
                                         samps.ZSMSNBH[paste0("alpha7")]*(log(2)-0.277)
))),0.025),3)

round(quantile(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha4")]+
                                         samps.ZSMSNBH[paste0("alpha7")]*(log(2)-0.277)
))),0.975),3)


# reemergence
#summary(ZSMSNB.samples.dynamic[,c("alpha1","alpha2","alpha3","gamma1", "delta1")])

round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha2")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha2")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha2")]))), 0.975),3)

round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha3")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha3")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("alpha3")]))), 0.975),3)

round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma1")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma1")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("gamma1")]))), 0.975),3)

round(mean(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta1")])))),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta1")]))), 0.025),3)
round(quantile(as.numeric(unlist(exp(samps.ZSMSNBH[paste0("delta1")]))), 0.975),3)

# for alpha1 prob: 

round(mean(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha1")])))),3)
round(quantile(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha1")]))), 0.025),3)
round(quantile(as.numeric(unlist(expit(samps.ZSMSNBH[paste0("alpha1")]))), 0.975),3)



# this is the end
