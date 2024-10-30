###################################################################
# RPS for Markov Switching models: ZS-MSNB
# large reporting rate: p0=1
###################################################################

args=commandArgs(trailingOnly = TRUE)

library(nimble)
library(readxl)
library(coda)
library(utils)
library(scoringRules)

if (length(args)==0){
  stop("At least one argument must be supplied", call.=FALSE)
}


# the hyper parameter
T0<- as.numeric(args[1])
print(T0)
K<-4 
#p0<- 1
Tsim<- 84
# set seed T0
set.seed(T0)

# read in the data files 

cases<-read.csv('/lustre03/project/6003552/mingchi/cases_full.csv')
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


# rps cluster for ZS-MSNB model



#log(cases) and mean of log(cases)
lpsi <- log(cases+1)
mlpsi <- mean(log(cases+1))


# define X, 384/ T0
X<- matrix(rep(NA,159*T0),nrow=159,ncol=T0)
X[cases[1:159,1:T0]>0] <- 1


# counts 
count<- rep(NA, 159)
num<- as.numeric(as.carAdjacency(N_matrix)$num)
count[1]<- 1
for (i in 1:158){
  count[i+1]<- count[i] + num[i]
}



#X.star<-X[which(is.na(X))]



dengeeConsts <- list(N=159,mT=T0,
                     HDI=HDI, 
                     temp=temp,
                     adj=as.numeric(as.carAdjacency(N_matrix)$adj),
                     num=as.numeric(as.carAdjacency(N_matrix)$num),
                     count=count)

dengeeData <- list(y=cases[1:159,1:T0],X=X) 
mT <- T0


parallel_nimble <- function( dengeeData, dengeeConsts,
                             n_iterations,
                             n_chains,
                             warm_up,
                             n_thin,N,mT) {
  
  require(nimble)
  
  
  # we have changes in the model defining
  dengeeCode <- nimbleCode({
    #priors
    beta0 ~ dnorm(0,sd=100)
    beta1 ~ dnorm(0, sd=100)
    beta2 ~ dnorm(0, sd=100)
    
    alpha1 ~ dnorm(0,sd=5)
    alpha2 ~ dnorm(0,sd=5)
    alpha3 ~ dnorm(0,sd=5)
    alpha4 ~ dnorm(0,sd=5)
    alpha5 ~ dnorm(0,sd=5)
    alpha6 ~ dnorm(0,sd=5)
    
    gamma1 ~ dnorm(0, sd=1)
    gamma2 ~ dnorm(0, sd=1)
    
    r ~ dunif(0,50)# add this
    
    
    #likelihood
    for (i in 1:N){
      X[i,1] ~ dbern(0.5)
      for (t in 2:mT){
        # the neighbor sum 
        Snei[num[i],i,(t-1)] <- X[adj[count[i]],(t-1)]
        for(j in 2:num[i]){
          Snei[num[i]-j+1,i,(t-1)] <- Snei[num[i]-j+2,i,(t-1)] + X[adj[count[i]+j-1],(t-1)]
        }
        
        # logistics regs: p01, p11 and r
        logit(p01[i,t])<- alpha1 + alpha2*(HDI[i]-mean(HDI[1:159]))+
          alpha3*(temp[t]-mean(temp[1:mT]))+
          gamma1*Snei[1, i, (t-1)]
        
        logit(p11[i,t])<- alpha4 + alpha5*(HDI[i]-mean(HDI[1:159]))+ 
          alpha6*(temp[t]-mean(temp[1:mT]))+
          gamma2*Snei[1, i, (t-1)]
        
        # by markov chain
        p1[i,t] <- X[i,t-1]*p11[i,t]+(1-X[i,t-1])*p01[i,t]
        X[i,t]~ dbern(p1[i,t])
      }
    }
    # X_final
    X_final[1:159] <- X[1:159,mT]
    ###################################
    # for mean 
    for(i in 1:N) {
      for(t in 2:mT){
        # mu
        mup[i,t]<- exp(beta0+
                         beta1*(HDI[i]-mean(HDI[1:159]))+
                         beta2*(temp[t]-mean(temp[1:mT])))
        
        p_nbb[i,t] <- r/(r+(X[i,t])*mup[i,t])-1e-10*(1-X[i,t])
        y[i,t] ~ dnegbin(prob=p_nbb[i,t],size=r)
        
        
      }
    }
  })
  inits <- list(beta0 = rnorm(n=1,mean=0,sd=1),
                beta1 = rnorm(n=1,mean=0,sd=1),
                beta2 = rnorm(n=1,mean=0,sd=1),
                alpha1 = rnorm(n=1,mean=0,sd=1),
                alpha2 = rnorm(n=1,mean=0,sd=1),
                alpha3 = rnorm(n=1,mean=0,sd=1),
                alpha4 = rnorm(n=1,mean=0,sd=1),
                alpha5 = rnorm(n=1,mean=0,sd=1),
                alpha6 = rnorm(n=1,mean=0,sd=1),
                gamma1 = rnorm(n=1,mean=0,sd=1),
                gamma2 = rnorm(n=1,mean=0,sd=1),
                r= runif(n=1,min=0,max=10),
                X=matrix(rbinom(n=159*mT,size=1,prob=.5),nrow=159))
  
  dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)
  
  Cdengee <- compileNimble(dengeemodel)
  
  dengeeConf <- configureMCMC(dengeemodel, print = TRUE)
  
  print(dengeeConf)
  #dengeeConf$resetMonitors()
  
  
  dengeeConf$addMonitors(c("beta0", "beta1", 'beta2',
                           "alpha1", 'alpha2', 'alpha3' ,"alpha4", 'alpha5', 'alpha6', 'gamma1', 'gamma2',
                           "r","X_final"))
  
  dengeeMCMC <- buildMCMC(dengeeConf)
  
  CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)
  
  
  initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),
                                   beta1 = rnorm(n=1,mean=0,sd=1),
                                   beta2 = rnorm(n=1,mean=0,sd=1),
                                   alpha1 = rnorm(n=1,mean=0,sd=1),
                                   alpha2 = rnorm(n=1,mean=0,sd=1),
                                   alpha3 = rnorm(n=1,mean=0,sd=1),
                                   alpha4 = rnorm(n=1,mean=0,sd=1),
                                   alpha5 = rnorm(n=1,mean=0,sd=1),
                                   alpha6 = rnorm(n=1,mean=0,sd=1),
                                   gamma1 = rnorm(n=1,mean=0,sd=1),
                                   gamma2 = rnorm(n=1,mean=0,sd=1),
                                   r= runif(n=1,min=0,max=10),
                                   X=matrix(rbinom(n=159*mT,size=1,prob=.5),nrow=159))
  
  samples.ZSMSNB <- runMCMC(CdengeeMCMC,  niter =n_iterations,nchains = n_chains,nburnin=warm_up
                            ,samplesAsCodaMCMC = TRUE,thin=n_thin,inits = initsFunction)
  return(samples.ZSMSNB)
}

library(parallel)
library(doParallel)
parallel::detectCores()
cl <- parallel::makeCluster(3,
                            outfile=paste0("/lustre03/project/6003552/mingchi/ZSMSNB_sim_large_p_",T0,".txt"))
doParallel::registerDoParallel(cl)
clusterExport(cl, list( "parallel_nimble"),
              envir = globalenv())
fit_par <- foreach::foreach(i = 1:3,
                            .packages = c("nimble", "Rcpp")) %dopar% {
                              parallel_nimble(dengeeData, dengeeConsts,
                                              n_iterations = 80000,
                                              n_chains = 1,
                                              warm_up =30000,
                                              n_thin = 5,N=159,mT=T0)
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

ZSMSNB.samples<- combine_par_mcmc(fit_par)


gelman.matrix<-gelman.diag(ZSMSNB.samples[,c("beta0","beta1","beta2", "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6",
                                             "gamma1","gamma2","r")],multivariate = FALSE)


write.csv(gelman.matrix$psrf,paste0("/lustre03/project/6003552/mingchi/rps_result/ZSMSNB_gelman_sim_large_p_",as.character(T0),".csv"))



######################################################
#---------- calculate rps for ZI model ----
#####################################################


samps.zsmsnb <- data.frame(rbind(ZSMSNB.samples[[1]],
                                 ZSMSNB.samples[[2]],
                                 ZSMSNB.samples[[3]]))


# 1. For ZS-MSNB model 
M <- nrow(samps.zsmsnb)
# for M samples
simulated.pred<- array(NA,dim=c(159,K,M))
# pred cases
simulated.pred.zsmsnb<- matrix(NA,nrow=159,ncol=K)

# rps and logs
rps_zsmsnb<- matrix(NA, nrow=159, ncol=K)
logs_zsmsnb<- matrix(NA, nrow=159, ncol=K)

# X matrix (updating)
X.simulated<- array(NA,dim=c(159,K+1,M)) 


# spatial constants 
adj = as.numeric(as.carAdjacency(N_matrix)$adj)
num = as.numeric(as.carAdjacency(N_matrix)$num)
count= count



# simulated X[t-1]
for(i in 1:159){
  X.simulated[i,1,]<- as.numeric(unlist(samps.zsmsnb[paste0("X_final.",i,".")]))
}

# iterations 
for (k in 1:K){
  print(k)
  for (i in 1:159){
    #neighborhood 
    nei <- which(N_matrix[i,]==1)
    #neighborhood sum
    nei_sum <- rep(0,M)
    for(j in nei){
      nei_sum <- nei_sum+X.simulated[j,k,]
    }
    
    # for t=T0+1,T0+2,...T0+K 
    beta0<- as.numeric(unlist(samps.zsmsnb[paste0("beta0")]))
    beta1<- as.numeric(unlist(samps.zsmsnb[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zsmsnb[paste0("beta2")]))
    
    alpha1<- as.numeric(unlist(samps.zsmsnb[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.zsmsnb[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.zsmsnb[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.zsmsnb[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zsmsnb[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zsmsnb[paste0("alpha6")]))
    
    gamma1<- as.numeric(unlist(samps.zsmsnb[paste0("gamma1")]))
    gamma2<- as.numeric(unlist(samps.zsmsnb[paste0("gamma2")]))
    
    r<- as.numeric(unlist(samps.zsmsnb[paste0("r")]))
    
    # mup
    mup.pred.zsmsnb<- exp(beta0+
                            beta1*(HDI[i]-mean(HDI[1:159]))+
                            beta2*(temp[T0+k]-mean(temp[1:Tsim])))
    # p01
    p01<- expit(alpha1+ alpha2*(HDI[i]-mean(HDI[1:159]))+
                  alpha3*(temp[T0+k]-mean(temp[1:Tsim]))+gamma1*nei_sum )
    
    # p11: change the parameter: 
    p11<- expit(alpha4+ +alpha5*(HDI[i]-mean(HDI[1:159]))+
                  alpha6*(temp[T0+k]-mean(temp[1:Tsim]))+gamma2*nei_sum )
    
    # r
    r_nb<- r
    
    # simulated X, calculate and save 
    Xk<- rbinom(n= rep(1,M), 
                size=1,
                p= X.simulated[i,k,]*p11 + (1-X.simulated[i,k,])*p01)
    X.simulated[i,k+1,]<- Xk
    
    p_nbbb <- (r_nb/(r_nb+(X.simulated[i,k+1,]*mup.pred.zsmsnb)))-.00000000001*(1-X.simulated[i,k+1,])
    
    simulated.pred[i,k,]<- rnbinom(n = M,size = r_nb,prob = p_nbbb)
    
    logs_zsmsnb[i,k] <- -log(mean(dnbinom(x=cases[i,T0+k] ,size = r_nb,prob = p_nbbb)))
    
    P_pred <- ecdf(simulated.pred[i,k,])(seq(0,10000,by=1))
    rps_zsmsnb[i,k]<- sum((P_pred -ifelse(cases[i,T0+k]<=seq(0,10000,by=1),1,0))^2)
    
  }
}

# rps 
K1<-cbind(mean(rps_zsmsnb[,1]))
K2<-cbind(mean(rps_zsmsnb[,2]))
K3<-cbind(mean(rps_zsmsnb[,3]))
K4<-cbind(mean(rps_zsmsnb[,4]))

L1<-cbind(mean(logs_zsmsnb[,1]))
L2<-cbind(mean(logs_zsmsnb[,2]))
L3<-cbind(mean(logs_zsmsnb[,3]))
L4<-cbind(mean(logs_zsmsnb[,4]))

dat.rps.zi<-data.frame(rbind(K1,K2,K3,K4))
dat.logs.zi<-data.frame(rbind(L1,L2,L3,L4))


row.names(dat.rps.zi)<- c("K=1","K=2","K=3","K=4")
row.names(dat.logs.zi)<- c("D=1","D=2","D=3","D=4")


write.csv(dat.rps.zi,paste0("/lustre03/project/6003552/mingchi/rps_result/ZSMSNB_rps_sim_large_p",as.character(T0),".csv"))
write.csv(dat.logs.zi,paste0("/lustre03/project/6003552/mingchi/dss_result/ZSMSNB_logs_sim_large_p",as.character(T0),".csv"))


