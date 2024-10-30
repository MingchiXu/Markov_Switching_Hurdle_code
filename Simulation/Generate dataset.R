######################################
# A summary for the final revisions
######################################


#==========================================================
# Simulation for Revise Version 1
#==========================================================
rm(list = ls())
# load packages
library(nimble)
library(microbenchmark)
library(ggplot2)
library(gridExtra)
library(compareMCMCs)
library(boot)
library(latex2exp)
library(ggpubr)
library(spdep)
library(leaflet)
library(MASS)
library(sf)
library(terra)

# import the map matrix 


shape_rj=maptools::readShapePoly("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/map/shape-RJ2-2.shp") 
shape_rj$oo <- seq(1:160)

nb2mat(poly2nb(shape_rj), style = "B",zero.policy = TRUE)# B is the basic binary coding
# this is the default order:
shape_rj$NOME


#==========================================================================
# Step 1: Simulate true counts from a hurdle model  -----------------------
#==========================================================================



# parameter setting 
set.seed(1998)
Tsim<- 84
beta0<- 0.5
alpha0_0<- -3
alpha0_1<- 1.5
r<- 1.5 # over-dispersion parameter 

beta1<- 0.1
beta2<- 0.4

alpha_HDI_p01<- 1.15 # 1
alpha_temp_01<- 1.1 # 1.1

alpha_HDI_p11<- 1.18 # 1
alpha_temp_11<- 1.2 # 1.2

# 0.8/ 0.5 // 2/1.4 // 0.6/0.3
gamma1<- 0.6
gamma2<- 0.3


#Read in the covariates: HDI, pop  
area.level<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/DistrictCovariates.csv")
load("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/final_data_cleaned.Rdata")

# a spatial covariate: HDI (159 districts)
HDI<- area.level$HDI
# a temporal covariate: Temperature (84 weeks)
temp <- d$MeanMaxTemp[1:Tsim]

# read in the matrix 
N_matrix<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/spatial_matrix.csv")
NM<-N_matrix

# neighbor's number
num <- as.numeric(as.carAdjacency(N_matrix)$num)
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
count <- rep(NA,159)
count[1] <- 1 
for(i in 1:158){
  count[i+1]  <- count[i]+num[i]
}

# create matrices to save the simulated results 

case_sim<- matrix(NA, nrow=159, ncol=Tsim)
X.sim<- matrix(NA, nrow=159, ncol=Tsim)
mu_cases<- matrix(NA, nrow=159, ncol=Tsim)
NeiXsum<- matrix(rep(0, 159*Tsim), nrow = 159)

p01<- matrix(NA, nrow=159, ncol=Tsim)
p11<- matrix(NA, nrow=159, ncol=Tsim)

# Simulate X_it, then we can simulate y_it 
# generate the initial stats 
X.sim[,1]<- rbinom(159, 1, 0.5)

# the initial status for p01 and p11
for (i in 1:159){
  p01[i,1]<- inv.logit(alpha0_0+
                         alpha_HDI_p01*(HDI[i]-mean(HDI))+
                         alpha_temp_01*(temp[1]-mean(temp))+
                         gamma1* NeiXsum[i,1])
  
  p11[i,1]<- inv.logit(alpha0_1+
                         alpha_HDI_p11*(HDI[i]-mean(HDI))+
                         alpha_temp_11*(temp[1]-mean(temp))+
                         gamma2* NeiXsum[i,1])
}

# generate the Markov Chain
for (t in 2: Tsim){
  for (i in 1:159){  
    for (j in which(NM[i,]==1)){
      NeiXsum[i,t]<-  NeiXsum[i,t]+X.sim[j,t-1]
    }
    # model the transition probs: p01 and p11
    p01[i,t]<- inv.logit(alpha0_0+
                           alpha_HDI_p01*(HDI[i]-mean(HDI))+
                           alpha_temp_01*(temp[t]-mean(temp))+
                           gamma1* NeiXsum[i,t])
    p11[i,t]<- inv.logit(alpha0_1+
                           alpha_HDI_p11*(HDI[i]-mean(HDI))+
                           alpha_temp_11*(temp[t]-mean(temp))+ 
                           gamma2* NeiXsum[i,t])
    
    # states 
    X.sim[i,t]<- rbinom(n=1, size=1, 
                        prob=p01[i,t]*(1-X.sim[i,t-1])+p11[i,t]*X.sim[i,t-1])
  }
}



#-------- Generate the simulated (true) cases-------------------------------

# change to Z~ ZTNB(p, r) & 0
for (i in 1:159){
  for (t in 1: Tsim){
    # generate the mean and the true cases 
    mu_cases[i,t]<- exp(beta0+
                          beta1*(HDI[i]-mean(HDI))+
                          beta2*(temp[t]-mean(temp)))
    p_nb<- r/(r+ mu_cases[i,t])
    
    # simulate Z~ hurdle model 
    if (X.sim[i,t]==0){
      case_sim[i,t]<- 0
    }else{
      xt<- rnbinom(1, prob=p_nb, size=r)
      while(xt==0){
        xt<- rnbinom(1, prob=p_nb, size=r)
      }
      case_sim[i,t]<- xt
    }
  }
}


#==========================================================================
# Step 2: Simulate the reported counts by different report rates ----------
#==========================================================================


# create a matrix to save the reported cases 
report_cases<- matrix(NA, nrow=159, ncol=Tsim)
# under a high reporting rate: p0=1/ p0=0.8/ p0=0.6/ p0=0.1: 
  p0<- 1     # p0^(1)
# p0<- 0.8   # p0^(2)
# p0<- 0.6   # p0^(3)
# p0<- 0.1   # p0^(4)  


# the reported cases
for (i in 1:159){
  for (t in 1: Tsim){
    report_cases[i,t]<- rbinom(1,size=case_sim[i,t], prob=p0)
  }
}

# save the simulated reported rates 
write.csv(report_cases, "C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_full.csv") # 100% reporting 
# write.csv(report_cases, "C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_80.csv") # 80% reporting 
# write.csv(report_cases, "C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_60.csv") # 60% reporting 
# write.csv(report_cases, "C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_10.csv") # 10% reporting 

report_cases<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_full.csv")
#report_cases<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_80.csv")
#report_cases<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_60.csv")
#report_cases<-read.csv("C:/Users/xumc7/OneDrive/Desktop/SiM_v1/cases_10.csv")


# remove the first column (heading)
cases<-report_cases[,-1]


# visualize the cases
little.f<- data.frame(seq(1:Tsim), as.vector(t(cases[103,])), X.sim[103,], (case_sim[103,]))
colnames(little.f)<- c("week","case", "state", "true_case")
library(ggplot2)
library(cowplot)

# reported cases:
p00<- ggplot(data=little.f)+
  geom_point(aes(x=week, y=case))
p11<- ggplot(data=little.f)+
  geom_point(aes(x=week, y=state))

plot_grid(p00, p11, align="v", nrow=2)

# actual cases 
p22<- ggplot(data=little.f)+
  geom_point(aes(x=week, y=true_case))
p33<- ggplot(data=little.f)+
  geom_point(aes(x=week, y=state))

plot_grid(p22, p33, align="v", nrow=2)

# check the mean disease present rate: about 50% 
mean(X.sim)
# check the mean cases
mean(as.matrix(cases))
# check the max cases
max(cases) 




#======================================================================
# Step 3: Fit ZS-MSNB and ZS-MSNBH model, and compare the rps ---------
#======================================================================
# refer to the other files in this folder.
