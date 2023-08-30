# Codes for all Figures in the paper

library(stringr)
library(ggplot2)
library("cowplot")
library(boot)
library('nimble')
library(grDevices)
library("scales")

# Read and clean the data
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
dataCHIKV<-datCHIKV[-c(14),]
cases<-matrix(seq(1,159*384,by=1),nrow=159)
cases<-as.data.frame(cases)
for (i in 1:159){
  cases[i,]<-as.numeric(dataCHIKV[i,][-1])
}
cases<-as.matrix(cases)




###################################################################################
####  Fig 1. Visualization of reported counts in two areas (No.1 and No.150) ######       
###################################################################################

# timeline
ts.dates.l<-list()
for (i in 1:159){
  ts.dates.l[[i]]<-ts(as.numeric(cases[i,]),
                      start=c(2015,1),end=c(2022,20),frequency=52)
}
districtname<-dataCHIKV$District 

zeros.in.area.150<-which(ts.dates.l[[150]]==0)

# plot for district 1 in introduction
# create a factor to distinguish 0 and non-zero cases
b0<- data.frame(Time=c(time(ts.dates.l[[1]])),
                district=c(ts.dates.l[[1]]))
b0$zero.index<- ifelse(b0$district==0,1,0)
b0$zero.index<- factor(b0$zero.index)

p_area1<- ggplot(b0,aes(x=Time,y=district,color=zero.index))+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  theme_classic()+
  labs(title="(a) Saúde district, population: 2,749", x="Year", y="Cases")+
  guides(color="none")+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))

# plot for district 150 
a0<- data.frame(Time=c(time(ts.dates.l[[150]])),
                district=c(ts.dates.l[[150]]))
a0$zero.index<- ifelse(a0$district==0,1,0)
a0$zero.index<- factor(a0$zero.index)
p_area150<- ggplot(a0,aes(x=Time,y=district,color=zero.index))+
  geom_point()+
  scale_color_manual(values=c("black","red"))+
  theme_classic()+
  xlab('Year')+
  ylab('Cases')+
  guides(color="none")+
  labs(title="(b) Campo Grande district, population: 328,370", x="Year", y="Cases")+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))

figure.1<-plot_grid(p_area1, p_area150, 
                    ncol = 1, nrow = 2,  align="v")

figure.1






####################################################################################
#### Fig 2. NB approximation: derive a NB to approximate the count distribution  ###
####################################################################################


p<- 0.1
r<- 2
lambda=8


# def exact distribution for Z
d.con.reported.cases<- function (z,r,p,lambda){
  if (z>0){
    return ( -( r/(lambda+r) )^r*p^z*factorial(r+z-1)*(lambda/(lambda+r))^z*( (lambda*p+r)/(lambda+r) )^(-r-z)/(
      factorial(r-1)*factorial(z)*(  (r/(lambda+r))^r-1 ))  
    )
  }
  if (z==0){
    return (  (r/lambda)^r*(1- ((lambda*p+r)/(lambda+r))^(-r))/((r/lambda)^r-( (lambda+r)/lambda )^r) )
  }
  else{
    print("density should be larger or equal to 0")
  }
}
d.con.reported.cases.V<-Vectorize(d.con.reported.cases)
seq.z<- seq(0,15,by=1)


#  plots
den.Z<-d.con.reported.cases.V(seq.z, r, p, lambda)

# r' and lambda' for a NB:
lambda.new<- p*lambda/(1- (r/(r+lambda) )^r) 
r.new<- 1/((1-(r/(r+lambda) )^r )*(1+1/r)-1)


den.NB<- dnbinom(seq.z, size= r.new, prob= (r.new/(r.new+lambda.new) ) )
df<- data.frame(seq.z,den.Z,den.NB)

# plot 3 cases: 

# Case 1: small p
library(ggplot2)

counts<- rep(seq(0, 15, by=1),2)
num<- c( df$den.Z, df$den.NB)
Distribution<- c(rep("True distribution of\nreported counts",16), rep("NB approximation",16))
df.small.p<- data.frame(counts,num, Distribution)

p1<- ggplot(df.small.p,aes(x=counts, y=num, group=Distribution))+
  geom_line(aes(linetype=Distribution))+
  geom_point(aes(), shape=1)+
  theme_classic()+
  labs(title="(a) 10% Reporting Rate", 
       x="Counts", 
       y="Probability")+
  theme(legend.position=c(0.656,0.94))+
  theme(legend.text = element_text(size=20))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  theme(legend.title = element_blank())+
  theme(legend.key.height=unit(1, "cm"))+
  theme(plot.title=element_text(face='bold'))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=20))+
  theme(legend.background = element_rect(fill='transparent', colour=NA))

# Case 2: moderate p
p.2<- 0.3
r.2<- 2
lambda.2=8

den.Z.2<-d.con.reported.cases.V(seq.z, r.2, p.2, lambda.2)

# r' and lambda' for a NB:
lambda.new.2<- p.2*lambda.2/(1- (r.2/(r.2+lambda.2) )^r.2)
r.new.2<- 1/((1-(r.2/(r.2+lambda.2) )^r.2 )*(1+1/r.2)-1)
den.NB.2<- dnbinom(seq.z, size= r.new.2, prob= (r.new.2/(r.new.2+lambda.new.2) ) )

counts<- rep(seq(0, 15, by=1),2)
num<- c( den.Z.2, den.NB.2)
Variable<- c(rep("True distribution of\nreported counts",16), rep("NB approximation",16))
df.mod.p<- data.frame(counts,num, Variable)

p2<- ggplot(df.mod.p,aes(x=counts, y=num, group=Variable))+
  geom_line(aes(linetype=Variable))+
  geom_point(aes(), shape=1)+
  theme_classic()+
  labs(title="(b) 30% Reporting Rate", 
       x="Counts", 
       y="Probability")+
  theme(legend.position=c(0.656,0.94))+
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_blank())+
  theme(legend.key.height=unit(1, "cm"))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  theme(plot.title=element_text(face='bold'))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=20))+
  theme(legend.background = element_rect(fill='transparent', colour=NA))

# Case 3: large p
p.3<- 0.8
r.3<- 2
lambda.3=8
den.Z.3<-d.con.reported.cases.V(seq.z, r.3, p.3, lambda.3)

# r' and lambda' for a NB:
lambda.new.3<- p.3*lambda.3/(1- (r.3/(r.3+lambda.3) )^r.3)
r.new.3<- 1/((1-(r.3/(r.3+lambda.3) )^r.3 )*(1+1/r.3)-1)
den.NB.3<- dnbinom(seq.z, size= r.new.3, prob= (r.new.3/(r.new.3+lambda.new.3) ) )

counts<- rep(seq(0, 15, by=1),2)
num<- c( den.Z.3, den.NB.3)
Variable<- c(rep("True distribution of\nreported counts",16), rep("NB approximation",16))
df.large.p<- data.frame(counts,num, Variable)

p3<- ggplot(df.large.p,aes(x=counts, y=num, group=Variable))+
  geom_line(aes(linetype=Variable))+
  geom_point(aes(), shape=1)+
  theme_classic()+
  labs(title="(c) 80% Reporting Rate", 
       x="Counts", 
       y="Probability")+
  theme(legend.position=c(0.656,0.94))+
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_blank())+
  theme(legend.key.height=unit(1, "cm"))+
  theme(plot.title=element_text(face='bold'))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size=20))+
  theme(panel.background = element_rect(fill = 'transparent', colour=NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = 'transparent', colour=NA))+
  theme(legend.background = element_rect(fill='transparent', colour=NA))


figure.2<-plot_grid(p1, p2, p3, 
                    ncol = 3, nrow = 1, align = "h")

figure.2



#####################################################
#            Fig  3. The endemic rates.             #
#####################################################


# for ZS-MSNB model: 

load("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNB_spatial.RData")# load the samples


samps.ZSMSNB<- data.frame(rbind(ZSMSNB.samples.dynamic[[1]],
                                ZSMSNB.samples.dynamic[[2]],
                                ZSMSNB.samples.dynamic[[3]]))



time.seq<-seq(1,105,by=1)
sin.time<-sin(time.seq*2*pi/52)
cos.time<-cos(time.seq*2*pi/52)


posterior.mean<- rep(NA, 105)
posterior.upper<- rep(NA, 105)
posterior.lower<- rep(NA, 105)

for (t in 1:105){
  print(t)
  # posterior mean
  posterior.mean[t]<-  mean(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",1,".")]))+
                                   as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                   as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ))
  # posterior upper 95% CI.
  posterior.upper[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",1,".")]))+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ),0.975)
  # posterior lower 95% CI. 
  posterior.lower[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",1,".")]))+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ),0.025)
}

ts.dates<-ts(time.seq, start=c(2015,1),end=c(2017,1),frequency=52)
df.zsmsnb<- data.frame(ts.dates, 
                       posterior.mean, 
                       posterior.lower, 
                       posterior.upper)

ts.posterior.mean<-ts(posterior.mean, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.upper<-ts(posterior.upper, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.lower<-ts(posterior.lower, start=c(2015,1),end=c(2017,1),frequency=52)



df.zsmsnb<- data.frame(Time=c(time(ts.posterior.mean)), 
                       post.mean=c(ts.posterior.mean), 
                       post.lower=c(ts.posterior.lower), 
                       post.upper=c(ts.posterior.upper))




ggplot(aes(x=Time), data=df.zsmsnb)+
  geom_line(aes(y=post.mean) )+
  geom_line(aes(y=post.upper), linetype='dashed')+
  geom_line(aes(y=post.lower), linetype='dashed')+
  theme_classic()+
  labs(title="(a) Saúde district", x="Time (weeks)", y="Endemic Rate")+
  theme(plot.title = element_text(hjust = 0.5))


# for ZS-MSNBH model: 

load("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_spatial.RData")# load the samples

samps.ZSMSNBH<- data.frame(rbind(ZSMSNBH.samples[[1]],
                                 ZSMSNBH.samples[[2]],
                                 ZSMSNBH.samples[[3]]))

time.seq<-seq(1,105,by=1)
sin.time<-sin(time.seq*2*pi/52)
cos.time<-cos(time.seq*2*pi/52)


posterior.mean.h<- rep(NA, 105)
posterior.upper.h<- rep(NA, 105)
posterior.lower.h<- rep(NA, 105)

for (t in 1:105){
  print(t)
  # posterior mean
  posterior.mean.h[t]<-  mean(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",1,".")]))+
                                     as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                     as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ))
  # posterior upper 95% CI.
  posterior.upper.h[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",1,".")]))+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ),0.975)
  # posterior lower 95% CI. 
  posterior.lower.h[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",1,".")]))+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ),0.025)
}

ts.dates<-ts(time.seq, start=c(2015,1),end=c(2017,1),frequency=52)

ts.posterior.mean.h<-ts(posterior.mean.h, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.upper.h<-ts(posterior.upper.h, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.lower.h<-ts(posterior.lower.h, start=c(2015,1),end=c(2017,1),frequency=52)


df.ms.two.mods<- data.frame(Time=c(time(ts.posterior.mean)),
                            zsmsnbh.mean= c(ts.posterior.mean),
                            zsmsnbh.mean= c(ts.posterior.mean.h))

colnames(df.ms.two.mods)<- c("time","ZS-MSNB","ZS-MSNBH")
df.ms.1<- data.frame(Time=c(time(ts.posterior.mean)), 
                     post.mean=c(ts.posterior.mean), 
                     post.lower=c(ts.posterior.lower), 
                     post.upper=c(ts.posterior.upper), 
                     post.mean.h=c(ts.posterior.mean.h), 
                     post.lower.h=c(ts.posterior.lower.h), 
                     post.upper.h=c(ts.posterior.upper.h))


cols=c("ZS-MSNB"='black',"ZS-MSNBH"='red')

summer.start<- c(time(ts.posterior.mean))[c(1,  49,  101)]
summer.end<-  c(time(ts.posterior.mean))[c(13,  65,  104)]
summer.len<- summer.end- summer.start

winter.start<- c(time(ts.posterior.mean))[c( 23,  75)]
winter.end<-  c(time(ts.posterior.mean))[c( 39,  91)]
winter.len<- winter.end- winter.start


bar.len.summer<- rep(max(df.ms.1$post.upper), length(summer.len))
summer.df<- data.frame(summer.start, summer.len, bar.len.summer)

bar.len.winter<- rep(max(df.ms.1$post.upper), length(winter.len))
winter.df<- data.frame(winter.start, winter.len, bar.len.winter)

x.lab.time<- c("2015", "",
               "2016", "",
               "2017")

p.endemic.1<- ggplot(aes(x=Time), data=df.ms.two.mods)+
  geom_bar(data=summer.df, aes(x=summer.start+summer.len, y=bar.len.summer), 
           stat='identity', width=summer.len, fill='red', alpha=.1)+
  geom_bar(data=winter.df, aes(x=winter.start+winter.len, y=bar.len.winter), 
           stat='identity', width=winter.len, fill='darkblue', alpha=.1)+
  geom_line(aes(y=post.upper), linetype='dotted', data=df.ms.1)+
  geom_line(aes(y=post.lower), linetype='dotted', data=df.ms.1)+
  geom_line(aes(y=post.upper.h), linetype='dotted', data=df.ms.1, color='red')+
  geom_line(aes(y=post.lower.h), linetype='dotted', data=df.ms.1, color='red')+
  geom_line(data=df.ms.1,aes(y=post.mean, color='ZS-MSNB'), size=0.7)+
  geom_line(data=df.ms.1, aes(y=post.mean.h,  color="ZS-MSNBH"), size=0.7)+
  theme_classic()+
  labs(title="(a) Saúde district", x="Year", y="Endemic Rate")+
  scale_color_manual(name='Models',
                     values=cols)+
  theme(legend.title.align=0.5)+
  theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position=c(0.82,0.96), legend.justification=c(0.01,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(labels=x.lab.time)+
  theme(legend.text = element_text(size=20))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))


# For district=150

posterior.mean.150<- rep(NA, 105)
posterior.upper.150<- rep(NA, 105)
posterior.lower.150<- rep(NA, 105)

for (t in 1:105){
  print(t)
  # posterior mean
  posterior.mean.150[t]<-  mean(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",150,".")]))+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                       as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ))
  # posterior upper 95% CI.
  posterior.upper.150[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",150,".")]))+
                                           as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                           as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ),0.975)
  # posterior lower 95% CI. 
  posterior.lower.150[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNB[paste0("b0.",150,".")]))+
                                           as.numeric(unlist(samps.ZSMSNB[paste0("beta1")]))*sin.time[t]+
                                           as.numeric(unlist(samps.ZSMSNB[paste0("beta2")]))*cos.time[t] ),0.025)
}



ts.posterior.mean.150<-ts(posterior.mean.150, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.upper.150<-ts(posterior.upper.150, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.lower.150<-ts(posterior.lower.150, start=c(2015,1),end=c(2017,1),frequency=52)

posterior.mean.h.150<- rep(NA, 105)
posterior.upper.h.150<- rep(NA, 105)
posterior.lower.h.150<- rep(NA, 105)

for (t in 1:105){
  print(t)
  # posterior mean
  posterior.mean.h.150[t]<-  mean(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",150,".")]))+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                         as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ))
  # posterior upper 95% CI.
  posterior.upper.h.150[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",150,".")]))+
                                             as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                             as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ),0.975)
  # posterior lower 95% CI. 
  posterior.lower.h.150[t]<- quantile(exp( as.numeric(unlist(samps.ZSMSNBH[paste0("b0.",150,".")]))+
                                             as.numeric(unlist(samps.ZSMSNBH[paste0("beta1")]))*sin.time[t]+
                                             as.numeric(unlist(samps.ZSMSNBH[paste0("beta2")]))*cos.time[t] ),0.025)
}

ts.posterior.mean.h.150<-ts(posterior.mean.h.150, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.upper.h.150<-ts(posterior.upper.h.150, start=c(2015,1),end=c(2017,1),frequency=52)
ts.posterior.lower.h.150<-ts(posterior.lower.h.150, start=c(2015,1),end=c(2017,1),frequency=52)




df.ms.two.mods.150<- data.frame(Time=c(time(ts.posterior.mean.150)),
                                zsmsnbh.mean= c(ts.posterior.mean.150),
                                zsmsnbh.mean= c(ts.posterior.mean.h.150))

colnames(df.ms.two.mods.150)<- c("time","ZS-MSNB","ZS-MSNBH")
df.ms.150<- data.frame(Time=c(time(ts.posterior.mean.150)), 
                       post.mean=c(ts.posterior.mean.150), 
                       post.lower=c(ts.posterior.lower.150), 
                       post.upper=c(ts.posterior.upper.150), 
                       post.mean.h=c(ts.posterior.mean.h.150), 
                       post.lower.h=c(ts.posterior.lower.h.150), 
                       post.upper.h=c(ts.posterior.upper.h.150))



cols=c("ZS-MSNB"='black',"ZS-MSNBH"='red')

summer.start<- c(time(ts.posterior.mean))[c(1,  49,  101)]
summer.end<-  c(time(ts.posterior.mean))[c(13,  65,  104)]
summer.len<- summer.end- summer.start

winter.start<- c(time(ts.posterior.mean))[c( 23,  75)]
winter.end<-  c(time(ts.posterior.mean))[c( 39,  91)]
winter.len<- winter.end- winter.start


bar.len.summer<- rep(max(df.ms.150$post.upper), length(summer.len))
summer.df<- data.frame(summer.start, summer.len, bar.len.summer)

bar.len.winter<- rep(max(df.ms.150$post.upper), length(winter.len))
winter.df<- data.frame(winter.start, winter.len, bar.len.winter)

x.lab.time<- c("2015", "",
               "2016", "",
               "2017")




cols=c("ZS-MSNB"='black',"ZS-MSNBH"='red')
p.endemic.150<-ggplot(aes(x=Time), data=df.ms.two.mods.150)+
  geom_bar(data=summer.df, aes(x=summer.start+summer.len, y=bar.len.summer), 
           stat='identity', width=summer.len, fill='red', alpha=.1)+
  geom_bar(data=winter.df, aes(x=winter.start+winter.len, y=bar.len.winter), 
           stat='identity', width=winter.len, fill='darkblue', alpha=.1)+
  
  geom_line(aes(y=post.upper), linetype='dotted', data=df.ms.150)+
  geom_line(aes(y=post.lower), linetype='dotted', data=df.ms.150)+
  geom_line(aes(y=post.upper.h), linetype='dotted', data=df.ms.150, color='red')+
  geom_line(aes(y=post.lower.h), linetype='dotted', data=df.ms.150, color='red')+
  geom_line(data=df.ms.150,aes(y=post.mean, color='ZS-MSNB'), size=0.7)+
  geom_line(data=df.ms.150, aes(y=post.mean.h,  color="ZS-MSNBH"), size=0.7)+
  theme_classic()+
  labs(title="(b) Campo Grande district", x="Year", y="Endemic Rate")+
  scale_color_manual(name='Models',
                     values=cols)+
  theme(legend.title.align=0.5)+
  theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position=c(0.82,0.96), legend.justification=c(0.01,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(labels=x.lab.time)+
  theme(legend.text = element_text(size=20))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))

# combine the two plots

figure.3<-plot_grid(p.endemic.1, p.endemic.150, 
                    ncol = 1, nrow = 2 , align="v")
figure.3




################################################################
#### Fig 7. Forecasting for two Markov switching models  #######
################################################################


# prediction up to 231, mainly focus on two Markov Switching models 

# ------------------ (A) Fit the model up to time point T0=231 
# (See the other file called prediction_231.R )


load("/Users/xumingchi/Desktop/hurdle-codes/NB_231.RData") # NB
load("/Users/xumingchi/Desktop/hurdle-codes/NBH_231.RData") # NBH
load("/Users/xumingchi/Desktop/hurdle-codes/ZINB_231.RData") # ZINB
load("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNB_231.RData") # ZSMSNB
load("/Users/xumingchi/Desktop/hurdle-codes/ZSMSNBH_231.RData") # ZSMSNBH


# five models: 
samps.nb<-data.frame(rbind(NB.samples[[1]],
                           NB.samples[[2]],
                           NB.samples[[3]]))# NB


samps.nbh<-data.frame(rbind(NBH.samples[[1]],
                            NBH.samples[[2]],
                            NBH.samples[[3]]))# NBH


samps.zinb<-data.frame(rbind(ZINB.samples[[1]],
                             ZINB.samples[[2]],
                             ZINB.samples[[3]]))# ZINB



samps.zsmsnb<-data.frame(rbind(ZSMSNB.samples.dynamic[[1]],
                               ZSMSNB.samples.dynamic[[2]],
                               ZSMSNB.samples.dynamic[[3]]))#ZSMSNB


samps.zsmsnbh<-data.frame(rbind(ZSMSNBH.samples[[1]],
                                ZSMSNBH.samples[[2]],
                                ZSMSNBH.samples[[3]]))# ZSMSNBH



# ------------------- (B) Fig 7: Plot the prediction values for the two Markov switching models 

K0<- 8

mcos.time <- mean(cos.time[1:T0])
msin.time <- mean(sin.time[1:T0])



#(B.1) ZS-MSNB: 

M <- nrow(samps.zsmsnb)
# for M samples
simulated.pred.k0<- array(0,dim=c(159,K0,M))
# pred cases
simulated.pred.zsmsnb.k0 <- matrix(0,nrow=159,ncol=K0)
# 95% upper and lower (within T0, for simulated prediction)? 
upper.simulated.pred.zsmsnb.k0 <-matrix(0,nrow=159,ncol=K0)
lower.simulated.pred.zsmsnb.k0 <-matrix(0,nrow=159,ncol=K0)


# X matrix 
X.simulated.k0<- array(0,dim=c(159,K0+1,M)) 
X.simulated.mean.k0<- array(0,dim=c(159,K0)) 

# spatial constants
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
num=as.numeric(as.carAdjacency(N_matrix)$num)
count=count


for (i in 1:159){
  print(i)
  X.simulated.k0[i,1,]<- as.numeric(unlist(samps.zsmsnb[paste0("X_final.",i,".")]))
  nei <- which(N_matrix[i,]==1)
  
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K0){
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.zsmsnb[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.zsmsnb[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.zsmsnb[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zsmsnb[paste0("beta2")]))
    alpha1<- as.numeric(unlist(samps.zsmsnb[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.zsmsnb[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.zsmsnb[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.zsmsnb[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zsmsnb[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zsmsnb[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.zsmsnb[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.zsmsnb[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.zsmsnb[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.zsmsnb[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.zsmsnb[paste0("alpha11")]))
    
    
    gamma1<- as.numeric(unlist(samps.zsmsnb[paste0("gamma1")]))
    gamma2<- as.numeric(unlist(samps.zsmsnb[paste0("gamma2")]))
    gamma3<- as.numeric(unlist(samps.zsmsnb[paste0("gamma3")]))
    
    delta1<- as.numeric(unlist(samps.zsmsnb[paste0("delta1")]))
    delta2<- as.numeric(unlist(samps.zsmsnb[paste0("delta2")]))
    
    #neighborhood sum
    nei_sum <- rep(0,M)
    if(k==1){
      for(j in nei){
        nei_sum <- nei_sum+samps.zsmsnb[,paste0("X_final.",j,".")]
      }
    } else{
      for(j in nei){
        nei_sum <- nei_sum+X.simulated.k0[j,k,]
      }
    }
    
    # mup
    if (k== 1){
      mup.pred.zsmsnb<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred- msin.time )+beta2*(cos.time.pred- mcos.time ))
    }else {
      mup.pred.zsmsnb <-  exp(b0)*simulated.pred.k0[i,k-1,]+exp(b+beta1*(sin.time.pred- msin.time )+beta2*(cos.time.pred- mcos.time ))
    }
    # p01 and p11 should be changed into these, and r: 
    
    # p01
    p01<- expit(alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
                  alpha3*(pop[i]-mean(pop[1:159]))+gamma1*(greenarea[i]-mean(greenarea[1:159]))+
                  delta1*nei_sum)
    
    # p11: change the parameter: 
    if(k==1){
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(cases[i,T0]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*nei_sum)
    }else{
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(simulated.pred.k0[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*nei_sum)
    }
    # r: change the parameters: 
    
    if (k==1){
      r_nb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(cases[i,T0]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      r_nb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                   alpha10*(pop[i]-mean(pop[1:159]))+
                   alpha11*(log(simulated.pred.k0[i,k-1,]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }
    #X
    Xk<- rbinom(n= rep(1,M), size=1,
                p= X.simulated.k0[i,k,]*p11 + (1-X.simulated.k0[i,k,])*p01)
    X.simulated.k0[i,k+1,]<- Xk
    # prediction matrix
    
    p_nb <- r_nb/(r_nb+(Xk)*mup.pred.zsmsnb) - 1e-10*(1-Xk)
    
    
    simulated.pred.k0[i,k,]<-rnbinom(M, prob=p_nb, size=r_nb) 
    
    
    # the mean and 95% CI
    #simulated.pred.zip[i,k]<-mean(simulated.pred[i,k,])
    simulated.pred.zsmsnb.k0[i,k]<- mean( simulated.pred.k0[i,k,])
    upper.simulated.pred.zsmsnb.k0[i,k]<- quantile(simulated.pred.k0[i,k,],probs = c(.975))
    lower.simulated.pred.zsmsnb.k0[i,k]<- quantile(simulated.pred.k0[i,k,],probs = c(.025))
    # state indicator: 
    X.simulated.mean.k0[i,k]<-mean(X.simulated.k0[i,k+1,])
  }
}



#---(B.2) ZSMSNBH: 



# define X
X<- matrix(rep(0,159*384),nrow=159,ncol=384)
X[cases>0] <- 1


NM<-N_matrix

NeiXsum<- matrix(rep(0, 159*384), nrow = 159)
for (i in 1:159){
  for (t in 2:384){
    for (j in which(NM[i,]==1)){
      NeiXsum[i,t]<-  NeiXsum[i,t]+X[j,t-1]
    }
  }
}




M <- nrow(samps.zsmsnbh)
# for M samples
simulated.pred.k0<- array(0,dim=c(159,K0,M))
# pred cases
simulated.pred.zsmsnbh.k0<- matrix(0,nrow=159,ncol=K0)
upper.simulated.pred.zsmsnbh.k0<- matrix(0,nrow=159,ncol=K0)
lower.simulated.pred.zsmsnbh.k0<- matrix(0,nrow=159,ncol=K0)


# spatial constants
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
num=as.numeric(as.carAdjacency(N_matrix)$num)
count=count
# X
X.h<- array(0,dim=c(159,K0,M))
X.simulated.mean.k0.h<- array(0,dim=c(159,K0)) 

for (i in 1:159){
  print(i)
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K0){
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
      mup.pred.zsmsnbh<-  exp(b0)*simulated.pred.k0[i,k-1,]+exp(b+beta1*(sin.time.pred-msin.time )+beta2*(cos.time.pred-mcos.time ))
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
                      alpha7*(log(simulated.pred.k0[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
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
                   alpha11*(log(simulated.pred.k0[i,k-1,]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    # X: 
    if (k==1){
      X.h[i,1,]<- ifelse(rep(cases[i,T0],M)>0,1,0)# X of "T0"
    }else{
      X.h[i,k,]<-ifelse(simulated.pred.k0[i,k-1,]>0,1,0)# X of "T0+1",...,"TO+K-1"
    }
    
    p_nb<- p01.h*(1-X.h[i,k,])+p11.h*X.h[i,k,]
    
    # simulated y
    simulated.pred.k0[i,k,]<-rhurdleNegBinMargV(n=rep(1,M), 
                                                p =p_nb, 
                                                lambda= mup.pred.zsmsnbh,
                                                r=r_nb) 
    
    # mean and 95% CI 
    simulated.pred.zsmsnbh.k0[i,k]<-mean(simulated.pred.k0[i,k,])
    upper.simulated.pred.zsmsnbh.k0[i,k]<- quantile(simulated.pred.k0[i,k,],probs = c(.975))
    lower.simulated.pred.zsmsnbh.k0[i,k]<- quantile(simulated.pred.k0[i,k,],probs = c(.025))
    # state indicator: 
    X.simulated.mean.k0.h[i,k]<-mean(X.h[i,k,])
  }
}


# plots the prediction and state indicator figures for the two Markov switching models 


# --- cases: 
ts.dates.l<-list()
for (i in 1:159){
  ts.dates.l[[i]]<-ts(as.numeric(cases[i,]),
                      start=c(2015,1),end=c(2022,20),frequency=52)
}
districtname<-dataCHIKV$District 


# 235: the end of prediction point
# 209: the start time, which is the first week of 2019

temp.series <- seq.Date(from = as.Date("2019/01/01",format = "%Y/%m/%d"), by = "week", length.out = 23)
temp.series2 <- seq.Date(from = as.Date("2019/06/11" ,format = "%Y/%m/%d"), by = "week", length.out = 8)

a1<- data.frame(Time=temp.series,
                district=c(ts.dates.l[[135]][209:231]))
datebreaks <- seq(as.Date("2018/12/31"), as.Date("2019-07-30"),by = "4 week")

a2<-data.frame(Time2=temp.series2,
               district2=c(ts.dates.l[[135]][232:239]))


# --- state indicators: 
ts.dates.states<-list()
for (i in 1:159){
  ts.dates.states[[i]]<-ts(as.numeric( ifelse(cases[i,]>0,1,0) ),
                           start=c(2015,1),end=c(2022,20),frequency=52)
}


a1.states<- data.frame(Time=temp.series,
                       states=c(ts.dates.states[[135]][209:231]))

a2.states<-data.frame(Time2=temp.series2,
                      states2=c(ts.dates.states[[135]][232:239]))





pred.compare.135.k0<-data.frame(Time2=temp.series2,
                                # ZSMSNB    
                                simulated.pred.zsmsnb.k0[135,],
                                upper.simulated.pred.zsmsnb.k0[135,],
                                lower.simulated.pred.zsmsnb.k0[135,],
                                # state indicator
                                state.indicator=X.simulated.mean.k0[135,],
                                state.indicator.h=X.simulated.mean.k0.h[135,],
                                # ZSMSNBH
                                simulated.pred.zsmsnbh.k0[135,],
                                upper.simulated.pred.zsmsnbh.k0[135,],
                                lower.simulated.pred.zsmsnbh.k0[135,])

colnames(pred.compare.135.k0)<-c("time2","pred.zsmsnb","pred.zsmsnb.upper",
                                 "pred.zsmsnb.lower",
                                 "state.indicator","state.indicator.h","pred.zsmsnbh",
                                 "pred.zsmsnbh.upper","pred.zsmsnbh.lower")



# ZSMSNB: 
pred.val.135<-ggplot(a1)+
  geom_point(aes(x=Time,y=district))+
  geom_point(data=a2,aes(x=Time2,y=district2),shape=21)+
  # ZS-MSNB, pred and 95% CI
  geom_line(data=pred.compare.135.k0,aes(x=time2,y=pred.zsmsnb),size=.8, linetype='solid')+
  # ZS-MSNBH, pred and 95% CI
  geom_line(data=pred.compare.135.k0,aes(x=time2,y=pred.zsmsnbh),size=.8, linetype='longdash')+
  geom_ribbon(data=pred.compare.135.k0,
              aes(x=time2, ymin=pred.zsmsnb.lower, ymax=pred.zsmsnb.upper),color = "grey",
              alpha=0.3)+
  geom_ribbon(data=pred.compare.135.k0,
              aes(x=time2, ymin=pred.zsmsnbh.lower, ymax=pred.zsmsnbh.upper),color = "grey36",
              linetype='dotted', alpha=0.5, lwd=.6)+
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15))+
  theme(plot.title = element_text(vjust=2.2))+
  ggtitle("(a) Vargem Pequena district")+
  xlab('')+
  ylab('Cases')



# states: 
state.ind.135<-ggplot(a1.states,aes(x=Time,y=states))+
  geom_point()+
  geom_point(data=a2.states,aes(x=Time2,y=states2),shape=21)+
  # for ZS-MSNB
  geom_line(data=pred.compare.135.k0,aes(x=time2,y=state.indicator),size=.8, linetype='solid')+
  # for ZS-MSNBH
  geom_line(data=pred.compare.135.k0,aes(x=time2,y=state.indicator.h),size=.8, linetype='longdash')+
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15))+
  scale_y_continuous(limits = c(0, 1))+
  xlab('')+
  ylab('Presence')



ms.district.135<-plot_grid(pred.val.135, 
                  state.ind.135, 
                        ncol = 1, nrow = 2, align="v")

# For district No.150: 


b1<- data.frame(Time=temp.series,
                district=c(ts.dates.l[[150]][209:231]))

b2<-data.frame(Time2=temp.series2,
               district2=c(ts.dates.l[[150]][232:239]))


# --- state indicators: 

b1.states<- data.frame(Time=temp.series,
                       states=c(ts.dates.states[[150]][209:231]))

b2.states<-data.frame(Time2=temp.series2,
                      states2=c(ts.dates.states[[150]][232:239]))


pred.compare.150.k0<-data.frame(Time2=temp.series2,
                                # ZSMSNB    
                                simulated.pred.zsmsnb.k0[150,],
                                upper.simulated.pred.zsmsnb.k0[150,],
                                lower.simulated.pred.zsmsnb.k0[150,],
                                # state indicator
                                state.indicator=X.simulated.mean.k0[150,],
                                state.indicator.h=X.simulated.mean.k0.h[150,],
                                # ZSMSNBH
                                simulated.pred.zsmsnbh.k0[150,],
                                upper.simulated.pred.zsmsnbh.k0[150,],
                                lower.simulated.pred.zsmsnbh.k0[150,])

colnames(pred.compare.150.k0)<-c("time2","pred.zsmsnb","pred.zsmsnb.upper",
                                 "pred.zsmsnb.lower",
                                 "state.indicator","state.indicator.h","pred.zsmsnbh",
                                 "pred.zsmsnbh.upper","pred.zsmsnbh.lower")



pred.val.150<-ggplot(b1)+
  geom_point(aes(x=Time,y=district))+
  geom_point(data=b2,aes(x=Time2,y=district2),shape=21)+
  # ZS-MSNB, pred and 95% CI
  geom_line(data=pred.compare.150.k0,aes(x=time2,y=pred.zsmsnb),size=.8, linetype='solid')+
  # ZS-MSNBH, pred and 95% CI
  geom_line(data=pred.compare.150.k0,aes(x=time2,y=pred.zsmsnbh),size=.8, linetype='longdash')+
  # ribbon 
  geom_ribbon(data=pred.compare.150.k0,
              aes(x=time2, ymin=pred.zsmsnb.lower, ymax=pred.zsmsnb.upper),color = "grey",
              alpha=0.3)+
  geom_ribbon(data=pred.compare.150.k0,
              aes(x=time2, ymin=pred.zsmsnbh.lower, ymax=pred.zsmsnbh.upper),color = "grey36",
              linetype='dotted', alpha=0.5, lwd=.6)+
  
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15))+
  theme(plot.title = element_text(vjust=2.2))+
  ggtitle("(b) Campo Grande district")+
  xlab('')+
  ylab('Cases')



# states: 
state.ind.150<-ggplot(b1.states,aes(x=Time,y=states))+
  geom_point()+
  geom_point(data=b2.states,aes(x=Time2,y=states2),shape=21)+
  # for ZS-MSNB
  geom_line(data=pred.compare.150.k0,aes(x=time2,y=state.indicator),size=.8, linetype='solid')+
  # for ZS-MSNBH
  geom_line(data=pred.compare.150.k0,aes(x=time2,y=state.indicator.h),size=.8,linetype='longdash' )+
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(0, 1))+
  theme(text = element_text(size=15),
        strip.text.x = element_text(size=15),
        strip.text.y = element_text(size=15))+
  xlab('')+
  ylab('Presence')



ms.district.150<-plot_grid(pred.val.150, 
                        state.ind.150, 
                        ncol = 1, nrow = 2, align="v")

figure.7<-plot_grid(ms.district.135, ms.district.150, 
          ncol = 1, nrow = 2, align="v")




#######################################################################
###### Fig 4. Prediction vals for 5 mods vs. the actual cases  ######
#######################################################################
K<- 4

mcos.time <- mean(cos.time[1:T0])
msin.time <- mean(sin.time[1:T0])


#---(C.1) NB: 
M <- nrow(samps.nb)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred by NB
simulated.pred.nb<- matrix(0,nrow=159,ncol=K)
# 95% upper and lower (within T0, for simulated prediction)
upper.simulated.pred.nb<-matrix(0,nrow=159,ncol=K)
lower.simulated.pred.nb<-matrix(0,nrow=159,ncol=K)


for (i in 1:159){
  print(i)
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K){
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.nb[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.nb[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.nb[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.nb[paste0("beta2")]))
    
    alpha8<- as.numeric(unlist(samps.nb[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.nb[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.nb[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.nb[paste0("alpha11")]))
    gamma3<- as.numeric(unlist(samps.nb[paste0("gamma3")]))
    
    
    # mup
    if (k== 1){
      mup.pred.nb<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred-msin.time)+beta2*(cos.time.pred-mcos.time))
    }else {
      mup.pred.nb <-  exp(b0)*simulated.pred[i,k-1,]+exp(b+beta1*(sin.time.pred- msin.time)+beta2*(cos.time.pred- mcos.time))
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
    # prediction matrix
    p_nb <- r_nb/(r_nb+mup.pred.nb) 
    simulated.pred[i,k,]<-rnbinom(M, prob=p_nb, size=r_nb) 
    # the mean and 95% CI
    simulated.pred.nb[i,k]<- mean(simulated.pred[i,k,])
    upper.simulated.pred.nb[i,k]<- quantile(simulated.pred[i,k,], probs = c(.975))
    lower.simulated.pred.nb[i,k]<- quantile(simulated.pred[i,k,], probs = c(.025))
  }
}

#---(C.2) NBH: 

M <- nrow(samps.nbh)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.nbh<- matrix(0,nrow=159,ncol=K)
upper.simulated.pred.nbh<- matrix(0,nrow=159,ncol=K)
lower.simulated.pred.nbh<- matrix(0,nrow=159,ncol=K)

for (i in 1:159){
  print(i)
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K){
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.nbh[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.nbh[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.nbh[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.nbh[paste0("beta2")]))
    alpha4<- as.numeric(unlist(samps.nbh[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.nbh[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.nbh[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.nbh[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.nbh[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.nbh[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.nbh[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.nbh[paste0("alpha11")]))
    gamma2<- as.numeric(unlist(samps.nbh[paste0("gamma2")]))
    gamma3<- as.numeric(unlist(samps.nbh[paste0("gamma3")]))
    
    # mup
    if (k== 1){
      mup.pred.nbh<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred-msin.time)+beta2*(cos.time.pred-mcos.time))
    }else {
      mup.pred.nbh <-  exp(b0)*simulated.pred[i,k-1,]+exp(b+beta1*(sin.time.pred-msin.time)+beta2*(cos.time.pred-mcos.time))
    }
    
    
    # p11
    if (k==1){
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(cases[i,T0]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    
    # r: change the parameters: 
    
    if (k==1){
      r_nbh<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                    alpha10*(pop[i]-mean(pop[1:159]))+
                    alpha11*(log(cases[i,T0]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      r_nbh<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                    alpha10*(pop[i]-mean(pop[1:159]))+
                    alpha11*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    
    # prediction matrix
    p_nbh <- r_nbh/(r_nbh+mup.pred.nbh) 
    
    simulated.pred[i,k,]<- rhurdleNegBinMargV(rep(1,M), p=p11, 
                                              lambda= mup.pred.nbh, r=r_nbh)
    
    # the mean and 95% CI
    simulated.pred.nbh[i,k]<- mean(simulated.pred[i,k,])
    upper.simulated.pred.nbh[i,k]<- quantile(simulated.pred[i,k,], probs = c(.975))
    lower.simulated.pred.nbh[i,k]<- quantile(simulated.pred[i,k,], probs = c(.025))
  }
}

#---(C.3) ZINB: 

M <- nrow(samps.zinb)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.zinb<- matrix(0,nrow=159,ncol=K)
# 95% upper and lower (within T0, for simulated prediction)
upper.simulated.pred.zinb<-matrix(0,nrow=159,ncol=K)
lower.simulated.pred.zinb<-matrix(0,nrow=159,ncol=K)



# X matrix (updating)
X.simulated<- array(0,dim=c(159,K+1,M)) 


for (i in 1:159){
  print(i)
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K){
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.zinb[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.zinb[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.zinb[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zinb[paste0("beta2")]))
    
    
    alpha4<- as.numeric(unlist(samps.zinb[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zinb[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zinb[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.zinb[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.zinb[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.zinb[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.zinb[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.zinb[paste0("alpha11")]))
    gamma2<- as.numeric(unlist(samps.zinb[paste0("gamma2")]))
    gamma3<- as.numeric(unlist(samps.zinb[paste0("gamma3")]))
    
    
    # mup
    if (k== 1){
      mup.pred.zinb<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred-msin.time)+beta2*(cos.time.pred-mcos.time))
    }else {
      mup.pred.zinb <-  exp(b0)*simulated.pred[i,k-1,]+exp(b+beta1*(sin.time.pred-msin.time)+beta2*(cos.time.pred-mcos.time))
    }
    
    
    # p11
    if (k==1){
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(cases[i,T0]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+ alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    
    # r: change the parameters: 
    
    if (k==1){
      r_zinb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                     alpha10*(pop[i]-mean(pop[1:159]))+
                     alpha11*(log(cases[i,T0]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }else{
      r_zinb<- exp(alpha8+ alpha9*(HDI[i]-mean(HDI[1:159]))+
                     alpha10*(pop[i]-mean(pop[1:159]))+
                     alpha11*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma3*(greenarea[i]-mean(greenarea[1:159])))
    }
    
    
    #X
    
    Xk<- rbinom(n= rep(1,M), size=1,
                p=p11)
    
    X.simulated[i,k,]<- Xk
    
    
    p_zinb <- r_zinb/(r_zinb+(Xk)*mup.pred.zinb) - 1e-10*(1-Xk)
    
    simulated.pred[i,k,]<-rnbinom(M, prob=p_zinb, size=r_zinb) 
    # the mean and 95% CI
    simulated.pred.zinb[i,k]<-mean(simulated.pred[i,k,])
    upper.simulated.pred.zinb[i,k]<- quantile(simulated.pred[i,k,],probs = c(.975))
    lower.simulated.pred.zinb[i,k]<- quantile(simulated.pred[i,k,],probs = c(.025))
  }
}
#---(C.4) ZSMSNB: 


M <- nrow(samps.zsmsnb)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.zsmsnb<- matrix(0,nrow=159,ncol=K)
# 95% upper and lower (within T0, for simulated prediction)? 
upper.simulated.pred.zsmsnb<-matrix(0,nrow=159,ncol=K)
lower.simulated.pred.zsmsnb<-matrix(0,nrow=159,ncol=K)
# X matrix (updating)
X.simulated<- array(0,dim=c(159,K+1,M)) 
# spatial constants
adj=as.numeric(as.carAdjacency(N_matrix)$adj)
num=as.numeric(as.carAdjacency(N_matrix)$num)
count=count


for (i in 1:159){
  print(i)
  X.simulated[i,1,]<- as.numeric(unlist(samps.zsmsnb[paste0("X_final.",i,".")]))
  nei <- which(N_matrix[i,]==1)
  
  ###  for t=T0+1,T0+2,...T0+K values 
  for (k in 1:K){
    # pred sin. and cos. constant
    sin.time.pred<-sin((T0+k)*2*pi/52)
    cos.time.pred<-cos((T0+k)*2*pi/52)
    
    b0<- as.numeric(unlist(samps.zsmsnb[paste0("b0.",i,".")]))
    b<-  as.numeric(unlist(samps.zsmsnb[paste0("b.",i,".")]))
    beta1<- as.numeric(unlist(samps.zsmsnb[paste0("beta1")]))
    beta2<- as.numeric(unlist(samps.zsmsnb[paste0("beta2")]))
    alpha1<- as.numeric(unlist(samps.zsmsnb[paste0("alpha1")]))
    alpha2<- as.numeric(unlist(samps.zsmsnb[paste0("alpha2")]))
    alpha3<- as.numeric(unlist(samps.zsmsnb[paste0("alpha3")]))
    alpha4<- as.numeric(unlist(samps.zsmsnb[paste0("alpha4")]))
    alpha5<- as.numeric(unlist(samps.zsmsnb[paste0("alpha5")]))
    alpha6<- as.numeric(unlist(samps.zsmsnb[paste0("alpha6")]))
    alpha7<- as.numeric(unlist(samps.zsmsnb[paste0("alpha7")]))
    alpha8<- as.numeric(unlist(samps.zsmsnb[paste0("alpha8")]))
    alpha9<- as.numeric(unlist(samps.zsmsnb[paste0("alpha9")]))
    alpha10<- as.numeric(unlist(samps.zsmsnb[paste0("alpha10")]))
    alpha11<- as.numeric(unlist(samps.zsmsnb[paste0("alpha11")]))
    
    
    gamma1<- as.numeric(unlist(samps.zsmsnb[paste0("gamma1")]))
    gamma2<- as.numeric(unlist(samps.zsmsnb[paste0("gamma2")]))
    gamma3<- as.numeric(unlist(samps.zsmsnb[paste0("gamma3")]))
    
    delta1<- as.numeric(unlist(samps.zsmsnb[paste0("delta1")]))
    delta2<- as.numeric(unlist(samps.zsmsnb[paste0("delta2")]))
    
    #neighborhood sum
    nei_sum <- rep(0,M)
    if(k==1){
      for(j in nei){
        nei_sum <- nei_sum+samps.zsmsnb[,paste0("X_final.",j,".")]
      }
    } else{
      for(j in nei){
        nei_sum <- nei_sum+X.simulated[j,k,]
      }
    }
    
    # mup
    if (k== 1){
      mup.pred.zsmsnb<- exp(b0)*cases[i,T0]+exp(b+beta1*(sin.time.pred- msin.time )+beta2*(cos.time.pred- mcos.time ))
    }else {
      mup.pred.zsmsnb <-  exp(b0)*simulated.pred[i,k-1,]+exp(b+beta1*(sin.time.pred- msin.time )+beta2*(cos.time.pred- mcos.time ))
    }
    # p01 and p11 should be changed into these, and r: 
    
    # p01
    p01<- expit(alpha1+alpha2*(HDI[i]-mean(HDI[1:159]))+
                  alpha3*(pop[i]-mean(pop[1:159]))+gamma1*(greenarea[i]-mean(greenarea[1:159]))+
                  delta1*nei_sum)
    
    
    # p11: change the parameter: 
    if(k==1){
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(cases[i,T0]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*nei_sum)
    }else{
      p11<- expit(alpha4+alpha5*(HDI[i]-mean(HDI[1:159]))+alpha6*(pop[i]-mean(pop[1:159]))+
                    alpha7*(log(simulated.pred[i,k-1,]+1)-mlpsi)+gamma2*(greenarea[i]-mean(greenarea[1:159]))+
                    delta2*nei_sum)
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
    #X
    Xk<- rbinom(n= rep(1,M), size=1,
                p= X.simulated[i,k,]*p11 + (1-X.simulated[i,k,])*p01)
    X.simulated[i,k+1,]<- Xk
    # prediction matrix
    p_nb <- r_nb/(r_nb+(Xk)*mup.pred.zsmsnb) - 1e-10*(1-Xk)
    simulated.pred[i,k,]<-rnbinom(M, prob=p_nb, size=r_nb) 
    # the mean and 95% CI
    #simulated.pred.zip[i,k]<-mean(simulated.pred[i,k,])
    simulated.pred.zsmsnb[i,k]<- mean( simulated.pred[i,k,])
    upper.simulated.pred.zsmsnb[i,k]<- quantile(simulated.pred[i,k,],probs = c(.975))
    lower.simulated.pred.zsmsnb[i,k]<- quantile(simulated.pred[i,k,],probs = c(.025))
  }
}


#---(C.5) ZSMSNBH: 


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

M <- nrow(samps.zsmsnbh)
# for M samples
simulated.pred<- array(0,dim=c(159,K,M))
# pred cases
simulated.pred.zsmsnbh<- matrix(0,nrow=159,ncol=K)
upper.simulated.pred.zsmsnbh<- matrix(0,nrow=159,ncol=K)
lower.simulated.pred.zsmsnbh<- matrix(0,nrow=159,ncol=K)

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
    
    # mean and 95% CI 
    simulated.pred.zsmsnbh[i,k]<-mean(simulated.pred[i,k,])
    upper.simulated.pred.zsmsnbh[i,k]<- quantile(simulated.pred[i,k,],probs = c(.975))
    lower.simulated.pred.zsmsnbh[i,k]<- quantile(simulated.pred[i,k,],probs = c(.025))
  }
}





#------------ Plot of the prediction vs. actual cases----------#

library(ggplot2)
library("scales")
# timeline
ts.dates.l<-list()
for (i in 1:159){
  ts.dates.l[[i]]<-ts(as.numeric(cases[i,]),
                      start=c(2015,1),end=c(2022,20),frequency=52)
}
districtname<-dataCHIKV$District 


# 235: the end of prediction point
# 209: the start time, which is the first week of 2019

temp.series <- seq.Date(from = as.Date("2018/12/31",format = "%Y/%m/%d"), by = "week", length.out = 23)
temp.series2 <- seq.Date(from = as.Date("2019/06/11" ,format = "%Y/%m/%d"), by = "week", length.out = 4)

a1<- data.frame(Time=temp.series,
                district=c(ts.dates.l[[135]][209:231]))
datebreaks <- seq(as.Date("2018/12/31"), as.Date("2019-07-02"),by = "4 week")
# the 4 predition points: 232, 233, 234, 235

a2<-data.frame(Time2=temp.series2,
               district2=c(ts.dates.l[[135]][232:235]))

pred.compare.135<-data.frame(Time2=temp.series2,
                             simulated.pred.nb[135,],
                             simulated.pred.nbh[135,],
                             simulated.pred.zinb[135,],
                             simulated.pred.zsmsnb[135,],
                             simulated.pred.zsmsnbh[135,])
colnames(pred.compare.135)<-c("time2","pred.nb","pred.nbh","pred.zinb",
                              "pred.zsmsnb","pred.zsmsnbh")

dat.CI<-data.frame(time.ci=temp.series2,
                   # for nb
                   simulated.pred.nb[135,],
                   upper.simulated.pred.nb[135,],
                   lower.simulated.pred.nb[135,],
                   # for nbh
                   simulated.pred.nbh[135,],
                   upper.simulated.pred.nbh[135,],
                   lower.simulated.pred.nbh[135,],
                   # for zinb
                   simulated.pred.zinb[135,],
                   upper.simulated.pred.zinb[135,],
                   lower.simulated.pred.zinb[135,],
                   # for zsmsnb
                   simulated.pred.zsmsnb[135,],
                   upper.simulated.pred.zsmsnb[135,],
                   lower.simulated.pred.zsmsnb[135,],
                   # for zsmsnbh
                   simulated.pred.zsmsnbh[135,],
                   upper.simulated.pred.zsmsnbh[135,],
                   lower.simulated.pred.zsmsnbh[135,])
colnames(dat.CI)<-c("time.ci","nb","upper.pred.nb","lower.pred.nb",
                    "nbh","upper.pred.nbh","lower.pred.nbh",
                    "zinb","upper.pred.zinb","lower.pred.zinb",
                    "zsmsnb","upper.pred.zsmsnb","lower.pred.zsmsnb",
                    "zsmsnbh","upper.pred.zsmsnbh","lower.pred.zsmsnbh")


cols=c("NB"="black","NBH"="#62BC9B","ZINB"='#EBA133',"ZS-MSNBH"="#E1706E","ZS-MSNB"="#6F8CBD")

pred.p135<- ggplot(a1,aes(x=Time,y=district))+
  geom_point()+
  geom_point(data=a2,aes(x=Time2,y=district2),shape=21)+
  # for NB
  geom_line(data=pred.compare.135,aes(x=time2,y=pred.nb,color="NB"),size=1.1)+
  geom_ribbon(data=dat.CI,aes(x=time.ci,y=nb, ymin=lower.pred.nb,ymax=upper.pred.nb),
              fill="black",alpha=.2,color="black",linetype = "dotted")+
  # for NBH
  geom_line(data=pred.compare.135,aes(x=time2,y=pred.nbh,color="NBH"),size=1.1)+
  geom_ribbon(data=dat.CI,aes(x=time.ci,y=nbh, ymin=lower.pred.nbh,ymax=upper.pred.nbh),
              fill="#62BC9B",alpha=.1,color="#62BC9B",linetype = "dotted")+
  # for ZINB
  geom_line(data=pred.compare.135,aes(x=time2,y=pred.zinb,color="ZINB"),size=1.1)+
  geom_ribbon(data=dat.CI,aes(x=time.ci,y=zinb,ymin=lower.pred.zinb,ymax=upper.pred.zinb),
              fill="#EBA133",alpha=.1,color="#EBA133",linetype = "dotted")+
  # for ZS-MSNBH
  geom_line(data=pred.compare.135,aes(x=time2,y=pred.zsmsnbh,color="ZS-MSNBH"),size=1.1)+
  geom_ribbon(data=dat.CI,aes(x=time.ci,y=zsmsnbh,ymin=lower.pred.zsmsnbh,ymax=upper.pred.zsmsnbh),
              fill="#E1706E",alpha=.2,color="#E1706E",linetype = "dotted")+
  # for ZS-MSNB
  geom_line(data=pred.compare.135,aes(x=time2,y=pred.zsmsnb,color="ZS-MSNB"),size=1.1)+
  geom_ribbon(data=dat.CI,aes(x=time.ci,y=zsmsnb,ymin=lower.pred.zsmsnb,ymax=upper.pred.zsmsnb),
              fill="#6F8CBD",alpha=.2,color="#6F8CBD",linetype = "dotted")+
  
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  scale_color_manual(name='Models',
                     values=cols)+
  theme_classic()+
  theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position=c(0.01,1), legend.justification=c(0.01,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  ggtitle("(a) Vargem Pequena district")+
  xlab('')+
  ylab('Cases')
pred.p135




# for district 150:

b1<- data.frame(Time=temp.series,
                district=c(ts.dates.l[[150]][209:231]))
b2<-data.frame(Time2=temp.series2,
               district2=c(ts.dates.l[[150]][232:235]))


pred.compare.150<-data.frame(Time2=temp.series2,
                             simulated.pred.nb[150,],
                             simulated.pred.nbh[150,],
                             simulated.pred.zinb[150,],
                             simulated.pred.zsmsnb[150,],
                             simulated.pred.zsmsnbh[150,])

colnames(pred.compare.150)<-c("time2","pred.nb","pred.nbh","pred.zinb",
                              "pred.zsmsnb","pred.zsmsnbh")

dat.CI.150<- data.frame(time.ci=temp.series2,
                        # for nb
                        simulated.pred.nb[150,],
                        upper.simulated.pred.nb[150,],
                        lower.simulated.pred.nb[150,],
                        # for nbh
                        simulated.pred.nbh[150,],
                        upper.simulated.pred.nbh[150,],
                        lower.simulated.pred.nbh[150,],
                        # for zinb
                        simulated.pred.zinb[150,],
                        upper.simulated.pred.zinb[150,],
                        lower.simulated.pred.zinb[150,],
                        # for zsmsnb
                        simulated.pred.zsmsnb[150,],
                        upper.simulated.pred.zsmsnb[150,],
                        lower.simulated.pred.zsmsnb[150,],
                        # for zsmsnbh
                        simulated.pred.zsmsnbh[150,],
                        upper.simulated.pred.zsmsnbh[150,],
                        lower.simulated.pred.zsmsnbh[150,])


colnames(dat.CI.150)<-c("time.ci","nb","upper.pred.nb","lower.pred.nb",
                        "nbh","upper.pred.nbh","lower.pred.nbh",
                        "zinb","upper.pred.zinb","lower.pred.zinb",
                        "zsmsnb","upper.pred.zsmsnb","lower.pred.zsmsnb",
                        "zsmsnbh","upper.pred.zsmsnbh","lower.pred.zsmsnbh")



cols=c("NB"="black","NBH"="#62BC9B","ZINB"='#EBA133',"ZS-MSNBH"="#E1706E","ZS-MSNB"="#6F8CBD")
pred.p150<- ggplot(b1,aes(x=Time,y=district))+
  geom_point()+
  geom_point(data=b2,aes(x=Time2,y=district2),shape=21)+
  # for NB
  geom_line(data=pred.compare.150,aes(x=time2,y=pred.nb,color="NB"),size=1.1)+
  geom_ribbon(data=dat.CI.150,aes(x=time.ci,y=nb,ymin=lower.pred.nb,ymax=upper.pred.nb),
              fill="black",alpha=.2,color="black",linetype = "dotted")+
  
  # for ZINB
  geom_line(data=pred.compare.150,aes(x=time2,y=pred.zinb,color="ZINB"),size=1.1)+
  geom_ribbon(data=dat.CI.150,aes(x=time.ci,y=zinb,ymin=lower.pred.zinb,ymax=upper.pred.zinb),
              fill="#EBA133",alpha=.1,color="#EBA133",linetype = "dotted")+
  # for ZSMSNBH
  geom_line(data=pred.compare.150,aes(x=time2,y=pred.zsmsnbh,color="ZS-MSNBH"),size=1.1)+
  geom_ribbon(data=dat.CI.150,aes(x=time.ci,y=zsmsnbh,ymin=lower.pred.zsmsnbh,ymax=upper.pred.zsmsnbh),
              fill="#E1706E",alpha=.2,color="#E1706E",linetype = "dotted")+
  # for ZSMSNB
  geom_line(data=pred.compare.150,aes(x=time2,y=pred.zsmsnb,color="ZS-MSNB"),size=1.1)+
  geom_ribbon(data=dat.CI.150,aes(x=time.ci,y=zsmsnb,ymin=lower.pred.zsmsnb,ymax=upper.pred.zsmsnb),
              fill="#6F8CBD",alpha=.2,color="#6F8CBD",linetype = "dotted")+
  # for NBH
  geom_line(data=pred.compare.150,aes(x=time2,y=pred.nbh,color="NBH"),size=1.1)+
  geom_ribbon(data=dat.CI.150,aes(x=time.ci,y=nbh,ymin=lower.pred.nbh,ymax=upper.pred.nbh),
              fill="#62BC9B",alpha=.1,color="#62BC9B",linetype = "dotted")+
  scale_x_date(breaks = datebreaks,labels=date_format("%b %Y"))+
  scale_color_manual(name='Models',
                     values=cols)+
  theme_classic()+
  theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position=c(0.01,1), legend.justification=c(0.01,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))+
  ggtitle("(b) Campo Grande district")+
  xlab('')+
  ylab('Cases')
pred.p150


figure.4<-plot_grid(pred.p135, pred.p150, 
                    ncol = 1, nrow = 2)




###########################################################
#####    Fig 5.  For the rps comparison plot       ########   
###########################################################

# 1. For NB model 
NB.dta<- data.frame()

i<- 280:380
NB.dta<- NULL
for (i in 280:380){
  path<- paste0("/Users/xumingchi/Desktop/paper/files/rps_results/NB/","NB_T0_",i,".csv")
  NB.dta<- cbind(NB.dta, read.csv(file=path, header=T)[,2])
}

NB.K1<-mean(NB.dta[1,])
NB.K2<-mean(NB.dta[2,])
NB.K3<-mean(NB.dta[3,])
NB.K4<-mean(NB.dta[4,])

# 2. For NBH model 
NBH.dta<- data.frame()

i<- 280:380
NBH.dta<- NULL
for (i in 280:380){
  path<- paste0("/Users/xumingchi/Desktop/paper/files/rps_results/NBH/","NBH_T0_",i,".csv")
  NBH.dta<- cbind(NBH.dta, read.csv(file=path, header=T)[,2])
}

NBH.K1<-mean(NBH.dta[1,])
NBH.K2<-mean(NBH.dta[2,])
NBH.K3<-mean(NBH.dta[3,])
NBH.K4<-mean(NBH.dta[4,])


# 3. For ZINB model 
ZINB.dta<- data.frame()

i<- 280:380
ZINB.dta<- NULL
for (i in 280:380){
  path<- paste0("/Users/xumingchi/Desktop/paper/files/rps_results/ZINB/","ZINB_T0_",i,".csv")
  ZINB.dta<- cbind(ZINB.dta, read.csv(file=path, header=T)[,2])
}

ZINB.K1<-mean(ZINB.dta[1,])
ZINB.K2<-mean(ZINB.dta[2,])
ZINB.K3<-mean(ZINB.dta[3,])
ZINB.K4<-mean(ZINB.dta[4,])



# 4. For ZS-MSNB model 
ZSMSNB.dta<- data.frame()


ZSMSNB.dta<- NULL
for (i in c(280:380)){
  path<- paste0("/Users/xumingchi/Desktop/paper/files/rps_results/ZS-MSNB/","ZSMSNB_T0_",i,".csv")
  ZSMSNB.dta<- cbind(ZSMSNB.dta, read.csv(file=path, header=T)[,2])
}

ZSMSNB.K1<-mean(ZSMSNB.dta[1,])
ZSMSNB.K2<-mean(ZSMSNB.dta[2,])
ZSMSNB.K3<-mean(ZSMSNB.dta[3,])
ZSMSNB.K4<-mean(ZSMSNB.dta[4,])



# 5. For ZS-MSNBH model 
ZSMSNBH.dta<- data.frame()


ZSMSNBH.dta<- NULL
for (i in c(280:380)){
  path<- paste0("/Users/xumingchi/Desktop/paper/files/rps_results/ZS-MSNBH/","ZSMSNBH_T0_",i,".csv")
  ZSMSNBH.dta<- cbind(ZSMSNBH.dta, read.csv(file=path, header=T)[,2])
}

ZSMSNBH.K1<-mean(ZSMSNBH.dta[1,])
ZSMSNBH.K2<-mean(ZSMSNBH.dta[2,])
ZSMSNBH.K3<-mean(ZSMSNBH.dta[3,])
ZSMSNBH.K4<-mean(ZSMSNBH.dta[4,])

# create a data frame

rps.df<-round(cbind(c(NB.K1, NB.K2, NB.K3, NB.K4),
                    c(NBH.K1, NBH.K2, NBH.K3, NBH.K4),
                    c(ZINB.K1, ZINB.K2, ZINB.K3, ZINB.K4), 
                    c(ZSMSNB.K1, ZSMSNB.K2, ZSMSNB.K3, ZSMSNB.K4),
                    c(ZSMSNBH.K1, ZSMSNBH.K2, ZSMSNBH.K3, ZSMSNBH.K4)),4)
colnames(rps.df)<- c('NB',"NBH","ZINB","ZS-MSNB","ZS-MSNBH")
rps.df<- as.data.frame(rps.df)
rps.df$time<- seq(1, 4, by=1)

Models<- c(rep("NB",4),
           rep("NBH", 4),
           rep("ZINB", 4),
           rep("ZS-MSNB",4),
           rep("ZS-MSNBH", 4))

rps.val<- c( rps.df$NB, 
             rps.df$NBH,
             rps.df$ZINB,
             rps.df$`ZS-MSNB`,
             rps.df$`ZS-MSNBH`)

Time<- rep(seq(1,4,by=1),5)

rps.df.plt<- data.frame(Models, rps.val, Time)


rps_p<- ggplot(data = rps.df.plt, aes(x=Time, y=rps.val, colour=Models))+
  geom_line()+
  geom_point()+
  scale_colour_brewer(palette = "Set1")+
  theme_classic()+
  labs(x="Forecast horizon (weeks)", y="Averaged ranked probability score")+
  theme(legend.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position=c(0.08, 0.95), legend.justification=c(0.01,1))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text = element_text(size=20))+
  theme(text = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20))
figure.5<-rps_p



# permutation t tests: 

x <- ZSMSNB.dta[1,]
y <- ZSMSNBH.dta[1,] 

B <- 99999
d <- x-y
m0 <- mean(d)
rndmdist <- replicate(B,mean((rbinom(length(d),1,.5)*2-1)*d))

# two tailed p-value:
sum( abs(rndmdist) >= abs(m0))/length(rndmdist)
t.test(ZSMSNB.dta[1,],ZSMSNBH.dta[1,],paired = TRUE)





