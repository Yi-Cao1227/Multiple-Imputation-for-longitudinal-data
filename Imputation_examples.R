
############################################################
##   Multiple imputation for longitudindal data          ###
############################################################


####################### Methods using wide-format data  ######################
load("toy_data.Rdata")
summary(sim.wide.data)

######### Fully conditional specification (FCS-Standard) ###
library(mice)

#incomplete variables: "hc_hosptstay","num.paydev","num.disease","BMI"


#initiate imputation
imp.0<-mice(sim.wide.data[,-1],maxit=0)
methods<-imp.0$method

#define univariate conditional models for incomplete variables
methods  #default options
# hc_hosptstay_2   num.paydev_2  num.disease_2     BMI_2    hc_hosptstay_3   num.paydev_3    num.disease_3    BMI_3   hc_hosptstay_4 
# "logreg"          "pmm"          "pmm"          "pmm"       "logreg"          "pmm"          "pmm"          "pmm"       "logreg" 
# num.paydev_4  num.disease_4    BMI_4           hc_hosptstay_1   num.paydev_1  num.disease_1  BMI_1        mc_meds         Age 
# "pmm"          "pmm"          "pmm"             ""             ""             ""             ""             ""             "" 
# hc_health   ss_painbothr   ip_covmedcad      r1dgender   ew_medpaovtm   rl1dracehisp 
# ""             ""             ""             ""             ""             "" 

## customize univariate models (Example:change predictive mean matching (pmm) to linear regressions)
methods[which(methods=="pmm")]<-"norm"

#check covariates used for imputing each incomplete variable 
imp.0$predictorMatrix
 
# perform imputation 
imputation<-mice(sim.wide.data[,-1],method = methods,print=T,m=5)
imp.data.mice<-imputation$imp



######### Joint modeling: General location model ########
library(mix)

## prepare for data--> re-arrange data frame into a matrix; place categorical variables into first p columns,then continuous variables 
#binary variables: factor->numeric (code 1,0 --> 1,2)
mat<-apply(sim.wide.data[,c("hc_hosptstay_1","hc_hosptstay_2","hc_hosptstay_3","hc_hosptstay_4","rl1dracehisp")],2,function(x)ifelse(x=="1",1,2))
#
mat<-as.matrix(cbind(mat,sim.wide.data[,c("num.paydev_2","num.disease_2","BMI_2","num.paydev_3","num.disease_3","BMI_3","num.paydev_4",
                                          "num.disease_4","BMI_4")]))

# first p columns are categorical variables (p=5);
pre.data<-prelim.mix(mat,p=5)
theta.hat<-em.mix(pre.data)
rngseed(123)
new.theta<-da.mix(pre.data,theta.hat,steps=1000,showits = F)

# multiple imputation
imp.data.temp<-list()
for(q in 1:5){
  imp.data<-imp.mix(pre.data,new.theta,x=mat)
  imp.data.temp[[q]]<-imp.data
}

# imputed datasets:
imp.data.all<-list()
for(j in 1:5){
  
  imp.data.j<-imp.data.temp[[j]]
  imp.data.j<-cbind(spid=sim.wide.data$spid,imp.data.j)
  
  ## rounding to integer
  for(g in c("num.paydev_2","num.disease_2","num.paydev_3","num.disease_3","num.paydev_4","num.disease_4")){
    
    imp.data.j[,g]<-ifelse(imp.data.j[,g]<0,0,round(imp.data.j[,g]))
    
  }
  imp.data.all[[j]]<-imp.data.j
}

## check imputed datasets
summary(imp.data.all[[1]])


####################### Methods using long-format data  ######################

######### Fully conditional specification (FCS-GLMM, FCS-LMM) ###
library(mice)
library(micemd)

#Find default imputation models for each incomplete variable, set subject id as cluster indicator
imp.mice.2l.methods<-find.defaultMethod(sim.long.data,ind.clust = 1)
imp.mice.2l.methods
# spid       mc_meds        redsid           Age     hc_health  ss_painbothr  ip_covmedcad  ew_medpaovtm     r1dgender 
# ""            ""            ""            ""            ""            ""            ""            ""            "" 
# rl1dracehisp  Round   hospital.stay    num.paydev      BMI            num.disease 
# ""            ""       "2l.jomo"       "2l.glm.norm"  "2l.glm.norm"   "2l.glm.norm" 

# Model reference:
#2l.jomo: linear mixed model with latent variable model
#2l.glm.norm: linear mixed model
#2l.glm.bin: logistic mixed model
#2l.glm.pois: Poisson mixed model

#For example: set customize models
imp.mice.2l.methods[c(12,13,15)]<-c("2l.glm.bin","2l.glm.pois","2l.glm.pois")
imp.2l.0<-mice(sim.long.data,method=imp.mice.2l.methods,maxit=0)

#set covariate to be included in the model 
pred.mat<-imp.2l.0$predictorMatrix

# set class variable:"spid" column-->-2  #
# set fixed effects variables -->1 #
# set random effects variables-->2 #
pred.mat[c("hospital.stay","num.paydev","BMI","num.disease"),"spid"]<--2
pred.mat[c("hospital.stay","num.paydev","BMI","num.disease"),-1]<-1
pred.mat[c("hospital.stay","num.paydev","BMI","num.disease"),"Round"]<-2

#set the diagonal value of incomplete variables to be 0 
pred.mat[c("hospital.stay"),c("hospital.stay")]<-0
pred.mat[c("num.paydev"),c("num.paydev")]<-0
pred.mat[c("num.disease"),c("num.disease")]<-0
pred.mat[c("BMI"),c("BMI")]<-0

sim.long.data$spid<-as.numeric(as.factor(sim.long.data$spid))
# multiple imputation (This method is very slow)
imputation_long_fcs<-mice(sim.long.data,method=imp.mice.2l.methods,pred=pred.mat,print=T)

imp.data<-imputation_long_fcs$imp


######### Joint modeling (JM-multivariate linear mixed model) ########
library(jomo)

## set incomplete variables
Y<-sim.long.data[,c("hospital.stay","num.paydev","BMI","num.disease")]

## set clusters (subject-level random effects)
clus<-sim.long.data[,c("spid")]

# complete covariates
X<-sim.long.data[,c("mc_meds","redsid","Age","hc_health","ss_painbothr","ip_covmedcad","ew_medpaovtm","r1dgender","rl1dracehisp","Round")]

# multiple imputation: linear mixed model 
imp.jomo<-jomo(Y,clus=clus,X=X, nimp = 5, meth="random")


########## Joint modeling MGLMM (shared-parameter model) ###
library(rjags)
library(reshape2)
library(tidyr)
library(dplyr)

Y.var<-c("hc_hosptstay","num.disease","BMI", "num.paydev")
X.var<-c("Age","r1dgender","rl1dracehisp","hh_martlstat","hh_dhshldnum", 
           "hc_health","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
 
           
## impute function ##           
cat(" 
    
    data{
    # separate random intercept
    zero.b[1]<-0
    zero.b[2]<-0
    
    }
    
    model
    {
    ## by subject random effect (J: total number of subjects)
    
    for(j in 1:J){
    
    b[j,1:2]~dmnorm(zero.b,invSigma.b)

    }
    
    ## define model for each observational unit
    
    for(i in 1:N){

    ## hospital stay
    logit(p[i])<- beta.1[1]+b[subj[i],1]+inprod(beta.1[2:np],X[i,])
    Y1[i] ~ dbern(p[i])

    ##num of disease   
   log(lambda.Y2[i])<-beta.2[1]+b[subj[i],2]+inprod(beta.2[2:np],X[i,])
   Y2[i]~dpois(lambda.Y2[i])

    }
    
    
    # Fixed intercept and slope (uninformative)
    for(q in 1:np){
    beta.1[q] ~ dnorm(0.0,1.0E-2)
    beta.2[q] ~ dnorm(0.0,1.0E-2)

    }
    
    # prior on the by-subjects random effects
    R.b[1,1] <- pow(sigma.1,2)
    R.b[2,2] <- pow(sigma.2,2)
    R.b[1,2] <- rho.b*sigma.1*sigma.2
    R.b[2,1] <- rho.b*sigma.1*sigma.2
    
    invSigma.b ~ dwish(R.b, 2)
    Sigma.b <- inverse(invSigma.b)
    
    # varying intercepts
    tau.1 ~ dgamma(1.5, pow(1.0E-4))
    tau.2 ~ dgamma(1.5, pow(1.0E-4))
    sigma.1 <- pow(tau.1,-1/2)
    sigma.2 <- pow(tau.2,-1/2)

    # correlation
    rho.b ~ dunif(-1,1)

    #rho.b ~ dnorm(0,tau_rho.b)T(-1,1)
    # rho.b ~ dnorm(mu_rho.b,tau_rho.b)T(-1,1)
    # mu_rho.b ~ dunif(-1,1)
    # tau_rho.b ~ dgamma(1.5,10E-2)
    
    }"
,file="imp_JM_SP.jag")


imp.mglmm<-function(sim.long){
  data.long<-data.frame(spid=sim.long$spid)
  data.long$Age_2<-ifelse(sim.long$Age=="2",1,0)
  data.long$Age_3<-ifelse(sim.long$Age=="3",1,0)
  data.long$Age_4<-ifelse(sim.long$Age=="4",1,0)
  data.long$Age_5<-ifelse(sim.long$Age=="5",1,0)
  data.long$Age_6<-ifelse(sim.long$Age=="6",1,0)
  data.long$gender_male<-ifelse(sim.long$r1dgender=="1",1,0)
  data.long$rl1dracehisp_1<-ifelse(sim.long$rl1dracehisp=="1",1,0)
  data.long$hh_martlstat_1<-ifelse(sim.long$hh_martlstat=="1",1,0)
  data.long$hh_dhshldnum<-sim.long$hh_dhshldnum
  data.long$hc_health_2<-ifelse(sim.long$hc_health=="2",1,0)
  data.long$hc_health_3<-ifelse(sim.long$hc_health=="3",1,0)
  data.long$hc_health_4<-ifelse(sim.long$hc_health=="4",1,0)
  data.long$hc_health_5<-ifelse(sim.long$hc_health=="5",1,0)
  data.long$ss_painbothr_1<-ifelse(sim.long$ss_painbothr=="1",1,0)
  data.long$ip_covmedcad_1<-ifelse(sim.long$ip_covmedcad=="1",1,0)
  data.long$ew_medpaovtm_1<-ifelse(sim.long$ew_medpaovtm=="1",1,0)
  data.long$mc_meds_1<-ifelse(sim.long$mc_meds=="1",1,0)
  
  data.long$num.paydev<-sim.long$num.paydev
  data.long$BMI<-sim.long$BMI
  data.long$Round<-sim.long$Round

  mat.X=as.matrix(data.long)

  p<-ncol(mat.X)
  n<-length(unique(sim.long$spid))
  data.list<-list(subj=rep(1:n,each=4),
                  Y1=as.numeric(as.character(sim.long$hospital.stay)),
                  Y2=sim.long$num.disease,
                  X=mat.X[,-1],
                  N=nrow(data.long),
                  J=length(unique(data.long$spid)),
                  np=p)
  
  ### Initial values ###
  sim.long0<-sim.long
  sim.long0$hospital.stay<-sim.long$hospital.stay
  sim.long0$num.disease<-sim.long$num.disease
  
  
  sim.long0$r1dgender<-as.factor(ifelse(sim.long$r1dgender=="1",1,0))
  sim.long0$ss_painbothr<-as.factor(ifelse(sim.long$ss_painbothr=="1",1,0))
  sim.long0$ip_covmedcad<-as.factor(ifelse(sim.long$ip_covmedcad=="1",1,0))
  sim.long0$ew_medpaovtm<-as.factor(ifelse(sim.long$ew_medpaovtm=="1",1,0))
  sim.long0$mc_meds<-as.factor(ifelse(sim.long$mc_meds=="1",1,0))
  
  
  fit0.binary<-glmer(hospital.stay~Age+r1dgender+rl1dracehisp+hh_martlstat+
                       hh_dhshldnum+hc_health+ss_painbothr+ip_covmedcad+
                       ew_medpaovtm+mc_meds+num.paydev+BMI+Round+(1|spid),
                     family = binomial("logit"),data=sim.long0,
                     nAGQ=0)
  
  res.Y1<-summary(fit0.binary)
  beta.1<-res.Y1$coefficients[,1]
  sigma.1<-res.Y1$varcor$spid[1,1] #residuals
  start.values.Y1<-list(beta=beta.1,sigma.1=sigma.1)
  
  fit0.poisson<-glmer(num.disease~Age+r1dgender+rl1dracehisp+hh_martlstat+
                       hh_dhshldnum+hc_health+ss_painbothr+ip_covmedcad+
                       ew_medpaovtm+mc_meds+num.paydev+BMI+Round+(1|spid),
                     family = poisson,data=sim.long0,
                     nAGQ=0)
  
  res.Y2<-summary(fit0.poisson)
  beta.2<-res.Y2$coefficients[,1]
  sigma.2<-res.Y2$varcor$spid[1,1] #residuals
  start.values.Y2<-list(beta=beta.2,sigma.2=sigma.2)
  
  start.values<-list(beta.1=start.values.Y1$beta,beta.2=start.values.Y2$beta)
  
  fit.imp <- jags.model(file ="imp_JM_SP.jag",data = data.list,n.chains = 1,quiet=F)
  
  id.na<-which(is.na(sim.long$hospital.stay))
  
  track.var<-c("beta1","beta2","rho.b",
               "sigma.1","sigma.2",paste0("b[1]","[",id.na,"]"),paste0("b[2]","[",id.na,"]"),
               paste0("p","[",id.na,"]"),paste0("lambda.Y2","[",id.na,"]"),
               paste0("Y1","[",id.na,"]"),paste0("Y2","[",id.na,"]"))
  
  fit.imp.sample <- coda.samples(fit.imp,var = track.var,n.iter = 5000)
  
  id.imp<-which(is.na(sim.long$hospital.stay))     
  id.Y1<-paste0("pred1[",id.imp,"]")
  id.Y2<-paste0("pred2[",id.imp,"]")
  id.Y3<-paste0("pred3[",id.imp,"]")
  imp.Y1<-as.matrix(fit.imp.sample[[1]][,id.Y1])
  imp.Y2<-as.matrix(fit.imp.sample[[1]][,id.Y2])
  imp.Y3<-as.matrix(fit.imp.sample[[1]][,id.Y3])
  
  id.multi<-seq(3000,5000,length.out = 5)
  imp.pay<-imp.Y1[id.multi,]
  imp.hos<-imp.Y2[id.multi,]
  imp.dis<-imp.Y3[id.multi,]
  
  result<-list()
  
  for(k in 1:5){
    
    imp.data<-sim.long
    imp.data$num.paydev[id.imp]<-imp.pay[k,]
    imp.data$hospital.stay[id.imp]<-ifelse(imp.hos[k,]==1,"1","2")
    imp.data$num.disease[id.imp]<-imp.dis[k,]
    
    imp.data$hospital.stay<-as.factor(imp.data$hospital.stay)
    
    result[[k]]<-imp.data
  }
  
  return(result)
  
}

## Perform imputation
imp.data.all<-imp.mglmm(sim.long.data)

