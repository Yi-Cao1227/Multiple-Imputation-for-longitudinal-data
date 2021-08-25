
# The following codes is used to conduct a simulation study to 
# compare currently available imputation methods for multivariate 
# longitudinal data.

# The simulation study supports three missing data configurations 
# Configuration 1 contains two incomplete variables
# Configuration 2 containsthree incomplete variables
# Configuration 3 contains four incomplete variables


## Before running the codes, one needs to download necessary data files 
## and save them in a working directory. Then, run the following R codes within that working dorectory.  

## codes for setting working directory
# uncomment to run: setwd("folder_location_path") #replace folder_path with actual folder location path. 


#####################################################################################
##                      Step 1: Data generating process                            ##
##                          Generate Simulations (MAR)                             ##
#####################################################################################


## loading required packages ## if packges not installed, installed them first ##

## to installing packages run: 
#  install.packages(c("reshape2","tidyr","dplyr"))

library(reshape2)
library(tidyr)
library(dplyr)
library(mice)

## Read in NHATS Rounds 1-4 data  ##
Sample_data_wide <- read.csv("Sample_data_complete_wide.csv")
Sample_data_long <- read.csv("Sample_data_complete_long.csv")

## Define functions to generate simulations  ##

## Function1: Convert long format to wide format ##
# function inputs:
# data: a long-formatted longitudinal data
# var.list: a vector of variable names
convert.wide<-function(data,var.list){
  
  temp.data<-data[,c("spid","Round",var.list)]
  wide.data<-data.frame(spid=unique(temp.data$spid))
  
  for(k in 1:length(var.list)){
    
    temp.wide<-dcast(temp.data,spid~Round,value.var = var.list[k])
    colnames(temp.wide)<-c("spid",paste(var.list[k],1:4,sep="_"))
    wide.data<-cbind(wide.data,temp.wide[,2:5])
    
  }
  
  return(wide.data)
  
}

## Function2: generate simulated data under MAR assumption
# Function inputs:
# longdata: long formatted original data 
# baseline.data: wide formatted original data 
# Y.var: a vector of characters including the name of incomplete variables
# X.var: a vector of characters including the name of baseline covariates used for modeling the drop-out process 
# p=c(a,b,c): a vector of three probabilities indicating the proportion of new drop-outs at round 2, 3, 4

simulate.data.mar<-function(longdata,baseline.data,Y.var, X.var, p=c(0.2,0.15,0.1)){
  
  #### prepare dataset: subset data to include pre-specified incomplete and covariates (wide data) #####
  temp.long<-longdata[,c("spid","Round",Y.var)]
  temp.baseline<-baseline.data[,c("spid",X.var)] #baseline information
  temp.wide<-convert.wide(temp.long,var.list=c(Y.var))
  
  temp.wide<-apply(temp.wide,2,function(x) ifelse(x==-7|x==-8,NA,x)) #code -7,-8 to NA
  temp.wide<-as.data.frame(temp.wide)
  
  ## check data missing pattern and locate drop-out people ##
  miss.pattern<-NULL
  for(v in seq(2,ncol(temp.wide),4)){
    
    temp.var<-temp.wide[,c(1,v:(v+3))]
    temp.miss.ind<-is.na(temp.var[,-1])
    miss.p<-apply(temp.miss.ind,1,function(x) paste(as.numeric(x),sep="-",collapse = "-"))
    
    miss.pattern<-cbind(miss.pattern,miss.p)
    
  }
  
  id.all.dropout<-list()
  for(l in 1:ncol(miss.pattern)){
    id.dropout<-temp.wide$spid[which(miss.pattern[,l]%in%c("0-0-0-0","0-0-0-1","0-0-1-1","0-1-1-1"))]
    id.all.dropout[[l]]<-id.dropout
  }
  
  id.dropouts<-Reduce(intersect,id.all.dropout) #5290
  id.inter<-temp.wide$spid[-which(temp.wide$spid%in%id.dropouts)] #Only 19 people had intermittent missing; Remove them)
  
  ### data only keep dropouts ###
  temp.wide3<-temp.wide[which(temp.wide$spid%in%id.dropouts),]  
  miss.pattern<-miss.pattern[which(temp.wide$spid%in%id.dropouts),] # check missing pattern (should be a block monotone misisng pattern )
  
  ## Create drop-out indicators at each round: miss.at2, miss.at3, miss.at4, for different configuration##
  if(length(Y.var)==2){
    comp<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-0"&x[2]=="0-0-0-0",1,0))
    miss.at2<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-1-1-1"&x[2]=="0-1-1-1",1,0))
    miss.at3<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-1-1"&x[2]=="0-0-1-1",1,0))
    miss.at4<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-1"&x[2]=="0-0-0-1",1,0))
  }else if(length(Y.var)==3){
    comp<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-0"&x[2]=="0-0-0-0"&x[3]=="0-0-0-0",1,0))
    miss.at2<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-1-1-1"&x[2]=="0-1-1-1"&x[3]=="0-1-1-1",1,0))
    miss.at3<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-1-1"&x[2]=="0-0-1-1"&x[3]=="0-0-1-1",1,0))
    miss.at4<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-1"&x[2]=="0-0-0-1"&x[3]=="0-0-0-1",1,0))  
  }else if(length(Y.var==4)){
    comp<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-0"&x[2]=="0-0-0-0"&x[3]=="0-0-0-0"&x[4]=="0-0-0-0",1,0))
    miss.at2<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-1-1-1"&x[2]=="0-1-1-1"&x[3]=="0-1-1-1"&x[4]=="0-1-1-1",1,0))
    miss.at3<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-1-1"&x[2]=="0-0-1-1"&x[3]=="0-0-1-1"&x[4]=="0-0-1-1",1,0))
    miss.at4<-apply(miss.pattern, 1, function(x) ifelse(x[1]=="0-0-0-1"&x[2]=="0-0-0-1"&x[3]=="0-0-0-1"&x[4]=="0-0-0-1",1,0))  
  }
  
  temp.wide3$comp<-comp
  temp.wide3$miss.at2<-miss.at2
  temp.wide3$miss.at3<-miss.at3
  temp.wide3$miss.at4<-miss.at4
  
  id.keep<-which(comp==1|miss.at2==1|miss.at3==1|miss.at4==1)
  temp.wide4<-temp.wide3[id.keep,] 
  
  ## Get covariate matrix consisting of (baseline information)
  temp.X<-apply(temp.baseline,2,function(x)ifelse(x==-7|x==-8,NA,x))
  temp.X<-na.omit(temp.X) #remove people with missing baseline information
  temp.X<-as.data.frame(temp.X)
  
  temp.wide.final<-merge(temp.wide4,temp.X,by="spid") #final dataset contains 4839 people
  
  
  ###### Modeling the drop-out process ########
  #### Fit a logistic model for the drop-out indicator at round 2, using data from round 1 #######
  fit.data.r2<-temp.wide.final
  
  cat.names<-c( "hc_hosptstay_1", "mc_meds","Age","hc_health", "ss_painbothr",  "ip_covmedcad",   
                "ew_medpaovtm" ,  "r1dgender", "rl1dracehisp","hh_martlstat")
  for(j in cat.names){
    
    fit.data.r2[,j]<-as.factor(fit.data.r2[,j])
    
  }
  
  ## fit logistic model ##
  covar.list<-c(paste(Y.var,1,sep="_"),X.var)
  fit.data.r2<-fit.data.r2[,c("spid","miss.at2",covar.list)]
  
  fit.mat.2<-glm(miss.at2~.,data=fit.data.r2,family = binomial)# for matrix (not the real model)
  
  fit.r2<-glm(miss.at2~.,data=fit.data.r2[,-1],family = binomial)
  
  #### Fit a logistic model for the drop-out indicator at round 3 using data from round 1&2 ###
  fit.data.r3<-temp.wide.final
  
  cat.names<-c("hc_hosptstay_1", "hc_hosptstay_2", "mc_meds","Age","hc_health", "ss_painbothr",  "ip_covmedcad",   
               "ew_medpaovtm" ,  "r1dgender", "rl1dracehisp","hh_martlstat")
  for(j in cat.names){
    
    fit.data.r3[,j]<-as.factor(fit.data.r3[,j])
    
  }
  
  covar.list<-c(paste(Y.var,1,sep="_"),paste(Y.var,2,sep="_"),X.var)
  fit.data.r3<-fit.data.r3[,c("spid","miss.at3",covar.list)]
  
  ## fit logistic model ##
  fit.mat.3<-glm(miss.at3~.,data=fit.data.r3,family = binomial)
  fit.r3<-glm(miss.at3~.,data=fit.data.r3[,-1],family = binomial)  
  
  ### Fit a logistic model for the drop-out indicator at round 4 using data from round 1, 2, 3 ###
  fit.data.r4<-temp.wide.final  
  
  cat.names<-c( "hc_hosptstay_1", "hc_hosptstay_2", "hc_hosptstay_3", 
                "mc_meds","Age","hc_health", "ss_painbothr",  "ip_covmedcad",   
                "ew_medpaovtm" ,  "r1dgender", "rl1dracehisp","hh_martlstat")
  for(j in cat.names){
    
    fit.data.r4[,j]<-as.factor(fit.data.r4[,j])
    
  }
  
  covar.list<-c(paste(Y.var,1,sep="_"),paste(Y.var,2,sep="_"),paste(Y.var,3,sep="_"),X.var)
  fit.data.r4<-fit.data.r4[,c("spid","miss.at4",covar.list)]
  
  ## fit logistic model ##
  fit.mat.4<-glm(miss.at4~.,data=fit.data.r4[,],family = binomial)
  fit.r4<-glm(miss.at4~.,data=fit.data.r4[,-1],family = binomial)
  
  ## obtain the coefficients estimates of each dro-out model
  coef.lp.r2<-coef(fit.r2)
  coef.lp.r3<-coef(fit.r3)
  coef.lp.r4<-coef(fit.r4)
  
  
  ###########################################################
  ###    Predict drop-out status on the complete data     ###
  ###########################################################
  
  comp.data<-temp.wide.final[temp.wide.final$comp==1,-which(colnames(temp.wide.final)%in%c("comp","miss.at2","miss.at3","miss.at4"))] 
  
  N<-length(unique(comp.data$spid))
  spid.sample<-unique(comp.data$spid)
  complete.data<-comp.data
  
  simulate.data<-complete.data
  
  ## Predict drop-out indicator at round 2 
  comp.data.r2<-model.matrix(fit.mat.2)
  comp.data.r2<-as.data.frame(comp.data.r2)
  comp.data.r2<-comp.data.r2[which(comp.data.r2$spid%in%spid.sample),-2]
  
  mat.X<-as.matrix(comp.data.r2)[,-1]
  lp<-mat.X%*%coef.lp.r2[-1]
  inv.logit<-function(x){
    exp(x)/(1+exp(x))
  }
  
  # calculate intercept parameter alpha0 to gurantee prespecified proportion of missing at each round
  f1<-function(x,p){
    sum(inv.logit((rep(x)+lp)))-p*nrow(mat.X)
  }
  alpha0<-uniroot(f1,p[1],interval=c(-10,10))$root
  probs<-exp(alpha0+lp)/(1+exp(alpha0+lp))
  R.mar.2<-rbinom(length(spid.sample),1,prob=probs)
  id.dp2<-which(R.mar.2==1)
  id.spid.r2<-spid.sample[id.dp2]
  spid.sample<-spid.sample[-id.dp2]
  
  ###############################################################  
  
  ## Predict drop-out indicator at round 3
  comp.data.r3<-model.matrix(fit.mat.3)
  comp.data.r3<-as.data.frame(comp.data.r3)
  comp.data.r3<-comp.data.r3[which(comp.data.r3$spid%in%spid.sample),-2]
  
  mat.X<-as.matrix(comp.data.r3)[,-1]
  lp<-mat.X%*%coef.lp.r3[-1]
  alpha0<-uniroot(f1,p[2],interval=c(-10,10))$root
  probs<-exp(alpha0+lp)/(1+exp(alpha0+lp))
  R.mar.3<-rbinom(length(spid.sample),1,prob=probs)
  id.dp3<-which(R.mar.3==1)
  id.spid.r3<- spid.sample[id.dp3]
  spid.sample<-spid.sample[-id.dp3]
  
  ################################################################
  
  ## Predict drop-out indicator at round 4 
  comp.data.r4<-model.matrix(fit.mat.4)
  comp.data.r4<-as.data.frame(comp.data.r4)
  comp.data.r4<-comp.data.r4[which(comp.data.r4$spid%in%spid.sample),-2]
  
  mat.X<-as.matrix(comp.data.r4)[,-1]
  lp<-mat.X%*%coef.lp.r4[-1]
  alpha0<-uniroot(f1,p[3],interval=c(-10,10))$root
  probs<-exp(alpha0+lp)/(1+exp(alpha0+lp))
  R.mar.4<-rbinom(length(spid.sample),1,prob=probs)
  id.dp4<-which(R.mar.4==1)
  id.spid.r4<-spid.sample[id.dp4]
  
  
  #### Remove observed values  #### 
  id.y1<-seq(2,ncol(temp.wide),4)
  for(k in id.y1){
    var.names<-colnames(temp.wide)[k:(k+3)]
    simulate.data[which(simulate.data$spid%in%id.spid.r2),c(var.names[2:4])]<-rep(NA)
    simulate.data[which(simulate.data$spid%in%id.spid.r3),c(var.names[3:4])]<-rep(NA)
    simulate.data[which(simulate.data$spid%in%id.spid.r4),c(var.names[4])]<-NA
  }
  
  ## final simulated dataset ##
  simulate.data.MAR<-simulate.data
  
  # return 2 objects: 1. the true complete data, 2. the simulated data
  return(list(comp=complete.data,sims=simulate.data.MAR))
  
}

##  Example for creating a simulated dataset with four incomplete variables 
# Specify incomplete variables and complete covariates
Y.var<-c("hc_hosptstay","num.disease","BMI", "num.paydev")
X.var<-c("Age","r1dgender","rl1dracehisp","hh_martlstat","hh_dhshldnum", 
         "hc_health","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
## generate data ##
sim.data.mar<-simulate.data.mar(longdata=Sample_data_long,baseline.data = Sample_data_wide ,Y.var=Y.var, X.var=X.var, p=c(0.2,0.15,0.1))

## save simulated data for imputation step ##
save(sim.data.mar,file="sim_mar_data.Rdata")

# Codes for generating 500 replications
# for(i in 1:500){
#   set.seed(i)
#   sim1.data.mar<-simulate.data.mar(longdata=Sample_data_long,baseline.data = Sample_data_wide ,Y.var=Y.var, X.var=X.var, p=c(0.2,0.15,0.1))
#   save(sim1.data.mar,file=paste0("sim_mar_data_study3_",i,".Rdata"))
# }


#####################################################################################
##                      Step 2: Multiple Imputation                                ##
#####################################################################################

################################################################
## Methods using the fully conditional specification strategy ##
################################################################

## Method 1. FCS-Standard: impute data arranged in a wide format, which treats the repeated measurements of the same variable at different time points as distinct variables                        #
## Imputation models for each incomplete variable:
## BMI: linear regression
## Hospital stay: logistic regression
## Assistive devices: linear regression, then rounding
## Comorbidity index: linear regression, then rounding

#load package 
library(mice)

# prepare datasets
# or run load("sim.data.mar.Rdata"); sim.data<-sim.data.mar$sims;true.data<-sim.data.mar$comp
true.data<-sim.data.mar$comp
sim.data<-sim.data.mar$sims

# make binary variables as factors
id.cont<-c("num.paydev_1","num.paydev_2","num.paydev_3","num.paydev_4",
             "num.disease_1","num.disease_2","num.disease_3","num.disease_4",
             "BMI_1", "BMI_2", "BMI_3", "BMI_4","hh_dhshldnum")
id.factor<-which(!colnames(sim.data)%in%id.cont)
for(l in id.factor){
  sim.data[,l]<-as.factor(sim.data[,l])
}


#initiate the multiple imputation to get default imputation models for each incomplete variables
imp.0<-mice(sim.data[,-1],maxit=0)
methods<-imp.0$method 

#check default methods 
methods
# hc_hosptstay_1 hc_hosptstay_2 hc_hosptstay_3 hc_hosptstay_4  num.disease_1  num.disease_2  num.disease_3  num.disease_4 
# ""       "logreg"       "logreg"       "logreg"             ""          "pmm"          "pmm"          "pmm" 
# BMI_1          BMI_2          BMI_3          BMI_4   num.paydev_1   num.paydev_2   num.paydev_3   num.paydev_4 
# ""          "pmm"          "pmm"          "pmm"             ""          "pmm"          "pmm"          "pmm" 
# Age      r1dgender   rl1dracehisp   hh_martlstat   hh_dhshldnum      hc_health   ss_painbothr   ip_covmedcad 
# ""             ""             ""             ""             ""             ""             ""             "" 
# ew_medpaovtm        mc_meds 
# ""             "" 

#change predictive mean matching (pmm) to linear regression (norm) for imputing continuous variables
methods[which(methods=="pmm")]<-"norm"
# one can customize this step by specifying different univariate imputation models for each variable

# Run imputation step
imp.step<-mice(sim.data[,-1],method = methods,print=F,m=5)
imp.data.mice<-imp.step$imp

# Replace missing values with imputed values (here, we generate 5 multiple imputations)
MI.imputed.data<-NULL
for(j in 1:5){
  
  imp.data.new<-sim.data
  imp.data.new$hc_hosptstay_2[which(is.na(imp.data.new$hc_hosptstay_2))]<-imp.data.mice$hc_hosptstay_2[[j]]
  imp.data.new$hc_hosptstay_3[which(is.na(imp.data.new$hc_hosptstay_3))]<-imp.data.mice$hc_hosptstay_3[[j]]
  imp.data.new$hc_hosptstay_4[which(is.na(imp.data.new$hc_hosptstay_4))]<-imp.data.mice$hc_hosptstay_4[[j]]
  
  imp.data.new$num.paydev_2[which(is.na(imp.data.new$num.paydev_2))]<-imp.data.mice$num.paydev_2[[j]]
  imp.data.new$num.paydev_3[which(is.na(imp.data.new$num.paydev_3))]<-imp.data.mice$num.paydev_3[[j]]
  imp.data.new$num.paydev_4[which(is.na(imp.data.new$num.paydev_4))]<-imp.data.mice$num.paydev_4[[j]]
  
  imp.data.new$num.disease_2[which(is.na(imp.data.new$num.disease_2))]<-imp.data.mice$num.disease_2[[j]]
  imp.data.new$num.disease_3[which(is.na(imp.data.new$num.disease_3))]<-imp.data.mice$num.disease_3[[j]]
  imp.data.new$num.disease_4[which(is.na(imp.data.new$num.disease_4))]<-imp.data.mice$num.disease_4[[j]]
  
  imp.data.new$BMI_2[which(is.na(imp.data.new$BMI_2))]<-imp.data.mice$BMI_2[[j]]
  imp.data.new$BMI_3[which(is.na(imp.data.new$BMI_3))]<-imp.data.mice$BMI_3[[j]]
  imp.data.new$BMI_4[which(is.na(imp.data.new$BMI_4))]<-imp.data.mice$BMI_4[[j]]
  
  MI.imputed.data[[j]]<-imp.data.new
  
}

## Put multiple completed datasets into a list
imp.data.all<-list()

for(j in 1:5){
  
  imp.data.j<-MI.imputed.data[[j]]
  ## rounding to integer
    for(g in c("num.paydev_2","num.disease_2","num.paydev_3","num.disease_3","num.paydev_4","num.disease_4")){
      
      imp.data.j[,g]<-ifelse(imp.data.j[,g]<0,0,round(imp.data.j[,g]))
      
  }
  
  imp.data.all[[j]]<-imp.data.j
}

# save imputed datasets
save(imp.data.all,file = "Imputed_data_FCS_standard.Rdata")


## Method 2-FCS-Standard (PMM): Follow the codes for FCS-Standard, and remove the code: methods[which(methods=="pmm")]<-"norm"
## Imputation models for each incomplete variable:
## BMI: linear regression
## Hospital stay: logistic regression
## Assistive devices: predictive mean matching
## Comorbidity index: predictive mean matching


## Method 3- FCS-Standard (Poisson) 
## Imputation models for each incomplete variable:
## BMI: linear regression
## Hospital stay: logistic regression
## Assistive devices: Poisson regression
## Comorbidity index: Poisson regression

## additional package required:
# To install package
devtools::install_git(url = "https://github.com/kkleinke/countimp", branch = "master")
# load package
library(countimp)

## To run this method: Follow the codes for FCS-Standard, and change the code: methods[which(methods=="pmm")]<-"pmm" to the following:

# methods[c("BMI_2","BMI_3","BMI_4")]<-"norm"
# methods[c("num.disease_2","num.disease_3","num.disease_4")]<-"pois"
# methods[c("num.paydev_2","num.paydev_3","num.paydev_4")]<-"pois"


## Method 4: FCS-LMM-Latent: impute data in its long-format
## Imputation models for each incomplete variable:
## BMI: multilevel linear regression
## Hospital stay: multilevel regression on latent variables
## Assistive devices: multilevel linear regression+rounding
## Comorbidity index: multilevel linear regression+rounding

## Additional packages required:
library(micemd) #install first if not installed
library(reshape2)
library(tidyr)
library(dplyr)

## convert to long-format ##
## hospital stay
sim.long.Y1<-gather(sim.data[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(sim.data)[2:5],factor_key = T)
## num disease ##
sim.long.Y2<-gather(sim.data[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(sim.data)[6:9],factor_key = T)
## bmi ##
sim.long.Y3<-gather(sim.data[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(sim.data)[10:13],factor_key = T)
## num.paydev ##
sim.long.Y4<-gather(sim.data[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(sim.data)[14:17],factor_key = T)

## combine together ##
sim.long<-cbind(sim.long.Y1,num.disease=sim.long.Y2$num.disease,BMI=sim.long.Y3$BMI,num.paydev=sim.long.Y4$num.paydev)
sim.long$Round<-as.numeric(sim.long$Round)
sim.long<-arrange(sim.long,spid)
sim.long$hospital.stay<-as.factor(sim.long$hospital.stay)
  

### Imputation step ###
imp.mice.2l.methods<-find.defaultMethod(sim.long,ind.clust = 1) #find default imputation models

#set imputation models for each incomplete variable 
imp.mice.2l.methods[c("hospital.stay")]<-"2l.jomo" #linear mixed model with latent variables
imp.mice.2l.methods["num.disease"]<-"2l.glm.norm"  #linear mixed model
imp.mice.2l.methods["BMI"]<-"2l.glm.norm"          #linear mixed model
imp.mice.2l.methods["num.paydev"]<-"2l.glm.norm"   #linear mixed model

imp.2l.0<-mice(sim.long,method=imp.mice.2l.methods,maxit=0) #initiate imputation, 

#customize predictor matrix
pred.mat<-imp.2l.0$predictorMatrix
pred.mat[c("hospital.stay","num.disease","BMI","num.paydev"),1]<--2 #set cluster indicator: "spid" to -2 
pred.mat[c("hospital.stay","num.disease","BMI","num.paydev"),-1]<-1  #set fixed effects

pred.mat[c("hospital.stay"),c("hospital.stay")]<-0 # set diagonals as zeros
pred.mat[c("num.paydev"),c("num.paydev")]<-0
pred.mat[c("num.disease"),c("num.disease")]<-0
pred.mat[c("BMI"),c("BMI")]<-0

sim.long$spid<-as.numeric(as.factor(sim.long$spid))

#imputation step
imp.step<-mice(sim.long,method=imp.mice.2l.methods,pred=pred.mat,print=T)
# imputed data
imp.data.mice<-imp.step$imp

# form multiple imputed datasets into a list
imp.data.all<-list()

#find missing value positions and replace them by imputed values
id.na<-which(is.na(sim.long$hospital.stay))
for(j in 1:5){
  
  imp.data.j<-sim.long
  imp.data.j$hospital.stay[id.na]<-imp.data.temp$hospital.stay[[j]]
  imp.data.j$num.paydev[id.na]<-ifelse(imp.data.temp$num.paydev[[j]]<0,0,round(imp.data.temp$num.paydev[[j]]))
  imp.data.j$num.disease[id.na]<-ifelse(imp.data.temp$num.disease[[j]]<0,0,round(imp.data.temp$num.disease[[j]]))
  imp.data.j$BMI[id.na]<-imp.data.temp$BMI[[j]]
  imp.data.all[[j]]<-imp.data.j
}

#save imputed datasets
save(imp.data.all,file = "Imputed_data_FCS_LMM_latent.Rdata")

## Method 5: FCS-GLMM (Gaussian): following the codes for FCS-LMM-latent; only need to change the models for binary variables and count variables
## Imputation models for each incomplete variable:
## BMI: multilevel linear regression
## Hospital stay: multilevel logistic regression 
## Assistive devices: multilevel linear regression
## Comorbidity index: multilevel linear regression

# Follow the above FCS-LMM-latent code, and change the code:imp.mice.2l.methods[c("hospital.stay")]<-"2l.jomo" to:
# imp.mice.2l.methods[c("hospital.stay")]<-"2l.glm.bin"


## Method 6: FCS-GLMM (Poisson): following the codes for FCS-LMM-latent; only need to change the models for binary variables and count variables
## Imputation models for each incomplete variable:
## BMI: multilevel linear regression
## Hospital stay: multilevel logistic regression 
## Assistive devices: multilevel Poisson regression
## Comorbidity index: multilevel Poisson regression

# Follow the above FCS-LMM-latent code, and change model specification to:
# imp.mice.2l.methods[c("hospital.stay")]<-"2l.glm.bin"
# imp.mice.2l.methods["num.disease"]<-"2l.glm.pois"
# imp.mice.2l.methods["BMI"]<-"2l.glm.norm"
# imp.mice.2l.methods["num.paydev"]<-"2l.glm.pois"

################################################################
###   Methods using the joint modeling strategy              ###
################################################################

## Method 7. JM-General location model (impute with wide format data)
library(mix) #install first if not installed

# Prepare datasets; set binary variables as factors 
id.cont<-c("num.paydev_1","num.paydev_2","num.paydev_3","num.paydev_4",
           "num.disease_1","num.disease_2","num.disease_3","num.disease_4",
           "BMI_1", "BMI_2", "BMI_3", "BMI_4","hh_dhshldnum")
id.factor<-which(!colnames(sim.data)%in%id.cont)
for(l in id.factor){
  
  sim.data[,l]<-as.factor(sim.data[,l])
  
}

# prepare data sets
sim.data<-sim.data[,c("spid", "hc_hosptstay_1","hc_hosptstay_2","hc_hosptstay_3","hc_hosptstay_4",
                      "num.disease_1", "num.disease_2","num.disease_3","num.disease_4",
                      "BMI_1","BMI_2","BMI_3","BMI_4",
                      "num.paydev_1","num.paydev_2","num.paydev_3","num.paydev_4",
                      X.var)]
# Code binary variables ("1","0")-->(1,2) 
id.factor<-c(paste("hc_hosptstay",1:4,sep="_"))
temp.hospital<-apply(sim.data[,id.factor],2,function(x)ifelse(x=="1",1,2))

# Code complete binary covarites to dummy variables,treat them as continuous variables; 
id.dummy<-c("r1dgender","rl1dracehisp","hh_martlstat","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
temp.predictor<-apply(sim.data[,id.dummy],2,function(x)ifelse(x=="1",1,0))
# reorder variables: incomplete factors takes the first p columns, then continuous variables
mat<-as.matrix(cbind(temp.hospital,sim.data[,c(paste("num.disease",1:4,sep="_"))],sim.data[,c(paste("BMI",1:4,sep="_"))],sim.data[,c(paste("num.paydev",1:4,sep="_"))],
                     sim.data[,c("hh_dhshldnum","Age","hc_health")],
                     temp.predictor))
#preliminary step
pre.data<-prelim.mix(mat,p=4)
theta.hat<-em.mix(pre.data)
rngseed(K)
new.theta<-da.mix(pre.data,theta.hat,steps=1000,showits = F)

#multiple imputation step
imp.data.all<-list()

for(j in 1:5){
  imp.data<-imp.mix(pre.data,new.theta,x=mat)
  
  #format the output data back to a data frame 
  imp.data.j<-imp.data[[j]]
  imp.data.j<-cbind(spid=true.data$spid,imp.data.j)
  imp.data.j<-apply(imp.data.j,2,function(x)as.numeric(as.character(x)))
  imp.data.j<-as.data.frame(imp.data.j)
  
  #code the binary variable from (1,2) back to (1,0)
  for(l in c("hc_hosptstay_1","hc_hosptstay_2","hc_hosptstay_3","hc_hosptstay_4")){
    imp.data.j[,l]<-ifelse(imp.data.j[,l]==1,1,0)
  }
  # round the count variables to integers
  for(g in c("num.paydev_2","num.disease_2","num.paydev_3","num.disease_3","num.paydev_4","num.disease_4")){
    imp.data.j[,g]<-ifelse(imp.data.j[,g]<0,0,round(imp.data.j[,g]))
  }
  
  imp.data.all[[j]]<-imp.data.j
  
  
}

# save multiple imputed data
save(imp.data.all,file = "Imputed_data_JM_GL.Rdata")


## Method 8,9. JM-MLMM_latent (common) or JM-MLMM-latent (random)  (impute with long format data)
#loading additional required package
library(jomo)
library(tidyr)
library(dplyr)

##  prepare long-format data-set

# set binary variable to factor
id.cont<-c("num.paydev_1","num.paydev_2","num.paydev_3","num.paydev_4",
             "num.disease_1","num.disease_2","num.disease_3","num.disease_4",
             "BMI_1", "BMI_2", "BMI_3", "BMI_4","hh_dhshldnum","Age","hc_health")
id.factor<-which(!colnames(sim.data)%in%id.cont)
sim.data[,l]<-as.factor(sim.data[,l])

# convert data to long-format
## hospital stay
sim.long.Y1<-gather(sim.data[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(sim.data)[2:5],factor_key = T)
## num disease ##
sim.long.Y2<-gather(sim.data[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(sim.data)[6:9],factor_key = T)
## bmi ##
sim.long.Y3<-gather(sim.data[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(sim.data)[10:13],factor_key = T)
## num.paydev ##
sim.long.Y4<-gather(sim.data[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(sim.data)[14:17],factor_key = T)

## combine together ##
sim.long<-cbind(sim.long.Y1,num.disease=sim.long.Y2$num.disease,BMI=sim.long.Y3$BMI,num.paydev=sim.long.Y4$num.paydev)
sim.long$Round<-as.numeric(sim.long$Round)
sim.long<-arrange(sim.long,spid)
sim.long$hospital.stay<-as.factor(sim.long$hospital.stay)


## Imputation step ##
Y<-sim.long[,c("hospital.stay","num.disease","BMI","num.paydev")] #get incomplete variables
clus<-sim.long[,c("spid")]                                        #cluster indicator 
X<-sim.long[,c(X.var,"Round")]                                    # complete variables

## Imputation step for JM-MLMM-latent with a common across subjects variance
imp.step<-jomo(Y,clus=clus,X=X,meth="common")

## Imputation step for JM-MLMM-latent with a random across subjects variance
# imp.jomo.1<-jomo(Y,clus=clus,X=X,meth="random")

## Get multiple imputed datasets
imp.data.all<-list()

for(j in 1:5){

  imp.data.j<-imp.step[which(imp.step$Imputation==j),!colnames(imp.step)%in%c("Z1","id","Imputation")]
  colnames(imp.data.j)[which(colnames(imp.data.j)=="clus")]<-"spid"

  imp.long<-imp.data.j
  
  # rounding count variables to integers
  for(g in c("num.paydev","num.disease")){
    imp.long[,g]<-ifelse(imp.long[,g]<0,0,round(imp.long[,g]))
  }
  imp.data.all[[j]]<-imp.long
}

save(imp.data.all,file = "Imputed_data_JM_MLMM_common.Rdata")


#####################################################################################
##                      Step 3: Statistical analysis                               ##
#####################################################################################

### An univariate multilevel logistic model ### 
# Hospital.stay ~ BMI+ Comorbidity+ Assistive devices + X.var + random.intercept

#####################################################################################
######      Study3: Statistical analysis  (univariate GLMM)                       ###
#####################################################################################

library(lme4)
library(tidyr)
library(dplyr)
library(glmmsr)

method<-"FCS_Standard" #(specify an imutation method)
## get statistical results using multplie imputed data 
temp.fixed.coef<-temp.fixed.sd<-temp.rand.coef<-temp.rand.sd<- temp.rand.corr<- temp.rand.corr.sd<-NULL
temp.fixed.coef.mode<-temp.rand.coef.mode<-temp.rand.corr.mode<-NULL

for(j in 1:5){
  
  imp.data.j<-imp.data.all[[j]]
  
  ## Prepare dataset ##
  if(method%in%c("FCS-Standard","FCS-Standard(Guassian)","FCS-Standard(Poisson)")){
    
    
    id.bin<-c("r1dgender","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]=="1","1","0")  # recode binary covariates (1,2)->(1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    imp.data.j$Age<-as.numeric(as.character(imp.data.j$Age))
    imp.data.j$hc_health<-as.numeric(as.character(imp.data.j$hc_health))
    
    ## convert data to long format ##
    ## hospital stay
    imp.long.Y1<-gather(imp.data.j[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(imp.data.j)[2:5],factor_key = T)
    ## num disease ##
    imp.long.Y2<-gather(imp.data.j[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(imp.data.j)[6:9],factor_key = T)
    ## bmi ##
    imp.long.Y3<-gather(imp.data.j[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(imp.data.j)[10:13],factor_key = T)
    ## num.paydev ##
    imp.long.Y4<-gather(imp.data.j[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(imp.data.j)[14:17],factor_key = T)
    
    ## combine together ##
    imp.long<-cbind(imp.long.Y1,num.disease=imp.long.Y2$num.disease,BMI=imp.long.Y3$BMI,num.paydev=imp.long.Y4$num.paydev)
    imp.long$Round<-as.numeric(imp.long$Round)
    imp.long<-arrange(imp.long,spid)
    imp.long$hospital.stay<-as.factor(imp.long$hospital.stay)
    
  }
  
  if(method %in% c("FCS-LMM-latent","FCS-GLMM-Gauss","FCS-GLMM-Poisson")){
    
    id.bin<-c("r1dgender","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]=="1","1","0")
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    imp.data.j$Age<-as.numeric(as.character(imp.data.j$Age))
    imp.data.j$hc_health<-as.numeric(as.character(imp.data.j$hc_health))
    
    imp.long<-imp.data.j
  }
  
  if(method=="JM_GL"){
    
    imp.data.j$spid<-as.factor(imp.data.j$spid)
    id.bin<-c(paste0("hc_hosptstay","_",1:4),"r1dgender","rl1dracehisp","hh_martlstat","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]==1,1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    ## convert to long format ##
    ## hospital stay
    imp.long.Y1<-gather(imp.data.j[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(imp.data.j)[2:5],factor_key = T)
    ## num disease ##
    imp.long.Y2<-gather(imp.data.j[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(imp.data.j)[6:9],factor_key = T)
    ## bmi ##
    imp.long.Y3<-gather(imp.data.j[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(imp.data.j)[10:13],factor_key = T)
    ## num.paydev ##
    imp.long.Y4<-gather(imp.data.j[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(imp.data.j)[14:17],factor_key = T)
    
    ## combine together ##
    imp.long<-cbind(imp.long.Y1,num.disease=imp.long.Y2$num.disease,BMI=imp.long.Y3$BMI,num.paydev=imp.long.Y4$num.paydev)
    imp.long$Round<-as.numeric(imp.long$Round)
    imp.long<-arrange(imp.long,spid)
    imp.long$hospital.stay<-as.factor(imp.long$hospital.stay)
    
    imp.long<-imp.data.j
  }
  
  if(method%in% c("JM_MLMM_common","JM_MLMM_random")){
    
    id.bin<-c("r1dgender","rl1dracehisp","hh_martlstat","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]==1,1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    
    imp.long<-imp.data.j
    
  }

  
  ####### Fit univariate model #######
  fit.glmm<-glmm(hospital.stay~num.disease+Age+r1dgender+hc_health+ew_medpaovtm+
                   mc_meds+num.paydev+Round+BMI+num.paydev+(1|spid),
                 data=imp.long, family = binomial,method="AGQ")
  
  
  ## Collect fixed effect
  beta.fix<-fit.glmm$estim[-1]
  sigma.rand<-fit.glmm$estim[1]
  
  res<-summary(fit.glmm)
  sd.beta.fix<-res$se[-1]
  sd.sigma.rand<-res$se[1]
  
  ## Combine results ##
  temp.fixed.coef<-rbind(temp.fixed.coef,beta.fix)
  temp.fixed.sd<-rbind(temp.fixed.sd,sd.beta.fix)
  
  temp.rand.coef<-rbind(temp.rand.coef,sigma.rand)
  temp.rand.sd<-rbind(temp.rand.sd,sd.sigma.rand)
  
  print(paste("j",j,sep="="))
}

##### Pool (Rubin's rule) #####
est.fixed.coef<-apply(temp.fixed.coef,2,mean)
est.fixed.between.var<-apply(temp.fixed.coef,2,var)
est.fixed.within.var<-apply(temp.fixed.sd,2,function(x)mean(x^2))
est.fixed.sd<-sqrt(est.fixed.within.var+(1+1/5)*est.fixed.between.var)
est.fixed.CI.lower<-apply(cbind(est.fixed.coef,est.fixed.sd),1,function(x) x[1]+qnorm(0.025)*x[2])
est.fixed.CI.upper<-apply(cbind(est.fixed.coef,est.fixed.sd),1,function(x) x[1]+qnorm(0.975)*x[2])

est.rand.coef<-apply(temp.rand.coef,2,mean)
est.rand.between.var<-apply(temp.rand.coef,2,var)
est.rand.within.var<-apply(temp.rand.sd,2,function(x)mean(x^2))
est.rand.sd<-sqrt(est.rand.within.var+(1+1/5)*est.rand.between.var)
est.rand.CI.lower<-apply(cbind(est.rand.coef,est.rand.sd),1,function(x) x[1]+qnorm(0.025)*x[2])
est.rand.CI.upper<-apply(cbind(est.rand.coef,est.rand.sd),1,function(x) x[1]+qnorm(0.975)*x[2])


### An bivariate multilevel model ### 
# Hospital.stay ~ BMI+ Comorbidity + X.var +random.intercept1
# Assistive devices ~ BMI+ Comorbidity + X.var +random.intercept2
# (random1,random2)~N(0, V_b)

# Loading additional package
library(MCMCglmm)
    
temp.fixed.coef<-temp.fixed.sd<-temp.rand.coef<-temp.rand.sd<- temp.rand.corr<- temp.rand.corr.sd<-NULL
temp.fixed.coef.mode<-temp.rand.coef.mode<-temp.rand.corr.mode<-NULL

for(j in 1:5){
  
  imp.data.j<-imp.data.all[[j]]
  
  ## Prepare dataset ##
  if(method%in%c("FCS-Standard","FCS-Standard(Guassian)","FCS-Standard(Poisson)")){
    
    
    id.bin<-c("r1dgender","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]=="1","1","0")  # recode binary covariates (1,2)->(1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    imp.data.j$Age<-as.numeric(as.character(imp.data.j$Age))
    imp.data.j$hc_health<-as.numeric(as.character(imp.data.j$hc_health))
    
    ## convert data to long format ##
    ## hospital stay
    imp.long.Y1<-gather(imp.data.j[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(imp.data.j)[2:5],factor_key = T)
    ## num disease ##
    imp.long.Y2<-gather(imp.data.j[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(imp.data.j)[6:9],factor_key = T)
    ## bmi ##
    imp.long.Y3<-gather(imp.data.j[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(imp.data.j)[10:13],factor_key = T)
    ## num.paydev ##
    imp.long.Y4<-gather(imp.data.j[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(imp.data.j)[14:17],factor_key = T)
    
    ## combine together ##
    imp.long<-cbind(imp.long.Y1,num.disease=imp.long.Y2$num.disease,BMI=imp.long.Y3$BMI,num.paydev=imp.long.Y4$num.paydev)
    imp.long$Round<-as.numeric(imp.long$Round)
    imp.long<-arrange(imp.long,spid)
    imp.long$hospital.stay<-as.factor(imp.long$hospital.stay)
    
  }
  
  if(method %in% c("FCS-LMM-latent","FCS-GLMM-Gauss","FCS-GLMM-Poisson")){
    
    id.bin<-c("r1dgender","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]=="1","1","0")
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    imp.data.j$Age<-as.numeric(as.character(imp.data.j$Age))
    imp.data.j$hc_health<-as.numeric(as.character(imp.data.j$hc_health))
    
    imp.long<-imp.data.j
  }
  
  if(method=="JM_GL"){
    
    imp.data.j$spid<-as.factor(imp.data.j$spid)
    id.bin<-c(paste0("hc_hosptstay","_",1:4),"r1dgender","rl1dracehisp","hh_martlstat","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]==1,1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    ## convert to long format ##
    ## hospital stay
    imp.long.Y1<-gather(imp.data.j[,c("spid",paste("hc_hosptstay",1:4,sep="_"),X.var)],Round,"hospital.stay",colnames(imp.data.j)[2:5],factor_key = T)
    ## num disease ##
    imp.long.Y2<-gather(imp.data.j[,c("spid",paste("num.disease",1:4,sep="_"))],Round,"num.disease",colnames(imp.data.j)[6:9],factor_key = T)
    ## bmi ##
    imp.long.Y3<-gather(imp.data.j[,c("spid",paste("BMI",1:4,sep="_"))],Round,"BMI",colnames(imp.data.j)[10:13],factor_key = T)
    ## num.paydev ##
    imp.long.Y4<-gather(imp.data.j[,c("spid",paste("num.paydev",1:4,sep="_"))],Round,"num.paydev",colnames(imp.data.j)[14:17],factor_key = T)
    
    ## combine together ##
    imp.long<-cbind(imp.long.Y1,num.disease=imp.long.Y2$num.disease,BMI=imp.long.Y3$BMI,num.paydev=imp.long.Y4$num.paydev)
    imp.long$Round<-as.numeric(imp.long$Round)
    imp.long<-arrange(imp.long,spid)
    imp.long$hospital.stay<-as.factor(imp.long$hospital.stay)
    
    imp.long<-imp.data.j
  }
  
  if(method%in% c("JM_MLMM_common","JM_MLMM_random")){
    
    id.bin<-c("r1dgender","rl1dracehisp","hh_martlstat","ss_painbothr","ip_covmedcad","ew_medpaovtm","mc_meds")
    for(v in id.bin){
      imp.data.j[,v]<-ifelse(imp.data.j[,v]==1,1,0)
      imp.data.j[,v]<-as.factor(imp.data.j[,v])
    }
    
    
    imp.long<-imp.data.j
    
  }
  
  ####### Fit bivariate model #######
  fit.glmm<-MCMCglmm(cbind(hospital.stay, num.paydev)~trait+trait:(num.disease+Age+r1dgender+hc_health+ew_medpaovtm+
                                                                         mc_meds+Round+BMI)-1, random = ~us(trait):spid, rcov=~idh(trait):units,nitt=30000,thin=10,burnin=10000,
                         data=imp.long, family=c("categorical","poisson"))
  sum.res<-summary(fit.glmm)
  
   
  ## Collect fixed effect
  beta.fix<-sum.res$solutions[,1]
  beta.fix.mode<-posterior.mode(fit.glmm$Sol)
  
  coef.sigma.rand<-sum.res$Gcovariances[,1]
  coef.sigma.rand.mode<-posterior.mode(fit.glmm$VCV[,1:4])
  
  corr.rand<-mean(posterior.cor(fit.glmm$VCV[,1:4])[,2])
  corr.rand.mode<-posterior.mode(posterior.cor(fit.glmm$VCV[,1:4]))[2]
  
  sd.beta.fix<-apply(fit.glmm$Sol,2,sd)
  sd.sigma.rand<-apply(fit.glmm$VCV[,1:4],2,sd)
  sd.corr.rand<-sd(posterior.cor(fit.glmm$VCV[,1:4])[,2])
  
  ## Combine results ##
  temp.fixed.coef<-rbind(temp.fixed.coef,beta.fix)
  temp.fixed.sd<-rbind(temp.fixed.sd,sd.beta.fix)
  
  temp.rand.coef<-rbind(temp.rand.coef,sigma.rand)
  temp.rand.sd<-rbind(temp.rand.sd,sd.sigma.rand)
  
  temp.rand.corr<-rbind(temp.rand.corr,corr.rand)
  temp.rand.corr.sd<-rbind(temp.rand.corr.sd,sd.corr.rand)
  
  temp.fixed.coef.mode<-rbind(temp.fixed.coef.mode,beta.fix.mode)
  temp.rand.coef.mode<-rbind(temp.rand.coef.mode,coef.sigma.rand.mode)
  
  temp.rand.corr.mode<-rbind(temp.rand.corr.mode,corr.rand.mode)
  
  print(paste("j",j,sep="="))

}

##### Pool (Rubin's rule) #####
est.fixed.coef<-apply(temp.fixed.coef,2,mean)
est.fixed.between.var<-apply(temp.fixed.coef,2,var)
est.fixed.within.var<-apply(temp.fixed.sd,2,function(x)mean(x^2))
est.fixed.sd<-sqrt(est.fixed.within.var+(1+1/5)*est.fixed.between.var)
est.fixed.CI.lower<-apply(cbind(est.fixed.coef,est.fixed.sd),1,function(x) x[1]+qnorm(0.025)*x[2])
est.fixed.CI.upper<-apply(cbind(est.fixed.coef,est.fixed.sd),1,function(x) x[1]+qnorm(0.975)*x[2])

est.rand.coef<-apply(temp.rand.coef,2,mean)
est.rand.between.var<-apply(temp.rand.coef,2,var)
est.rand.within.var<-apply(temp.rand.sd,2,function(x)mean(x^2))
est.rand.sd<-sqrt(est.rand.within.var+(1+1/5)*est.rand.between.var)
est.rand.CI.lower<-apply(cbind(est.rand.coef,est.rand.sd),1,function(x) x[1]+qnorm(0.025)*x[2])
est.rand.CI.upper<-apply(cbind(est.rand.coef,est.rand.sd),1,function(x) x[1]+qnorm(0.975)*x[2])

est.rand.corr<-apply(temp.rand.corr,2,mean)
est.rand.corr.between.var<-apply(temp.rand.corr,2,var)
est.rand.corr.within.var<-apply(temp.rand.corr.sd,2,function(x)mean(x^2))
est.rand.corr.sd<-sqrt(est.rand.corr.within.var+(1+1/5)*est.rand.corr.between.var)
est.rand.corr.CI.lower<-apply(cbind(est.rand.corr,est.rand.corr.sd),1,function(x) x[1]+qnorm(0.025)*x[2])
est.rand.corr.CI.upper<-apply(cbind(est.rand.corr,est.rand.corr.sd),1,function(x) x[1]+qnorm(0.975)*x[2])

est.fixed.coef.mode<-apply(temp.fixed.coef.mode,2,mean)
est.rand.coef.mode<-apply(temp.rand.coef.mode,2,mean)
est.rand.corr.mode<-apply(temp.rand.corr.mode,2,mean)
  




