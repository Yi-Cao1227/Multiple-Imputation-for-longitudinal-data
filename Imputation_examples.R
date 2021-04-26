
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


####################### Methods using wide-format data  ######################

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



