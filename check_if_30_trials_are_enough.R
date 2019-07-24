## check whether a sample of 30 trials is enough to reliably estimate the lambda parameter

## preparation
# clear workspace
rm(list=ls())
root_wd  = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(root_wd)

## load libraries and functions
setwd('..')
setwd('R')
source('agk_library.R')

## set the pwd needed
path_data = 'C:/Users/genaucka/Google Drive/09_Diplom/LA/daten_behav_test_finale_SP_Diplom/Data'

## options
# currently centered-only
use.z    <- 1 
# use absolute value of loss as predictor
loss.abs <- 1
# calculate lmlist?
do_lmlist <- 1
# use only 30 trials?
use_sample <- 1

# if there are more than x% percent missings in response then drop-out
missing.cutoff <- 0.07

## prep the final matrix
subject <- c()
p       <- c()
r       <- c()

## prep some variables
missing.check  <- 1 # if set to 1 then people with too many missings will be dropped 
missing.people <- list()

## get all the data in long format
setwd(root_wd)
source("get_data_la.R")
# set the pwd needed
setwd(path_data)

## calc lmlist lambda
## calculate lmlist (easier than mixed model)
if (do_lmlist == 1){
  lmlist_ed         <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = data.la, na.action = na.omit, family = "binomial")
}

## calculate Lambda
if(do_lmlist == 1){
  crm_lmlist <- c()
  for (ii in 1:length(lmlist_ed)) {
    crm_lmlist <- rbind(crm_lmlist,t(as.matrix(as.numeric(lmlist_ed[[ii]][[1]]))))
  }
  crm_lmlist <- as.data.frame(crm_lmlist)
  names(crm_lmlist) <- c("intercept", "ed.abs", "Gewinn", "Verlust")
  crm_lmlist$id       <- as.factor(as.numeric(as.matrix(names(lmlist_ed))))
}

crm <- crm_lmlist

lambda_ed_bl <- calc_lambda(crm$Gewinn,crm$Verlust)

## get only a subsample of 30 trials per subject
# run this a couple of times to see how the sampling leads to shakes in estimation of lambda
small.lambdas <- c()

for (jj in 1:100) {
  print(paste("now running round",jj))
  cur.data.la <- make_sample_of_trials("subject",data.la,45,verbosity = 1, same=FALSE,check=FALSE)
  ## calculate lmlist (easier than mixed model)
  if (do_lmlist == 1) {
    lmlist_ed         <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = cur.data.la, na.action = na.omit, family = "binomial")
  }
  
  ## calculate Lambda
  if(do_lmlist == 1){
    crm_lmlist <- c()
    for (ii in 1:length(lmlist_ed)) {
      crm_lmlist <- rbind(crm_lmlist,t(as.matrix(as.numeric(lmlist_ed[[ii]][[1]]))))
    }
    crm_lmlist <- as.data.frame(crm_lmlist)
    names(crm_lmlist) <- c("intercept", "ed.abs", "Gewinn", "Verlust")
    crm_lmlist$id       <- as.factor(as.numeric(as.matrix(names(lmlist_ed))))
  }
  
  crm <- crm_lmlist
  
  crm$lambda_30ed <- calc_lambda(crm$Gewinn,crm$Verlust)
  small.lambdas <- rbind(small.lambdas,crm$lambda_30ed)
}


## now compare to original Lambda (from lmlist)
small.lambdas <- t(small.lambdas)

# corr
corr_with_orig_lambda <- c()
p_of_corr_with_orig_lam <- c()
for (kk in 1: length(small.lambdas[1,])){
  cur.cor <- cor.test(lambda_ed_bl,small.lambdas[,kk],method="spearman")
  corr_with_orig_lambda <- rbind(corr_with_orig_lambda,as.numeric(cur.cor$estimate))
  p_of_corr_with_orig_lam <- rbind(p_of_corr_with_orig_lam,as.numeric(cur.cor$p.value))
}

# assessment
hist(corr_with_orig_lambda)
print(p_of_corr_with_orig_lam < 0.01)
