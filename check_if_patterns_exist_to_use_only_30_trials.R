## check whether a sample of 30 trials is enough to reliably estimate the lambda parameter

## preparation
# clear workspace
rm(list = ls())

# path
pfad <- "C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"
#pfad <- "C:\\Users\\Alexander\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"

setwd(pfad)

## load libraries and functions
setwd(paste(pfad, "\\Scripts", sep=""))
source ('LA_functions.R')

## set the pwd needed
setwd(paste(pfad, "\\Data", sep=""))

## options
# currently centered-only
use.z    <- 1 
# use absolute value of loss as predictor
loss.abs <- 1
# calculate lmlist?
do_lmlist <- 1
# use only n_trials trials?
use_sample <- 1
n_trials   <- 35

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
setwd(paste(pfad, "\\Scripts", sep=""))
source("get_data_la.R")
# set the pwd needed
setwd(paste(pfad, "\\Data", sep=""))

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
r.lambdas <- c()
p.lambdas <- c()
cur.r     <- 0
count <- 0
while (cur.r < 0.8) {
  if (count > 200) {break}
  count <- count + 1
  print(paste("now running round",count))
  
  tmp <- make.gamble.sample("subject",data.la,40,verbosity = 0)
  tmp[[2]]
  cur.data.la <- tmp[[1]]
  
  ## calculate lmlist (easier than mixed model)
  if (do_lmlist == 1){
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
  cur.cor <- cor.test(lambda_ed_bl,crm$lambda_30ed,method="spearman")
  r.lambdas[count] <- as.numeric(cur.cor$estimate)
  p.lambdas[count] <- as.numeric(cur.cor$p.value)
  cur.r <- as.numeric(cur.cor$estimate)
}


# darstellung der gamble matrix

# assessment
hist(corr_with_orig_lambda)
p_of_corr_with_orig_lam < 0.01
