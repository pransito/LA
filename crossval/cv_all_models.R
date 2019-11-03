## cv of all candidate (120) models (or in case of the sample of 51 we use less)

# in case of sample of 51
setwd(pfad_scripts)
source("all_models_51_no_ed.R")

# evaluate all models in parallel
cl<-makeCluster(6) #change to your number of CPU cores  
registerDoSNOW(cl)
modc <- foreach(ff=1:length(all_modc), .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
  eval(parse(text=all_modc[[ff]]))
} 
stopCluster(cl)
setwd(pfad_results)
save(file="modc51_noed_majrev.RData",list = c("modc"))

setwd(pfad_results)
load("modc51_noed_majrev.RData")

## cv of modc
flds <- agk.getflds(10,"subject",modc[[1]])

# clusters: number of CPU cores
# CV of modc not based on subjects but on subjects (only using group fixed effects)
cl<-makeCluster(5)  
registerDoSNOW(cl)
cvsetmodc <- foreach(ff=1:length(modc), .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
  agk.setcv.glmer(modc[[ff]],"subject",flds)
} 
stopCluster(cl)
save(file="cvsetmodc_3.RData",list=c("cvsetmodc"))

## CV of modc not based on subjects but on trials ("observations") 
flds <- agk.getflds(10,byvar = NULL,modc[[1]])

# clusters: number of CPU cores 
cl<-makeCluster(5)  
registerDoSNOW(cl)
cvsetmodc_obs <- foreach(ff=1:length(modc), .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
  agk.setcv.glmer(modc[[ff]],"observations",flds)
} 
stopCluster(cl)
save(file="cvsetmodc2_obs.RData",list=c("cvsetmodc_obs"))

## plot the cv mean and 95% CI
setwd(pfad_results)
load("modc.RData")
load("cvsetmodc_3.RData")
cvsetmodc <- cvsetmodc

modc        <- modc
cvsetmodc   <- cvsetmodc
cur_pattern <- "ed.abs"
modc_pruned <- list()
cvsetmodc_p <- list()
all_mod_p   <- list()
count       <- 0
for (ii in 1:length(modc)) {
  cur_str <- as.character(modc[[ii]]@call)
  if (length(grep(cur_pattern,cur_str)) == 0) {
    next
  } else {
    count <- count + 1
    modc_pruned[count] <- modc[ii]
    cvsetmodc_p[count] <- cvsetmodc[ii]
    all_mod_p[[count]] <- all_modc[[ii]]
  }
}

modc      <- modc_pruned
cvsetmodc <- cvsetmodc_p

cvci <- matrix(NA,nrow=length(cvsetmodc),ncol=5)
for (kk in 1:length(cvsetmodc)) {
  b1     <- one.boot(cvsetmodc[[kk]]$pred_score, mean, R=2000, tr=0)
  cur_ci <- boot.ci(b1, type=c("perc"))
  cvci[kk,1] <- b1$t0
  cvci[kk,2] <- cur_ci$percent[[4]]
  cvci[kk,3] <- cur_ci$percent[[5]]
  cvci[kk,4] <- AIC(modc[[kk]])
  cvci[kk,5] <- BIC(modc[[kk]])
}
cvci <- as.data.frame(cvci)
names(cvci) <- c("cvmean","cvcilo", "cvcihi","AIC","BIC")
cvci$model <- as.factor(1:length(modc))

# which models to consider?
chosen_formulas <- c()
for (ii in 1:length(modc)) {
  cur_form <- as.character(formula(modc[[ii]]@call))
  chosen_formulas[ii] <- cur_form[3]
}

model_family <- c()
pos_families <- c("RiskMar","ratio","diff","EV","Gewinn")
for (ii in 1:length(chosen_formulas)) {
  for (jj in 1:length(pos_families)) {
    if (length(grep(chosen_formulas[[ii]],pattern = pos_families[jj]))>0) {model_family[ii] <- jj}
  }
}

cvci$model_family <- model_family
cvci$model_family <- as.factor(cvci$model_family)
levels(cvci$model_family) <- pos_families

# plot
p2 <- ggplot(data = cvci, aes(x=model,y=cvmean,fill=model_family))
p2 <- p2+geom_bar(stat="identity")
p2 <- p2 + geom_errorbar(aes(ymin=cvcilo,ymax=cvcihi), size=.2, width=0.5,col="black",fill=NULL) + ylab("mean (95% boots. CI)")
p2 <- p2 + ggtitle("EG estimates folded by subjects of models with CIs\n")
p2 

## only using AIC
setwd(pfad_results)
load("modc.RData")
load("modcc.RData")
modc <- c(modc,modcc)

all_AIC <- lapply(modc,FUN = AIC)
all_BIC <- lapply(modc,FUN = BIC)
modc[[which(min(unlist(all_AIC))==all_AIC)]]
modc[[which(min(unlist(all_BIC))==all_BIC)]]

boot_la <-bootMer(modc[[which(min(unlist(all_AIC))==all_AIC)]],use.u = F,
                  FUN = get_la_fixef,
                  nsim = 50,
                  type = "parametric",
                  verbose = TRUE,
                  parallel = "snow",
                  ncpus = 6)

## bootstrap p-values of group effects on LA
# use external script
