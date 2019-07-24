# function to compute a p-value for a fixed effect of interest
# in a glmer model
# mod_bl ist das baseline model gefittet
# mod ist ein string der das fitting des complete models angbit
# fun_extract ist die funktion, die fiexedeffect und contrasts of interests extrahiert
# mod_compl_fit ist das komplette model fit on the original data

# in case of sample of 51
setwd(pfad_scripts)
source("all_models_51_no_ed.R")
setwd(pfad_results)
load("modc51_noed_majrev.RData")

## parameters
num_cpus      = 5
num           = 1000
fun_extract   = get_la_fixef
mod_bl        = modc[[2]]
mod_compl_fit = modc[[3]]
mod           = all_modc[[3]]

## functions
agk.boot.p.lmer.subfun = function(mod_bl,mod,fun_extract,mod_compl_fit) {
  data.la = mod_compl_fit@frame
  # simulate under 0 hyothesis
  data.la$accept.reject <- (simulate(mod_bl))[[1]]
  # fit under alternative hypothesis
  cur_fit <- eval(parse(text=mod))
  return(fun_extract(cur_fit))
}
  
## running the bootstrap (no cov)
cl<-makeCluster(num_cpus)  
registerDoSNOW(cl)
effects_under_0 <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
  agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
} 
stopCluster(cl)

setwd(pfad_results)
save(file="effects_under_0_51_noed_majrev.RData",list=c("effects_under_0"))
load("effects_under_0_51_noed_majrev.RData")


cur_vars    <- c("bg_HC","bg_PG","bg_AD","bl_HC", "bl_PG", "bl_AD","x_la_HCgrAD", "x_la_HCgrPG", "x_la_PGgrAD")
for (ii in 1:length(cur_vars)) {
    cur_fun     <- function(l) {return(l[cur_vars[ii]])}
  eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0,FUN = cur_fun))",sep="")))
}
con_bl_PG_AD <- bl_PG-bl_AD 
con_bg_PG_AD <- bg_PG-bg_AD 
con_bl_PG_HC <- bl_PG-bl_HC 
con_bg_PG_HC <- bg_PG-bg_HC 
con_bl_AD_HC <- bl_AD-bl_HC 
con_bg_AD_HC <- bg_AD-bg_HC 

## p-values for HC gr PG
obs_fixef <- get_la_fixef(modc[[3]]) # old study f 48: modc[[20]]
HCgrPG_density <- density(x_la_HCgrPG,bw = c("ucv"))
HCgrPG_density$y <- cumsum(HCgrPG_density$y/sum(HCgrPG_density$y))
f <- approxfun(HCgrPG_density, rule=2)
1-f(obs_fixef[14])

# p-value for PG gr HC: beta gain
PGgrHC_density_bg <- density(con_bg_PG_HC,bw = c("ucv"))
PGgrHC_density_bg$y <- cumsum(PGgrHC_density_bg$y/sum(PGgrHC_density_bg$y))
f <- approxfun(PGgrHC_density_bg, rule=2)
1-f(obs_fixef[6]-obs_fixef[5])

# p-value for PG gr HC: beta loss
PGgrHC_density_bl <- density(con_bl_PG_HC,bw = c("ucv"))
PGgrHC_density_bl$y <- cumsum(PGgrHC_density_bl$y/sum(PGgrHC_density_bl$y))
f <- approxfun(PGgrHC_density_bl, rule=2)
1-f(obs_fixef[9]-obs_fixef[8])

## p-value for HC gr AD
obs_fixef <- get_la_fixef(modc[[3]])
HCgrAD_density <- density(x_la_HCgrAD,bw = c("ucv"))
HCgrAD_density$y <- cumsum(HCgrAD_density$y/sum(HCgrAD_density$y))
f <- approxfun(HCgrAD_density, rule=2)
1-f(obs_fixef[13])

# p-value for HC gr AD: beta gain
ADgrHC_density_bg <- density(con_bg_AD_HC,bw = c("ucv"))
ADgrHC_density_bg$y <- cumsum(ADgrHC_density_bg$y/sum(ADgrHC_density_bg$y))
f <- approxfun(ADgrHC_density_bg, rule=2)
f(obs_fixef[4]-obs_fixef[5])

# p-value for HC gr AD: beta loss
ADgrHC_density_bl <- density(con_bl_AD_HC,bw = c("ucv"))
ADgrHC_density_bl$y <- cumsum(ADgrHC_density_bl$y/sum(ADgrHC_density_bl$y))
f <- approxfun(ADgrHC_density_bl, rule=2)
1-f(obs_fixef[7]-obs_fixef[8])

## p-value for PG gr AD
obs_fixef <- get_la_fixef(modc[[3]])
PGgrAD_density <- density(x_la_PGgrAD,bw = c("ucv"))
PGgrAD_density$y <- cumsum(PGgrAD_density$y/sum(PGgrAD_density$y))
f <- approxfun(PGgrAD_density, rule=2)
1-f(obs_fixef[15])

# p-value for PG gr AD: beta loss
PGgrAD_density_bl <- density(con_bl_PG_AD,bw = c("ucv"))
PGgrAD_density_bl$y <- cumsum(PGgrAD_density_bl$y/sum(PGgrAD_density_bl$y))
f <- approxfun(PGgrAD_density_bl, rule=2)
f(obs_fixef[9]-obs_fixef[7])

# p-value for PG gr AD: beta gain
PGgrAD_density_bg <- density(con_bg_PG_AD,bw = c("ucv"))
PGgrAD_density_bg$y <- cumsum(PGgrAD_density_bg$y/sum(PGgrAD_density_bg$y))
f <- approxfun(PGgrAD_density_bg, rule=2)
1-f(obs_fixef[6]-obs_fixef[4])

## bootstrap CIs for actual effects;
mod_compl_fit = modc[[3]]
boot_la <-bootMer(mod_compl_fit,use.u = F,
                  FUN = get_la_fixef,
                  nsim = num,
                  type = "parametric",
                  verbose = TRUE,
                  parallel = "snow",
                  ncpus = num_cpus)

setwd(pfad_results)
save(file="boot_la_compl_mod_no_cov51_noed_majrev.RData", list=c("boot_la"))

# ## modeling the covariates to check if results stay the same; first: EDU
# mod_compl_fit <- eval(parse(text=all_modc[[120]]))     # fit complete model
# mod           <- all_modc[[120]] # unfitted model formular of complete model
# mod_bl        <- eval(parse(text=all_modc[[118]]))
# 
# 
# ## running the bootstrap
# cl<-makeCluster(4)  
# registerDoSNOW(cl)
# effects_under_0_edu <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
#   agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
# } 
# stopCluster(cl)
# 
# setwd(pfad_results)
# save(file="effects_under_0_edu.RData",list=c("effects_under_0_edu"))
# load(file="effects_under_0_edu.RData")
# 
# cur_vars    <- c("bg_HC","bg_PG","bg_AD","bl_HC", "bl_PG", "bl_AD","x_la_HCgrAD", "x_la_HCgrPG", "x_la_PGgrAD")
# for (ii in 1:length(cur_vars)) {
#   cur_fun     <- function(l) {return(l[cur_vars[ii]])}
#   eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0_edu,FUN = cur_fun))",sep="")))
# }
# 
# con_bl_PG_AD <- bl_PG-bl_AD 
# con_bg_PG_AD <- bg_PG-bg_AD 
# con_bl_PG_HC <- bl_PG-bl_HC 
# con_bg_PG_HC <- bg_PG-bg_HC 
# con_bl_AD_HC <- bl_AD-bl_HC 
# con_bg_AD_HC <- bg_AD-bg_HC 
# 
# # p-value for HC gr PG
# obs_fixef <- get_la_fixef(modc[[120]])
# HCgrPG_density <- density(x_la_HCgrPG,bw = c("ucv"))
# HCgrPG_density$y <- cumsum(HCgrPG_density$y/sum(HCgrPG_density$y))
# f <- approxfun(HCgrPG_density, rule=2)
# 1-f(obs_fixef[14])
# 
# # p-value for HC gr AD
# obs_fixef <- get_la_fixef(modc[[120]])
# HCgrAD_density <- density(x_la_HCgrAD,bw = c("ucv"))
# HCgrAD_density$y <- cumsum(HCgrAD_density$y/sum(HCgrAD_density$y))
# f <- approxfun(HCgrAD_density, rule=2)
# 1-f(obs_fixef[13])
# 
# # p-value for PG gr AD
# obs_fixef <- get_la_fixef(modc[[120]])
# PGgrAD_density <- density(x_la_PGgrAD,bw = c("ucv"))
# PGgrAD_density$y <- cumsum(PGgrAD_density$y/sum(PGgrAD_density$y))
# f <- approxfun(PGgrAD_density, rule=2)
# 1-f(obs_fixef[15])

## modeling the covariates to check if results stay the same;AGE
mod_compl_fit <- modc[[7]]
mod           <- all_modc[[7]]
mod_bl        <- modc[[6]]


## running the bootstrap (AGE)
cl<-makeCluster(num_cpus)  
registerDoSNOW(cl)
effects_under_0_age <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
  agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
} 
stopCluster(cl)

setwd(pfad_results)
save(file="effects_under_0_age51_noed_majrev.RData",list=c("effects_under_0_age"))
load("effects_under_0_age51_noed_majrev.RData")

cur_vars    <- c("bg_HC","bg_PG","bg_AD","bl_HC", "bl_PG", "bl_AD","x_la_HCgrAD", "x_la_HCgrPG", "x_la_PGgrAD")
for (ii in 1:length(cur_vars)) {
  cur_fun     <- function(l) {return(l[cur_vars[ii]])}
  eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0_age,FUN = cur_fun))",sep="")))
}

con_bl_PG_AD <- bl_PG-bl_AD 
con_bg_PG_AD <- bg_PG-bg_AD 
con_bl_PG_HC <- bl_PG-bl_HC 
con_bg_PG_HC <- bg_PG-bg_HC 
con_bl_AD_HC <- bl_AD-bl_HC 
con_bg_AD_HC <- bg_AD-bg_HC 

# p-value for HC gr PG
obs_fixef <- get_la_fixef(modc[[7]])
HCgrPG_density <- density(x_la_HCgrPG,bw = c("ucv"))
HCgrPG_density$y <- cumsum(HCgrPG_density$y/sum(HCgrPG_density$y))
f <- approxfun(HCgrPG_density, rule=2)
1-f(obs_fixef[14])

# p-value for HC gr AD
obs_fixef <- get_la_fixef(modc[[7]])
HCgrAD_density <- density(x_la_HCgrAD,bw = c("ucv"))
HCgrAD_density$y <- cumsum(HCgrAD_density$y/sum(HCgrAD_density$y))
f <- approxfun(HCgrAD_density, rule=2)
1-f(obs_fixef[13])

# p-value for PG gr AD
obs_fixef <- get_la_fixef(modc[[7]])
PGgrAD_density <- density(x_la_PGgrAD,bw = c("ucv"))
PGgrAD_density$y <- cumsum(PGgrAD_density$y/sum(PGgrAD_density$y))
f <- approxfun(PGgrAD_density, rule=2)
1-f(obs_fixef[15])

## beta loss and beta gain ##
# p-value for PG gr AD: beta loss
obs_fixef <- get_la_fixef(modc[[7]])
PGgrAD_density_bl <- density(con_bl_PG_AD,bw = c("ucv"))
PGgrAD_density_bl$y <- cumsum(PGgrAD_density_bl$y/sum(PGgrAD_density_bl$y))
f <- approxfun(PGgrAD_density_bl, rule=2)
f(obs_fixef[9]-obs_fixef[7])

# p-value for PG gr AD: beta gain
PGgrAD_density_bg <- density(con_bg_PG_AD,bw = c("ucv"))
PGgrAD_density_bg$y <- cumsum(PGgrAD_density_bg$y/sum(PGgrAD_density_bg$y))
f <- approxfun(PGgrAD_density_bg, rule=2)
1-f(obs_fixef[6]-obs_fixef[4])

# p-value for PG gr HC: beta gain
PGgrHC_density_bg <- density(con_bg_PG_HC,bw = c("ucv"))
PGgrHC_density_bg$y <- cumsum(PGgrHC_density_bg$y/sum(PGgrHC_density_bg$y))
f <- approxfun(PGgrHC_density_bg, rule=2)
1-f(obs_fixef[6]-obs_fixef[5])

# p-value for PG gr HC: beta loss
PGgrHC_density_bl <- density(con_bl_PG_HC,bw = c("ucv"))
PGgrHC_density_bl$y <- cumsum(PGgrHC_density_bl$y/sum(PGgrHC_density_bl$y))
f <- approxfun(PGgrHC_density_bl, rule=2)
1-f(obs_fixef[9]-obs_fixef[8])

# p-value for HC gr AD: beta gain
ADgrHC_density_bg <- density(con_bg_AD_HC,bw = c("ucv"))
ADgrHC_density_bg$y <- cumsum(ADgrHC_density_bg$y/sum(ADgrHC_density_bg$y))
f <- approxfun(ADgrHC_density_bg, rule=2)
f(obs_fixef[4]-obs_fixef[5])

# p-value for HC gr AD: beta loss
ADgrHC_density_bl <- density(con_bl_AD_HC,bw = c("ucv"))
ADgrHC_density_bl$y <- cumsum(ADgrHC_density_bl$y/sum(ADgrHC_density_bl$y))
f <- approxfun(ADgrHC_density_bl, rule=2)
1-f(obs_fixef[7]-obs_fixef[8])

# ## running the bootstrap (EDU AND AGE BOTH AS COVARIATES; STUDY OF 51)
# # models
# mod_compl_fit <- modc[[4]]
# mod           <- all_modc[[4]]
# mod_bl        <- modc[[5]]
# 
# # running...
# cl<-makeCluster(6)  
# registerDoSNOW(cl)
# effects_under_0_edu_age <- foreach(ff=1:num, .packages='lme4',.verbose=T,.export = c("data.la")) %dopar% {
#   agk.boot.p.lmer.subfun(mod_bl,mod,fun_extract,mod_compl_fit)
# } 
# stopCluster(cl)
# 
# setwd(pfad_results)
# save(file="effects_under_0_eduage51.RData",list=c("effects_under_0_edu_age"))
# load("effects_under_0_eduage51.RData")
# 
# cur_vars    <- c("bg_HC","bg_PG","bg_AD","bl_HC", "bl_PG", "bl_AD","x_la_HCgrAD", "x_la_HCgrPG", "x_la_PGgrAD")
# for (ii in 1:length(cur_vars)) {
#   cur_fun     <- function(l) {return(l[cur_vars[ii]])}
#   eval(parse(text=paste(cur_vars[ii], "<- unlist(lapply(effects_under_0_edu_age,FUN = cur_fun))",sep="")))
# }
# 
# con_bl_PG_AD <- bl_PG-bl_AD 
# con_bg_PG_AD <- bg_PG-bg_AD 
# con_bl_PG_HC <- bl_PG-bl_HC 
# con_bg_PG_HC <- bg_PG-bg_HC 
# con_bl_AD_HC <- bl_AD-bl_HC 
# con_bg_AD_HC <- bg_AD-bg_HC 
# 
# # p-value for HC gr PG
# obs_fixef <- get_la_fixef(modc[[119]])
# HCgrPG_density <- density(x_la_HCgrPG,bw = c("ucv"))
# HCgrPG_density$y <- cumsum(HCgrPG_density$y/sum(HCgrPG_density$y))
# f <- approxfun(HCgrPG_density, rule=2)
# 1-f(obs_fixef[14])
# 
# # p-value for HC gr AD
# obs_fixef <- get_la_fixef(modc[[119]])
# HCgrAD_density <- density(x_la_HCgrAD,bw = c("ucv"))
# HCgrAD_density$y <- cumsum(HCgrAD_density$y/sum(HCgrAD_density$y))
# f <- approxfun(HCgrAD_density, rule=2)
# 1-f(obs_fixef[13])
# 
# # p-value for PG gr AD
# obs_fixef <- get_la_fixef(modc[[119]])
# PGgrAD_density <- density(x_la_PGgrAD,bw = c("ucv"))
# PGgrAD_density$y <- cumsum(PGgrAD_density$y/sum(PGgrAD_density$y))
# f <- approxfun(PGgrAD_density, rule=2)
# 1-f(obs_fixef[15])