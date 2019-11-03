# MODEL comparison, LA study

# compare candidate models
# which one is best
# and are there differences in model fit between groups?
# there should not be because, if so, then reviewers could say:
# "Your group differences are only due to differences in model fit."

data_pdt = data.la
data_pdt$gain = data_pdt$Gewinn
data_pdt$loss = data_pdt$Verlust
data_pdt$ed_abs = data_pdt$ed.abs
data_pdt$accept_reject = data_pdt$accept.reject

## get the models and their BICs/AICs per subject
# 1st model: our favorite model:
la_list      = lmList(accept_reject ~ gain + loss + ed_abs | subject,data=data_pdt,family=binomial, pool =F,na.action=NULL)
BICs_lae     = unlist(lapply(la_list,FUN=BIC))
AICs_lae     = unlist(lapply(la_list,FUN=AIC))
c_lae        = coef(la_list)
c_lae$lambda = get_log(-c_lae$loss/c_lae$gain)

# 2nd model: our favorite model, but without ed
la_list      = lmList(accept_reject ~ gain + loss | subject,data=data_pdt,family=binomial, pool =F,na.action=NULL)
BICs_la      = unlist(lapply(la_list,FUN=BIC))
AICs_la      = unlist(lapply(la_list,FUN=AIC))
c_la         = coef(la_list)
c_la$lambda  = get_log(-c_la$loss/c_la$gain)

# 3rd model: the TOM model
la_list      = lmList(accept_reject ~ 0 + gain + loss | subject,data=data_pdt,family=binomial, pool =F,na.action=NULL)
BICs_lat     = unlist(lapply(la_list,FUN=BIC))
AICs_lat     = unlist(lapply(la_list,FUN=AIC))
c_lat        = coef(la_list)
c_lat$lambda = get_log(-c_lat$loss/c_lat$gain)

# 4th model: the ratio model
la_list      = lmList(accept_reject ~ 1 + ratio | subject,data=data_pdt,family=binomial, pool =F,na.action=NULL)
BICs_lar     = unlist(lapply(la_list,FUN=BIC))
AICs_lar     = unlist(lapply(la_list,FUN=AIC))
c_lar        = coef(la_list)
c_lar$lambda = get_log(c_lar$ratio)


# 5th model the Charpentier model
setwd('C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom\\Scripts\\model_calculation')
source("MLE_Charpentier_updated.R")
c_lac_lambda = get_log(all_params$lambda*(-1))


# get the BICs and AICs of the Charpentier model
BICs_lac     = c()
AICs_lac     = c()
k = 2
for (ii in 1:length(all_params[,1])) {
  cur_params    = all_params[ii,(2:3)]
  cur_params[3] = 0
  cur_sub       = all_params[ii,1]
  # get all the trials
  cur_data = subset(data_pdt, subject == as.character(cur_sub))
  x        = cur_data
  x$accnum=as.numeric(as.character(x$accept_reject))
  cur_lh       = exp(fn.c(theta = cur_params,x =x))
  BICs_lac[ii] = -((-2)*log(cur_lh)+k*log(length(cur_data[,1])))
  AICs_lac[ii] = -(2*k-2*log(cur_lh))
}

# BIC comparison
BICs_lac = data.frame(all_params$VPNR,BICs_lac)
BIC_df = data.frame(names(la_list),BICs_lae,BICs_la,BICs_lat,BICs_lar)
BIC_df = merge(BIC_df,BICs_lac,by.x='names.la_list.',by.y = 'all_params.VPNR')
names(BIC_df)[1] = "subject"
BIC_df$group = agk.recode.c(as.character(BIC_df$subject),as.character(data_pdt$subject),as.character(data_pdt$group))

BIC_dfs = subset(BIC_df, (group == "HC" | group == "PG"))
BIC_dfs = melt(BIC_dfs)

# test between models first
describeBy(BIC_dfs$value,BIC_dfs$variable)
BIC_dfs$modelgroup = paste0(BIC_dfs$variable,BIC_dfs$group)
describe.by(BIC_dfs$value,BIC_dfs$modelgroup)
aggregate(BIC_dfs$value,by = list(BIC_dfs$group,BIC_dfs$variable),FUN=mean)

mod_00 = lm(value ~ variable,data=BIC_dfs)
mod_01 = lm(value ~ variable+group,data=BIC_dfs)
mod_02 = lm(value ~ variable*group,data=BIC_dfs)
anova(mod_00,mod_01)
anova(mod_01,mod_02)

# AIC comparison
AICs_lac = data.frame(all_params$VPNR,AICs_lac)
AIC_df   = data.frame(names(la_list),AICs_lae,AICs_la,AICs_lat,AICs_lar)
AIC_df   = merge(AIC_df,AICs_lac,by.x='names.la_list.',by.y = 'all_params.VPNR')
names(AIC_df)[1] = "subject"
AIC_df$group = agk.recode.c(as.character(AIC_df$subject),as.character(data_pdt$subject),as.character(data_pdt$group))

AIC_dfs = subset(AIC_df, (group == "HC" | group == "PG" | group == "AD"))
AIC_dfs = melt(AIC_dfs)

# test between models first
describeBy(AIC_dfs$value,AIC_dfs$variable)
AIC_dfs$modelgroup = paste0(AIC_dfs$variable,AIC_dfs$group)
describe.by(AIC_dfs$value,AIC_dfs$modelgroup)
m = aggregate(AIC_dfs$value,by = list(AIC_dfs$group,AIC_dfs$variable),FUN=mean)

mod_00 = lm(value ~ variable,data=AIC_dfs)
mod_01 = lm(value ~ variable+group,data=AIC_dfs)
mod_02 = lm(value ~ variable*group,data=AIC_dfs)
anova(mod_00,mod_01)
anova(mod_01,mod_02)

summary(mod_00)
summary(mod_01)
summary(mod_02)

AIC_dfs_HC = subset(AIC_dfs, group == "HC")
AIC_dfs_PG = subset(AIC_dfs, group == "PG")
AIC_dfs_AD = subset(AIC_dfs, group == "AD")
# HC < PG
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lae"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lae"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lat"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lar"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lar"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lac"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lac"])
# HC < AD
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lae"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lae"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_la"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lat"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lar"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lar"])
t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lac"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lac"])
# PG < AD
t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lae"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lae"])
t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_la"])
t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lat"])
t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lar"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lar"])
t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lac"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lac"])

## Correlation
# getting the lambda of all the models
all_lambdas = data.frame(c_lae$lambda,c_la$lambda,c_lat$lambda,c_lar$lambda,c_lac_lambda)
corr.test(all_lambdas)
tmp = corr.test(all_lambdas)

setwd(pfad_results)
write.table(tmp$r,file="lambda_corr")

