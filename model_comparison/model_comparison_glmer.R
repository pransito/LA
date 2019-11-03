# MODEL comparison, LA study

# compare candidate models
# which one is best
# and are there differences in model fit between groups?
# there should not be because, if so, then reviewers could say:
# "Your group differences are only due to differences in model fit."

# params
# if 0 then already estimated; will only load
estimate_models = 0

# data preparation
data_pdt = data.la
data_pdt$gain = data_pdt$Gewinn
data_pdt$loss = data_pdt$Verlust
data_pdt$ed_abs = data_pdt$ed.abs
data_pdt$accept_reject = data_pdt$accept.reject
data_pdt_HC = subset(data_pdt,group == "HC")
data_pdt_PG = subset(data_pdt,group == "PG")
data_pdt_AD = subset(data_pdt,group == "AD")

# path
setwd(pfad_model_calc)

# get the models and their BICs/AICs per subject
## 1st model: our favorite model:
if (estimate_models) {
  #null_m           = glmer(accept_reject ~ 1 + (1|subject),data=data_pdt,family=binomial)
  lae_list         = glmer(accept_reject ~ (gain + loss + ed_abs)*group +(gain + loss + ed_abs | subject),data=data_pdt,family=binomial)
  lae_list_nofe    = glmer(accept_reject ~ 0 +(0+gain + loss + ed_abs | subject),data=data_pdt,family=binomial)
  lae_list_HC      = glmer(accept_reject ~ (gain + loss + ed_abs) +(gain + loss + ed_abs | subject),data=data_pdt_HC,family=binomial)
  lae_list_PG      = glmer(accept_reject ~ (gain + loss + ed_abs) +(gain + loss + ed_abs | subject),data=data_pdt_PG,family=binomial)
  lae_list_AD      = glmer(accept_reject ~ (gain + loss + ed_abs) +(gain + loss + ed_abs | subject),data=data_pdt_AD,family=binomial)
  save(file="lae_list.RData",list = c("lae_list","lae_list_nofe","lae_list_HC","lae_list_PG","lae_list_AD"))
} else {
  load("lae_list.RData")
}

# model fit
AIC_lae          = AIC(lae_list)/51
AIC_lae_HC       = AIC(lae_list_HC)/17
AIC_lae_PG       = AIC(lae_list_PG)/19
AIC_lae_AD       = AIC(lae_list_AD)/15
AIC_lae_splj     = mean(c(AIC_lae_HC,AIC_lae_PG,AIC_lae_AD))

BIC_lae          = BIC(lae_list)/51
BIC_lae_HC       = BIC(lae_list_HC)/17
BIC_lae_PG       = BIC(lae_list_PG)/19
BIC_lae_AD       = BIC(lae_list_AD)/15
BIC_lae_splj     = mean(c(BIC_lae_HC,BIC_lae_PG,BIC_lae_AD))

# lambda
c_lae            = agk.get.compl.coef(lae_list,int="group")
c_lae$lambda     = -c_lae$loss/c_lae$gain

## 2nd model: our favorite model, but without ed
if (estimate_models) {
  la_list         = glmer(accept_reject ~ (gain + loss)*group +(gain + loss | subject),data=data_pdt,family=binomial)
  la_list_HC      = glmer(accept_reject ~ (gain + loss) +(gain + loss | subject),data=data_pdt_HC,family=binomial)
  la_list_PG      = glmer(accept_reject ~ (gain + loss) +(gain + loss | subject),data=data_pdt_PG,family=binomial)
  la_list_AD      = glmer(accept_reject ~ (gain + loss) +(gain + loss | subject),data=data_pdt_AD,family=binomial)
  save(file="la_list.RData",list = c("la_list","la_list_HC","la_list_PG","la_list_AD"))
} else {
  load("la_list.RData")
}

# model fit
AIC_la          = AIC(la_list)/51
AIC_la_HC       = AIC(la_list_HC)/17
AIC_la_PG       = AIC(la_list_PG)/19
AIC_la_AD       = AIC(la_list_AD)/15
AIC_la_splj     = mean(c(AIC_la_HC,AIC_la_PG,AIC_la_AD))

BIC_la          = BIC(la_list)/51
BIC_la_HC       = BIC(la_list_HC)/17
BIC_la_PG       = BIC(la_list_PG)/19
BIC_la_AD       = BIC(la_list_AD)/15
BIC_la_splj     = mean(c(BIC_la_HC,BIC_la_PG,BIC_la_AD))

# lambda
c_la            = agk.get.compl.coef(la_list,int="group")
c_la$lambda     = -c_la$loss/c_la$gain

# 4th model: the ratio model
# use orig variable to get interpretable results! (the break even ratio...)
# no difference in model fit to centered variable
if (estimate_models) {
  lar_list         = glmer(accept_reject ~ (ratio_bcp)*group +(ratio_bcp | subject),data=data_pdt,family=binomial)
  lar_list_HC      = glmer(accept_reject ~ (ratio_bcp) +(ratio_bcp | subject),data=data_pdt_HC,family=binomial)
  lar_list_PG      = glmer(accept_reject ~ (ratio_bcp) +(ratio_bcp | subject),data=data_pdt_PG,family=binomial)
  lar_list_AD      = glmer(accept_reject ~ (ratio_bcp) +(ratio_bcp | subject),data=data_pdt_AD,family=binomial)
  save(file="lar_list.RData",list = c("lar_list","lar_list_HC","lar_list_PG","lar_list_AD"))
} else {
  load("lar_list.RData")
}


# the lambda
c_lar      = agk.get.compl.coef(lar_list, int="group")
lambda_lar = c()
for (ii in 1:length(c_lar[,1])) {
  cur_par        = c_lar[ii,]
  lambda_lar[ii] = as.numeric(-cur_par[1]/cur_par[2])  
}
c_lar$lambda = lambda_lar

# model fit
AIC_lar          = AIC(lar_list)/51
AIC_lar_HC       = AIC(lar_list_HC)/17
AIC_lar_PG       = AIC(lar_list_PG)/19
AIC_lar_AD       = AIC(lar_list_AD)/15
AIC_lar_splj     = mean(c(AIC_lar_HC,AIC_lar_PG,AIC_lar_AD))

BIC_lar          = BIC(lar_list)/51
BIC_lar_HC       = BIC(lar_list_HC)/17
BIC_lar_PG       = BIC(lar_list_PG)/19
BIC_lar_AD       = BIC(lar_list_AD)/15
BIC_lar_splj     = mean(c(BIC_lar_HC,BIC_lar_PG,BIC_lar_AD))

# the models just as lmlist
lae_lmlist = lmList(accept_reject ~ gain + loss + ed_abs|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)
la_lmlist = lmList(accept_reject ~ gain + loss|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)
lar_lmlist = lmList(accept_reject ~ ratio_bcp|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)

# coefs lmlist
c_lae_lmlist = coef(lae_lmlist)
c_lae_lmlist$lambda = -c_lae_lmlist$loss/c_lae_lmlist$gain
c_lae_lmlist$group = agk.recode.c(row.names(c_lae_lmlist),data_pdt$subject,data_pdt$group)

c_la_lmlist = coef(la_lmlist)
c_la_lmlist$lambda = -c_la_lmlist$loss/c_la_lmlist$gain
c_la_lmlist$group = agk.recode.c(row.names(c_la_lmlist),data_pdt$subject,data_pdt$group)

c_lar_lmlist = coef(lar_lmlist)
c_lar_lmlist$group = agk.recode.c(row.names(c_lar_lmlist),data_pdt$subject,data_pdt$group)
lambda_lar_lmlist = c()
for (ii in 1:length(c_lar_lmlist[,1])) {
  cur_par               = c_lar_lmlist[ii,]
  lambda_lar_lmlist[ii] = as.numeric(-cur_par$`(Intercept)`/cur_par$ratio_bcp)  
}
c_lar_lmlist$lambda = lambda_lar_lmlist
lmlists = list(c_lae_lmlist,c_la_lmlist,c_lar_lmlist)
names(lmlists) = c("c_lae_lmlist","c_la_lmlist","c_lar_lmlist")

# group testing lmlists
for (ii in 1:length(lmlists)) {
  disp(names(lmlists)[ii])
  cur_dat = lmlists[[ii]]
  cur_mod = lm(lambda ~ group, data=cur_dat)
  print(summary(cur_mod))
  print(describeBy(cur_dat$lambda,cur_dat$group))
  plot(cur_dat$lambda~as.factor(cur_dat$group))
}


# 5th model the Charpentier model and all other models with mu
# scripts will estimate if estimate_models == 1, else just load params and function for likelihood estimation
source("MLE_Charpentier_updated.R")
c_lacm_lambda  = all_params$lambda*(-1)
source("MLE_Charpentier_updated_no_mu.R")
c_lac_lambda  = all_params$lambda*(-1)
source("MLE_la_mu.R")
c_lam_lambda = all_params$beta_loss*(-1)/all_params$beta_gain
source("MLE_lar_mu.R")
crm = all_params
lambda_lar = c()
for (ii in 1:length(crm[,1])) {
  cur_par        = crm[ii,]
  lambda_lar[ii] = as.numeric(-cur_par$beta_int/cur_par$beta_ra)  
}
crm$lambda    = lambda_lar
c_larm_lambda = lambda_lar
source("MLE_lae_mu.R")
c_laem_lambda = all_params$beta_loss*(-1)/all_params$beta_gain

crm = all_params

crm$lambda = crm$lambda*(-1)
crm$lambda = crm$beta_loss*(-1)/crm$beta_gain
crm$group = as.factor(agk.recode.c(as.character(crm$VPNR),as.character(data_pdt$subject),as.character(data_pdt$group)))
crm_HC = subset(crm,group == "HC")
#crm_HC$lambda = winsor(crm_HC$lambda,trim=0.1)
crm_PG = subset(crm,group == "PG")
crm_AD = subset(crm,group == "AD")
#crm_AD$lambda = winsor(crm_AD$lambda,trim=0.1)
crm = rbind(crm_HC,crm_PG,crm_AD)

c_lac = crm

# lambda group testing of glmer model
cur_mod = la_list
crm = agk.get.compl.coef(cur_mod,int = "group")
crm = c_lar

# testing
contrasts(crm$group) = contr.treatment(levels(crm$group),base=1)
contrasts(crm$group) = contr.treatment(levels(crm$group),base=2)
crm$lambda = -crm$loss/crm$gain
cur_model = lm(lambda~group,crm)
cur_model = lm(mu~group,crm)
#cur_model = lm(lambda_winsor~group,crm)
summary(cur_model)
kruskal.test(crm$lambda,crm$group)
describeBy(crm$lambda,crm$group)
plot(crm$lambda~crm$group)
describeBy(crm$mu,crm$group)
plot(crm$mu~crm$group)

# get the BICs and AICs of the Charpentier model/models with mu
BICs_lac     = c()
AICs_lac     = c()
LH_lac       = c()
k = 3 # numbers of params in optim model
for (ii in 1:length(all_params[,1])) {
  cur_params    = all_params[ii,c(2:(k+1))]
  cur_sub       = all_params$VPNR[ii]
  # get all the trials
  cur_data = subset(data_pdt, subject == as.character(cur_sub))
  x        = cur_data
  x$accnum=as.numeric(as.character(x$accept_reject))
  cur_lh       = exp(fn.c(theta = cur_params,x =x))
  LH_lac[ii]   = cur_lh
  BICs_lac[ii] = -((-2)*log(cur_lh)+k*log(length(cur_data[,1])))
  AICs_lac[ii] = -(2*k-2*log(cur_lh))
}

# log-likelihood ratio test
agk.load.ifnot.install("lmtest")
lrtest(lae_list,lar_list)
lrtest(lae_list,la_list)
m = anova(lae_list,la_list)
lh_lae = m$logLik[2]
lh_lac = -sum(log(LH_lac))
df_lae = m$Df[2]
df_lac = k
df_dif = df_lae-df_lac
D = 2*(lh_lae - lh_lac)
D = 2*(lh_lac - lh_lae)
dchisq(D,df_dif)


# BIC comparison
BICs_lac = data.frame(all_params$VPNR,BICs_lac)
# BIC_df = data.frame(names(la_list),BICs_lae,BICs_la,BICs_lat,BICs_lar)
# BIC_df = merge(BIC_df,BICs_lac,by.x='names.la_list.',by.y = 'all_params.VPNR')
# names(BIC_df)[1] = "subject"
# BIC_df$group = agk.recode.c(as.character(BIC_df$subject),as.character(data_pdt$subject),as.character(data_pdt$group))
BICs_lac$group = agk.recode.c(as.character(BICs_lac$all_params.VPNR),as.character(data_pdt$subject),as.character(data_pdt$group))
desc = describeBy(BICs_lac$BICs_lac,BICs_lac$group)
print(desc)
mean(c(as.numeric(desc$AD$mean),as.numeric(desc$HC$mean),as.numeric(desc$PG$mean)))
sum(BICs_lac$BICs_lac/51)

# BIC_dfs = subset(BIC_df, (group == "HC" | group == "PG"))
# BIC_dfs = melt(BIC_dfs)

# test between models first
# describeBy(BIC_dfs$value,BIC_dfs$variable)
# BIC_dfs$modelgroup = paste0(BIC_dfs$variable,BIC_dfs$group)
# describe.by(BIC_dfs$value,BIC_dfs$modelgroup)
# aggregate(BIC_dfs$value,by = list(BIC_dfs$group,BIC_dfs$variable),FUN=mean)
# 
# mod_00 = lm(value ~ variable,data=BIC_dfs)
# mod_01 = lm(value ~ variable+group,data=BIC_dfs)
# mod_02 = lm(value ~ variable*group,data=BIC_dfs)
# anova(mod_00,mod_01)
# anova(mod_01,mod_02)

# AIC comparison
AICs_lac = data.frame(all_params$VPNR,AICs_lac)
#AIC_df   = data.frame(names(la_list),AICs_lae,AICs_la,AICs_lat,AICs_lar)
#AIC_df   = merge(AIC_df,AICs_lac,by.x='names.la_list.',by.y = 'all_params.VPNR')
#names(AIC_df)[1] = "subject"
#AIC_df$group = agk.recode.c(as.character(AIC_df$subject),as.character(data_pdt$subject),as.character(data_pdt$group))
AICs_lac$group = agk.recode.c(as.character(all_params$VPNR),as.character(data_pdt$subject),as.character(data_pdt$group))
desc = describeBy(AICs_lac$AICs_lac,AICs_lac$group)
print(desc)
mean(c(as.numeric(desc$AD$mean),as.numeric(desc$HC$mean),as.numeric(desc$PG$mean)))
sum(AICs_lac$AICs_lac/51)

#AIC_dfs = subset(AIC_df, (group == "HC" | group == "PG" | group == "AD"))
#AIC_dfs = melt(AIC_dfs)

# test between models first
# describeBy(AIC_dfs$value,AIC_dfs$variable)
# AIC_dfs$modelgroup = paste0(AIC_dfs$variable,AIC_dfs$group)
# describe.by(AIC_dfs$value,AIC_dfs$modelgroup)
# m = aggregate(AIC_dfs$value,by = list(AIC_dfs$group,AIC_dfs$variable),FUN=mean)
# 
# mod_00 = lm(value ~ variable,data=AIC_dfs)
# mod_01 = lm(value ~ variable+group,data=AIC_dfs)
# mod_02 = lm(value ~ variable*group,data=AIC_dfs)
# anova(mod_00,mod_01)
# anova(mod_01,mod_02)
# 
# summary(mod_00)
# summary(mod_01)
# summary(mod_02)

# AIC_dfs_HC = subset(AIC_dfs, group == "HC")
# AIC_dfs_PG = subset(AIC_dfs, group == "PG")
# AIC_dfs_AD = subset(AIC_dfs, group == "AD")
# # HC < PG
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lae"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lae"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lat"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lar"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lar"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lac"],AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lac"])
# # HC < AD
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lae"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lae"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_la"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lat"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lar"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lar"])
# t.test(AIC_dfs_HC$value[AIC_dfs_HC$variable=="AICs_lac"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lac"])
# # PG < AD
# t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lae"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lae"])
# t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_la"])
# t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_la"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lat"])
# t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lar"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lar"])
# t.test(AIC_dfs_PG$value[AIC_dfs_PG$variable=="AICs_lac"],AIC_dfs_AD$value[AIC_dfs_AD$variable=="AICs_lac"])

## Correlation
# getting the lambda of all the models
all_lambdas = data.frame(c_lae$lambda,c_la$lambda,c_lar$lambda,c_lac$lambda,
                         c_laem_lambda,c_lam_lambda,c_larm_lambda,c_lacm_lambda)
names(all_lambdas) = c('lae','la','lar','lac','laem','lam','larm','lacm')
tmp_s = corr.test(all_lambdas,method = "spearman")
tmp_p = corr.test(all_lambdas,method = "pearson")

# edit for printing
cur_r = tmp_s
cur_r$r[cur_r$p>=0.05] = 0
cur_r$r = round(cur_r$r,2)
setwd(pfad_results)
write.table(cur_r$r,file="lambda_corr.txt",quote = FALSE)

# group comparisons
all_lambdas$subject = row.names(c_lar)
all_lambdas$group = agk.recode.c(as.character(all_lambdas$subject),as.character(data_pdt$subject),as.character(data_pdt$group))
all_lamm = melt(all_lambdas)

m = aggregate(all_lamm$value,by = list(all_lamm$group,all_lamm$variable),FUN=winsor.mean,trim = 0.1)
aggregate(all_lamm$value,by = list(all_lamm$group,all_lamm$variable),FUN=median)

all_lamm$group = as.factor(all_lamm$group)
all_lamm_HCPG  = subset(all_lamm, (group == "HC" | group == "PG"))
all_lamm_HCAD  = subset(all_lamm, (group == "HC" | group == "AD"))
all_lamm_PGAD  = subset(all_lamm, (group == "PG" | group == "AD"))

# HC < PG
kruskal.test(all_lamm_HCPG$value[all_lamm_HCPG$variable == "c_lae.lambda"],all_lamm_HCPG$group[all_lamm_HCPG$variable == "c_lae.lambda"])
kruskal.test(all_lamm_HCAD$value[all_lamm_HCAD$variable == "c_lae.lambda"],all_lamm_HCAD$group[all_lamm_HCAD$variable == "c_lae.lambda"])
kruskal.test(all_lamm_PGAD$value[all_lamm_PGAD$variable == "c_lae.lambda"],all_lamm_PGAD$group[all_lamm_PGAD$variable == "c_lae.lambda"])


