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

## 3rd model: lac (charpentier, like Tom but without intercept, no ed)
if (estimate_models) {
  lac_list         = glmer(accept_reject ~ (0+gain + loss)*group +(0+gain + loss | subject),data=data_pdt,family=binomial)
  lac_list_HC      = glmer(accept_reject ~ (0+gain + loss) +(0+gain + loss | subject),data=data_pdt_HC,family=binomial)
  lac_list_PG      = glmer(accept_reject ~ (0+gain + loss) +(0+gain + loss | subject),data=data_pdt_PG,family=binomial)
  lac_list_AD      = glmer(accept_reject ~ (0+gain + loss) +(0+gain + loss | subject),data=data_pdt_AD,family=binomial)
  save(file="lac_list.RData",list = c("lac_list","lac_list_HC","lac_list_PG","lac_list_AD"))
} else {
  load("lac_list.RData")
}

# model fit
AIC_lac          = AIC(lac_list)/51
AIC_lac_HC       = AIC(lac_list_HC)/17
AIC_lac_PG       = AIC(lac_list_PG)/19
AIC_lac_AD       = AIC(lac_list_AD)/15
AIC_lac_splj     = mean(c(AIC_lac_HC,AIC_lac_PG,AIC_lac_AD))

BIC_lac          = BIC(lac_list)/51
BIC_lac_HC       = BIC(lac_list_HC)/17
BIC_lac_PG       = BIC(lac_list_PG)/19
BIC_lac_AD       = BIC(lac_list_AD)/15
BIC_lac_splj     = mean(c(BIC_lac_HC,BIC_lac_PG,BIC_lac_AD))

# lambda
c_lac            = agk.get.compl.coef(lac_list,int="group")
c_lac$lambda     = -c_lac$loss/c_lac$gain

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
all_lambdas = data.frame(c_lae$lambda,c_la$lambda,c_lar$lambda,c_lac$lambda)
names(all_lambdas) = c('lae','la','lar','lac')
tmp_s = corr.test(all_lambdas,method = "spearman")
tmp_p = corr.test(all_lambdas,method = "pearson")

# edit for printing
cur_r = tmp_s
cur_r$r[cur_r$p>=0.05] = 0
cur_r$r = round(cur_r$r,2)
setwd(pfad_results)
write.table(cur_r$r,file="lambda_corr.txt",quote = FALSE)

# group comparisons
all_lambdas$subject = unique(lae_list@frame$subject)
all_lambdas$group = agk.recode.c(as.character(all_lambdas$subject),as.character(data_pdt$subject),as.character(data_pdt$group))
all_lamm = melt(all_lambdas)

m = aggregate(all_lamm$value,by = list(all_lamm$group,all_lamm$variable),FUN=winsor.mean,trim = 0.1)
aggregate(all_lamm$value,by = list(all_lamm$group,all_lamm$variable),FUN=median)

all_lamm$group = as.factor(all_lamm$group)
all_lamm_HCPG  = subset(all_lamm, (group == "HC" | group == "PG"))
all_lamm_HCAD  = subset(all_lamm, (group == "HC" | group == "AD"))
all_lamm_PGAD  = subset(all_lamm, (group == "PG" | group == "AD"))

# HC < PG
kruskal.test(all_lamm_HCPG$value[all_lamm_HCPG$variable == "lae"],all_lamm_HCPG$group[all_lamm_HCPG$variable == "lae"])
kruskal.test(all_lamm_HCAD$value[all_lamm_HCAD$variable == "lae"],all_lamm_HCAD$group[all_lamm_HCAD$variable == "lae"])
kruskal.test(all_lamm_PGAD$value[all_lamm_PGAD$variable == "lae"],all_lamm_PGAD$group[all_lamm_PGAD$variable == "lae"])

cur_mod = 'lac'
ttest.group(all_lamm_HCPG$value[all_lamm_HCPG$variable == cur_mod],as.factor(as.character(all_lamm_HCPG$group[all_lamm_HCPG$variable == cur_mod])),cur_alt = 'greater')
ttest.group(all_lamm_HCAD$value[all_lamm_HCAD$variable == cur_mod],as.factor(as.character(all_lamm_HCAD$group[all_lamm_HCAD$variable == cur_mod])),cur_alt = 'less')
ttest.group(all_lamm_PGAD$value[all_lamm_PGAD$variable == cur_mod],as.factor(as.character(all_lamm_PGAD$group[all_lamm_PGAD$variable == cur_mod])),cur_alt = 'two.sided')

aggregate(all_lamm$value,FUN=mean,by=list(all_lamm$group,all_lamm$variable))

