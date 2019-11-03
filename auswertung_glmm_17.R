## VD-ANALYSEN LOSS AVERSION ==================================================

# 02.11.2012
# path update: 03.11.2019

# fitting a Generalized Mixed Model to the LA data
# DV: accept.reject
# IV (fixed): Gewinn, Verlust, Group
# IV (randm): subject random

## clear workspace

rm(list = ls())

## YET CHANGE TO NEW MODEL AFTER ED FIX ##

## preparation
# path
path_gd = 'E:/Google Drive'
path = paste0(path_gd, '/09_Diplom/LA/daten_behav_test_finale_SP_Diplom')
path_data = paste0(path,'/Data')
path_scripts = 'E:/GitHub/LA'
path_scripts_data = paste0(path_scripts,'/get_data')
path_scripts_components = paste0(path_scripts,'/get_components') 
path_scripts_components = paste0(path_scripts,'/model_calculation')
path_results = paste0(path,'/Results')
path_library = 'E:/GitHub/R'
path_ransim = paste0(path_results, "/rand_sim")
path_R_cit = paste0(path_gd, '/Diplom/Manuskript/Lit to be impl in Zotero')
path_model_calc = paste0(path_scripts, '/model_calculation')

# parameters
run_ana          = 0 # run analysis? 
des_agg          = 3
# exporting model paramters to a .mat file for matlab (for model based fmri)
export_to_matlab = 0
# crm for a coefficients from a random effects model (see below which one)
# cfm for a fixed effects model (see below which one) (mu models)
# cfml for fixed effects model based on lmlist (no mu)
param_df         = 'crm'

# load libraries and functions, options (can be modified)
setwd(path_library)
source ('agk_library.R')
setwd(path_scripts_data)
source ('LA_options.R')

# get all the data in long format
setwd(path_scripts_data)
source("get_data_la_2.R")

# only select uncertain gambles
if (only.some == 1){
  setwd(path_scripts)
  source('only_uncertain_gambles.R')  
  # set the pwd needed
  setwd(path_data)
}

setwd(path_scripts_data)
source("LA_add_demograpics.R")

# my contrasts
contrasts(data.la$group) <- cbind(c(1, 0, 0), c(0, 0, 1))
colnames(contrasts(data.la$group)) <- c("AD>HC", "PG>HC")

# data preparation for model extraction (see model comparison script)
estimate_models = 0
data_pdt = data.la
data_pdt$gain = data_pdt$Gewinn
data_pdt$loss = data_pdt$Verlust
data_pdt$ed_abs = data_pdt$ed.abs
data_pdt$accept_reject = data_pdt$accept.reject
data_pdt_HC = subset(data_pdt,group == "HC")
data_pdt_PG = subset(data_pdt,group == "PG")
data_pdt_AD = subset(data_pdt,group == "AD")

## get the crm/cfm
if (param_df == 'crm') {
  setwd(path_results)
  load("modc51_noed_majrev.RData")
  crm <- agk.get.compl.coef(modc[[3]],"group")
  names(crm) <- c("Intercept","betagain","betaloss","group","subject")
  crm <- crm[,-4]
  crm$lambda <- crm$betaloss*(-1)/crm$betagain
  data.la <- merge(data.la, crm, by.x = "subject", by.y = "subject", all.x = TRUE)
} else if (param_df == 'cfm'){
  setwd(path_model_calc)
  source("MLE_la_mu.R")
  crm        = all_params
  crm$lambda = crm$beta_loss*(-1)/crm$beta_gain
  crm$group  = as.factor(agk.recode.c(as.character(crm$VPNR),as.character(data_pdt$subject),as.character(data_pdt$group)))
  names(crm) = c("subject","mu","Intercept","betagain","betaloss","lambda","group")
  crm <- crm[,-7]
  data.la <- merge(data.la, crm, by.x = "subject", by.y = "subject", all.x = TRUE)
} else if (param_df == 'cfml') {
  # case if lmlist models used
  lae_lmlist = lmList(accept_reject ~ gain + loss + ed_abs|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)
  la_lmlist  = lmList(accept_reject ~ gain + loss|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)
  lar_lmlist = lmList(accept_reject ~ ratio_bcp|subject,data=data_pdt,family=binomial,na.action = NA,pool = F)
  cur_mod    = la_lmlist
  crm = coef(cur_mod)
  crm$lambda  = -crm$loss/crm$gain
  crm$group   = agk.recode.c(row.names(crm),data_pdt$subject,data_pdt$group)
  crm$subject = row.names(crm)
  names(crm) = c("Intercept","betagain","betaloss","lambda","group","subject")
  crm <- crm[,-5]
  data.la <- merge(data.la, crm, by.x = "subject", by.y = "subject", all.x = TRUE)
}

## saving the data for export to matlab
setwd(path_data)
data.la.aggr           <- aggregate(data.la, by = list(data.la$subject),FUN="first")
data.la.aggr           <- data.la.aggr[,-1]
data.la.aggr$loglambda <- get.log(data.la.aggr$lambda)

# REPLACE GBQ MISSINGS
# overall imputation function (also used later in severity)
imput_fun = function(x) {return(mean(x,na.rm=T))}
data.la.aggr_PG = subset(data.la.aggr,group == "PG")
data.la.aggr_HC = subset(data.la.aggr,group == "HC")
data.la.aggr_AD = subset(data.la.aggr,group == "AD")

data.la.aggr_AD$GBQ_mean_rec[is.nan(data.la.aggr_AD$GBQ_mean_rec)] = imput_fun(data.la.aggr_AD$GBQ_mean_rec)
data.la.aggr_AD$GBQ_persi[is.nan(data.la.aggr_AD$GBQ_persi)] = imput_fun(data.la.aggr_AD$GBQ_persi)
data.la.aggr_AD$GBQ_illus[is.nan(data.la.aggr_AD$GBQ_illus)] = imput_fun(data.la.aggr_AD$GBQ_illus)

data.la.aggr = rbind(data.la.aggr_HC,data.la.aggr_PG,data.la.aggr_AD)

#write.table(data.la.aggr,file = "la_ed.txt",sep = "\t",quote = F,row.names = F)

# get the models
setwd(path_model_calc)
source('LA_comp_mod.R')
setwd(path_model_calc)
source('LA_comp_mod_cov.R')
source('LA_comp_mod_covc.R')
all_modc <- c(all_mod,all_modc,all_modcc)

# prep data
data.la$Age_bcp = data.la$Age
data.la$Bildungsjahre_ges_bcp = data.la$Bildungsjahre_ges
data.la$Age <- scale(data.la$Age,center=T,scale=F)
data.la$Bildungsjahre_ges <- scale(data.la$Bildungsjahre_ges,center=T,scale=F)

if (export_to_matlab) {
  path_model_exp <- paste0("C:\\Users\\", user, "\\Google Drive\\Library\\MATLAB\\LA\\SS_Analysis\\LA_model")
  setwd(path_model_exp)
  writeMat(x=crm,con = "la_model.mat")
  # for value-based baseline based analysis
  # feHC = fixef(modc[[3]])
  # feHC = feHC[c(1,2,3,4)]
  # writeMat(x=feHC,con = "la_model_HC.mat")
}


## DEMOGRAPHICS FOR DEMOGR TABLE ==============================================
# # some demographics
# #contrasts(data.la.aggr$group) <- matrix(c(0,-1,1,1,-1,0),ncol=2,dimnames = list(c("AD","HC","PG"),c("PG>HC","AD>HC")))
# contrasts(data.la.aggr$group) <- contr.treatment(n=3,base=2)
# cur_var    = data.la.aggr$BIS_own_mean_BIS11
# 
# #replacing missings by median
# cur_var_HC = cur_var[data.la.aggr$group =="HC"]
# cur_var_PG = cur_var[data.la.aggr$group =="PG"]
# cur_var_AD = cur_var[data.la.aggr$group =="AD"]
# cur_med    = median(cur_var_HC,na.rm = T)
# cur_var_HC[is.na(cur_var_HC)] = cur_med
# cur_med    = median(cur_var_HC,na.rm = T)
# cur_var_PG[is.na(cur_var_PG)] = cur_med
# cur_med    = median(cur_var_HC,na.rm = T)
# cur_var_AD[is.na(cur_var_AD)] = cur_med
# cur_var    = c(cur_var_HC,cur_var_PG,cur_var_AD)
# 
# # modeling
# mod = lm(cur_var~data.la.aggr$group)
# shapiro.test(as.numeric(resid(mod)))
# summary(mod)
# describeBy(cur_var,data.la.aggr$group)
# xtabs(is.na(cur_var) ~ data.la.aggr$group)
# xtabs(is.nan(cur_var) ~ data.la.aggr$group)
# bartlett.test(cur_var,data.la.aggr$group)
# kruskal.test(cur_var,data.la.aggr$group)
# 
# # fisher
# cur_tab= xtabs( ~ cur_var+data.la.aggr$group)
# fisher.test(cur_tab)
# 
# # modeling two-sample tests
# data.la.aggr_HCAD = subset(data.la.aggr, (group == "HC" | group == "AD"))
# cur_var = data.la.aggr_HCAD$Schulden_amount
# cur_med = median(cur_var,na.rm = T)
# cur_var[is.na(cur_var)] = cur_med
# data.la.aggr_HCAD$group = as.factor(as.character(data.la.aggr_HCAD$group))  
# mod = lm(cur_var~data.la.aggr_HCAD$group)
# shapiro.test(as.numeric(resid(mod)))
# summary(mod)
# kruskal.test(cur_var,data.la.aggr_HCAD$group)
# cur_tab= xtabs( ~ cur_var+data.la.aggr_HCAD$group)
# fisher.test(cur_tab)
# 
# data.la.aggr_HCPG = subset(data.la.aggr, (group == "HC" | group == "PG"))
# cur_var = data.la.aggr_HCPG$Schulden_amount
# cur_med = median(cur_var,na.rm = T)
# cur_var[is.na(cur_var)] = cur_med
# data.la.aggr_HCPG$group = as.factor(as.character(data.la.aggr_HCPG$group))  
# mod = lm(cur_var~data.la.aggr_HCPG$group)
# shapiro.test(as.numeric(resid(mod)))
# summary(mod)
# kruskal.test(cur_var,data.la.aggr_HCPG$group)
# cur_tab= xtabs( ~ cur_var+data.la.aggr_HCPG$group)
# fisher.test(cur_tab)
# 
# 
# data.la.aggr_ADPG = subset(data.la.aggr, (group == "PG" | group == "AD"))
# cur_var = data.la.aggr_ADPG$Schulden_amount
# cur_med = median(cur_var,na.rm = T)
# cur_var[is.na(cur_var)] = cur_med
# data.la.aggr_ADPG$group = as.factor(as.character(data.la.aggr_ADPG$group))  
# mod = lm(cur_var~data.la.aggr_ADPG$group)
# shapiro.test(as.numeric(resid(mod)))
# summary(mod)
# kruskal.test(cur_var,data.la.aggr_ADPG$group)
# cur_tab= xtabs( ~ cur_var+data.la.aggr_ADPG$group)
# fisher.test(cur_tab)

## for Guillaume ==============================================================
# 
# # BGG study
# setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/OFC_project_Guillaume/data_sent_to_Guillaume_Sescousse/Data/BGG/Behav')
# cur_dat = xlsx::read.xlsx("BGG_01_all_behav_ed_2016_12_05.xlsx",1)
# desired_vars=c('subject','BIS_own_sum_BIS11','BIS_own_mean_BIS11')
# cur_dat_match = data.la.aggr[desired_vars]
# cur_dat=merge(cur_dat,cur_dat_match,by.x="ID",by.y="subject",all.x = T,all.y = F)
# 
# xlsx::write.xlsx(cur_dat,file="BGG_01_all_behav_ed_2016_12_21.xlsx")
# 
# # compare BGG study and VPPG study
# setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/OFC_project_Guillaume/data_sent_to_Guillaume_Sescousse/Data/BGG/Behav')
# cur_dat_BGG  = xlsx::read.xlsx("BGG_01_all_behav_ed_2016_11_25.xlsx",1)
# setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/OFC_project_Guillaume/data_sent_to_Guillaume_Sescousse/Data/VPPG/Behav')
# cur_dat_VPPG = xlsx::read.xlsx("VPPG_01_all_behav_AG_ed_2016_11_25.xlsx",1)
# cur_dat_VPPG$study = "VPPG"
# cur_dat_BGG$study = "BGG"
# desired_vars = c("ID","Age","Sex","KFG_Score","GBQ_sum","study","Group")
# cur_dat_BGG  = cur_dat_BGG[desired_vars]
# cur_dat_VPPG = cur_dat_VPPG[desired_vars]
# cur_dat = rbind(cur_dat_BGG,cur_dat_VPPG)
# 
# summary(lm(GBQ_sum ~ study,data=cur_dat))
# 
# ## HC
# cur_dat_HC = subset(cur_dat, Group == "control")
# summary(lm(GBQ_sum ~ study,data=cur_dat_HC))
# 
# cur_dat_HC$KFG_Score = as.numeric(cur_dat_HC$KFG_Score)
# mod_00 = lm(GBQ_sum ~ KFG_Score,data=cur_dat_HC)
# mod_01 = lm(GBQ_sum ~ KFG_Score + study,data=cur_dat_HC)
# anova(mod_00,mod_01)
# summary(mod_01)
# 
# summary(lm(KFG_Score ~ study,data=cur_dat_HC))
# kruskal.test(cur_dat_HC$KFG_Score,as.factor(cur_dat_HC$study))
# 
# ## PG
# cur_dat_PG = subset(cur_dat, Group == "PG")
# summary(lm(GBQ_sum ~ study,data=cur_dat_PG))
# 
# cur_dat_PG$KFG_Score = as.numeric(cur_dat_PG$KFG_Score)
# mod_00 = lm(GBQ_sum ~ KFG_Score,data=cur_dat_PG)
# mod_01 = lm(GBQ_sum ~ KFG_Score + study,data=cur_dat_PG)
# anova(mod_00,mod_01)
# summary(mod_01)
# 
# summary(lm(KFG_Score ~ study,data=cur_dat_PG))
# kruskal.test(cur_dat_PG$KFG_Score,as.factor(cur_dat_PG$study))
# 
# # compare BGG study and VPPG study
# setwd('S:/AG/AG-Spielsucht2/Daten/VPPG_Daten/OFC_project_Guillaume/data_sent_to_Guillaume_Sescousse/Data/VPPG/Behav')
# cur_dat_VPPG = xlsx::read.xlsx("VPPG_01_all_behav_AG_ed_2016_11_25.xlsx",1)
# cur_dat_VPPG$KFG_Score = as.numeric(as.character(cur_dat_VPPG$KFG_Score))
# cor.test(cur_dat_VPPG$KFG_Score,cur_dat_VPPG$SOGS)
# plot(cur_dat_VPPG$KFG_Score,cur_dat_VPPG$SOGS)
# 
# cur_dat_VPPG_PG = subset(cur_dat_VPPG, Group == "PG")
# cor.test(cur_dat_VPPG_PG$KFG_Score,cur_dat_VPPG_PG$SOGS)
# cor.test(cur_dat_VPPG_PG$KFG_Score,cur_dat_VPPG_PG$SOGS,method="spearman")
# plot(cur_dat_VPPG_PG$KFG_Score,cur_dat_VPPG_PG$SOGS)


## ANALYSIS ===================================================================

if (run_ana) {
  # exclude too short trials
  data.la_pr = data.la[data.la$RT >=0.5,]
  # glmer model
  bl_mod           = glmer(accept.reject ~ 1 + (1 | subject), data = data.la,family="binomial") # bl
  la_mod           = glmer(accept.reject ~ (Gewinn+Verlust) + (1| subject),data = data.la,family = "binomial") # la model
  anova(bl_mod,la_mod)
  lae_mod          = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs) + (1| subject),data = data.la,family = "binomial") # la model
  
  lae_mod          = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs)*group + (Gewinn+Verlust+ed.abs| subject),data = data.la,family = "binomial")    # la model
  AIC(lae_mod)/51
  la_mod           = glmer(accept.reject ~ (Gewinn+Verlust)*group + (Gewinn+Verlust| subject),data = data.la,family = "binomial")    # la model
  AIC(la_mod)/51
  lar_mod          = glmer(accept.reject ~ (ratio_bcp)*group + (ratio_bcp| subject),data = data.la,family = "binomial")    # la model
  AIC(lar_mod)/51
  lac_mod          = glmer(accept.reject ~ (0+Gewinn+Verlust)*group + (0+Gewinn+Verlust| subject),data = data.la,family = "binomial")    # la model
  
  lac_mod_AD       = glmer(accept.reject ~ (0+Gewinn+Verlust) + (0+Gewinn+Verlust| subject),data = data.la.AD,family = "binomial") # la model
  lac_mod_PG       = glmer(accept.reject ~ (0+Gewinn+Verlust) + (0+Gewinn+Verlust| subject),data = data.la.PG,family = "binomial") # la model
  lac_mod_HC       = glmer(accept.reject ~ (0+Gewinn+Verlust) + (0+Gewinn+Verlust| subject),data = data.la.HC,family = "binomial") # la model
  AIC(lac_mod)/51
  AIC(lac_mod_AD)/15
  AIC(lac_mod_HC)/17
  AIC(lac_mod_PG)/19
  mean(c(AIC(lac_mod_AD)/15,AIC(lac_mod_HC)/17,AIC(lac_mod_PG)/19))
  
  anova(lae_mod,la_mod) # hardly any effect
  la_mod_ranef     = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs) + (Gewinn+Verlust+ed.abs| subject),data = data.la,family = "binomial") # la model
  anova(la_mod,la_mod_ranef)
  la_modg_ranefg   = glmer(accept.reject ~ (Gewinn+Verlust) + Verlust*group + (Gewinn| subject),data = data.la_pr,family = "binomial") # la model
  la_modg_ranefv   = glmer(accept.reject ~ (Gewinn+Verlust) + Verlust*group + (Verlust| subject),data = data.la_pr,family = "binomial") # la model
  la_modg_ranef    = glmer(accept.reject ~ (Gewinn+Verlust) + Verlust*group + (Gewinn+Verlust| subject),data = data.la_pr,family = "binomial") # la model
  la_mod_ranef     = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs) + (Gewinn+Verlust+ed.abs| subject),data = data.la,family = "binomial") # la model
  la_modg_ranef    = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs)*group + (Gewinn+Verlust+ed.abs| subject),data = data.la,family = "binomial") # la model
  anova(la_mod_ranef,la_modg_ranef)
  anova(la_mod_ranef,la_modg_ranefg)
  anova(la_mod_ranef,la_modg_ranefv)
  summary(la_modg_ranef)
  
  cur_fe = as.data.frame(t(as.matrix(fixef(la_modg_ranef))))
  cur_LA_HC = -cur_fe$Verlust/cur_fe$Gewinn
  cur_LA_PG = -((cur_fe$Verlust+cur_fe$'Verlust:groupPG>HC')/(cur_fe$Gewinn))
  cur_LA_AD = -((cur_fe$Verlust+cur_fe$'Verlust:groupAD>HC')/(cur_fe$Gewinn))
  
  cur_cc    = agk.get.compl.coef(la_modg_ranef,"group")
  
  # risk model
  risk_modg_ranef    = glmer(accept.reject ~ (RiskMar+ed.abs)*group + (RiskMar+ed.abs| subject),data = data.la,family = "binomial",
                             nAGQ = 0,control=glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=300)))
  evrisk_modg_ranef  = glmer(accept.reject ~ (EV+RiskMar)*group + (EV+RiskMar| subject),data = data.la,family = "binomial",
                             nAGQ = 0,control=glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=300)))
  evriske_modg_ranef = glmer(accept.reject ~ (EV+RiskMar+ed.abs)*group + (EV+RiskMar+ed.abs| subject),data = data.la,family = "binomial",
                             nAGQ = 0,control=glmerControl(check.conv.grad='ignore',check.conv.singular='ignore',check.conv.hess='ignore',optCtrl=list(optimizer = 'nloptwrap',maxfun=300)))
  anova(risk_modg_ranef,modc[[3]])
  anova(evrisk_modg_ranef,modc[[3]])
  anova(evriske_modg_ranef,modc[[3]])
  
  # how much did they win?
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  lmp = cmpfun(lmp)
  all_ps = rep(NA,150)
  for (kk in 1:150) {
    all_subs = levels(data.la_unagg$subject)
    win      = c()
    for(ii in 1:length(all_subs)) {
      cur_dat = subset(data.la_unagg, subject == all_subs[ii])
      cur_dat = subset(cur_dat, accept.reject == 1)
      cur_sam = sample(1:length(cur_dat[,1]),5)
      cur_dat = cur_dat[cur_sam,]
      cur_out = c()
      for (jj in 1:length(cur_dat[,1])) {
        cur_out[jj] = ifelse(rand()>0.5,cur_dat$Gewinn_bcp[jj],cur_dat$Verlust_bcp[jj])
      }
      win[ii] = mean(cur_out)
    }
    
    data.la.aggr$win = win
    tmp        = lm(win~group, data.la.aggr)
    all_ps[kk] = lmp(tmp)
    disp(kk)
  }
  
  
  # estimate models: external script (cv)
  # estimate p-values: external script (bootstrapping p-values)
  # get bootstr CIs: external script (bootstrapping p-values)
  
  # plotting the fixed effects with bootstrapped CI's
  # prep
  cur_name <- paste(path_results,"boot_la_compl_mod_no_cov51_noed_majrev.RData",sep="\\")
  ne <- new.env()
  load(file=cur_name, env=ne)
  cur_boot <- eval(parse(text=ls(env=ne,pattern = "*boot*")),envir = ne)
  boot_pdt_df <- cur_boot
  boot_pdt_df <- data.frame(boot_pdt_df$t)
  boot_pdt_df <- na.omit(boot_pdt_df)
  names(boot_pdt_df)[1:3] <- c("in_AD","in_HC","in_PG")
  boot_pdt_df = boot_pdt_df[,c(2,3,1,5,6,4,8,9,7,11,12,10,13,14,15)]
  
  p1 <- agk.barplot.boot(boot_pdt_df,cur_n = 20,cur_cat = "group",nsim = "20",study_name = "BGG")
  p1 <- p1 + geom_segment(aes(x=1.3,xend=2.3,y=3.5,yend=3.5))
  p1 <- p1 + geom_text(aes(label="*",x=1.75,y=3.6),size=10)
  p1 <- p1 + scale_fill_manual(labels=c(expression(beta[gain]), expression(beta[loss]), expression(lambda)),
                               values=cbbPalette,
                               name="fixed effect")
  p1 <- p1 + ggtitle("Differences in loss aversion by group\n")
  p1 + theme_la()
  
  # get the bootstrapped p-values (# check outside script for this)
  
  # compare ratio/... and la model loss aversion estimates
  la_ratio_model  <- glmer(formula = accept.reject ~ (ed.abs + ratio) * group + ((ed.abs + ratio) | subject), data = data.la, family = "binomial", control = glmerControl(optimizer = "bobyqa"))
  crm_ratio      <- agk.get.compl.coef(la_ratio_model,"group")
  setwd(path_results)
  load("cur_mod.Rdata")
  crm <- agk.get.compl.coef(cur_mod,"group") # no covs
  crm$la_ratio <- crm_ratio$ratio
  crm$lambda <- crm$Verlust*(-1)/crm$Gewinn
  describeBy(crm$lambda,crm$group)
  plot(crm$la_ratio,crm$lambda)
  plot(crm$la_ratio,crm$lambda)
  
  # ratio lm list
  la_ratio_lmlist <- lmList(accept.reject ~ ed.abs + ratio | subject, data = data.la, family="binomial",na.action=NULL,pool = F)
  coef_ratio <- coef(la_ratio_lmlist)
  crm$la_ratio_lml <- coef_ratio$ratio
  plot(crm$la_ratio_lml,log(crm$lambda-min(crm$lambda)+0.01))
  
  # diff lm list
  la_diff_lmlist <- lmList(accept.reject ~ ed.abs + diff | subject, data = data.la, family="binomial",na.action=NULL,pool = F)
  coef_diff <- coef(la_diff_lmlist)
  crm$la_diff_lml <- coef_diff$diff
  plot(crm$la_diff_lml,log(crm$lambda-min(crm$lambda)+0.01))
  cor.test(crm$la_diff_lml,log(crm$lambda-min(crm$lambda)+0.01),method = "spearman")
  
  # ev lm list
  la_EV_lmlist <- lmList(accept.reject ~ ed.abs + EV | subject, data = data.la, family="binomial",na.action=NULL,pool = F)
  coef_EV <- coef(la_EV_lmlist)
  crm$la_EV_lml <- coef_EV$EV
  plot(crm$la_EV_lml,log(crm$lambda-min(crm$lambda)+0.01))
  cor.test(crm$la_EV_lml,log(crm$lambda-min(crm$lambda)+0.01),method = "spearman")
  
  # gain+loss lm list
  la_list = lmList(accept.reject ~ (Gewinn+Verlust) | subject,data = data.la,family = "binomial",pool = F, na.action=NULL,method="ML")
  bl_list = lmList(accept.reject ~ 1 | subject,data = data.la,family = "binomial",pool = F, na.action=NULL,method="ML")
  bad_subs = c()
  all_subs = names(la_list)
  for (ii in 1:length(all_subs)) {
    cur_mod = la_list[[ii]]
    if (is.null(cur_mod)) {next}
    cur_bl  = bl_list[[ii]]
    cur_sum = summary(cur_mod)$coefficients
    cur_anova = anova(cur_bl,cur_mod,test = "Chisq")
    if (cur_sum[2,4] >= 0.05 & cur_sum[2,4] >= 0.05) {
      bad_subs = c(bad_subs,all_subs[ii])
      next
    }
    if (((cur_anova$'Pr(>Chi)')[2] <0.05) == F) {
      bad_subs = c(bad_subs,all_subs[ii])
      next
    }
    if ((AIC(cur_bl) - AIC(cur_mod))<=0) {
      bad_subs = c(bad_subs,all_subs[ii])
      next
    }
    if ((BIC(cur_bl) - BIC(cur_mod))<=0) {
      bad_subs = c(bad_subs,all_subs[ii])
      next
    }
  }
  
  ## reaction times ===========================================================
  
  data.la$group = factor(data.la$group, levels=c("AD","HC","PG"))
  data.la$group = factor(data.la$group, levels=c("HC","PG","AD"))
  mod00 = lmer(RT ~ Gewinn + Verlust + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML = FALSE)
  
  mod01 = lmer(RT ~ (Gewinn + Verlust)*group + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML=FALSE)
  mod02 = lmer(RT ~ (Gewinn + Verlust + ed.abs)*group + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML=FALSE)
  
  mod03 = lmer(RT ~ Gewinn + Verlust + Age + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML = FALSE)
  mod04 = lmer(RT ~ (Gewinn + Verlust + Age)*group + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML=FALSE)
  
  mod05 = lmer(RT ~ (Gewinn + Verlust)*Age + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML = FALSE)
  mod06 = lmer(RT ~ (Gewinn + Verlust)*Age*group + ed.abs + (Gewinn + Verlust + ed.abs|subject), data = data.la,REML=FALSE)
  
  mod00p = glmer(RT ~ Gewinn + Verlust + ed.abs + (Gewinn + Verlust + ed.abs|subject), family = "poisson", data = data.la,REML = FALSE,
                 control=glmerControl(check.conv.grad="ignore",check.conv.singular="ignore",
                                      check.conv.hess="ignore",optCtrl=list(optimizer = "nloptwrap",maxfun=10)))
  mod01p = glmer(RT ~ (Gewinn + Verlust)*group + ed.abs + (Gewinn + Verlust + ed.abs|subject), family = "poisson",  data = data.la,REML=FALSE)
  
  anova(mod00,mod01)
  anova(mod00,mod02)
  anova(mod03,mod04)
  anova(mod05,mod06)
  summary(mod01)
  summary(mod04)
  summary(mod06)
  fe = fixef(mod01)
  agk.normality.check(mod01)
  
  # make a plot of this
  mean_RT_HC  = fe['(Intercept)']
  mean_RT_PG  = mean_RT_HC + fe['groupPG>HC']
  mean_RT_AD  = mean_RT_HC + fe['groupAD>HC']
  incr_RT_HCg = mean_RT_HC + fe['Gewinn'] 
  incr_RT_PGg = mean_RT_PG + fe['Gewinn:groupPG>HC'] 
  incr_RT_ADg = mean_RT_AD + fe['Gewinn:groupAD>HC'] 
  incr_RT_HCl = mean_RT_HC + fe['Verlust'] 
  incr_RT_PGl = mean_RT_PG + fe['Verlust:groupPG>HC'] 
  incr_RT_ADl = mean_RT_AD + fe['Verlust:groupAD>HC'] 
  
  cur_b = c(incr_RT_HCl,mean_RT_HC,incr_RT_HCg,0,incr_RT_PGl,mean_RT_PG,incr_RT_PGg,0,incr_RT_ADl,mean_RT_AD,incr_RT_ADg)
  
  # prepare data_frame
  group = c("HC", "HC", "HC","PG", "PG", "PG","AD", "AD", "AD")
  cond  = c("byloss","mean","bygain","byloss","mean","bygain","byloss","mean","bygain")
  value = cur_b
  
  # using lmlist for RT
  rtlist = lmList(RT ~ Gewinn + Verlust + ed.abs | subject, data = data.la,pool =F, na.action=F)
  rtcoef = coef(rtlist)
  rtcoef$subject = row.names(rtcoef)
  rtcoef$group = agk.recode(as.character(rtcoef$subject),as.character(data.la.aggr$subject),as.character(data.la.aggr$group))
  names(rtcoef)[1] = "Intercept"
  rtcoef$Gewinn = rtcoef$Gewinn*5 + rtcoef$Intercept # 5 Euros up
  rtcoef$Verlust = rtcoef$Verlust*5 + rtcoef$Intercept # 5 Euro down
  rtcoef$ed.abs = NULL # controlled for ed.abs; but now taken out
  names(rtcoef)[1] = "Mean"
  
  # stats tests
  rtcoef$group = factor(rtcoef$group,levels=c("HC","PG","AD"))
  summary(lm(rtcoef$Mean ~ rtcoef$group))
  
  rtcoef_PG = subset(rtcoef, group == "PG")
  t.test(rtcoef_PG$Mean,rtcoef_PG$Gewinn,paired = T,alternative = "greater")
  t.test(rtcoef_PG$Mean,rtcoef_PG$Verlust,paired = T,alternative = "less")
  
  rtcoef_AD = subset(rtcoef, group == "AD")
  t.test(rtcoef_AD$Mean,rtcoef_AD$Gewinn,paired = T,alternative = "greater")
  t.test(rtcoef_AD$Mean,rtcoef_AD$Verlust,paired = T,alternative = "less")
  
  
  # aggregate
  meanRT   = aggregate(rtcoef$Mean,by=list(rtcoef$group),FUN=agk.boot.ci,lower=0.025,upper=0.975,R=1000,cur_fun=median)
  meanRT   = data.frame(meanRT$Group.1,as.data.frame(meanRT$x))
  names(meanRT) = c("group","mean","lower","upper")
  bygain   = aggregate(rtcoef$Gewinn,by=list(rtcoef$group),FUN=agk.boot.ci,lower=0.025,upper=0.975,R=1000,cur_fun=median)
  bygain   = data.frame(bygain$Group.1,as.data.frame(bygain$x))
  names(bygain) = c("group","mean","lower","upper")
  byloss   = aggregate(rtcoef$Verlust,by=list(rtcoef$group),FUN=agk.boot.ci,lower=0.025,upper=0.975,R=1000,cur_fun=median)
  byloss   = data.frame(byloss$Group.1,as.data.frame(byloss$x))
  names(byloss) = c("group","mean","lower","upper")
  cur_bar = rbind(meanRT,bygain,byloss)
  cur_bar$RT = factor(rep(c("meanRT","bygain","byloss"),each=3),levels=c("byloss","meanRT","bygain"))
  cur_bar$group  = factor(cur_bar$group,levels=c("HC","PG","AD"))
  
  # plot RT
  p = ggplot(data = cur_bar, aes(group,mean,fill=RT))
  p = p+geom_bar(position="dodge",stat="identity")
  p = p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                        position=position_dodge(width=0.9), color=cbbPalette[4],
                        width=0.1) + ylab("mean (s) (95% boots. CI)\n")
  p = p + scale_fill_manual(labels=c("byloss", "meanRT", "bygain"),
                            values=cbbPalette,
                            name="Reaction Time")
  p = p + theme_la()
  p = p + coord_cartesian(ylim = c(0, 1.95)) 
  p = p + ggtitle("Reaction Times\n")
  p
  
  
  
  query_dat     = rtlist[[1]]$model[1,]
  query_dat$Gewinn  = query_dat$Gewinn-query_dat$Gewinn
  query_dat$Verlust = query_dat$Verlust-query_dat$Verlust
  query_dat$ed.abs  = query_dat$ed.abs-query_dat$ed.abs
  query_dat_int = query_dat
  query_dat_gin = query_dat
  query_dat_lss = query_dat
  
  predict(rtlist[[1]],newdata = query_dat_int,interval="confidence")
  
  ## correlations with severity ===============================================
  setwd(path_results)
  load("modc51_noed_majrev.RData")
  cur_mod <- modc[[3]]
  crm <- agk.get.compl.coef(cur_mod,"group") # no covs
  tmp <- names(crm)
  tmp[1] <- "Intercept"
  names(crm) <- tmp
  
  crm$lambda <- crm$Verlust*(-1)/crm$Gewinn
  describeBy(crm$lambda,crm$group)
  
  #desired_vars <- c("GBQ_mean_rec","PGYBOCS_Summe","KFG_Summe","GSAS_Summe")
  desired_vars <- c("GBQ_mean_rec","PGYBOCS_Summe","KFG_Summe","GSAS_Summe","GBQ_persi","GBQ_illus","ADS","OCDS_Summe")
  #desired_vars <- c("ADS","OCDS_Summe")
  
  #des_vars_df = data.la.aggr
  des_vars_df <- aggregate(data.la[desired_vars],by=list(data.la$subject),FUN = "mean")
  names(des_vars_df) <- c("subject",desired_vars)
  des_vars_df$group <- crm$group
  
  ## impute values instead of missing values
  des_vars_df_HC = subset(des_vars_df,group == "HC")
  des_vars_df_PG = subset(des_vars_df,group == "PG")
  des_vars_df_AD = subset(des_vars_df,group == "AD")
  all_groups_df = list()
  all_groups_df[[1]]=des_vars_df_HC
  all_groups_df[[2]]=des_vars_df_PG
  all_groups_df[[3]]=des_vars_df_AD
  
  for (ii in 1:length(all_groups_df)) {
    cur_df = all_groups_df[[ii]]
    for (jj in 1:length(cur_df)) {
      if(sum(is.na(cur_df[[jj]]))) {
        cur_df[[jj]][is.na(cur_df[[jj]])] = imput_fun(cur_df[[jj]])
      } else if (sum(is.nan(cur_df[[jj]]))) {
        cur_df[[jj]][is.nan(cur_df[[jj]])] = imput_fun(cur_df[[jj]])
      }
    }
    all_groups_df[[ii]] = cur_df
  }
  des_vars_df = rbind(all_groups_df[[1]],all_groups_df[[2]],all_groups_df[[3]])
  
  
  cur_cor_info <- list()
  cor_plots    <- list()
  cur_levels   <- levels(data.la$group)
  for (ii in 1:length(cur_levels)) {
    cur_data <- subset(crm,group == (levels(data.la$group))[ii])
    #cur_data <- agk.scale.ifpossible(subset(crm,group == (levels(data.la$group))[ii]))
    #cur_des  <- agk.scale.ifpossible(subset(des_vars_df,group == (levels(data.la$group))[ii]))
    cur_des  <- subset(des_vars_df,group == (levels(data.la$group))[ii])
    tmp_2 <- get.log(cur_data$lambda)
    cur_data$lambda <- tmp_2
    tmp <-  corr.test(cur_data[c("lambda")],cur_des[desired_vars],method = "spearman",adjust = "none")
    cur_cor_info[[ii]] <- tmp
    cor_plots[[ii]] <- agk.extract.sig.corr(tmp,cur_data[c("Gewinn","Verlust","lambda","Intercept")],cur_des[desired_vars],thresh=0.1)  
  }
  
  # bootstrapping rho-p-value
  x=cur_data$lambda
  y=cur_des$GBQ_illus
  rho <- function(x, y, indices){
    rho <- cor.test(x[indices], y[indices],  method = c("spearman"))
    return(rho$estimate)
  }
  library(boot)    
  boot.rho <- boot(x ,y=y, rho, R=20000)
  boot.ci(boot.rho,type = 'basic',conf = c(0.999999))
  
  rhodens   = density(boot.rho$t,bw = c("ucv"))
  rhodens$y = cumsum(rhodens$y/sum(rhodens$y))
  f <- approxfun(rhodens, rule=2)
  1-f(0)
  
  # plot
  all_plots <- list()
  cur_count <- 0
  for (ii in 1:length(cor_plots)) {
    if (length(cor_plots[[ii]]) == 0) {next}
    for (jj in 1:length(cor_plots[[ii]])) {
      cur_count <- cur_count + 1
      tmp <- cor_plots[[ii]][[jj]]
      cur_strings <- names(tmp)
      all_plots[[cur_count]] <- ggplot(tmp, aes_string(x=cur_strings[1],y=cur_strings[2])) + geom_point(size =5) + theme_la() + xlab(label = paste0(cur_strings[1],'\n'))+ylab(label = paste0(cur_strings[2],'\n'))
      all_plots[[cur_count]] <- all_plots[[cur_count]] + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) + ggtitle(paste(cur_levels[ii],"group\n"))
      all_plots[[cur_count]]
    }
  }
  all_plots
  
  # export plots for publication
  all_plots[[1]] <- all_plots[[1]] + ylab(label = "GBQ (z-score)\n")
  all_plots[[1]] <- all_plots[[1]] + xlab(label = expression(paste("loss aversion (log(",lambda,"))",sep="")))
  all_plots[[1]] 
  #all_plots[[3]] <- all_plots[[3]] + xlab(label = expression(paste("acceptance (int.)")))
  #multiplot(plotlist = all_plots,cols = 3)
  
  # plot the PG gPPI GBQ_persi relationship
  setwd(path_results)
  tmp = read.table(file="GBQ_persi_gPPI_PG.csv",sep=";",header=T,dec=",")
  tmp$gPPI <- as.numeric(tmp$gPPI)
  tmp$GBQpersi <- as.numeric(tmp$GBQpersi)
  p <- ggplot(tmp, aes(x=GBQpersi,y=gPPI)) + geom_point(size =5) + theme_la() + xlab(label = "GBQ persistence") + ylab(label = "functional connectivity\nL DLPFC to L caudate")
  p <- p + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) + ggtitle("PG group")
  p
  
  ## cor GBQ whole data set
  crm$lambda = get.log(crm$lambda)
  cor.test(crm$lambda,des_vars_df$GBQ_persi)
  
  ## CORRELATION LIKE IN PAPER (ONE-SIDED, THE LOWER LA THE HIGHER SEVERITY)
  library(boot)  
  
  myrho <- function(cur_d, indices){
    # boot sample
    cur_d = cur_d[indices,]
    x = cur_d$x
    y = cur_d$y
    # corr function
    rho <- cor.test(x, y,  method = 'pearson')
    cur_est = rho$estimate
    if (is.na(cur_est)) {
      # NA will be 0 (conservative)
      cur_est = 0
    }
    return(cur_est)
  }
  
  agk.corr.severity.la.single.p = function(x,y,R=10000) {
    if (any(is.na(y))) {
      return(NaN)
    }
    # bootstrap the p-value for r, one-sided, data is x,y
    
    boot.rho <- boot(data.frame(x,y), statistic = myrho, R=R)
    
    cur_t = boot.rho$t
    
    rhodens   = density(cur_t,bw = c("ucv"))
    rhodens$y = cumsum(rhodens$y/sum(rhodens$y))
    f <- approxfun(rhodens, rule=2)
    return(1-f(0))
  }
  
  # use the bootstrap function in all groups for all desired vars
  cor_boot = list()
  cur_levels   <- levels(data.la$group)
  ct = 0
  for (ii in 1:length(cur_levels)) {
    cur_data <- subset(crm,group == (levels(data.la$group))[ii])
    cur_des  <- subset(des_vars_df,group == (levels(data.la$group))[ii])
    tmp_2 <- get.log(cur_data$lambda)
    cur_data$lambda <- tmp_2
    for (dv in desired_vars) {
      # get the bootstr corr
      cur_p = agk.corr.severity.la.single.p(cur_data$lambda,cur_des[[dv]],R=1000)
      res = c((levels(data.la$group))[ii],dv,cur_p)
      print(res)
      ct = ct + 1
      cor_boot[[ct]] = res
    }
  }
  
  # reduce results to excluded the NAN
  cur_fun = function(x) {return(x[3])}
  all_ps = as.numeric(lapply(cor_boot,cur_fun))
  cor_boot = cor_boot[all_ps != 'NaN']
  
  # get the fdr p-value:
  all_ps = as.numeric(lapply(cor_boot,cur_fun))
  all_ps_fdr = p.adjust(all_ps,'bonferroni')
  
  # what corr won?
  print(cor_boot[all_ps_fdr < 0.05])
  
  ## heatmap; standardize! ####################################################
  # per group per subject; 2 step approach
  # get the unaggr data
  data.la.bcp = data.la
  rm(data.la)
  des_agg     = 1
  # get all the data in long format
  setwd(path_scripts)
  source("get_data_la_2.R")
  data.la_unagg = data.la
  data.la       = data.la.bcp
  des_agg       = 3
  hmdf <- aggregate(as.numeric(as.character(data.la_unagg$accept.reject)),by=list(data.la_unagg$subject, data.la_unagg$group,data.la_unagg$Gewinn,data.la_unagg$Verlust),FUN="mean.rmna")
  names(hmdf) <- c("subject", "group","Gewinn_centered","Verlust_centered","PoA")
  
  # # entropy per subject
  # f <- function(x) {
  #   cur_entr <- agk.entr.st(x)
  #   return(cur_entr)
  # }
  # 
  # tmp <- aggregate(hmdf$PoA,by = list(hmdf$subject),FUN=f)
  # pred_perc$entr <- tmp[,2]
  
  # get a group specific heatmap by agregating over subs
  tmp <- aggregate(hmdf$PoA,by = list(hmdf$group,hmdf$Gewinn_centered,hmdf$Verlust_centered),FUN=mean)
  names(tmp) <- c("group","gain", "loss","PoA")
  tmp$loss <- as.factor((tmp$loss + mean(seq(7,18)))*(-1))
  tmp$gain <-  as.factor(tmp$gain  + mean(seq(14,36,by = 2)))
  tmp$group = factor(tmp$group,levels=c("HC","PG","AD"),labels=c("HC","PG","AD")) 
  
  # heatmap
  p_hm <- qplot(x=gain, y=loss, data=tmp, fill=PoA, geom="tile")
  #p_hm <- p_hm + scale_fill_gradientn(colours=rainbow(4))
  p_hm <- p_hm + scale_fill_gradientn(colours=c(Palette_heat[1], Palette_heat[3]))
  p_hm <- p_hm + facet_grid(.~group)
  p_hm <- p_hm + scale_x_discrete(expand=c(0,0),breaks=c(14,36))+ scale_y_discrete(expand=c(0,0),breaks=c(-18,-7))
  p_hm <- p_hm + labs(x = "possible gain (euros)", y = "possible loss (euros)")
  p_hm  = p_hm + theme_la()
  p_neu = (p_hm + theme_update(strip.text.x = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                               strip.text.y = element_text(colour = "black",size=25,face="bold",vjust=2,margin=margin(5,5,5,5)),
                               strip.background = element_rect(fill="white", color = "white"),
                               strip.switch.pad.grid = unit(c(3),"cm"),
                               strip.switch.pad.wrap = unit(c(3),"cm"),
                               panel.margin.x =unit(0.5,"cm"),
                               panel.margin.y =unit(1,"cm"),
                               plot.title = element_text(colour = "black",size=25),
                               axis.text = element_text(colour = "black",size=25,face="bold"),
                               axis.text = element_text(colour = "black",size=25,face="bold"),
                               axis.ticks = element_line(colour ="black", size=1.5),
                               axis.title.x = element_text(colour = "black",size=25,face="bold",vjust=0.2),
                               axis.title.y = element_text(colour = "black",size=25,face="bold",vjust=2.0,angle = 90),
                               plot.margin = unit(c(1,0.5,1,1),"cm"),
                               legend.background = element_rect(fill="white", color = "white"),
                               legend.text = element_text(colour = "black",size=20,face="bold"),
                               legend.title = element_text(colour = "black",size=20,face="bold"),
                               #legend.margin=unit(.85,"cm"),
                               legend.text.align=0,
                               legend.margin=unit(0.2,"cm")))
  p_neu + guides(fill = guide_legend(title.theme = element_text(size=25, face="bold", colour = "black", angle = 0)))
  
  # heatmap only HC
  tmp <- subset(tmp, group = "HC")
  p_hm <- qplot(x=gain, y=loss, data=tmp, fill=PoA, geom="tile")
  #p_hm <- p_hm + scale_fill_gradientn(colours=rainbow(4))
  p_hm <- p_hm + scale_fill_gradientn(colours=c(cbbPalette[2], cbbPalette[3]))
  p_hm <- p_hm +scale_x_discrete(expand=c(0,0),breaks=c(14,32))+ scale_y_discrete(expand=c(0,0),breaks=c(-16,-7))
  p_hm <- p_hm + labs(x = "possible gain (???)", y = "possible loss (???)")
  p_hm + theme_la() + ggtitle("HC")
  
  # pred perc per group
  data.la$pred <- round(predict(modc[[78]],type="response",newdata = data.la))
  data.la$pred_perc <- as.numeric(data.la$pred == data.la$accept.reject)
  
  pred_perc <- aggregate(data.la$pred_perc,by=list(data.la$subject, data.la$group),FUN="mean.rmna")
  names(pred_perc) <- c("group", "subject","pred_perc")
  summary(lm(pred_perc ~ group,data=pred_perc))
  summary(lm(entr ~ group,data=pred_perc))
  crm$pred_perc <- pred_perc$pred_perc
  
  contrasts(crm$group) <- contr.treatment(n = 3,base = 2)
  colnames(contrasts(crm$group)) <- c("AD","PG")
  summary(lm(lambda ~ group, data = crm))
  summary(lm(lambda ~ pred_perc + group, data = crm))
  ## also possible to run it as covariate in glmer model!
  
  f <- function(x) {return(mean.rmna(as.logical(x)))}
  tmp <- unlist(lapply(pred_perc_list,FUN = f))
  
  ## entr
  
  x <- agk.normalize(x)
  y <- agk.normalize(y)
  x <- x/sum(x)
  y <- y/(sum(y))
  hist(x)
  hist(y)
  
  agk.entr(x)
  agk.entr(y)
  
  
  f <- function(x) {print(x@call)}
  
  ####### MODEL AND PREDICTING KFG #####
  fm  = glmer(accept.reject ~ (Gewinn+Verlust+ed.abs)*group + (Gewinn+Verlust+ed.abs| subject),data = data.la,family = "binomial") # la model
  tmp = coef(modc[[20]])$subject
  names(tmp)[1] = "Intercept"
  tmp$subject = row.names(tmp)
  tmp$KFG     = as.numeric(agk.recode.c(as.character(tmp$subject),as.character(data.la$subject),as.character(data.la$KFG_Summe)))
  tmp$group   = as.factor(agk.recode.c(as.character(tmp$subject),as.character(data.la$subject),as.character(data.la$group)))
  tmpPG       = subset(tmp,group == "PG")
  
  tmpPGs      = agk.scale.ifpossible(tmpPG)
  lmlnet      = cvAlpha.glmnet(KFG ~ Intercept + Gewinn + Verlust + ed.abs,data=tmpPGs,nfolds=5)
  KFG_pred    = predict(lmlnet,alpha=0,s='lambda.min',newdata=tmpPGs)
  
  plot(KFG_pred,tmpPGs$KFG)
  cor.test(KFG_pred,tmpPGs_MRI$KFG)
  
  coefs = as.matrix(coef(lmlnet,alpha=1,s='lambda.min'))
  
  #######
  # OLD #
  #######
  
  # evaluate model with interaction effect
  mod_int <- paste('glmer(formula = accept.reject ~ (ed.abs + Gewinn + Verlust + Gewinn*Verlust) *', 
                   'group + ((ed.abs + Gewinn + Verlust + Gewinn*Verlust) | subject), data = data.la,', 
                   'family = "binomial", control = glmerControl(optimizer = "bobyqa"))',sep='')
  mod_int_est <- eval(parse(text=mod_int))
  
  
  ###################################################################
  
  ## estimate glmer model with pred_perc as predictor
  cur_mod   <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = data.la,family = "binomial",pool = F,na.action =NULL)
  
  # check whether all subjects have at least one beta (gain or loss) significant
  cur_est <- summary(cur_mod)
  cur_est <- cur_est$coefficients
  table(cur_est[,4,3] < 0.05 | cur_est[,4,4] < 0.05) # p-values Gewinn and p-values Verlust
  
  pred_perc <- agk.predperc.lmList.logit(cur_mod)
  data.la   <- merge(data.la, pred_perc, by.x = "subject", by.y = "subs", all.x = TRUE)
}
