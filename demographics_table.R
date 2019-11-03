# demographics table
des_vars = c("Age_bcp","Rauchen_upd","WIE_Score","Bildungsjahre_ges_bcp","Schuljahre_vp","BIS_gesamt_Mean","BDI_Summe",
             "Schulden","Schulden_amount","Einkommen_VP_merge","EHI_lat","KFG_Summe","GBQ_mean_rec","GBQ_illus","GBQ_persi",
             "PGYBOCS_Summe","OCDS_Summe","ADS")


cur_dat = data.la[c("subject","group", des_vars)]
cur_dat = aggregate(data.la[des_vars],by=list(cur_dat$group,cur_dat$subject),FUN=first)

contrasts(cur_dat$Group.1)

tmp = matrix(c(-1,1,0,0,1,-1,-1,0,1),nrow = 3,byrow=F)
colnames(tmp) = c('AD<HC','PG<HC','AD<PG')
rownames(tmp) = c('AD','HC','PG')

contrasts(cur_dat$Group.1) = tmp
cur_dat$Group.3 = cur_dat$Group.1
contrasts(cur_dat$Group.3) = tmp[1:3,2:3]

cur_dat_HCAD = cur_dat[(cur_dat$Group.1=="HC" | cur_dat$Group.1=="AD"),] 
cur_dat_HCPG = cur_dat[(cur_dat$Group.1=="HC" | cur_dat$Group.1=="PG"),] 
cur_dat_PGAD = cur_dat[(cur_dat$Group.1=="PG" | cur_dat$Group.1=="AD"),] 

cur_dat_HCAD$Group.1 = droplevels(cur_dat_HCAD$Group.1) 
cur_dat_HCPG$Group.1 = droplevels(cur_dat_HCPG$Group.1) 
cur_dat_PGAD$Group.1 = droplevels(cur_dat_PGAD$Group.1) 

####### START HERE ####################
# cur_var
cur_var      = cur_dat$Bildungsjahre_ges_bcp
cur_var_HCAD = cur_dat_HCAD$Bildungsjahre_ges_bcp
cur_var_HCPG = cur_dat_HCPG$Bildungsjahre_ges_bcp
cur_var_PGAD = cur_dat_PGAD$Bildungsjahre_ges_bcp

# getting the descriptives
median.rmna = function(x) {return(median(x,na.rm = T))}
tmp = aggregate(cur_var,by=list(cur_dat$Group.1),FUN=agk.boot.ci,lower=0.025,upper=0.975,R=1000,fun=median.rmna)
describeBy(cur_var,cur_dat$Group.1)

xtabs(~Schulden+Group.1,data=cur_dat)

# getting the statistics
mod1= lm(cur_var ~ cur_dat$Group.1)
agk.normality.check(mod1)
shapiro.test(resid(mod1))
summary(mod1)

t.test(cur_var[cur_dat$Group.1=="AD"],cur_var[cur_dat$Group.1=="HC"],var.equal = F)
t.test(cur_var[cur_dat$Group.1=="PG"],cur_var[cur_dat$Group.1=="HC"],var.equal = F)
t.test(cur_var[cur_dat$Group.1=="PG"],cur_var[cur_dat$Group.1=="AD"],var.equal = F)

bartlett.test(cur_var,cur_dat$Group.1)

kruskal.test(cur_var,cur_dat$Group.1)
kruskal.test(cur_var_HCAD,cur_dat_HCAD$Group.1)
kruskal.test(cur_var_HCPG,cur_dat_HCPG$Group.1)
kruskal.test(cur_var_PGAD,cur_dat_PGAD$Group.1)

wilcox.test(cur_var_PGAD[cur_dat_PGAD$Group.1 == "PG"],cur_var_PGAD[cur_dat_PGAD$Group.1 == "AD"],paired = F,correct = F)

fisher.test((xtabs(~Schulden+Group.1,data=cur_dat)))
fisher.test(xtabs(~Schulden+Group.1,data=cur_dat_HCAD))
fisher.test(xtabs(~Schulden+Group.1,data=cur_dat_HCPG))
fisher.test(xtabs(~Schulden+Group.1,data=cur_dat_PGAD))


# bootstrapping p-value
cur_var_df = data.frame(cur_var,cur_dat$Group.2)
cur_fun = function(d,indices) {
  cur_d = d[indices,]
  fit9              = lm(cur_var ~ 1, data=cur_d)
  coef(fit9)
}

tmp = boot(cur_var_df,cur_fun,R=1000)
obs_fixef  <- as.numeric(mod1$coefficients[1]) + as.numeric(mod1$coefficients[3])
nulldens   <- density(tmp$t[,1],bw = c("ucv"),n = 1000)
nulldens$y <- nulldens$y/sum(nulldens$y)
nulldens$y <- cumsum(nulldens$y)
f <- approxfun(nulldens, rule=2)
1-f(obs_fixef)
