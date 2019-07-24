# test craving gambling/alcohol pre-post;

# functions
mean.rmna.rtna = function(x) {
  if(sum(is.na(x)) == length(x)) {return(NA)} else {
    return(mean(x,na.rm=T))
  }
}

# for alcohol
cur_names = names(data.la)[grep(pattern = "^VAS_A_P", names(data.la),perl = F)]
cur_dat   = data.la[c("subject", "group", cur_names)]


agg_dat       = aggregate(. ~ subject,data = cur_dat,mean.rmna.rtna)

ln_dat        = as.data.frame(rbind(as.matrix(agg_dat$subject),as.matrix(agg_dat$subject)))
names(ln_dat) = c("subject")
ln_dat$group  = as.factor(rbind(as.matrix(agg_dat$group),as.matrix(agg_dat$group)))
contrasts(ln_dat$group) = contr.treatment(3,base=2)
ln_dat$crav_A = rbind(as.matrix(agg_dat$VAS_A_PRE_mean),as.matrix(agg_dat$VAS_A_POST_mean))
ln_dat$time   = as.factor(rbind(as.matrix(rep(1,47)),as.matrix(rep(2,47))))

require(nlme)
crav_A_mod_nlme_0 = lme(crav_A ~ 1, random = ~1|subject/time, data=ln_dat,method = "ML")
crav_A_mod_nlme_1 = lme(crav_A ~ group, random = ~1|subject/time, data=ln_dat,method = "ML")
crav_A_mod_nlme_2 = lme(crav_A ~ group*time, random = ~1|subject/time, data=ln_dat,method = "ML")
anova(crav_A_mod_nlme_0,crav_A_mod_nlme_1,crav_A_mod_nlme_2)
summary(crav_A_mod_nlme_0)

# for gambling
cur_names = names(data.la)[grep(pattern = "^VAS_G_P", names(data.la),perl = F)]
cur_dat   = data.la[c("subject", "group", cur_names)]


agg_dat       = aggregate(. ~ subject,data = cur_dat,mean.rmna.rtna)

ln_dat        = as.data.frame(rbind(as.matrix(agg_dat$subject),as.matrix(agg_dat$subject)))
names(ln_dat) = c("subject")
ln_dat$group  = as.factor(rbind(as.matrix(agg_dat$group),as.matrix(agg_dat$group)))
contrasts(ln_dat$group) = contr.treatment(3,base=2)
ln_dat$crav_G = rbind(as.matrix(agg_dat$VAS_G_PRE_mean),as.matrix(agg_dat$VAS_G_POST_mean))
ln_dat$time   = as.factor(rbind(as.matrix(rep(1,47)),as.matrix(rep(2,47))))

require(nlme)
crav_G_mod_nlme_0 = lme(crav_G ~ 1, random = ~1|subject/time, data=ln_dat,method = "ML")
crav_G_mod_nlme_1 = lme(crav_G ~ group, random = ~1|subject/time, data=ln_dat,method = "ML")
crav_G_mod_nlme_2 = lme(crav_G ~ group*time, random = ~1|subject/time, data=ln_dat,method = "ML")
anova(crav_G_mod_nlme_0,crav_G_mod_nlme_1,crav_G_mod_nlme_2)
summary(crav_G_mod_nlme_1)

