## simulating random answering and its effect on lambda
# draws from data.la (the complete data frame in long fromat)
# also needs the cdf_HC aggregated data

load("cdf_ed.Rdata")
cdf <- cdf_ed
cdf_HC <- cdf[cdf$group == "HC",]

# first find a subject with especially good model fit and middle acceptance rate (over 40%)
# and lambda under 5
data.la.HC        <- data.la[which(data.la$group   == "HC"),]
cdf_HC            <- cdf_HC[cdf_HC$perc.accept > 0.4,]
cdf_HC            <- cdf_HC[cdf_HC$Lambda_ed < 5,]
cur_sub = cdf_HC$subject[cdf_HC$pred_perc == max(cdf_HC$pred_perc)]
cur_sub = as.numeric(as.matrix(cur_sub))
data.la.HC_best <- data.la.HC[which(as.numeric(as.matrix(data.la.HC$subject))   == cur_sub),]
# remove NA trials
data.la.HC_best <- data.la.HC_best[!is.na(data.la.HC_best$accept.reject),]


# ony by one replace choices with random choices
n_replaced = seq(0,100,length.out = 2000)
lambdas_ed = c()
betas_g_ed = c()
betas_l_ed = c()
pred_perc  = c()
save_warn  = list()

for (jj in 1:length(n_replaced)) {
  cur.data   = data.la.HC_best
  cur.data$accept.reject <- as.factor(v_ran(cur.data$accept.reject, n_replaced[jj]))
  mylogit    = glm(accept.reject ~ ed.abs + Gewinn + Verlust, data=cur.data, family="binomial")
  
  lambdas_ed[jj] = calc_lambda(as.numeric(coef(mylogit)[3]), as.numeric(coef(mylogit)[4]))
  betas_g_ed[jj] = as.numeric(coef(mylogit)[3])
  betas_l_ed[jj] = as.numeric(coef(mylogit)[4])
  tmp.pred       = round(predict(mylogit, type="response"))
  tmp.targ       = as.numeric(as.matrix(na.omit(data.la.HC_best$accept.reject)))
  pred_perc      = rbind(pred_perc,1-mean(abs(tmp.targ-tmp.pred),na.rm=T))
}

# plots
plot(lambdas_ed,ylim=c(-10,10))
plot(pred_perc,ylim=c(0,1))
plot(n_replaced,pred_perc)
plot(n_replaced,lambdas_ed)
plot(n_replaced,betas_g_ed)
plot(n_replaced,betas_l_ed)
plot(pred_perc,n_replaced)

## sampling now from single cases to make an "AD" group; reference is the HC group
# actually now we have to estimate the number of random responses from the empirical AD pred_perc
# use the pred_perc as an indicator to pull an instance of "random answering"
# draw a sample from the "AD" group (which here is the group made from one HC subject with different levels of random responses)
# i.e. the random responses group "RR"

# the AD's pred_perc empirically
load ("cdf_ed.Rdata")
cdf_HC <- cdf[cdf$group == "HC",]
cdf_AD <- cdf[cdf$group == "AD",]

tmp.prob = round(cdf_AD$pred_perc_ed,digits = 1)
# prep HC data
group_1 <- as.matrix(rep("HC",17))
group_2 <- as.matrix(rep("RR", 12))
group <- as.factor(rbind(group_1,group_2))
load ("cdf_ed.Rdata")
cdf_HC <- cdf[cdf$group == "HC",]
ks.test.result.la  <- c()
ks.test.result.bg  <- c()
ks.test.result.bl  <- c()
ks.test.result.con.la <- c()
ks.test.result.con.bg <- c()
ks.test.result.con.bl <- c()
cur.jitter         <- 0.02

for (jj in 1:100) { 
  # get the current sample
  cur.sample.la <- c()
  cur.sample.bg <- c()
  cur.sample.bl <- c()
  cur.pred.perc <- c()
  for (kk in 1:length(tmp.prob)){
    if (tmp.prob[kk] == 1) {tmp.prob[kk] <- 0.9}
    # estimate the n_replaced using the empirical pred_perc and the simulated data
    cur.n.replaced    <- n_replaced[pred_perc>(tmp.prob[kk]-cur.jitter) & pred_perc<(tmp.prob[kk]+cur.jitter)]
    cur.n.replaced    <- sample(cur.n.replaced,1)
    cur.sample.la[kk] <- lambdas_ed[n_replaced==cur.n.replaced]
    cur.sample.bg[kk] <- betas_g_ed[n_replaced==cur.n.replaced]
    cur.sample.bl[kk] <- betas_l_ed[n_replaced==cur.n.replaced]
    cur.pred.perc[kk] <- pred_perc[n_replaced == cur.n.replaced] 
  }

  # without controlling for pred_perc
  HC_RR_betas_l_ed <- rbind(as.matrix(cdf_HC$Verlust_ed),as.matrix(cur.sample.bl))
  HC_RR_betas_g_ed <- rbind(as.matrix(cdf_HC$Gewinn_ed),as.matrix(cur.sample.bg))
  HC_RR_lambdas_ed <- rbind(as.matrix(cdf_HC$Lambda_ed),as.matrix(cur.sample.la))
  plot(HC_RR_betas_g_ed ~ group)
  plot(HC_RR_betas_l_ed ~ group)
  plot(HC_RR_lambdas_ed ~ group)
  
  tmp <- kruskal.test(HC_RR_betas_g_ed,group)
  ks.test.result.bg[jj] <- tmp$p.value
  tmp <- kruskal.test(HC_RR_betas_l_ed,group)
  ks.test.result.bl[jj] <- tmp$p.value
  tmp <- kruskal.test(HC_RR_lambdas_ed,group)
  ks.test.result.la[jj] <- tmp$p.value
  
  # with controlling for pred_perc
  # tie together HC pred perc (empirical) and RR pred perc from simulation
  cur.pred.perc <- rbind(as.matrix(cdf_HC$pred_perc_ed), as.matrix(cur.pred.perc))
  # partial out
  m <- lm(HC_RR_lambdas_ed ~ cur.pred.perc)
  HC_RR_lambdas_ed_con <- resid(m)
  plot(HC_RR_lambdas_ed_con ~ group)
  tmp <- kruskal.test(HC_RR_lambdas_ed_con,group)
  ks.test.result.con.la[jj] <- tmp$p.value
  
  m <- lm(HC_RR_betas_l_ed ~ cur.pred.perc)
  HC_RR_betas_l_ed_con <- resid(m)
  plot(HC_RR_betas_l_ed_con ~ group)
  tmp <- kruskal.test(HC_RR_betas_l_ed_con,group)
  ks.test.result.con.bl[jj] <- tmp$p.value
  
  m <- lm(HC_RR_betas_g_ed ~ cur.pred.perc)
  HC_RR_betas_g_ed_con <- resid(m)
  plot(HC_RR_betas_g_ed_con ~ group)
  tmp <- kruskal.test(HC_RR_betas_g_ed_con,group)
  ks.test.result.con.bg[jj] <- tmp$p.value
  
}
# without controlling for pred perc
hist(ks.test.result.la)
table((ks.test.result.la)<0.05)

hist(ks.test.result.bg)
table((ks.test.result.bg)<0.05)

hist(ks.test.result.bl)
table((ks.test.result.bl)<0.05)

# with controlling for pred perc
hist(ks.test.result.con.la)
table((ks.test.result.con.la)<0.05)

hist(ks.test.result.con.bg)
table((ks.test.result.con.bg)<0.05)

hist(ks.test.result.con.bl)
table((ks.test.result.con.bl)<0.05)

# did it always help to control for pred_perc?; check whether p-value always got smaller
check_ks_la <- ks.test.result.la - ks.test.result.con.la
check_ks_bg <- ks.test.result.bg - ks.test.result.con.bg
check_ks_bl <- ks.test.result.bl - ks.test.result.con.bl
table(sign(check_ks_la))
table(sign(check_ks_bg))
table(sign(check_ks_bl))
