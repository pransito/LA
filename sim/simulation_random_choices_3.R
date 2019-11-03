## simulating random answering and its effect on lambda
# draws from data.la (the complete data frame in long fromat)

# data set generation at initial point (HC vs. RR)
data.la.HC         <- subset(data.la,group=="HC")
data.la.HC         <- na.omit(data.la.HC[c("subject","group","accept.reject","Gewinn","Verlust")])
data.la.RR         <- data.la.HC
data.la.RR$group   <- as.character(data.la.RR$group)
data.la.RR$group   <- "RR"
data.la.HC$group   <- as.character(data.la.HC$group)
data.la.HCRR       <- rbind(data.la.HC,data.la.RR)
data.la.HCRR$group <- as.factor(data.la.HCRR$group) 

# one by one replace choices with random choices
n_replaced = seq(0,100,length.out = 700)

# function for extracting fixefs of interest
get.la.fixef.RR <- function(.) {
  
  # intercepts
  an_RR <- as.numeric((fixef(.)["(Intercept)"] + fixef(.)["groupRR"]))
  an_HC <- as.numeric(fixef(.)["(Intercept)"])
  
  la_RR <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupRR"])*(-1)/(fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupRR"]))
  la_HC <- as.numeric((fixef(.)["Verlust"])*(-1)/(fixef(.)["Gewinn"]))
  
  bg_RR <- as.numeric((fixef(.)["Gewinn"] + fixef(.)["Gewinn:groupRR"]))
  bg_HC <- as.numeric(fixef(.)["Gewinn"])
  
  bl_RR <- as.numeric((fixef(.)["Verlust"] + fixef(.)["Verlust:groupRR"]))
  bl_HC <- as.numeric(fixef(.)["Verlust"])
  
  # contrasts
  x_la_HCgrRR <- la_HC - la_RR
  
  # return the coefficients
  tmp   <- ls()
  out_v <- c() 
  for (ii in 1:length(tmp)) {
    out_v[ii] <- eval(parse(text=tmp[ii]))
  }
  names(out_v) <- tmp
  return(out_v)
}

# the function to get the entr
f <- function(x) {
  cur_entr <- agk.entr.st(x,n = length(x))
  return(cur_entr)
}

# function to impute the random choices subject by subject given a distr sd
# df has accept.reject as variable with simulated (non-random) choices
# cur_ran is the n_replaced percentage of this randomness iteration
agk.impute.random <- function(df,cur_ran,cur_sd) {
  levs       <- levels(as.factor(as.character(df$subject)))
  if (cur_sd == 0) {
    n_repl_vec <- rep(cur_ran,times=length(levs))
  } else {
    n_repl_vec <- (randn(length(levs),1)*cur_sd)+cur_ran
  }
  n_repl_vec <- ifelse(n_repl_vec<0,0,identity(n_repl_vec))
  n_repl_vec <- ifelse(n_repl_vec>100,100,identity(n_repl_vec))
  for (ii in 1:length(levs)) {
    df$accept.reject[df$subject == levs[ii]] <- as.factor(v_ran(df$accept.reject[df$subject == levs[ii]], n_repl_vec[ii]))
  }
  return(df)
}

# decide which sd should be used to spread the levels of random choices in RR group
hmdf <- aggregate(as.numeric(as.character(data.la.HC$accept.reject)),by=list(data.la.HC$subject, data.la.HC$group,data.la.HC$Gewinn,data.la.HC$Verlust),FUN="mean.rmna")
names(hmdf) <- c("subject", "group","Gewinn_centered","Verlust_centered","PoA")
tmp <- aggregate(hmdf$PoA,by = list(hmdf$subject,hmdf$group),FUN=f)
names(tmp) <- c("subject", "group", "entropy")
rand_sd <- sd(tmp$entropy)*100

data.la.HCRRbl  <- data.la.HCRR
bl_mod          <-  glmer(accept.reject ~ (Gewinn+Verlust)*group + (Gewinn+Verlust|subject),data=data.la.HCRRbl,family="binomial")

## simulate function
agk.sim.rand <- function(cur_n_replaced, data.la.HCRR,data.la.HCRRbl,rand_sd,bl_mod) {
  
  # impute random choices
  #tmp <- simulate(bl_mod)
  #data.la.HCRR$accept.reject <- tmp$sim_1 
  #data.la.HCRR$accept.reject[data.la.HCRR$group == "HC"] <- data.la.HCRRbl$accept.reject[data.la.HCRRbl$group == "HC"]
  cur_HC <- subset(data.la.HCRR,group=="RR")
  cur_HC <- agk.impute.random(cur_HC,cur_ran = cur_n_replaced,cur_sd = rand_sd) 
  data.la.HCRR$accept.reject[data.la.HCRR$group == "RR"] <- cur_HC$accept.reject
  
  # get the current PoA matrix
  hmdf          <- aggregate(as.numeric(as.character(data.la.HCRR$accept.reject)),by=list(data.la.HCRR$subject, data.la.HCRR$group,data.la.HCRR$Gewinn,data.la.HCRR$Verlust),FUN="mean",na.rm=T)
  names(hmdf)   <- c("subject", "group","Gewinn_centered","Verlust_centered","PoA")
  tmp           <- aggregate(hmdf$PoA,by = list(hmdf$subject,hmdf$group),FUN=f)
  names(tmp)    <- c("subject", "group", "entropy")
  data.la.HCRR  <- merge(data.la.HCRR,tmp,by = c("subject","group"))
  entr_by_grp   <- describeBy(tmp$entropy,tmp$group)
  entr          <- data.la.HCRR$entropy
  data.la.HCRR$entropy <- scale(data.la.HCRR$entropy)
  data.la.HCRR$entropy_sqrd <- (data.la.HCRR$entropy)^2 
  
  # fit model with random answers inserted
  RR_mod <- glmer(accept.reject ~ (Gewinn+Verlust)*group + (Gewinn+Verlust|subject),data=data.la.HCRR,family="binomial")
  param_est <- get.la.fixef.RR(RR_mod) 
  # fit model with random answers inserted but controlled for entropy
  RR_modce <- glmer(accept.reject ~ (Gewinn+Verlust)*(group + entropy) + (Gewinn+Verlust|subject),data=data.la.HCRR,family="binomial")
  param_est_conte <- get.la.fixef.RR(RR_modce)
  # fit model with random answers inserted but controlled for entropy^2
  RR_modcesq <- glmer(accept.reject ~ (Gewinn+Verlust)*(group + entropy_sqrd) + (Gewinn+Verlust|subject),data=data.la.HCRR,family="binomial")
  param_est_contesq <- get.la.fixef.RR(RR_modcesq)
  
  # get the current pred_perc
  pred_perc_list <- round(predict(RR_mod,type="response")) == data.la.HCRR$accept.reject
  data.la.HCRR$pred <- pred_perc_list
  tmp <- aggregate(data.la.HCRR$pred,by=list(data.la.HCRR$subject, data.la.HCRR$group),FUN="mean",na.rm=T)
  names(tmp) <- c("subject","group","pred_perc")
  data.la.HCRR <- merge(data.la.HCRR,tmp,by = c("subject","group"))
  predperc_by_grp <- describeBy(tmp$pred_perc,tmp$group)
  data.la.HCRR$pred_perc <- scale(data.la.HCRR$pred_perc)
  
  # fit model with random answers inserted but controlled for pred_perc
  RR_modcp <- glmer(accept.reject ~ (Gewinn+Verlust)*(group + pred_perc) + (Gewinn+Verlust|subject),data=data.la.HCRR,family="binomial")
  param_est_contp <- get.la.fixef.RR(RR_modcp) 
  
  # get sig scores
  sig_scores <- list(RR_mod = coef(summary(RR_mod)),RR_modce = coef(summary(RR_modce)),RR_modcp = coef(summary(RR_modcp)))
  
  # saving...
  pred_perc_group_sim <- list(cur_perc_rand = cur_n_replaced,param_est=param_est,param_est_contp=param_est_contp,param_est_conte=param_est_conte,
                              param_est_contesq=param_est_contesq,
                              sig_scores=sig_scores,
                              entr=entr,entr_by_grp=entr_by_grp,predperc_by_grp=predperc_by_grp)
  return(pred_perc_group_sim)
}

## simulate
# evaluate all models in parallel
cl<-makeCluster(7) #change to your number of CPU cores  
registerDoSNOW(cl)
sim_rand <- foreach(ff=1:length(n_replaced), .packages=c('lme4','pracma','psych'),.verbose=T) %dopar% {
  agk.sim.rand(n_replaced[ff],data.la.HCRR,data.la.HCRRbl,rand_sd,bl_mod)
} 
stopCluster(cl)
setwd(pfad_ransim)
save(file="sim_rand.RData",list = c("sim_rand"))


## analyze the simulations
# pred_perc with growing randomness
f <- function (x) {tmp <- x[[9]]; return(tmp[[2]][[3]])} 
plot(unlist(lapply(sim_rand,f)),main="Randomness and Pred Perc",xlab="level of rand",ylab="group mean pred perc")

# entropy with growing randomness
f <- function (x) {tmp <- x[[8]]; return(tmp[[2]][[3]])} 
plot(unlist(lapply(sim_rand,f)),main="Entropy and Pred Perc", xlab="level of rand",ylab="group mean entropy")

# lambda con with growing randomness
f <- function (x) {tmp <- x[[2]]; return(tmp[9])} 
plot(unlist(lapply(sim_rand,f)),main="Lambda HC>RR with growing randomness in RR",ylim =c(-10,10),xlab="level of rand",ylab="group mean diff HC>RR")

# lambda con with growing randomness
f <- function (x) {tmp <- x[[3]]; return(tmp[9])} 
plot(unlist(lapply(sim_rand,f)),main="Lambda HC>RR with growing randomness in RR\n con for pred_perc",ylim =c(-10,10),xlab="level of rand",ylab="group mean diff HC>RR")

# lambda con with growing randomness
f <- function (x) {tmp <- x[[4]]; return(tmp[9])} 
plot(unlist(lapply(sim_rand,f)),main="Lambda HC>RR with growing randomness in RR\n con for entropy",ylim =c(-10,10),xlab="level of rand",ylab="group mean diff HC>RR")

# lambda with growing randomness
f <- function (x) {tmp <- x[[5]]; return(tmp[9])} 
plot(unlist(lapply(sim_rand,f)),main="Lambda HC>RR with growing randomness in RR\n con for entropy^2",ylim =c(-10,10),xlab="level of rand",ylab="group mean diff HC>RR")

# plot the local average smooth onto Lambda contrast
f <- function (x) {tmp <- x[[2]]; return(tmp[9])} 
x <- unlist(lapply(sim_rand,f))
dat <- data.frame(x=x)
dat$nr <- 1:nrow(as.matrix(x))

ggplot() + geom_point(aes(x=nr,y=x), data=dat)+ geom_hline(yintercept=0) + geom_smooth(aes(x=nr,y=x), data=dat, method="loess",level=0.9)  + coord_cartesian(ylim=c(-10,10))
