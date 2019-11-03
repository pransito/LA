###########################################################
#         Checking gain and loss by group                 #
###########################################################

# 09.04.2015
# using the gen mixed model with ed as additional predictor;

## preparation

rm(list = ls())

# path
#pfad <- "C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"
#pfad <- "C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"
pfad <- "C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"

setwd(pfad)

# load libraries and functions
setwd(paste(pfad, "\\Scripts", sep=""))
source ('LA_functions.R')

# set the pwd needed and load data
setwd(paste(pfad, "\\Data", sep=""))
load("cdf_ed.Rdata")
cdf <- cdf_ed

# plot bar
prep_ggplot_la()
# prep data
cur.cdf   <- data.frame(cdf$Gewinn_ed,cdf$Verlust_ed, cdf$Lambda_ed)
plot.data <- data.frame()
var.names <- c()
for (ii in 1:length(cur.cdf[1,])){
  # bootstrapped 95% CI
  bt.df     <- as.data.frame(as.list(aggregate(cur.cdf[ii],by = list(group = cdf$group), smean.cl.boot)))
  names(bt.df) <- c("group", "mean","lower", "upper")
  var.names <- c(var.names,rep(paste(ii,names(cur.cdf[ii]),sep="_"),length(levels(bt.df$group))))
  plot.data <- rbind(plot.data,bt.df)
}
plot.data$variable <- var.names

# plot
p <- ggplot(data = plot.data, aes(group,mean,fill=variable))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=.2, width=0.5,
              position=position_dodge(width=0.9)) + ylab("mean (95% boots. CI)")

p

# plot bar (geometric mean)
# get log, then mean, then retransform and subtract constant
prep_ggplot_la()
# prep data
cur.cdf   <- data.frame(cdf$Gewinn_ed_Aga,cdf$Verlust_ed_Aga, cdf$Lambda_ed_Aga)
# get log
dat.log <- list()
con.log <- list()
for (ii in 1:length(cur.cdf[1,])){
  cur.log <- get_log(cur.cdf[[ii]])
  dat.log[[ii]] <- cur.log$v1 
  con.log[[ii]] <- cur.log$const 
}

plot.data <- data.frame()
var.names <- c()
for (ii in 1:length(cur.cdf[1,])){
  # bootstrapped 95% CI
  bt.df     <- as.data.frame(as.list(aggregate(dat.log[[ii]],by = list(group = cdf$group), smean.cl.boot)))
  # retransform
  bt.df <- cbind(bt.df$group,exp(bt.df[,c(2:4)])+con.log[[ii]])
  names(bt.df) <- c("group", "mean","lower", "upper")
  var.names <- c(var.names,rep(paste(ii,names(cur.cdf[ii]),sep="_"),length(levels(bt.df$group))))
  plot.data <- rbind(plot.data,bt.df)
}
plot.data$variable <- var.names

# plot
p <- ggplot(data = plot.data, aes(group,mean,fill=variable))
p <- p+geom_bar(position="dodge",stat="identity")
p
p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=.2, width=0.5,
                       position=position_dodge(width=0.9)) + ylab("geometric mean with 95% CI")

p

## tests with controlling for covariates
## complete covariate model
dar.log <- list()
for (ii in 1:length(dat.log)) {
  la               <- lm(dat.log[[ii]] ~ cdf$group, data=cdf)
  la.cleaned       <- lm(dat.log[[ii]] ~ scale(cdf$Age,s=F) + scale(cdf$Bildungsjahre_ges,s=F)+scale(cdf$pred_perc_ed,s=F)+group, data=cdf)
  la.cleaned.nogrp <- lm(dat.log[[ii]] ~ scale(cdf$Age,s=F) + scale(cdf$Bildungsjahre_ges,s=F)+scale(cdf$pred_perc,s=F), data=cdf)
  dar.log[[ii]] <- resid(la.cleaned.nogrp)
}

