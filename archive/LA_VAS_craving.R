## LA: VAS craving effects?
# first load data with main script

vars = c("A", "G")
m <- list()

for (ii in 1:length(vars)) {
  cur_filter <- paste("^VAS_[", vars[ii],"]_P",sep="")
  
  names_VAS <- grep(cur_filter, names(data.la), value=TRUE)
  
  first <- function(x) {return(x[1])}
  data.la.VAS <- data.la[c("subject","group",names_VAS)]
  data.la.VAS <- aggregate(data.la.VAS,by = list(data.la.VAS$subject),FUN = "first")
  data.la.VAS <- data.la.VAS[,-1]
  names(data.la.VAS) <- c("subject","group","x.1","x.2")
  data.la.VAS <- na.omit(data.la.VAS)
  data.la.VAS$subject <- as.character(data.la.VAS$subject)
  data.la.VAS$subject <- as.factor(data.la.VAS$subject)
  data.la.VAS$diff <- scale(data.la.VAS$x.2-data.la.VAS$x.1) 
  contrasts(data.la.VAS$group) <- contr.treatment(n = 3,base = 2)
  colnames(contrasts(data.la.VAS$group)) <- c("AD","PG")
                  
  m[[ii]] <- lmer(diff ~ group + (1 | subject),data=data.la.VAS,REML=F,control=lmerControl(check.nobs.vs.nlev="ignore",
                                                                                      check.nobs.vs.rankZ="ignore",
                                                                                      check.nobs.vs.nRE="ignore"))
  
}

boot_m1 <- bootMer(m[[1]],FUN = fixef,nsim = 100,type = "parametric")
hist(boot_m1$t[,3],main ="PG > HC")
hist(boot_m1$t[,2], main ="AD > HC")
hist(boot_m1$t[,1], main ="HC")
boot_m2 <- bootMer(m[[2]],FUN = fixef,nsim = 1000,type = "parametric")
hist(boot_m2$t[,3],main ="PG > HC")
hist(boot_m2$t[,2], main ="AD > HC")
hist(boot_m2$t[,1], main ="HC")

cur_con = function(x,i) {
  x <- x[i,]
  return(mean(x$diff[x$group =="HC"]) - mean(x$diff[x$group =="PG"]))
}

HCvsPG_diff <- boot(data.la.VAS[data.la.VAS$group != "AD",],statistic = cur_con,R = 2000)

# wilcox
cur_dat <- subset(data.la.VAS,group=="HC") 
wilcox.test(x=cur_dat$x.1,y=cur_dat$x.2,mu = 0,paired = T,exact = T)
cur_dat <- subset(data.la.VAS,group=="PG") 
wilcox.test(x=cur_dat$x.1,y=cur_dat$x.2,mu = 0,paired = T,exact = T)
cur_dat <- subset(data.la.VAS,group=="AD") 
wilcox.test(x=cur_dat$x.1,y=cur_dat$x.2,mu = 0,paired = T,exact = T)




