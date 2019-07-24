## check whether subjects did not really pay attention

# via entropy
# the function to get the entr of a probability distribution
f <- function(x) {
  cur_entr <- agk.entr.st(x)
  return(cur_entr)
}

hmdf          <- aggregate(as.numeric(as.character(data.la$accept.reject)),by=list(data.la$subject, data.la$group,data.la$Gewinn,data.la$Verlust),FUN="mean",na.rm=T)
names(hmdf)   <- c("subject", "group","Gewinn_centered","Verlust_centered","PoA")
tmp           <- aggregate(hmdf$PoA,by = list(hmdf$subject,hmdf$group),FUN=f)
names(tmp)    <- c("subject", "group", "entropy")

tmp = merge(crm,tmp,by=c("subject","group"))
kruskal.test(tmp$entropy,tmp$group)


# via sig of beta loss / beta gain
cur_mod   <- lmList(accept.reject ~ ed.abs + Gewinn + Verlust | subject,data = data.la,family = "binomial",pool = F,na.action =NULL)
cur_est <- summary(cur_mod)
cur_est <- cur_est$coefficients
cur_sel <- data.frame(cur_est[,4,3] < 0.1 | cur_est[,4,4] < 0.1) # p-values Gewinn and p-values Verlust
cur_sel$sub <- row.names(cur_sel)
names(cur_sel) <- c("sel","sub")
