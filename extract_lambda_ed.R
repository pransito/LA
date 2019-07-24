## get the complete coefficient
## make group vector
tmp <- aggregate(data.la$RT, by=list(data.la$subject, data.la$group), FUN=mean)
HC<-rep("HC",as.numeric(table(tmp$Group.2)[2]))
AD<-rep("AD",as.numeric(table(tmp$Group.2)[1]))
PG<-rep("PG",as.numeric(table(tmp$Group.2)[3]))
group <- rbind(as.matrix(HC),as.matrix(PG),as.matrix(AD))
## get intercept, gain, loss of HC fixef
intercept <- matrix(rep(fixef(la05)["(Intercept)"], length.out=nrow(tmp)))
ed        <- matrix(rep(fixef(la05)["ed"], length.out=nrow(tmp)))
## get group specific int
int.AD    <- matrix(rep(fixef(la05)["groupAD>HC"], as.numeric(table(group)["AD"])))
int.PG    <- matrix(rep(fixef(la05)["groupPG>HC"], as.numeric(table(group)["PG"])))
int.HC    <- matrix(rep(0,as.numeric(table(group)["HC"])))
int.group <- rbind(int.HC, int.PG, int.AD)
intercept <- intercept + int.group
##get group specific gain
ed.HC   <- matrix(rep(0,as.numeric(table(group)["HC"])))
ed.PG   <- matrix(rep(fixef(la05)["ed:groupPG>HC"], as.numeric(table(group)["PG"])))
ed.AD   <- matrix(rep(fixef(la05)["ed:groupAD>HC"], as.numeric(table(group)["AD"])))
ed.group<- rbind(ed.HC,ed.PG,ed.AD)
ed      <- ed + ed.group
## fixef as a whole
fixef.la  <- cbind(intercept, ed)
## ad ranef
test      <- fixef.la + ranef(la05)$subject

crm       <- data.frame(id=rownames(ranef(la05)$subject),
                        group =group,
                        intercept=test["(Intercept)"],
                        gain=test["ed"])

## make predictions from the crm coefficients

## functions
prob.log.reg <-function(int,b1,x1){
  
  return (1/(1+exp(-(int+x1*b1))))
}
## predict function
mylogit_predict <- function (cur.data,cur.coeff) {
  cur.pred <- c()
  for (ii in 1:length(cur.data[,1])) {
    tmp <- prob.log.reg(cur.coeff[3], cur.coeff[4], cur.data$ed[ii])
    cur.pred <- rbind(cur.pred, tmp)    
  }
  return (cur.pred)
}


subjects <- as.numeric(row.names(as.matrix(table(data.la$subject))))
preds    <- c()
for (jj in 1:length(subjects)) {
  cur.data <- subset(data.la, subject == subjects[jj])
  cur.coef <- subset(crm, id == subjects[jj])
  preds    <- rbind(preds, mylogit_predict(cur.data,cur.coef))
  
}

preds <- round(as.numeric(preds[[1]]))
data.la <- data.frame(data.la, preds)
data.la$cor.pred <- abs(as.numeric(as.matrix(data.la$accept.reject))-data.la$preds)
data.la$cor.pred <- ifelse(data.la$cor.pred==0,1,0)

test <- aggregate(x=data.la$cor.pred, by=list(data.la$subject), FUN = "mean",na.rm=T)
pred_perc <- as.numeric(test[[2]])
crm <- data.frame(crm, pred_perc)
