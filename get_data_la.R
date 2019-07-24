## prep some variables
missing.check  <- 1 # if set to 1 then people with too many missings will be dropped 
missing.people <- list()

## prep the calculation of euclidean distance

Q1 <- as.matrix(c(14, -7, 0))
Q2 <- as.matrix(c(36, -18, 0))

## get all the necessary .txt files from the current directory
setwd(path_data)
all.behav.files <- list.files(all.files=T, pattern = "*_raw.*")

## start for loop to get all the data of the subjects and put them in one wide data frame
# i is the variable containing the current txt file

for (i in 1:length(all.behav.files)) {  
  
  
  ## read in data for one subject
  df.cur.person <- read.table(all.behav.files[i], header=T, skip=1, nrows=144, sep="\t")
  
  ## euclidean distance ## IS NOT USED ANYMORE BECAUSE WRONG; GETS OVERWRITTEN LATER IN SCRIPT ##
  ed <- c()
  library(pracma,warn.conflicts = F)
  
  for (g in 1:length(df.cur.person[,1])) {
    
    P       <- as.matrix(c(df.cur.person$Gewinn[g], df.cur.person$Verlust[g], 0))
    nenner  <- (pracma::cross(Q2-Q1, P-Q1))
    zaehler <- (Q2-Q1) 
    cur.ed  <- nenner[3,]/zaehler[1,]
    ed      <- rbind(ed, cur.ed)
  }
  # negative ed: further into the high loss low gain half
  # positive ed: further into the low loss high gain half
  
  ed.abs <- abs(ed)
  df.cur.person <- data.frame(df.cur.person, ed, ed.abs)
  detach("package:pracma", unload=TRUE)
  
  if (loss.abs == 1) {
    df.cur.person$Verlust <- abs(df.cur.person$Verlust) ## work as Tom et al. (2007)
  }
  
  ## die Formeln von Minati und Martino und mir sind hier gleich
  if (loss.abs == 1) {
    df.cur.person$EV      <- 0.5*df.cur.person$Gewinn + 0.5*df.cur.person$Verlust*(-1)
    df.cur.person$ratio   <- abs(df.cur.person$Gewinn/df.cur.person$Verlust*(-1))
    df.cur.person$diff    <- df.cur.person$Gewinn-df.cur.person$Verlust*(-1)
  } else {
    df.cur.person$EV      <- 0.5*df.cur.person$Gewinn + 0.5*df.cur.person$Verlust
    df.cur.person$ratio   <- abs(df.cur.person$Gewinn/df.cur.person$Verlust)
    df.cur.person$diff    <- df.cur.person$Gewinn-df.cur.person$Verlust
  }
  
  if (loss.abs == 1) {
    df.cur.person$RiskMin <- 0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2+0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2 # typo?? no becuase it is symmetrically far away
    df.cur.person$RiskAG  <- 0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2+0.5*(df.cur.person$Verlust*(-1)-df.cur.person$EV)^2
    df.cur.person$RiskMar <- (0.5*df.cur.person$Gewinn-0.5*df.cur.person$Verlust*(-1))^2
  } else {
    df.cur.person$RiskMin <- 0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2+0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2 # typo??; no because it is symmetric 
    df.cur.person$RiskAG  <- 0.5*(df.cur.person$Gewinn-df.cur.person$EV)^2+0.5*(df.cur.person$Verlust-df.cur.person$EV)^2
    df.cur.person$RiskMar <- (0.5*df.cur.person$Gewinn-0.5*df.cur.person$Verlust)^2
  }

  ## Skew aufbauend auf Formel von Minati
  df.cur.person$SkewMin  <- 0.5*(df.cur.person$Gewinn-df.cur.person$EV)^3+0.5*(df.cur.person$Gewinn-df.cur.person$EV)^3
  
  if (use.z == 1) {
    
    df.cur.person$Gewinn  <- scale(df.cur.person$Gewinn, center = T, scale = F)
    df.cur.person$Verlust <- scale(df.cur.person$Verlust, center = T, scale = F)
    df.cur.person$EV      <- scale(df.cur.person$EV, center = T, scale = F)
    df.cur.person$ed.abs  <- scale(df.cur.person$ed.abs, center = T, scale = F)
    df.cur.person$RiskMar <- scale(df.cur.person$RiskMar, center = T, scale = F)
    df.cur.person$ratio   <- scale(df.cur.person$ratio,center =T,scale=F)
    df.cur.person$diff    <- scale(df.cur.person$diff,center =T,scale=F)
  }
  
  ## collect per subject the p and r value;
  p <- rbind(p, as.numeric(cor.test(df.cur.person$Verlust, df.cur.person$RT, method="spearman")[3]))
  r <- rbind(r, as.numeric(cor.test(df.cur.person$Verlust, df.cur.person$RT, method="spearman")[4]))
  
  accept.reject <- df.cur.person$Button
  
  ## checking for missings; if there are too many
  
  if (missing.check == 1) {
    m.test <- which(accept.reject==5)
    if (length (m.test) > missing.cutoff*length(df.cur.person[,1])) {
      missing.people[length(missing.people)+1] <- all.behav.files[i]
      next}}
  
  ## create subject variable
  
  current.subject <- all.behav.files[i]
  k               <- strsplit(current.subject, "_")
  k               <- unlist(k)
  current.subject <- k[2]
  current.subject <- as.numeric(current.subject)
  subject <- rep(current.subject, length =144)
  
  ## it so happens that we have missing data in the categorical response variable,
  ## we have to correct for this by setting "5" in accept.reject to "NA"
  ## and we recode the ordinal response scale into a dichotomous response scale
  
  ## recode accept.reject "5" into NA
  ## recode 1,2 into 1 and 3,4 into 0
  
  for (ii in 1:length(accept.reject)) {
    
    if(accept.reject[ii]==1) {accept.reject[ii]<-1}
    if(accept.reject[ii]==2) {accept.reject[ii]<-1}
    if(accept.reject[ii]==3) {accept.reject[ii]<-0}
    if(accept.reject[ii]==4) {accept.reject[ii]<-0}
    if(accept.reject[ii]==5) {accept.reject[ii]<-NA}
  }
  
  df.cur.person <- data.frame(subject,df.cur.person, accept.reject)
  
  if (!exists ("data.la")) {data.la <- df.cur.person} else if (exists ("data.la")) {
    data.la <- rbind(data.la, df.cur.person)}
  
}

## group vector

subjects.numeric <- as.numeric(data.la$subject)

group <- c()

for (i in 1:length(subjects.numeric)) {
  
  if (subjects.numeric[i] > 999 & subjects.numeric[i] < 2000) {
    current.group <- "HC"}
  
  if (subjects.numeric[i] > 1999 & subjects.numeric[i] < 3000) {
    current.group <- "PG"}
  
  if (subjects.numeric[i] > 2999 & subjects.numeric[i] < 4000) {
    current.group <- "AD"}
  
  group <- cbind(group, current.group)
  
}

group   <- t(group)
data.la <- (cbind(data.la, group))
rm(group)

data.la$group <- as.factor(data.la$group)

## data handling

data.la$Gewinn    <- as.numeric(data.la$Gewinn)
data.la$Verlust   <- as.numeric(data.la$Verlust)
data.la$subject   <- as.factor(data.la$subject)
data.la$diff      <- as.numeric(data.la$diff)

data.la$accept.reject <- as.factor(data.la$accept.reject)

## here we select the HC subsample and then per subject only a sample of 30 trials
data.la_HC <- subset (data.la, data.la$group=="HC")

## here we select only the odd or even trials
data.la_odd  <- subset (data.la, mod(as.numeric(row.names(data.la)),2)==0)
data.la_even <- subset (data.la, mod(as.numeric(row.names(data.la)),2)!=0)

data.la.HC        <- data.la[which(data.la$group   == "HC"),]
data.la.PG        <- data.la[which(data.la$group   == "PG"),]
data.la.AD        <- data.la[which(data.la$group   == "AD"),]