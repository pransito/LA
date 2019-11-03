## check how randomness in answering influences single-person's lambda distribution
## for this we use the HC group's data and add different levels of randomness to there answers
## do this by replacing x amount of there answers with random 0's and 1's
## i.e. we pick a random HC subject, then we picka random amount of randomness and distort the data,
## then we calculate lambda and the goodness of prediction
## in the end we can plot on x-axis level of randomness and on y-axis level of prediction

###########################################################
#          VD-ANALYSEN LOSS AVERSION                      #
###########################################################

## A. Genauck
## 08.2.2014

## clear workspace
rm(list = ls())

## preparation

# zu hause
pfad <- "C:\\Users\\genaucka\\Google Drive\\Diplom\\LA\\daten_behav_test_finale_SP_Diplom"
setwd(pfad)

## options
use.z    <- 1 # currently centered-only
loss.abs <- 1 # use absolute value of loss as predictor

## prep the final matrix; will be a multiple column matrix
subject <- c()
p       <- c()
r       <- c()

## prep some variables

missing.check  <- 1 # if set to 1 then people with too many missings will be dropped 
missing.cutoff <- 0.07 # if there are more than x% percent missings in response then drop-out
missing.people <- list()
missing.count  <- c()

# load functions
setwd(paste(pfad, "\\Scripts", sep=""))
source("LA_functions.R")
# set the pwd needed
setwd(paste(pfad, "\\Data", sep=""))

# get the data.la
# get all the data in long format
setwd(paste(pfad, "\\Scripts", sep=""))
source("get_data_la.R")
# set the pwd needed
setwd(paste(pfad, "\\Data", sep=""))

### get all the necessary .txt files from the current directory

all.behav.files <- list.files(all.files=T, pattern = "*_raw.*")

### how many subjects shall be simulated?
rep_sub <- 100

## prep the output
step_col <-c()
Lambda <- c()
fits <- c()

### start for loop to run the logistic regression in one subject

for (jj in 1:rep_sub) {
  
  perc    <- round(runif(1,0,100)) # the current amount of trials which will be replaced
  cur.sub <- sample(1:17,1) # a random subject will be picked from the HC group
  
  ### read in data for one subject
  df.cur.person <- read.table(all.behav.files[cur.sub], header=T, skip=1, nrows=144, sep="\t")
  
  if (loss.abs == 1) {
    df.cur.person$Verlust <- abs(df.cur.person$Verlust) ## work as Tom et al. (2007)
  }
  
  if (use.z == 1) {
    
    df.cur.person$Gewinn  <- scale(df.cur.person$Gewinn, center = T, scale = F)
    df.cur.person$Verlust <- scale(df.cur.person$Verlust, center = T, scale = F)
  }
  
  accept.reject <- df.cur.person$Button
  
  ### checking for missings; if there are too many
  if (missing.check == 1) {
    m.test <- which(accept.reject==5)
    if (length (m.test) > missing.cutoff*length(df.cur.person[,1])) {
      missing.people[length(missing.people)+1] <- all.behav.files[i]
      next}
    missing.count <- rbind(missing.count, length(m.test))
  }
  
  ## recode accept.reject "5" into NA
  ## recode 1,2 into 1 and 3,4 into 0
  
  for (ii in 1:length(accept.reject)) {
    
    if(accept.reject[ii]==1) {accept.reject[ii]<-1}
    if(accept.reject[ii]==2) {accept.reject[ii]<-1}
    if(accept.reject[ii]==3) {accept.reject[ii]<-0}
    if(accept.reject[ii]==4) {accept.reject[ii]<-0}
    if(accept.reject[ii]==5) {accept.reject[ii]<-NA}
  }
  
  ## here randomize
  accept.reject <- v_ran(accept.reject, perc)
  
  df.cur.person <- data.frame(df.cur.person, accept.reject)
  
  ######################
  ### log regression ###
  ######################
  
  mylogit <- glm(accept.reject ~ Gewinn + Verlust, data = df.cur.person, family = binomial)
  
  ## get current betas/lambdas/goodness of fits
  current.beta.Gewinn     <- coef (mylogit)[[2]]
  current.beta.Verlust    <- coef (mylogit)[[3]]
  current.lambda          <- calc_lambda(current.beta.Gewinn,current.beta.Verlust)
  Lambda                  <- rbind(Lambda, as.numeric(current.lambda))
  
  predictions <- round(predict(mylogit, type=c("response")))
  disparity   <- predictions-(accept.reject[!is.na(accept.reject)]) ## no NA's!
  cur.fit     <- 1 - abs(disparity)
  cur.fit     <- mean(cur.fit)
  
  fits <- rbind(fits, as.numeric(cur.fit))
  
  step <- round (perc/100*144)
  step_col <- rbind(step_col, as.numeric(step))
  
}

## now we look into "group" Lambdas
intv <- (seq(1,144, 2))
spb  <- 5 # groups per bin
gsize<- 17 # group size
med.Lambda <- c()

for (jj in 1:length(intv)) {

  if (jj == length(intv)) {
    cur.intv.l <- intv[jj]
    cur.group  <- sample(Lambda[step_col>=cur.intv.l], 17)
    med.Lambda <- rbind(med.Lambda, as.numeric(median(cur.group)))
    next
  }
  cur.intv.l <- intv[jj]
  cur.intv.u <- intv[jj+1]
  cur.group  <- sample(Lambda[step_col>=cur.intv.l & step_col<cur.intv.u], 17)
  med.Lambda <- rbind(med.Lambda, as.numeric(median(cur.group)))
  
}

#########################
## glmer log regression #
#########################

## using the function to predict


sub_prd <- c()
sub_list <- as.numeric(names(table(data.la$subject)))

for (jj in 1: length(sub_list)) {
  tmp_data <- subset(data.la, subject == sub_list[jj])
  cur_coef <- subset(crm, subject==sub_list[jj])
  prob_tmp <- c()
  
  for (ii in 1:length(tmp_data[,1])) {
    tmp      <- prob.log.reg(cur_coef$rmIntercept,((cur_coef$rmGewinn)),cur_coef$rmVerlust,
                             tmp_data$Gewinn[ii],tmp_data$Verlust[ii])
    prob_tmp <- c(prob_tmp,as.numeric(tmp))
  }
  ## hier fehlt noch CODE!!! man muss die estimated probs überprüfen!
}