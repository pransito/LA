# cross validation function for a glmer model
# given a single random effect unit (subject)
# k-fold
# one 20% subjects will be left out and the data fit
# but during prediction we would only use the fixed effects; right?!
# what does the predict function do in such a case?

agk.cv.glmer <- function(model,byvar,k) {
  
  # get data
  dat <- model@frame
  
  # get the levels of the byvar
  lev <- eval(parse(text=paste("levels(dat$",byvar,")",sep="")))
  
  # create folds
  flds <- f_K_fold(length(lev),k)
  f = function(.,x) {any(.==x)}
  
  # get the column that has the byvar that the folds will be done over ("subject")
  bycol <- eval(parse(text=(paste("as.matrix(dat$",byvar,")",sep=""))))
  
  pred_score <- c()
  for (ii in 1:k){
    print(paste("fold...",ii))
    
    # create train data
    ind_train  <- flds[[ii]]$train
    ind_train  <- apply(bycol,1,f,x=lev[ind_train])
    data_train <- dat[ind_train,]
    
    # create test data
    ind_test   <- flds[[ii]]$test
    ind_test   <- apply(bycol,1,f,x=lev[ind_test])
    data_test  <- dat[ind_test,]
    
    # training
    curr_model <- update(model,data = data_train)
    
    # testing
    curr_targ  <- data_test$accept.reject
    preds      <- predict(curr_model,newdata = data_test,type="response",re.form=~0)
    preds      <- round(preds)
    cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
    cor.preds  <- ifelse(cor.preds==0,1,0)
    
    pred_score[ii] <- mean(cor.preds,na.rm = T)
  }
  EG <- mean(pred_score)
  return(EG)
}
