
## cross validation ##
## f_K_fold function by http://stackoverflow.com/questions/7402313/generate-sets-for-cross-validation-in-r
## Wojciech Sobala

f_K_fold <- function(Nobs,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d) 
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}

## how many folds?
k = 5

## create folds
folds <- f_K_fold(length(data.la[,1]), K=k)

## train on data and testing k-fold whole group
pred_score <- c()
for (ii in 1:k){
  
  ind_train  <- folds[[ii]]$train
  data_train <- data.la[ind_train,]
  ind_test   <- folds[[ii]]$test
  data_test  <- data.la[ind_test,]
  
  curr_model <- glmer(accept.reject ~ (Gewinn + Verlust)*group +(Gewinn + Verlust|subject), data=data_train, family="binomial",
                      control=glmerControl(optimizer="bobyqa"))
  
  curr_targ  <- data_test$accept.reject
  preds      <- predict(curr_model,newdata = data_test,type="response")
  preds      <- round(preds)
  cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
  cor.preds  <- ifelse(cor.preds==0,1,0)
  
  pred_score[ii] <- mean(cor.preds,na.rm = T)
  
}


## train on data and testing k-fold HC

## create folds
folds <- f_K_fold(length(data.la.HC[,1]), K=k)

pred_score.HC <- c()
for (ii in 1:k){
  
  ind_train  <- folds[[ii]]$train
  data_train <- data.la.HC[ind_train,]
  ind_test   <- folds[[ii]]$test
  data_test  <- data.la.HC[ind_test,]
  
  curr_model <- glmer(accept.reject ~ (Gewinn + Verlust) +(Gewinn + Verlust|subject), data=data_train, family="binomial",
                      control=glmerControl(optimizer="bobyqa"))
  
  curr_targ  <- data_test$accept.reject
  preds      <- predict(curr_model,newdata = data_test,type="response")
  preds      <- round(preds)
  cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
  cor.preds  <- ifelse(cor.preds==0,1,0)
  
  pred_score.HC[ii] <- mean(cor.preds,na.rm = T)
  
}

## train on data and testing k-fold PG

## create folds
folds <- f_K_fold(length(data.la.PG[,1]), K=k)

pred_score.PG <- c()
for (ii in 1:k){
  
  ind_train  <- folds[[ii]]$train
  data_train <- data.la.PG[ind_train,]
  ind_test   <- folds[[ii]]$test
  data_test  <- data.la.PG[ind_test,]
  
  curr_model <- glmer(accept.reject ~ (Gewinn + Verlust) +(Gewinn + Verlust|subject), data=data_train, family="binomial",
                      control=glmerControl(optimizer="bobyqa"))
  
  curr_targ  <- data_test$accept.reject
  preds      <- predict(curr_model,newdata = data_test,type="response")
  preds      <- round(preds)
  cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
  cor.preds  <- ifelse(cor.preds==0,1,0)
  
  pred_score.PG[ii] <- mean(cor.preds,na.rm = T)
  
}

## train on data and testing k-fold AD

## create folds
folds <- f_K_fold(length(data.la.AD[,1]), K=k)

pred_score.AD <- c()
for (ii in 1:k){
  
  ind_train  <- folds[[ii]]$train
  data_train <- data.la.AD[ind_train,]
  ind_test   <- folds[[ii]]$test
  data_test  <- data.la.AD[ind_test,]
  
  curr_model <- glmer(accept.reject ~ (Gewinn + Verlust) +(Gewinn + Verlust|subject), data=data_train, family="binomial",
                      control=glmerControl(optimizer="bobyqa"))
  
  curr_targ  <- data_test$accept.reject
  preds      <- predict(curr_model,newdata = data_test,type="response")
  preds      <- round(preds)
  cor.preds  <- abs((as.numeric(curr_targ)-1)-preds)
  cor.preds  <- ifelse(cor.preds==0,1,0)
  
  pred_score.AD[ii] <- mean(cor.preds,na.rm = T)
  
}