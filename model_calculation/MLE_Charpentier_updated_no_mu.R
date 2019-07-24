# Calculates Beta_0, Beta_Gain and Beta_Loss parameters -> see Charpentier et al. 2015 or Charp. et al. 2012
# but no mu!

# starting values for model params
lambdas  = 2
all_subs = unique(data_pdt$subject)

# function to compute log-likelihood given certain model params
fn<- function(theta,x){
  lambda    = theta[1];
  res1<-0
  for (k in 1:nrow(x)){ #go through all trials
    if (is.na(x$accnum[k])){#skip trials without reaction
    } 
    else{ 
      value    = 0.5*x$gain[k]+0.5*x$loss[k]*lambda;
      prob_acc = (1+exp(-1*value))^-1;
      if(prob_acc > 0.9999999999999999) {
        prob_acc = 0.9999999999999999
      }
      res<-x$accnum[k]*log(prob_acc)+(1-x$accnum[k])*log(1-prob_acc);
      
      res1<-res1-res;
    }
    
  }
  return(res1)
}
fn.c = cmpfun(fn)

# maximization function
do_maxim = function(cur_sub,lambdas,x,fn.c) {
  y3       = as.character(cur_sub)
  x$acc    = x$accept_reject
  x$accnum = as.numeric(as.character(x$acc))
  
  x1_list  = list()
  x1_value = c()
  for (kk in 1:length(lambdas)) {
    cur_theta = c(lambdas[kk])
    x1_list[[kk]] = optim(cur_theta,fn.c,x=x,method = "Brent",lower=-11,upper=11)     # maximize likelihood of parameters 
    x1_value[kk]  = x1_list[[kk]]$value
  }
  cur_winner = which(min(x1_value) == x1_value) # get global minimums of ML-estimates
  cur_winner = cur_winner[length(cur_winner)]   # if more than one have same min, take highest mu
  x1 = x1_list[[cur_winner]]
  x2 = list()
  x2$name   = y3
  x2$params = c(lambda = as.numeric(x1$par[1]))
  
  return(x2)
}

# parallel processing
if (estimate_models == 1) {
  cl<-makeCluster(5)  
  registerDoSNOW(cl)
  maxim_results <- foreach(ff=1:length(all_subs),.verbose=T) %dopar% {
    x = subset(data_pdt, subject == all_subs[ff])
    do_maxim(all_subs[ff],lambdas,x,fn.c)
  } 
  stopCluster(cl)
  
  # bind all params after all subs have been fit
  all_params = t(as.matrix(maxim_results[[1]]$params))
  VPNR       = c()
  VPNR[1]    = maxim_results[[1]]$name
  for (ii in 2:length(maxim_results)) {
    VPNR[ii]   = maxim_results[[ii]]$name
    all_params = rbind(all_params,maxim_results[[ii]]$params)
  }
  all_params = as.data.frame(all_params)
  all_params = cbind(VPNR,all_params)
  
  # saving
  save(file="Charpentier_model_no_mu.RData",list = c("all_params","fn.c"))
} else {
  load("Charpentier_model_no_mu.RData")
}


# put results in dataframe and save as csv
#write.table(all_params, "C:/Users/genaucka/Google Drive/Promotion/VPPG/VPPG_Exchange/Library/MLE+Model/LA_Charpentier/results_CH.txt", sep="\t")
