# Calculates Beta_0, Beta_Gain and Beta_Loss parameters -> see Charpentier et al. 2015 or Charp. et al. 2012

# starting values for model params
mus         = seq(0.2,30,by = 6)
beta_ints   = rep(0,5)
beta_ras    = rep(2,5)
all_subs    = unique(data_pdt$subject)

# function to compute log-likelihood given certain model params
fn<- function(theta,x){
  mu        = theta[1];
  beta_int  = theta[2];
  beta_ra   = theta[3];
  res1 = 0
  for (k in 1:nrow(x)){ #go through all trials
    if (is.na(x$accnum[k])){#skip trials without reaction
    } 
    else{ 
      value    = beta_int+x$ratio_bcp[k]*beta_ra;
      prob_acc = (1+exp(-mu*value))^-1;
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
do_maxim = function(cur_sub,mus,beta_ints,beta_ras,x,fn.c) {
  y3       = as.character(cur_sub)
  x$acc    = x$accept_reject
  x$accnum = as.numeric(as.character(x$acc))
  
  x1_list  = list()
  x1_value = c()
  for (kk in 1:length(mus)) {
    cur_theta = c(mus[kk],beta_ints[kk],beta_ras[kk])
    x1_list[[kk]] = optim(cur_theta,fn.c,x=x)     # maximize likelihood of parameters 
    x1_value[kk]  = x1_list[[kk]]$value
  }
  cur_winner = which(min(x1_value) == x1_value) # get global minimums of ML-estimates
  cur_winner = cur_winner[length(cur_winner)]   # if more than one have same min, take highest mu
  x1 = x1_list[[cur_winner]]
  x2 = list()
  x2$name   = y3
  x2$params = c(mu        = as.numeric(x1$par[1]),
                beta_int  = as.numeric(x1$par[2]),
                beta_ra   = as.numeric(x1$par[3]))
  
  return(x2)
}

if (estimate_models == 1) {
  # parallel processing
  cl<-makeCluster(5)  
  registerDoSNOW(cl)
  maxim_results <- foreach(ff=1:length(all_subs),.verbose=T) %dopar% {
    x = subset(data_pdt, subject == all_subs[ff])
    do_maxim(all_subs[ff],mus,beta_ints,beta_ras,x,fn.c)
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
  save(file="larm.RData",list = c("all_params","fn.c"))
} else {
  load("larm.RData")
}