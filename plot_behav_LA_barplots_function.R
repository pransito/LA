# function to plot a boot object and save the plot (with CI)
agk.barplot.boot <- function(boot_pdt_df, cur_n, cur_cat, nsim,study_name) {
  
  temp        <- 1:length(boot_pdt_df[,1])
  boot_pdt_df <- data.frame(boot_pdt_df,temp)
  boot_pdt_df <- melt(boot_pdt_df,id.vars = c("temp"))
  boot_pdt_df$variable <- factor(boot_pdt_df$variable,ordered = T)
  
  # now make a grouping variable
  # first discard all the contrasts; won't plot them here
  # also intercepts
  to_discard  <- grep("gr",as.character(boot_pdt_df$variable))
  boot_pdt_df <- boot_pdt_df[-to_discard,]
  to_discard  <- grep("in_",as.character(boot_pdt_df$variable))
  boot_pdt_df <- boot_pdt_df[-to_discard,]
  #to_discard  <- grep("AD",as.character(boot_pdt_df$variable))
  #boot_pdt_df <- boot_pdt_df[-to_discard,]
  # make the grouping var
  ADgroup <- grep("AD",as.character(boot_pdt_df$variable))
  PGgroup <- grep("PG",as.character(boot_pdt_df$variable))
  HCgroup <- grep("HC",as.character(boot_pdt_df$variable))
  grouping_var <- c()
  for (kk in 1:length(boot_pdt_df$variable)) {
    if (any(ADgroup==kk)) {
      grouping_var[kk] <- "AD" 
    } else if (any(PGgroup==kk)) {
      grouping_var[kk] <- "PG" 
    } else if (any(HCgroup == kk)) {
      grouping_var[kk] = "HC"
    }
  }
  
  boot_pdt_df <- data.frame(boot_pdt_df,as.factor(grouping_var))
  names(boot_pdt_df)[4] <- "group"
  boot_pdt_df$group = factor(boot_pdt_df$group, levels=c("HC","PG","AD"), labels = c("HC","PG","AD"))
  
  # rename the variables
  boot_pdt_df$variable <- as.character(boot_pdt_df$variable)
  for (kk in 1:length(boot_pdt_df[,1])) {
    boot_pdt_df$variable[kk] <- strsplit(as.character(boot_pdt_df$variable[kk]),"_")[[1]][1]
  }
  
  boot_pdt_df <- as.data.frame(acast(boot_pdt_df, group + temp  ~ variable, value.var = "value"))
  boot_pdt_df$group <- row.names(boot_pdt_df)
  
  boot_pdt_df$group <- as.character(boot_pdt_df$group)
  for (kk in 1:length(boot_pdt_df$group)) {
    boot_pdt_df$group[kk] <- strsplit(as.character(boot_pdt_df$group[kk]),"_")[[1]][1]
    
  }
  
  # get the mean and errorbar values
  # reorder
  #boot_pdt_df <- boot_pdt_df[,c(3,1,2,4,5)]
  plot.dat <- data.frame()
  var.names <- c()
  for (ii in 1:(length(boot_pdt_df[1,])-1)){
    # bootstrapped 95% CI
    bt.df     <- as.data.frame(as.list(aggregate(boot_pdt_df[ii],by = list(group = boot_pdt_df$group), agk.mean.quantile,lower=0.025,upper=0.975)))
    names(bt.df) <- c("group", "mean", "lower","upper")
    var.names <- c(var.names,rep(names(boot_pdt_df[ii]),length(levels(bt.df$group))))
    plot.dat <- rbind(plot.dat,bt.df)
  }
  plot.dat$variable <- var.names
  plot.dat$variable <- as.factor(plot.dat$variable)
  #plot.dat$variable <- factor(plot.dat$variable,levels(plot.dat$variable)[c(3,1,2,4)])
  #plot.dat$group    <- factor(plot.dat$group,levels(plot.dat$group)[c(2,3,1)])
  plot.dat$group = factor(plot.dat$group, levels=c("HC","PG","AD"), labels = c("HC","PG","AD"))
  
  # plot
  p <- ggplot(data = plot.dat, aes(group,mean,fill=variable))
  p <- p+geom_bar(position="dodge",stat="identity")
  p
  p <- p + geom_errorbar(aes(ymin=lower,ymax=upper), size=1.3, width=0,
                         position=position_dodge(width=0.9), color=cbbPalette[4],
                         width=0.1) + ylab("mean (95% boots. CI)\n")
  
  p <- p + ggtitle("Fixed effects \n")
  
  return(p)
}