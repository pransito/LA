## complete cov models
all_modcc <- list()

ct         <- length(all_modcc)+1
vars       <- c("RiskMar", "ratio","diff", "EV", "Gewinn+Verlust", "RiskMar+ed.abs", "ratio+ed.abs","diff+ed.abs", "EV+ed.abs", "Gewinn+Verlust+ed.abs")

group      <- c()
group[1]   <- "*("
group[2]   <- "*(group"

incl       <- list()
incl[[1]]  <- c("Age)", "Bildungsjahre_ges)", "Age+Bildungsjahre_ges)")
incl[[2]]  <- c("+Age)", "+Bildungsjahre_ges)", "+Age+Bildungsjahre_ges)")

for (ii in 1:length(vars)) {
  text1    <- paste("glmer(accept.reject ~ (",vars[ii],")",sep="")
  text2    <- paste("+ (",vars[ii],"|subject), data=data.la, family='binomial',control=glmerControl(optimizer='bobyqa'))",sep="")
  for (jj in 1:length(group)) {
    for (kk in 1:length(incl)) {
      all_modcc[[ct]] <- paste(text1,group[jj],incl[[jj]][kk],text2,sep="")
      ct <- ct +1
    }
  }
}
