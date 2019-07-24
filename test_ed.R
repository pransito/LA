tmp = lmList(accept.reject ~ Gewinn + Verlust + ed.abs | subject, data = data.la, family = "binomial",na.action=NULL,pool=F)
tmp = coef(tmp)
tmp_lambda = (tmp$Verlust)*(-1)/tmp$Gewinn
cor.test(crm$lambda,tmp_lambda,method="spearman")
plot(crm$lambda,tmp_lambda)

tmp_data = subset(data.la,subject =="1001")
tmp_data = tmp_data[!is.na(tmp_data$accept.reject),]
write.table(tmp_data$ed.abs,file="tmp.csv",sep = ";",dec=",")
