data.la_HC <- data.la[data.la$group == "HC",]
data.la_PG <- data.la[data.la$group == "PG",]
data.la_AD <- data.la[data.la$group == "AD",]

gmat_HC <- aggregate(as.numeric(as.matrix(data.la_HC$accept.reject)), by = list(data.la_HC$Gewinn, data.la_HC$Verlust),FUN = mean.rmna)
gmat_PG <- aggregate(as.numeric(as.matrix(data.la_PG$accept.reject)), by = list(data.la_PG$Gewinn, data.la_PG$Verlust),FUN = mean.rmna)
gmat_AD <- aggregate(as.numeric(as.matrix(data.la_AD$accept.reject)), by = list(data.la_AD$Gewinn, data.la_AD$Verlust),FUN = mean.rmna)

gmat_HC$take <- ifelse(gmat_HC$x < uncertainty.up & gmat_HC$x > uncertainty.lo,1,0)
gmat_PG$take <- ifelse(gmat_PG$x < uncertainty.up & gmat_PG$x > uncertainty.lo,1,0)
gmat_AD$take <- ifelse(gmat_AD$x < uncertainty.up & gmat_AD$x > uncertainty.lo,1,0)

gmat_HC$take_allgroups <- gmat_HC$take + gmat_PG$take + gmat_AD$take
gmat_HC$take_allgroups <- ifelse(gmat_HC$take_allgroups >0, 1,0)

# now generate a take variable for data.la
take_data.la <- rep(0,times = length(data.la[,1]))
takes <- list()
for (ii in 1:length(gmat_HC[,1])) {
  if(gmat_HC$take_allgroups[ii] == 0) {next}
  takes[[ii]] <- ifelse(data.la$Gewinn == gmat_HC$Group.1[ii] & data.la$Verlust == gmat_HC$Group.2[ii],1,0)
}
for (ii in 1:length(takes)) {
  if(length(takes[[ii]])==0) {next
  } else {take_data.la <- take_data.la + takes[[ii]]}
}

data.la <- data.la[take_data.la == 1,]
