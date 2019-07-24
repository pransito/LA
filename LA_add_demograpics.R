############################################################
## get the covariates per subject and add them to data.la ##
############################################################

agk.load.ifnot.install("pracma")

## read in spss file with demographics and survey data
setwd(paste(pfad, "\\Data", sep = ""))
# only the included subjects
demographics <- read.spss("BGG_4Alex_original.sav", use.value.labels = TRUE, to.data.frame = TRUE,
                          max.value.labels = Inf, trim.factor.names = FALSE,
                          trim_values = TRUE, reencode = NA, use.missings = to.data.frame)
# all subjects
demogr.compl    <- read.spss("BGG_4Alex_original.sav", use.value.labels = TRUE, to.data.frame = F,
                          max.value.labels = Inf, trim.factor.names = FALSE,
                          trim_values = FALSE, reencode = NA, use.missings = to.data.frame)
demogr.compl.nl <- read.spss("BGG_4Alex_original.sav", use.value.labels = FALSE, to.data.frame = T,
                          max.value.labels = Inf, trim.factor.names = FALSE,
                          trim_values = FALSE, reencode = NA, use.missings = to.data.frame)

# recode GBQ_x items
names_GBQ     <- grep("^GBQ_[0-9]", names(demogr.compl), value=TRUE)
cur_y <- c(1:7,5.5)
cur_z <- c(7:1,2.5)
demogr.compl.nl[names_GBQ]   <- apply(demogr.compl.nl[names_GBQ],FUN = agk.recode,MARGIN = 2,y=cur_y,z=cur_z)
demogr.compl.nl$GBQ_mean_rec <- apply(demogr.compl.nl[names_GBQ],FUN = mean.rmna,MARGIN = 1)
tmp <- demogr.compl.nl[names_GBQ] 
demogr.compl.nl$GBQ_illus    <- apply(tmp[c(1,2,3,5,7,8,9,19)],FUN = mean.rmna,MARGIN = 1)
demogr.compl.nl$GBQ_persi    <- apply(tmp[c(4,6,10,11,12,13,14,15,16,17,18,20,21)],FUN = mean.rmna,MARGIN = 1)
demo_GBQ     <- demogr.compl.nl[c(c("VP_Code","GBQ_mean_rec","GBQ_persi","GBQ_illus"),names_GBQ)]
demographics <- merge(demographics, demo_GBQ, by.x="VP_Code", by.y="VP_Code", all.x =T)

# do my own BIS calculation
names_BIS     <- grep("^BIS_[0-9]", names(demogr.compl), value=TRUE)
rec_BIS       <- c()
for (kk in 1:34) {
  cur_name = names_BIS[kk]
  rec_name = paste0(cur_name,'_umkodiert')
  rec_id   = which(rec_name == names_BIS)
  if (length(rec_id) != 0) {
    rec_BIS[kk] = rec_name
  } else {
    rec_BIS[kk] = cur_name
  }
}
demogr.compl.nl$BIS_own_mean = apply(demogr.compl.nl[rec_BIS],FUN = mean.rmna,MARGIN = 1)
demo_BIS     = demogr.compl.nl[c("VP_Code","BIS_own_mean")]
demographics = merge(demographics, demo_BIS, by.x="VP_Code", by.y="VP_Code", all.x =T)
demographics$BIS_own_sum = demographics$BIS_own_mean * 34 # we used BIS-10, some missings in some items in some people, so extrapolate the sum; 34 items  

# do my own BIS calculation; BIS-11; it is the BIS-10 (which we used) but minus 4 items
items_to_drop = c(19,26,27,29)
names_BIS     <- grep("^BIS_[0-9]", names(demogr.compl), value=TRUE)
rec_BIS       <- c()
for (kk in 1:34) {
  cur_name = names_BIS[kk]
  rec_name = paste0(cur_name,'_umkodiert')
  rec_id   = which(rec_name == names_BIS)
  if (length(rec_id) != 0) {
    rec_BIS[kk] = rec_name
  } else {
    rec_BIS[kk] = cur_name
  }
}
rec_BIS = rec_BIS[-items_to_drop]
demogr.compl.nl$BIS_own_mean_BIS11 = apply(demogr.compl.nl[rec_BIS],FUN = mean.rmna,MARGIN = 1)
demo_BIS     = demogr.compl.nl[c("VP_Code","BIS_own_mean_BIS11")]
demographics = merge(demographics, demo_BIS, by.x="VP_Code", by.y="VP_Code", all.x =T)
demographics$BIS_own_sum_BIS11 = demographics$BIS_own_mean_BIS11 * 30 # some missings in some items in some people, so extrapolate the sum; 34 items  


## now get the components
setwd(pfad_scripts)
source("LA_get_components_2.R")
demographics <- merge(demographics, demogr.compl.nl[c("VP_Code","loc","noa")],by.x="VP_Code", by.y="VP_Code", all.x =T)

## get the merged Einkommen per subject
# best guess from two variables combined: household and personal income)
setwd(pfad_data)
VP_Einkommen_merged <- read.table("Einkommen_VP_merge.csv", sep = ";", header= T) ## CHANGE THIS TABLE!
demographics <- merge(demographics, VP_Einkommen_merged, by.x="VP_Code", by.y="Subject", all.y =T)

## one BDI score in the CT group which is missing has to be replaced by the mean
demographics_ct <- demographics[which(demographics$VP_Code < 2000),]
demographics$BDI_Summe[which(is.na(demographics$BDI_Summe))] <- mean(demographics_ct$BDI_Summe[!is.na(demographics_ct$BDI_Summe)])

## merger with data.la 
colnames(demographics)[1] <- c("subject")
# ACHTUNG!!! we change the name of subject 3030 into 2030 because it
# was a coding mistake in the scanner; in the behavioral data I changed this as well
# and also in the scanner data stored on my local hard drive
demographics$subject[demographics$subject == "3030"] <- c("2030")

data.la <- merge(data.la, demographics, by.x = "subject", by.y = "subject", all.x = TRUE)



