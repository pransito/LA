## creating loss of control & neglect of other areas in life components

# set dist fun for component vector vs. found component (cor or dist or none)
# if none then extracted factors are neglected, instead the loc and noa components
# are in fact the rule to compute the loc and noa scores
dist_fun = "none"
comp_exp = "Nina"

# get the file with component rules
#setwd(pfad_scripts)
#locnoaalexnina <- read.table(file="loc_noa_components_Nina_Alex_Nina.txt",header = T,sep="\t")

# function to select the component for loc and noa
# if the (first) component is picked for both then go to next component in both cases
# until they differ
agk.select.component <- function(tmp_loc,tmp_noa) {
  control_var <- 0
  while (control_var == 0) {
    PCloc <- which(max(abs(tmp_loc))==abs(tmp_loc))
    PCnoa <- which(max(abs(tmp_noa))==abs(tmp_noa))
    if (PCloc != PCnoa) {control_var <- 1} else {
      tmp_loc <- tmp_loc[-1]
      tmp_noa <- tmp_noa[-1]
    }
  }
  return(list(PCloc=PCloc,PCnoa=PCnoa))
}

# get the variables
var_labels <- attr(demogr.compl,which = "variable.labels")

names_KFG     <- grep("^KFG_[0-9]", names(demogr.compl), value=TRUE)
names_GBQ     <- grep("^GBQ_[0-9]", names(demogr.compl), value=TRUE)
names_BIS     <- grep("^BIS_[0-9]{,6}$", names(demogr.compl), value=TRUE)
names_BIS_BAS <- grep("^BIS_BAS_[0-9]{,10}$", names(demogr.compl), value=TRUE)
names_GSAS    <- grep("^GSAS_[0-9]{,7}$", names(demogr.compl), value=TRUE)
names_PGYBOCS <- grep("^PGYBOCS_[0-9]{,10}$", names(demogr.compl), value=TRUE)
names_debt    <- grep("^Schulden_amount",names(demogr.compl), value=TRUE)
names_smoking <- grep("^Rauchen_upd",names(demogr.compl), value=TRUE)
names_BDI     <- grep("^BDI2_",names(demogr.compl), value=TRUE)
names_OCDS    <- grep("^OCDS_[0-9]",names(demogr.compl), value=TRUE)
names_ADS     <- grep("^ADS",names(demogr.compl), value=TRUE)
names_detox   <- grep("^Entgiftung",names(demogr.compl), value=TRUE)
all_names_PG  <- c(names_GBQ,names_BIS,names_BIS_BAS,names_GSAS,names_KFG,names_PGYBOCS,names_debt,names_BDI)
all_names_HC  <- c(names_GBQ,names_BIS,names_BIS_BAS,names_debt,names_smoking,names_BDI)
all_names_AD  <- c(names_GBQ,names_BIS,names_BIS_BAS,names_debt,names_smoking,names_BDI,names_OCDS,names_ADS,names_detox)

# # drop BDI2_14 in HC (no variance)
# which_to_drop <- which(all_names_HC=="BDI2_14")
# all_names_HC  <- all_names_HC[-which_to_drop]

# drop GBQ 20 in AD (no variance)
which_to_drop <- which(all_names_HC=="GBQ_20")
all_names_AD  <- all_names_AD[-which_to_drop]

# # Alex or Ninas collection?
# if (comp_exp == "Nina") {
#   loc_PG <- locnoaalexnina$loc_nina[locnoaalexnina$item %in% all_names_PG]
#   noa_PG <- locnoaalexnina$noa_nina[locnoaalexnina$item %in% all_names_PG]
#   loc_AD <- locnoaalexnina$loc_nina[locnoaalexnina$item %in% all_names_AD]
#   noa_AD <- locnoaalexnina$noa_nina[locnoaalexnina$item %in% all_names_AD]
#   loc_HC <- locnoaalexnina$loc_nina[locnoaalexnina$item %in% all_names_HC]
#   noa_HC <- locnoaalexnina$noa_nina[locnoaalexnina$item %in% all_names_HC]  
# } else if (comp_exp == "Alex") {
#   loc_PG <- locnoaalexnina$loc_alex[locnoaalexnina$item %in% all_names_PG]
#   noa_PG <- locnoaalexnina$noa_alex[locnoaalexnina$item %in% all_names_PG]
#   loc_AD <- locnoaalexnina$loc_alex[locnoaalexnina$item %in% all_names_AD]
#   noa_AD <- locnoaalexnina$noa_alex[locnoaalexnina$item %in% all_names_AD]
#   loc_HC <- locnoaalexnina$loc_alex[locnoaalexnina$item %in% all_names_HC]
#   noa_HC <- locnoaalexnina$noa_alex[locnoaalexnina$item %in% all_names_HC]
# }

# # names BDI HC adjustment
# names_BDI_HC = names_BDI[-which(names_BDI == "BDI2_14")]


####################################################
### get the components #############################
####################################################
demogr.compl.nl$loc <- NA
demogr.compl.nl$noa <- NA
max_fac <- 20

## for PG
cur_df <- subset(demogr.compl.nl,Gruppe==2)
cur_df$Schulden_amount <- log(cur_df$Schulden_amount+0.01)
cur_df <- as.data.frame(cur_df[all_names_PG])
cur_df <- apply(cur_df,MARGIN = 2,FUN = agk.impute.mean)
cur_df <- scale(cur_df)

# get pr. components in data
pca     <- eigen(cov(cur_df))
cur_rot <- pca$vectors

# normalize
cur_rot <- apply(cur_rot,FUN=agk.norm.range,MARGIN = 2,x=-1,y=1)

# measure distance and pick component
if (dist_fun == "dist") {
  tmp_loc <- (distmat(t(as.matrix(loc_PG)),t(cur_rot)))[which(!pca$values<=1)]
  tmp_noa <- (distmat(t(as.matrix(noa_PG)),t(cur_rot)))[which(!pca$values<=1)]
  PCloc <- which(min(tmp_loc)==tmp_loc)
  PCnoa <- which(min(tmp_noa)==tmp_noa)
} else if (dist_fun == "cor") {
  # polyserial (# same result like normal correlation)
  tmp_loc <- (polyserial(as.data.frame(cur_rot),as.data.frame(loc_PG)))[which(!pca$values<=1)]
  tmp_noa <- (polyserial(as.data.frame(cur_rot),as.data.frame(noa_PG)))[which(!pca$values<=1)]
  # select the component
  locnoa <- agk.select.component(tmp_loc,tmp_noa)
  PCloc <- locnoa$PCloc
  PCnoa <- locnoa$PCnoa
} else if (dist_fun == "rbf") {
  tmp_loc <- (kernelMatrix(t(as.matrix(loc_PG)),t(cur_rot),kernel=rbf))[which(!pca$values<=1)]
  tmp_noa <- (kernelMatrix(t(as.matrix(noa_PG)),t(cur_rot),kernel=rbf))[which(!pca$values<=1)]
  PCloc <- which(max(abs(tmp_loc))==abs(tmp_loc))
  PCnoa <- which(max(abs(tmp_noa))==abs(tmp_noa)) 
}

# compute score and attach
if (dist_fun != "none") {
  demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 2] <- t(t(pca$vectors[,PCloc]) %*% t(cur_df))
  demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 2] <- t(t(pca$vectors[,PCnoa]) %*% t(cur_df))
} else {
  demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 2] <- t(t(loc_PG) %*% t(cur_df))
  demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 2] <- t(t(noa_PG) %*% t(cur_df))
}


## for HC
cur_df <- subset(demogr.compl.nl,Gruppe==1)
cur_df$Schulden_amount <- log(cur_df$Schulden_amount+0.01)
cur_df <- as.data.frame(cur_df[all_names_HC])
cur_df <- apply(cur_df,MARGIN = 2,FUN = agk.impute.mean)
cur_df <- as.data.frame(scale(cur_df))
# impute for BDI the min(BDI) at nan(BDI)
cur_min <- min(cur_df[names_BDI],na.rm = T)
f <- function(x) {ifelse(is.nan(x),cur_min,identity(x))}
cur_df[names_BDI] <- apply(cur_df[names_BDI],MARGIN = 2,FUN=f)

  
# get components in data
pca <- eigen(cov(cur_df))
cur_rot <- pca$vectors

# normalize
cur_rot <- apply(cur_rot,FUN=agk.norm.range,MARGIN = 2,x=-1,y=1)

# measure distance and pick component
if (dist_fun == "dist") {
  tmp_loc <- (distmat(t(as.matrix(loc_HC)),t(cur_rot)))[which(!pca$values<=1)]
  tmp_noa <- (distmat(t(as.matrix(noa_HC)),t(cur_rot)))[which(!pca$values<=1)]
  PCloc <- which(min(tmp_loc)==tmp_loc)
  PCnoa <- which(min(tmp_noa)==tmp_noa)
} else if (dist_fun == "cor") {
  # polyserial (# same result like normal correlation)
  tmp_loc <- (polyserial(as.data.frame(cur_rot),as.data.frame(loc_HC)))[which(!pca$values<=1)]
  tmp_noa <- (polyserial(as.data.frame(cur_rot),as.data.frame(noa_HC)))[which(!pca$values<=1)]
  # select the component
  locnoa <- agk.select.component(tmp_loc,tmp_noa)
  PCloc <- locnoa$PCloc
  PCnoa <- locnoa$PCnoa
} else if (dist_fun == "rbf") {
  tmp_loc <- (kernelMatrix(t(as.matrix(loc_HC)),t(cur_rot),kernel=rbf))[which(!pca$values<=1)]
  tmp_noa <- (kernelMatrix(t(as.matrix(noa_HC)),t(cur_rot),kernel=rbf))[which(!pca$values<=1)]
  PCloc <- which(max(abs(tmp_loc))==abs(tmp_loc))
  PCnoa <- which(max(abs(tmp_noa))==abs(tmp_noa)) 
}

if (dist_fun != "none") {
  demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 1] <- t(t(pca$vectors[,PCloc]) %*% t(cur_df))
  demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 1] <- t(t(pca$vectors[,PCnoa]) %*% t(cur_df))
} else {
  demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 1] <- t(t(loc_HC) %*% t(cur_df))
  demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 1] <- t(t(noa_HC) %*% t(cur_df))
}

## for AD
cur_df <- subset(demogr.compl.nl,Gruppe==3)
cur_df$Schulden_amount <- log(cur_df$Schulden_amount+0.01)
cur_df <- as.data.frame(cur_df[all_names_AD])
cur_df <- apply(cur_df,MARGIN = 2,FUN = agk.impute.mean)
cur_df <- as.data.frame(scale(cur_df))

# get components in data
pca <- eigen(cov(cur_df))
cur_rot <- pca$vectors

# normalize
cur_rot <- apply(cur_rot,FUN=agk.norm.range,MARGIN = 2,x=-1,y=1)

# measure distance and pick component
if (dist_fun == "dist") {
  tmp_loc <- (distmat(t(as.matrix(loc_AD)),t(cur_rot)))[which(!pca$values<=1)]
  tmp_noa <- (distmat(t(as.matrix(noa_AD)),t(cur_rot)))[which(!pca$values<=1)]
  PCloc <- which(min(tmp_loc)==tmp_loc)
  PCnoa <- which(min(tmp_noa)==tmp_noa)
} else if (dist_fun == "cor") {
  # polyserial (# same result like normal correlation)
  tmp_loc <- (polyserial(as.data.frame(cur_rot),as.data.frame(loc_AD)))[which(!pca$values<=1)]
  tmp_noa <- (polyserial(as.data.frame(cur_rot),as.data.frame(noa_AD)))[which(!pca$values<=1)]
  # select the component
  locnoa <- agk.select.component(tmp_loc,tmp_noa)
  PCloc <- locnoa$PCloc
  PCnoa <- locnoa$PCnoa
}

if (dist_fun != "none") {
demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 3] <- t(t(pca$vectors[,PCloc]) %*% t(cur_df))
demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 3] <- t(t(pca$vectors[,PCnoa]) %*% t(cur_df))
} else {
  demogr.compl.nl$loc[demogr.compl.nl$Gruppe == 3] <- t(t(loc_AD) %*% t(cur_df))
  demogr.compl.nl$noa[demogr.compl.nl$Gruppe == 3] <- t(t(noa_AD) %*% t(cur_df))
}







### OLD ###

# # plot
# cur_PC <- pca$rotation[,3]
# cur_name <- "PC3"
# cur_pl_df <- data.frame(all_names_PG,cur_PC)
# p1 <- ggplot(cur_pl_df, aes(all_names_PG, cur_PC)) + geom_point()
# p1 <- p1 + theme(axis.text.x=element_text(angle=-90, size =10, color = "black"))
# tmp_title <- paste("Loadings of", cur_name)
# p1 <- p1 + ggtitle(tmp_title)
# p1
# 
# # entropy
# barplot(apply(pca$rot,MARGIN = 2,FUN = agk.entr))

# fit <- principal(cur_df, nfactors=10,rotate = "oblimin")

# tmp_loc <-apply(cur_rot,FUN = agk.corboxmetric,MARGIN = 2,cur_box=loc_PG)
# tmp_noa <-apply(cur_rot,FUN = agk.corboxmetric,MARGIN = 2,cur_box=noa_PG)


#   tmp_loc <- (corr.test(as.data.frame(loc_PG),as.data.frame(cur_rot),adjust = "none"))[which(!pca$values<=1)]
#   tmp_noa <- (corr.test(as.data.frame(noa_PG),as.data.frame(cur_rot),adjust = "none"))[which(!pca$values<=1)]
#   PCloc <- which(max(abs(tmp_loc$r))==abs(tmp_loc$r))
#   PCnoa <- which(max(abs(tmp_noa$r))==abs(tmp_noa$r))
