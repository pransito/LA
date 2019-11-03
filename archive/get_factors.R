## CFA
# LOC
all_items = all_names_PG[as.logical(loc_PG)]

cur_str = paste('PG_loc_fac =~', all_items[1])

for(ii in 2:length(all_items)) {
  cur_str = paste(cur_str,all_items[ii],sep='+')  
}

# NOA
all_items_noa = all_names_PG[as.logical(noa_PG)]

cur_str_noa = paste('PG_noa_fac =~', all_items_noa[1])

for(ii in 2:length(all_items_noa)) {
  cur_str_noa = paste(cur_str_noa,all_items_noa[ii],sep='+')  
}

loc_noa_mod = (paste(cur_str,
                    cur_str_noa,sep='\n'))

fit = cfa(loc_noa_mod,cur_df)

## EFA
# PG
cor_PG             <- cor(cur_df)
cor_PG[cor_PG<0.2] <- 0 
tmp=fa(cor_PG,n.obs = length(cur_df[,1]),nfactors = 2,covar = F,rotate = "varimax")

cor.test(tmp$loadings[,1],noa_PG)
cor.test(tmp$loadings[,2],loc_PG)

# HC
cor_HC             <- cor(cur_df)
cor_HC[cor_HC<0.2] <- 0 
tmp=fa(cor_HC,n.obs = length(cur_df[,1]),nfactors = 2,covar = F,rotate = "varimax")

cor.test(tmp$loadings[,2],noa_HC)
cor.test(tmp$loadings[,1],loc_HC)

# AD
cor_AD             <- cor(cur_df)
cor_AD[cor_AD<0.2] <- 0 
tmp=fa(cor_AD,n.obs = length(cur_df[,1]),nfactors = 2,covar = F)

cor.test(tmp$loadings[,2],noa_AD)
cor.test(tmp$loadings[,1],loc_AD)

## FASTICA

# PG
pca   = eigen(cov(cur_df))
cur_k = sum((pca$values<1)==FALSE)
tmp   = fastICA(X = cur_df,n.comp = 6)

cri_loc = t(t(loc_PG) %*% t(cur_df))
cri_noa = t(t(noa_PG) %*% t(cur_df))

corr.test(as.data.frame(tmp$S),as.data.frame(cri_loc),method="pearson",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_noa),method="pearson",adjust = "none")

# LOC: 5 NOA: 1

# HC
pca   = eigen(cov(cur_df))
plot(pca$values)
tmp   = fastICA(X = cur_df,n.comp = 6)

cri_loc = t(t(loc_HC) %*% t(cur_df))
cri_noa = t(t(noa_HC) %*% t(cur_df))

corr.test(as.data.frame(tmp$S),as.data.frame(cri_loc),method="pearson",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_noa),method="pearson",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_loc),method="spearman",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_noa),method="spearman",adjust = "none")

# LOC: 6 NOA: 3

# AD
pca   = eigen(cov(cur_df))
plot(pca$values)
tmp   = fastICA(X = cur_df,n.comp = 7)

cri_loc = t(t(loc_AD) %*% t(cur_df))
cri_noa = t(t(noa_AD) %*% t(cur_df))

corr.test(as.data.frame(tmp$S),as.data.frame(cri_loc),method="pearson",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_noa),method="pearson",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_loc),method="spearman",adjust = "none")
corr.test(as.data.frame(tmp$S),as.data.frame(cri_noa),method="spearman",adjust = "none")

# LOC: 4 NOA: 6

## OLD ##
cur_unmix = tmp$K%*%tmp$W
zeroitems    = ((loc_PG == 0) == (noa_PG == 0))==FALSE
cur_unmix = cur_unmix[zeroitems,]
loc_PGp = loc_PG[zeroitems]
noa_PGp = noa_PG[zeroitems]

corr.test(as.data.frame(cur_unmix),as.data.frame(loc_PGp))
corr.test(as.data.frame(cur_unmix),as.data.frame(noa_PGp))

cor.test(loc_PGp,noa_PGp)


