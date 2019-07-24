# neural correlates

## neural correlates
# get the neural correlates
# fmri
path_eigv_fmri <- "F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\sig_results\\sig_res_3_noacc\\owancova\\eigenvariates\\export"
setwd(path_eigv_fmri)
fmri_dat <- read.table("eigv_export.dat",sep="\t",header=T)
nms_anly <- read.table("names_analyses.dat",sep="\t",header=F)

# get ME and EOI out
throw_out <- (-1)*c(grep(pattern = "ME", as.character(nms_anly$V1)),grep(pattern = "EOI", as.character(nms_anly$V1)))
if (length(throw_out) != 0) {
  fmri_dat <- fmri_dat[,throw_out]
  nms_anly <- as.data.frame(nms_anly[throw_out,])
  names(nms_anly) <- c("V1") }

fmri_dat <- fmri_dat[,-2]
names(fmri_dat) <- c("rmfg","lmfg")

# get crm
setwd(pfad_results)
# save(file="cur_mod.Rdata",list = c("cur_mod"))
load("cur_mod.Rdata")
crm <- agk.get.compl.coef(cur_mod,"group") # no covs
fmri_dat$group =crm$group
fmri_dat$group = factor(fmri_dat$group,levels(fmri_dat$group)[c(2,3,1)])

# get a subject variable
agk.get.lastdir = function(cur_str) {
  all_str = unlist(strsplit(cur_str,"[\\]"))
  n = length(all_str)
  return(all_str[n])
}
pfad_lapipeline = "C:\\Users\\genaucka\\Google Drive\\Library\\MATLAB\\LA"
setwd(pfad_lapipeline)
all_subs = readMat('base_dirs.mat')
all_subs = all_subs$P.tmp[[4]]
alll_subs = rbind(unlist(all_subs[[1]][[1]]),unlist(all_subs[[2]][[1]]),unlist(all_subs[[3]][[1]]))
all_subs = unlist(lapply(alll_subs,FUN=agk.get.lastdir))
fmri_dat$subjects = all_subs

lac_subs = readMat('base_dirs_control.mat')
lac_subs = lac_subs$P.tmp[[4]]
lacc_subs = rbind(unlist(lac_subs[[1]][[1]]),unlist(lac_subs[[2]][[1]]),unlist(lac_subs[[3]][[1]]))
lac_subs = unlist(lapply(lacc_subs,FUN=agk.get.lastdir))

# get the GMD
pfad_GMD = 'F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\GMD_extract'
setwd(pfad_GMD)
cur_name = 'VOI_neg39.5_38.5_22.5_LeftMFG.mat'
cur_extr = readMat(paste(getwd(),cur_name,sep="/"))
fmri_dat$gmd_lmfg = cur_extr$Y

cur_name = 'VOI_neg36_neg42_54_LeftSPLs.mat_ls.mat'
cur_extr = readMat(paste(getwd(),cur_name,sep="/"))
fmri_dat$gmd_lspl = cur_extr$Y

# cur_name = 'VOI_neg29_neg42_65.mat'
# cur_extr = readMat(paste(getwd(),cur_name,sep="/"))
# fmri_dat$gmd_lspl2 = cur_extr$Y

cur_name = 'VOI_55_neg42_38.mat'
cur_extr = readMat(paste(getwd(),cur_name,sep="/"))
fmri_dat$gmd_rsmg = cur_extr$Y

# get the eigenvariates of functional contrasts
pfad_SPL = 'F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\owancova\\ep2d_bold_LA_00\\nocov\\loss\\3_noacc'
setwd(pfad_SPL)
cur_name = 'VOI_neg36_neg42_54.mat'
cur_extr = readMat(paste(getwd(),cur_name,sep="/"))
fmri_dat$lspl = cur_extr$Y

cur_name = 'VOI_neg40_39_23.mat'
cur_extr = readMat(cur_name)
fmri_dat$lmfg = cur_extr$Y

cur_name = 'VOI_55_neg42_37.mat'
cur_extr = readMat(cur_name)
fmri_dat$rsmg = cur_extr$Y

# get SPL and MFG Lac
pfad_SPL = 'F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\owancova\\ep2d_bold_LA_control_00\\nocov\\loss\\3_noacc'
setwd(pfad_SPL)
cur_name = 'VOI_neg36_neg42_54.mat'
cur_extr = readMat(cur_name)
cur_lac  = data.frame(lac_subs,cur_extr$Y)
names(cur_lac) = c("subject","lspl_lac")

cur_name = 'VOI_neg40_39_23.mat'
cur_extr = readMat(cur_name)
cur_lac$lmfg_lac = cur_extr$Y

cur_name = 'VOI_55_neg42_37.mat'
cur_extr = readMat(cur_name)
cur_lac$rsmg_lac = cur_extr$Y

fmri_dat_lac = merge(fmri_dat,cur_lac,by.x=c("subjects"),by.y = c("subject"))

# get the lMFG and lSPL data controlled for GMD using robust regression
require(foreign)
require(MASS)

# LMFG
rr_bisquare <- rlm(lmfg ~ gmd_lmfg, data=fmri_dat, psi = psi.bisquare)
mod0        <-  lm(lmfg ~ gmd_lmfg, data=fmri_dat)
summary(rr_bisquare)
summary(mod0)
plot(fmri_dat$gmd_lmfg,fmri_dat$lmfg)
abline(rr_bisquare$coefficients[1],rr_bisquare$coefficients[2])
abline(mod0$coefficients[1],mod0$coefficients[2],col="red")

fmri_dat$lmfg_con_GMD = resid(rr_bisquare)

# LSPL
rr_bisquare <- rlm(lspl ~ gmd_lspl, data=fmri_dat, psi = psi.bisquare)
mod0        <-  lm(lspl ~ gmd_lspl, data=fmri_dat)
summary(rr_bisquare)
summary(mod0)
plot(fmri_dat$gmd_lspl,fmri_dat$lspl)
abline(rr_bisquare$coefficients[1],rr_bisquare$coefficients[2])
abline(mod0$coefficients[1],mod0$coefficients[2],col="red")

fmri_dat$lspl_con_GMD = resid(rr_bisquare)

# RSMG
rr_bisquare <- rlm(rsmg ~ gmd_rsmg, data=fmri_dat, psi = psi.bisquare)
mod0        <-  lm(rsmg ~ gmd_rsmg, data=fmri_dat)
summary(rr_bisquare)
summary(mod0)
plot(fmri_dat$gmd_rsmg,fmri_dat$rsmg)
abline(rr_bisquare$coefficients[1],rr_bisquare$coefficients[2])
abline(mod0$coefficients[1],mod0$coefficients[2],col="red")

fmri_dat$rsmg_con_GMD = resid(rr_bisquare)


# 1st result (lmfg)
p <- ggplot(data = fmri_dat, aes(x = group, y = lmfg))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.3,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-40, 39, 23]"))
p <- p + geom_segment(aes(x=1, y=0.24, xend=3, yend=0.24), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.27), label = "HC<AD pFWE = 0.015*"), size =10)
p <- p +theme_la()
p

# 1st result controlled for GMD (lmfg)
p <- ggplot(data = fmri_dat, aes(x = group, y = lmfg_con_GMD))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.3,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-40, 39, 23]"))
p <- p + geom_segment(aes(x=1, y=0.24, xend=3, yend=0.24), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.27), label = "HC<AD pFWE = 0.035*"), size = 10)
p <- p +theme_la()
p

# 1st result LAc (lmfg)
p <- ggplot(data = fmri_dat_lac, aes(x = group, y = lmfg_lac))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.3,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-40, 39, 23]"))
p <- p + geom_segment(aes(x=1, y=0.24, xend=3, yend=0.24), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.27), label = "HC<AD pFWE > 0.05"), size = 10)
p <- p +theme_la()
p

# 2nd result (lspl)
p <- ggplot(data = fmri_dat, aes(x = group, y = lspl))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.4,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-36, -42, 54]"))
p <- p + geom_segment(aes(x=1, y=0.26, xend=3, yend=0.26), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<AD pFWE = 0.004*"), size = 10)
p <- p +theme_la()
p

# 2nd result controlled for GMD 
p <- ggplot(data = fmri_dat, aes(x = group, y = lspl_con_GMD))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.4,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-36, -42, 54]"))
p <- p + geom_segment(aes(x=1, y=0.26, xend=3, yend=0.26), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<AD pFWE = 0.017*"), size = 10)
p <- p +theme_la()
p

# 2nd result controlled for lac (lspl)
p <- ggplot(data = fmri_dat_lac, aes(x = group, y = lspl_lac))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.4,0.3))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [-36, -42, 54]"))
p <- p + geom_segment(aes(x=1, y=0.26, xend=3, yend=0.26), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<AD p > 0.05"), size = 10)
p <- p +theme_la()
p

# 3rd result (rsmg)
p <- ggplot(data = fmri_dat, aes(x = group, y = rsmg))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.2,0.4))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [55, -42, 37]"))
p <- p + geom_segment(aes(x=1, y=0.25, xend=2, yend=0.25), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<PG pFWE = 0.004*"), size = 10)
p <- p +theme_la()
p

# 3rd result controlled for GMD 
p <- ggplot(data = fmri_dat, aes(x = group, y = rsmg_con_GMD))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.2,0.4))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [55, -42, 37]"))
p <- p + geom_segment(aes(x=1, y=0.25, xend=2, yend=0.25), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<PG pFWE = 0.050*"), size = 10)
p <- p +theme_la()
p

# 3rd result controlled for lac (rsmg)
p <- ggplot(data = fmri_dat_lac, aes(x = group, y = rsmg_lac))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.9, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.2,0.4))
p <- p + xlab("Group") + ylab(expression("parameter esimate at [55, -42, 37]"))
p <- p + geom_segment(aes(x=1, y=0.25, xend=2, yend=0.25), size =1.2)
p <- p + geom_text (mapping = aes(x=c(2.0), y=c(0.285), label = "HC<PG p > 0.05"), size = 10)
p <- p +theme_la()
p




all_plots_stl[[1]] = p

multiplot(plotlist = all_plots_stl,cols=1)

# 3rd result (gPPI)
# fmri
path_eigv_fmri <- "F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\sig_results\\sig_res_3_noacc\\gPPI\\eigenvariates\\export"
setwd(path_eigv_fmri)
fmri_dat <- read.table("eigv_export.dat",sep="\t",header=T)
nms_anly <- read.table("names_analyses.dat",sep="\t",header=F)

# get ME and EOI out
throw_out <- (-1)*c(grep(pattern = "ME", as.character(nms_anly$V1)),grep(pattern = "EOI", as.character(nms_anly$V1)))
if (length(throw_out) != 0) {
  fmri_dat <- fmri_dat[,throw_out]
  nms_anly <- as.data.frame(nms_anly[throw_out,])
  names(nms_anly) <- c("V1") }

fmri_dat <- as.data.frame(fmri_dat[,-2])
names(fmri_dat) <- c("gPPI")

# get crm
setwd(pfad_results)
# save(file="cur_mod.Rdata",list = c("cur_mod"))
load("cur_mod.Rdata")
crm <- agk.get.compl.coef(cur_mod,"group") # no covs
fmri_dat$group =crm$group
fmri_dat$group = factor(fmri_dat$group,levels(fmri_dat$group)[c(2,3,1)])

p <- ggplot(data = fmri_dat, aes(x = group, y = gPPI))
p <- p+stat_summary(fun.y = mean, geom = "bar", width=0.2, fill="#333333")
p <- p + geom_abline(intercept = 0, slope = 0, color="black", size=1)
p <- p + stat_summary(fun.data = mean_cl_boot, geom = "linerange",
                      color="#FF3333", size =2.5)
p <-  p + coord_cartesian(ylim=c( -0.06,0.17))
p <- p + xlab("Group") + ylab(expression("gPPI: l. MFG to l. Putamen (a.u.)"))
p <- p + geom_line (mapping = aes(x=c(1,2), y=c(0.13, 0.13)), size =1.2)
p <- p + geom_text (mapping = aes(x=c(1.5), y=c(0.15), label = "HC>PG pFWE = 0.048*"), size = 9)
p <- p +theme_la()
p

## CORRELATIONS WITH NEURAL
# GSAS with nLA
# get the GMD
pfad_neur = 'F:\\AG_Diplomarbeit\\Preprocessed_on_Windows_Machine_swuaf_VBM8\\results_2nd_level\\multreg\\ep2d_bold_LA_00\\GSAS_Summe\\nLA_PG_GSAS\\3_noacc'
setwd(pfad_neur)
cur_name  = 'VOI_VOI.mat'
cur_extr  = readMat(paste(getwd(),cur_name,sep="/"))
cur_neur  = cur_extr$Y
cur_behav = data.la.aggr$GSAS_Summe[data.la.aggr$group=="PG"]
cur_dat   = data.frame(cur_behav,cur_neur)
p = ggplot(cur_dat, aes(x=cur_neur,y=cur_behav)) + geom_point(size =5) + theme_la() + xlab(label = paste0('neural LA at right ant. insula','\n'))+ylab(label = paste0('GSAS in PG group','\n'))
p = p + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) + ggtitle(paste('PG',"group\n"))
p
