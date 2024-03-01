
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

load("RDATA/res_bone_SSPC_afterQC.RData")
load("RDATA/example_realData_bone/orig_cls_afterQC.RData")


#######################################
######### Trajectory Analysis #########
#######################################

# Take using Slingshot to do trajectory analysis as an example
library(slingshot)
library(DescTools)
library(vioplot)
library(scales)
library(mclust)
set.seed(123)
# kmeans_nstart= ifelse(ncol(norm_counts)>2000, 50, 1000)

# for subset 1
dimred <- subdata@reductions$umap@cell.embeddings
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_outorig <- ti_out
rawpseorig <- rawpse
pseorig <- pse
fitLineorig <- as.data.frame(fitLine)


source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")

SS_objorig <- info_traj(dimred, rawpse = rawpseorig, fitLine = fitLineorig)
ss_orig_eval <- TEvalobj(dimred, traj = SS_objorig)

# dir.create("RDATA/example_realData_bone")
save(ss_orig_eval, file = "RDATA/example_realData_bone/step0_ss_orig_afterQC.RData")


