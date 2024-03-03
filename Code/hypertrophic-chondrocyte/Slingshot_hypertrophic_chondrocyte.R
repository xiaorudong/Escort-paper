
library(Escort)

setwd("path/to/your_file")
load("RDATA/Manuscript/dr_bone_SSPC.RData")
load("RDATA/Manuscript/orig_bone_SSPC_cls_afterQC.RData")


#######################################
######### Trajectory Analysis #########
#######################################

# Take using Slingshot to do trajectory analysis as an example
library(slingshot)
set.seed(123)

ls_fitLines <- list()
rawpses <- list()
for (i in 1:13) {
  dimred <- embeddings[[i]]
  cl1 <- list.cluster
  ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
  rawpse <- slingPseudotime(ti_out, na=T)
  ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
  
  ls_fitLines[[i]] <- ls_fitLine
  rawpses[[i]] <- rawpse
}


library(parallel)
fitLines <- mclapply(ls_fitLines, segFormat)
resobjs <- mclapply(1:9, function(x) {
  prepTraj(embeddings[[x]], PT=rawpses[[x]], fitLine=fitLines[[x]])
})


save(fitLines, resobjs, file = "RDATA/Manuscript/slingshot_bone_SSPC.RData")







