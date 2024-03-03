
library(Escort)

setwd("path/to/your_file")
load("RDATA/Manuscript/clean_dyntoy_DC3.RData")
load("RDATA/Manuscript/dr_dyntoy_DC3.RData")


#######################################
######### Trajectory Analysis #########
#######################################

# Take using Slingshot to do trajectory analysis as an example
library(slingshot)
library(mclust)
set.seed(123)
startcell <- names(which.min(scdata_time))

ls_fitLines <- list()
rawpses <- list()
for (i in 1:9) {
  dimred <- embeddings[[i]]
  cl1 <- Mclust(dimred)$classification
  ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus=as.character(cl1[startcell]))
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


save(fitLines, resobjs, file = "RDATA/Manuscript/slingshot_dyntoyDC3.RData")







