
rm(list = ls())
setwd("/Users/ruby/Desktop/Graduate/Research/Escort/")
devtools::load_all()

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

load(file = "RDATA/example_realData_bone/step0_clean_afterQC.RData")

load(file = "RDATA/example_realData_bone/res_cls_afterQC.RData")
scoredf_ss <- scoredf
load(file = "RDATA/example_realData_bone/res_m3_cls_afterQC.RData")
scoredf_m3 <- scoredf



scoredf_ss$USHAPE <- scoredf_ss$USHAPE/ncol(norm_counts)
escort_ss <- calcScore(scoredf_ss)

scoredf_m3$USHAPE <- scoredf_m3$USHAPE/ncol(norm_counts)
escort_m3 <- calcScore(scoredf_m3)

rownames(scoredf_ss) <- paste("SS", rownames(scoredf_ss), sep = "_")
rownames(scoredf_m3) <- paste("M3", rownames(scoredf_m3), sep = "_")

scoredf_all <-rbind(scoredf_ss, scoredf_m3)

escort_all <- calcScore(scoredf_all)
escort_all <- escort_all[order(escort_all$score, decreasing = T),]

write.csv(escort_all, file="OUT/supp_bone_example_escort_results.csv", row.names = F)

c(2006, 4011, 10028)

library(stringr)
library(reshape2)
small_escort <- escort_all[,c(1, 8, 9, 10)]
split_col <- as.data.frame(str_split_fixed(small_escort$Row.names, "_", 3))
split_dr <- as.data.frame(colsplit(split_col$V2, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit")))
split_dr$HVGs <- ifelse(split_dr$digit%%3==0, "10028 (50%Genes)",
                        ifelse(split_dr$digit%%3==1, "2006 (10%Genes)", "4011 (20%Genes)"))
split_dr$HVGs[is.na(split_dr$HVGs)] <- "2000"
small_escort <- cbind(Trajectory=split_col$V1, split_dr[,c(1,3)], small_escort)
small_escort$char[small_escort$char=="orig"] <- "PCA+UMAP"

write.csv(small_escort, file="OUT/supp_bone_example_escort_results_v2.csv", row.names = F)

