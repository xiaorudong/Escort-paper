
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")
# dir.create("PLOTS/Script_dyntoy_DC3")


#################################
######### Load Data Set #########
#################################

load("RDATA/step0_clean_dyntoy_DC3.RData")
load("RDATA/step0_ss_dyntoy_DC3.RData")
load("RDATA/step1_dyntoy_DC3.RData")



######################## STEP 2:########################
######################## STEP 3:########################


scoredf <- data.frame(DCcheck=rep(NA, 9), SimiRetain=rep(NA, 9),
                      GOF=rep(NA, 9),USHAPE=rep(NA, 9))
rownames(scoredf) <- c("MDS1_eval", "MDS2_eval", "MDS3_eval",
                       "TSNE4_eval", "TSNE5_eval", "TSNE6_eval",
                       "UMAP7_eval", "UMAP8_eval", "UMAP9_eval")  

scoredf


###################################
######## Calculate Scores #########
###################################

final_df <- score_cal(scoredf)

save(scoredf, final_df, file = "RDATA/res_dyntoy_DC3.RData")
save.image(file = "RDATA/allres_dyntoy_DC3.RData")










