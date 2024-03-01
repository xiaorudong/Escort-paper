
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")
# dir.create("PLOTS/Script_sim_DEG0_500")


#################################
######### Load Data Set #########
#################################

load("RDATA/step0_clean_sim_DEG0_500.RData")
load("RDATA/step0_m3_sim_DEG0_500.RData")
load("RDATA/step1_sim_DEG0_500.RData")



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

save(scoredf, final_df, file = "RDATA/res_m3_sim_DEG0_500.RData")
save.image(file = "RDATA/allres_m3_sim_DEG0_500.RData")







