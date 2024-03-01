
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")
# dir.create("PLOTS/Script_RealCollection1")


#################################
######### Load Data Set #########
#################################

load("RDATA/step0_clean_RealCollection1.RData")
load("RDATA/step0_m3_RealCollection1.RData")
load("RDATA/step1_RealCollection1.RData")



######################## STEP 2:########################

load("RDATA/res_RealCollection1.RData")


######################## STEP 3:########################


#######################################################
######## Check the shape of trajectory curves #########
#######################################################

source("CODE/MethodScripts/Funtions/updated_ushape_cal.R")

pdf("PLOTS/Script_RealCollection1/ushape_m3_evaluation.pdf", onefile = TRUE, width = 8, height = 4.5)
par(mfrow=c(1, 2))
ushap_eval1 <- UshapeDetector(MDS1_eval)
ushap_eval2 <- UshapeDetector(MDS2_eval)
ushap_eval3 <- UshapeDetector(MDS3_eval)
ushap_eval4 <- UshapeDetector(TSNE4_eval, outlierdetect = "resistant")
ushap_eval5 <- UshapeDetector(TSNE5_eval, outlierdetect = "resistant")
ushap_eval6 <- UshapeDetector(TSNE6_eval, outlierdetect = "resistant")
ushap_eval7 <- UshapeDetector(UMAP7_eval, outlierdetect = "resistant")
ushap_eval8 <- UshapeDetector(UMAP8_eval, outlierdetect = "resistant")
ushap_eval9 <- UshapeDetector(UMAP9_eval, outlierdetect = "resistant")
dev.off()

goodmed <-  c(ushap_eval1$NoACells, ushap_eval2$NoACells, ushap_eval3$NoACells,
              ushap_eval4$NoACells, ushap_eval5$NoACells, ushap_eval6$NoACells,
              ushap_eval7$NoACells, ushap_eval8$NoACells, ushap_eval9$NoACells) 
scoredf$USHAPE <- goodmed



###################################
######## Calculate Scores #########
###################################

final_df <- score_cal(scoredf)

save(scoredf, final_df, file = "RDATA/res_m3_RealCollection1.RData")
save.image(file = "RDATA/allres_m3_RealCollection1.RData")







