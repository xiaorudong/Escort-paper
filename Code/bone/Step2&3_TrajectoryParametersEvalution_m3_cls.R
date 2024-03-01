
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")
# dir.create("PLOTS/example_realData_bone/")


#################################
######### Load Data Set #########
#################################

load("RDATA/example_realData_bone/step0_clean_afterQC.RData")
load("RDATA/example_realData_bone/step0_m3_cls_afterQC.RData")
load("RDATA/example_realData_bone/step1_afterQC.RData")
load("RDATA/example_realData_bone/step0_m3_orig_afterQC.RData")



######################## STEP 2:########################

load("RDATA/example_realData_bone/res_cls_afterQC.RData")


######################## STEP 3:########################


#######################################################
######## Check the shape of trajectory curves #########
#######################################################

source("CODE/MethodScripts/Funtions/updated_ushape_cal.R")

pdf("PLOTS/example_realData_bone/ushape_m3_cls_evaluation_afterQC.pdf", onefile = TRUE, width = 8, height = 4.5)
par(mfrow=c(1, 2))
ushap_eval1 <- UshapeDetector(MDS1_eval, outlierdetect = "resistant")
ushap_eval2 <- UshapeDetector(MDS2_eval, outlierdetect = "resistant")
ushap_eval3 <- UshapeDetector(MDS3_eval, outlierdetect = "resistant")
ushap_eval4 <- UshapeDetector(TSNE4_eval, outlierdetect = "resistant")
ushap_eval5 <- UshapeDetector(TSNE5_eval, outlierdetect = "resistant")
ushap_eval6 <- UshapeDetector(TSNE6_eval, outlierdetect = "resistant")
ushap_eval7 <- UshapeDetector(UMAP7_eval, outlierdetect = "resistant")
ushap_eval8 <- UshapeDetector(UMAP8_eval, outlierdetect = "resistant")
ushap_eval9 <- UshapeDetector(UMAP9_eval, outlierdetect = "resistant")
ushap_eval10 <- UshapeDetector(PCA10_eval, outlierdetect = "resistant")
ushap_eval11 <- UshapeDetector(PCA11_eval, outlierdetect = "resistant")
ushap_eval12 <- UshapeDetector(PCA12_eval, outlierdetect = "resistant")
ushap_evalorig <- UshapeDetector(orig_eval, outlierdetect = "resistant")
dev.off()

goodmed <-  c(ushap_eval1$NoACells, ushap_eval2$NoACells, ushap_eval3$NoACells,
              ushap_eval4$NoACells, ushap_eval5$NoACells, ushap_eval6$NoACells,
              ushap_eval7$NoACells, ushap_eval8$NoACells, ushap_eval9$NoACells,
              ushap_eval10$NoACells, ushap_eval11$NoACells, ushap_eval12$NoACells,
              ushap_evalorig$NoACells) 
scoredf$USHAPE <- goodmed



###################################
######## Calculate Scores #########
###################################

final_df <- score_cal(scoredf)

save(scoredf, final_df, file = "RDATA/example_realData_bone/res_m3_cls_afterQC.RData")
save.image(file = "RDATA/example_realData_bone/allres_m3_cls_afterQC.RData")







