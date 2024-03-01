
# results from escort
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
load("RDATA/allres_m3_dyntoy_M4.RData")
scoredf_m3 <- scoredf
load("RDATA/allres_dyntoy_M4.RData")
scoredf_ss <- scoredf
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")


##
all_scoredf <- rbind(scoredf_m3, scoredf_ss)
all_escort <- score_cal(all_scoredf)



# evaluation results
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORResEval_v4.R")
load(file = "RDATA/step0_clean_dyntoy_M4.RData")
load(file = "RDATA/step0_m3_dyntoy_M4.RData")
res_m3 <- rbind(round(DRMethodsEval(obj=UMAP1_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP2_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP3_eval, TruePT=scdata_time)$summary_df,3))
res_m3$Row.names <- rownames(scoredf_m3)
res_m3$TI <- rep("M3", 3)


load(file = "RDATA/step0_clean_dyntoy_M4.RData")
load(file = "RDATA/step0_ss_dyntoy_M4.RData")
res_ss <- rbind(round(DRMethodsEval(obj=MDS1_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=MDS2_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=MDS3_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE4_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE5_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE6_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP7_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP8_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP9_eval, TruePT=scdata_time)$summary_df,3))
res_ss$Row.names <- rownames(scoredf_ss)
res_ss$TI <- rep("SS", 9)

res_eval <- rbind(res_ss, res_m3)

# merge:
final_res_merge <- merge(all_escort, res_eval, by="Row.names")

par(mfrow = c(2, 2))
boxplot(final_res_merge$MeanMSE~final_res_merge$decision)
boxplot(final_res_merge$MeanCOR~final_res_merge$decision)


# for ss:
final_res_ss <- merge(score_cal(scoredf_ss), res_ss, by="Row.names")
boxplot(final_res_ss$MeanMSE~final_res_ss$decision)
boxplot(final_res_ss$MeanCOR~final_res_ss$decision)

# for m3:
final_res_m3 <- merge(score_cal(scoredf_m3), res_m3, by="Row.names")
boxplot(final_res_m3$MeanMSE~final_res_m3$decision)
boxplot(final_res_m3$MeanCOR~final_res_m3$decision)
