
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
escort_ss <- score_cal(scoredf_ss)

scoredf_m3$USHAPE <- scoredf_m3$USHAPE/ncol(norm_counts)
escort_m3 <- score_cal(scoredf_m3)

rownames(scoredf_ss) <- paste("SS", rownames(scoredf_ss), sep = "_")
rownames(scoredf_m3) <- paste("M3", rownames(scoredf_m3), sep = "_")

scoredf_all <-rbind(scoredf_ss, scoredf_m3)

escort_all <- score_cal(scoredf_all, SimiRetaincutoff = 0.5)
escort_all[order(escort_all$score, decreasing = T),]


# PT correlations:
PT_list <- list(MDS1_SS_PT, MDS2_SS_PT, MDS3_SS_PT,
                TSNE4_SS_PT, TSNE5_SS_PT, TSNE6_SS_PT,
                UMAP7_SS_PT, UMAP8_SS_PT, UMAP9_SS_PT,
                PCA10_SS_PT, PCA11_SS_PT, PCA12_SS_PT,
                MDS1_M3_PT, MDS2_M3_PT, MDS3_M3_PT,
                TSNE4_M3_PT, TSNE5_M3_PT, TSNE6_M3_PT,
                UMAP7_M3_PT, UMAP8_M3_PT, UMAP9_M3_PT,
                PCA10_M3_PT, PCA11_M3_PT, PCA12_M3_PT,
                orig_SS_PT,orig_M3_PT)

cor.mod <- cor(do.call(cbind, PT_list), method = "spearman")
colnames(cor.mod) <- c("MDS1_SS", "MDS2_SS", "MDS3_SS", 
                       "TSNE4_SS", "TSNE5_SS", "TSNE6_SS",
                       "UMAP7_SS", "UMAP8_SS", "UMAP9_SS",
                       "PCA10_SS", "PCA11_SS", "PCA12_SS",
                       "MDS1_M3", "MDS2_M3", "MDS3_M3", 
                       "TSNE4_M3", "TSNE5_M3", "TSNE6_M3",
                       "UMAP7_M3", "UMAP8_M3", "UMAP9_M3",
                       "PCA10_M3", "PCA11_M3", "PCA12_M3", 
                       "orig_SS", "orig_M3")
rownames(cor.mod) <- c("MDS1_SS", "MDS2_SS", "MDS3_SS", 
                       "TSNE4_SS", "TSNE5_SS", "TSNE6_SS",
                       "UMAP7_SS", "UMAP8_SS", "UMAP9_SS",
                       "PCA10_SS", "PCA11_SS", "PCA12_SS",
                       "MDS1_M3", "MDS2_M3", "MDS3_M3", 
                       "TSNE4_M3", "TSNE5_M3", "TSNE6_M3",
                       "UMAP7_M3", "UMAP8_M3", "UMAP9_M3",
                       "PCA10_M3", "PCA11_M3", "PCA12_M3", 
                       "orig_SS", "orig_M3")



pdf("PLOTS/example_realData_bone/method_evaluation_corr_cls_afterQC.pdf", onefile = TRUE, width = 16, height = 16)
par(mfrow=c(1,1))
library(corrplot)
corrplot(cor.mod, method = 'color', diag = F, type = 'lower',
         tl.cex=0.8, tl.col = "black", tl.srt = 0, col = COL2('RdYlBu'),
         addgrid.col = 'white', addCoef.col = 'black',number.cex = 0.8)
dev.off()


