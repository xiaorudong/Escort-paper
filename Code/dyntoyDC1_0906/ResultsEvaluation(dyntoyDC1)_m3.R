
rm(list = ls())

##############################################
######### Load All Related Data Sets #########
##############################################

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORResEval_v4.R")
load(file = "RDATA/step0_clean_dyntoy_DC1.RData")
load(file = "RDATA/step0_m3_dyntoy_DC1.RData")



###########################################
######### Estimated PT vs True PT #########
###########################################

pdf("PLOTS/Script_dyntoy_DC1/method_m3_evaluation.pdf", onefile = TRUE, width = 8, height = 8)
DRMethodsEval(obj=MDS1_eval, TruePT=scdata_time)
DRMethodsEval(obj=MDS2_eval, TruePT=scdata_time)
DRMethodsEval(obj=MDS3_eval, TruePT=scdata_time)
DRMethodsEval(obj=TSNE4_eval, TruePT=scdata_time)
DRMethodsEval(obj=TSNE5_eval, TruePT=scdata_time)
DRMethodsEval(obj=TSNE6_eval, TruePT=scdata_time)
DRMethodsEval(obj=UMAP7_eval, TruePT=scdata_time)
DRMethodsEval(obj=UMAP8_eval, TruePT=scdata_time)
DRMethodsEval(obj=UMAP9_eval, TruePT=scdata_time)
dev.off()

res <- rbind(round(DRMethodsEval(obj=MDS1_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=MDS2_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=MDS3_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE4_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE5_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=TSNE6_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP7_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP8_eval, TruePT=scdata_time)$summary_df,3),
             round(DRMethodsEval(obj=UMAP9_eval, TruePT=scdata_time)$summary_df,3))


MDS1_PT <- rowMeans(MDS1_eval$pse, na.rm = T)
MDS2_PT <- rowMeans(MDS2_eval$pse, na.rm = T)
MDS3_PT <- rowMeans(MDS3_eval$pse, na.rm = T)
TSNE4_PT <- rowMeans(TSNE4_eval$pse, na.rm = T)
TSNE5_PT <- rowMeans(TSNE5_eval$pse, na.rm = T)
TSNE6_PT <- rowMeans(TSNE6_eval$pse, na.rm = T)
UMAP7_PT <- rowMeans(UMAP7_eval$pse, na.rm = T)
UMAP8_PT <- rowMeans(UMAP8_eval$pse, na.rm = T)
UMAP9_PT <- rowMeans(UMAP9_eval$pse, na.rm = T)



# PT distributions
pdf("PLOTS/Script_dyntoy_DC1/method_m3_evaluation2.pdf", onefile = TRUE, width = 12, height = 12)
library(gridExtra)
# select genes
library(scuttle)
library(scran)
gene.var <- modelGeneVar(norm_counts)
topg <- getTopHVGs(gene.var, n=20)

for (gene in topg) {
  # MDS1_eval
  Gene_PT <- data.frame(Time=MDS1_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(MDS1_eval$pse[,1])))
  Gene_TT <- data.frame(Time=scdata_time/max(scdata_time), Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("TT", length(scdata_time)))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p1 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("MDS2000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # MDS2_eval
  Gene_PT <- data.frame(Time=MDS2_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(MDS2_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p2 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("MDS5000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # MDS3_eval
  Gene_PT <- data.frame(Time=MDS3_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(MDS3_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p3 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("MDS11000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # TSNE4_eval
  Gene_PT <- data.frame(Time=TSNE4_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(TSNE4_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p4 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("TSNE2000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # TSNE5_eval
  Gene_PT <- data.frame(Time=TSNE5_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(TSNE5_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p5 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("TSNE5000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # TSNE6_eval
  Gene_PT <- data.frame(Time=TSNE6_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(TSNE6_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p6 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("TSNE11000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # UMAP7_eval
  Gene_PT <- data.frame(Time=UMAP7_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(UMAP7_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p7 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("UMAP2000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # UMAP8_eval
  Gene_PT <- data.frame(Time=UMAP8_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(UMAP8_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p8 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("UMAP5000_eval_", gene, sep = "")) +
    theme_bw()
  
  
  # UMAP9_eval
  Gene_PT <- data.frame(Time=UMAP9_PT, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(UMAP9_eval$pse[,1])))
  Gene_df <- rbind(Gene_PT, Gene_TT)
  
  library(ggplot2)
  p9 <- ggplot(data=Gene_df, aes(x=Time, y=Expresssion, color = Group)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "loess", se = F) +
    ggtitle(paste("UMAP11000_eval_", gene, sep = "")) +
    theme_bw()
  
  grid.arrange(p1, p2, p3, 
               p4, p5, p6, 
               p7, p8, p9, nrow = 3)
  
}

dev.off()


# PT correlations:
PT_list <- list(MDS1_PT, MDS2_PT, MDS3_PT,
                TSNE4_PT, TSNE5_PT, TSNE6_PT,
                UMAP7_PT, UMAP8_PT, UMAP9_PT)

cor.mod <- cor(do.call(cbind, PT_list), method = "spearman")
colnames(cor.mod) <- c("MDS1", "MDS2", "MDS3", 
                       "TSNE4", "TSNE5", "TSNE6",
                       "UMAP7", "UMAP8", "UMAP9")
rownames(cor.mod) <- c("MDS1", "MDS2", "MDS3", 
                       "TSNE4", "TSNE5", "TSNE6",
                       "UMAP7", "UMAP8", "UMAP9")



pdf("PLOTS/Script_dyntoy_DC1/method_m3_evaluation_corr.pdf", onefile = TRUE, width = 5, height = 5)
par(mfrow=c(1,1))
library(corrplot)
corrplot(cor.mod, method = 'color', diag = F, type = 'lower',
         tl.cex=0.8, tl.col = "black", tl.srt = 0, col = COL2('RdYlBu'),
         addgrid.col = 'white', addCoef.col = 'black',number.cex = 0.8)
dev.off()





