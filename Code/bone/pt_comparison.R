
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORResEval_v4.R")


# plot PCA10_SS, orig_ss

par(mfrow=c(1, 2))
load(file="RDATA/example_realData_bone/step0_ss_cls_afterQC.RData")
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

pca_obj <- PCA10_eval
pca_pse <- pca_obj$pse
pca_dimred <- pca_obj$Embedding
pca_fitLine <- pca_obj$fitLine
pca_pse$Ave <- rowMeans(pca_pse, na.rm = T)
plotcol <- colors[cut(pca_pse$Ave, breaks=100)]
plot(pca_dimred, col = plotcol, pch=16, main="PCA_10%genes_SS",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(pca_fitLine)$x0, as.data.frame(pca_fitLine)$y0, 
         as.data.frame(pca_fitLine)$x1, as.data.frame(pca_fitLine)$y1, 
         lwd=2, col='black')


load(file = "RDATA/example_realData_bone/step0_ss_orig_afterQC.RData")
orig_obj <- ss_orig_eval
orig_pse <- orig_obj$pse
orig_dimred <- orig_obj$Embedding
orig_fitLine <- orig_obj$fitLine
orig_pse$Ave <- rowMeans(orig_pse, na.rm = T)
plotcol <- colors[cut(orig_pse$Ave, breaks=100)]
plot(orig_dimred, col = plotcol, pch=16, main="from paper_ss",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(orig_fitLine)$x0, as.data.frame(orig_fitLine)$y0, 
         as.data.frame(orig_fitLine)$x1, as.data.frame(orig_fitLine)$y1, 
         lwd=2, col='black')


# check each lineage:
# pca 
par(mfrow=c(1, 2))
plotcol <- colors[cut(pca_pse$Ave, breaks=100)]
plot(pca_dimred, col = "grey", pch=16, main="PCA_Lineage1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca_dimred[!is.na(pca_pse$Lineage1),1], pca_dimred[!is.na(pca_pse$Lineage1),2],
       col = plotcol[!is.na(pca_pse$Lineage1)], pch=16)

plot(pca_dimred, col = "grey", pch=16, main="PCA_Lineage2",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(pca_dimred[!is.na(pca_pse$Lineage2),1], pca_dimred[!is.na(pca_pse$Lineage2),2],
       col = plotcol[!is.na(pca_pse$Lineage2)], pch=16)

# orig
par(mfrow=c(1, 3))
plotcol <- colors[cut(orig_pse$Ave, breaks=100)]
plot(orig_dimred, col = "grey", pch=16, main="orig_Lineage1",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(orig_dimred[!is.na(orig_pse$Lineage1),1], orig_dimred[!is.na(orig_pse$Lineage1),2],
       col = plotcol[!is.na(orig_pse$Lineage1)], pch=16)

plot(orig_dimred, col = "grey", pch=16, main="orig_Lineage2",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(orig_dimred[!is.na(orig_pse$Lineage2),1], orig_dimred[!is.na(orig_pse$Lineage2),2],
       col = plotcol[!is.na(orig_pse$Lineage2)], pch=16)

plot(orig_dimred, col = "grey", pch=16, main="orig_Lineage3",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(orig_dimred[!is.na(orig_pse$Lineage3),1], orig_dimred[!is.na(orig_pse$Lineage3),2],
       col = plotcol[!is.na(orig_pse$Lineage3)], pch=16)



load(file = "RDATA/example_realData_bone/step0_clean_afterQC.RData")
library(gridExtra)
library(ggplot2)
library(tidyr)
# paper markers
paperg <- c("Ihh", "Col10a1", "Mmp13", "Sp7", "Clec11a", "Col1a1", "Bglap", 
            "Cre", "Sox9", "Acan", "Runx2","Spp1", "Alpl", "Ibsp", "Pdgfra", 
            "Ly6a", "Lepr", "Cxcl12")

# more genes from HVGs:
library(scuttle)
library(scran)
gene.var <- modelGeneVar(norm_counts)
genes.HVGs <- getTopHVGs(gene.var, n = 20)

topg <- c(paperg, genes.HVGs[! genes.HVGs %in% paperg])

pdf("PLOTS/example_realData_bone/pt_comparison_afterQC.pdf", onefile = TRUE, width = 7, height = 3)
for (gene in topg) {
  Gene_PT <- data.frame(pca_obj$pse, Expresssion=norm_counts[which(rownames(norm_counts)==gene),])
  Gene_PT_long_pca <- Gene_PT %>% pivot_longer(!Expresssion, names_to = "Lineage", values_to = "PT")
  Gene_PT_long_pca$From <- rep("PCA_SS", nrow(Gene_PT_long_pca))
  
  Gene_PT <- data.frame(orig_obj$pse[,1:2], Expresssion=norm_counts[which(rownames(norm_counts)==gene),])
  Gene_PT_long_orig <- Gene_PT %>% pivot_longer(!Expresssion, names_to = "Lineage", values_to = "PT")
  Gene_PT_long_orig$From <- rep("Orig_SS", nrow(Gene_PT_long_orig))
  
  Gene_PT_long <- as.data.frame(rbind(Gene_PT_long_pca, Gene_PT_long_orig))
  
  gene_intro <- ifelse(gene %in% paperg & gene %in% genes.HVGs, "Top 20 HVG; Paper Marker",
                       ifelse(gene %in% paperg, "Paper Marker", "Top 20 HVG"))
  p <- ggplot(data=Gene_PT_long, aes(x=PT, y=Expresssion, color = Lineage)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "gam", se = F) +
    facet_grid(cols = vars(From))+
    ggtitle(paste(gene, gene_intro, sep = "--")) +
    theme_bw()
  print(p)
}
dev.off()



