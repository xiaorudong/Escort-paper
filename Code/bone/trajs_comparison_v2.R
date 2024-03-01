
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORResEval_v4.R")


# plot PCA10_SS, MDS3_M3, original

# pdf("PLOTS/example_realData_bone/trajs_comparison_afterQC.pdf", onefile = TRUE, width = 12, height = 3)
par(mfrow=c(1, 4))
load(file="RDATA/example_realData_bone/step0_ss_cls_afterQC.RData")
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

obj <- PCA10_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[cut(pse$Ave, breaks=100)]
plot(dimred, col = plotcol, pch=16, main="TSNE_10%genes_SS",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')


load(file = "RDATA/example_realData_bone/step0_ss_orig_afterQC.RData")
obj <- ss_orig_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[cut(pse$Ave, breaks=100)]
plot(dimred, col = plotcol, pch=16, main="from paper",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')

load(file="RDATA/example_realData_bone/step0_m3_cls_afterQC.RData")
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

obj <- MDS3_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[cut(pse$Ave, breaks=100)]
plot(dimred, col = plotcol, pch=16, main="TSNE_50%genes_M3",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')

load(file = "RDATA/example_realData_bone/step0_m3_orig_afterQC.RData")
obj <- orig_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[cut(pse$Ave, breaks=100)]
plot(dimred, col = plotcol, pch=16, main="from paper",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')
# dev.off()



# make a similar figure 2
load("RDATA/res_bone_SSPC_afterQC.RData")
load(file="RDATA/example_realData_bone/step0_ss_cls_afterQC.RData")
update_data <- subdata

update_data@reductions$pca_new <- CreateDimReducObject(
  embeddings = as.matrix(PCA10_eval$Embedding),
  key = "PCA_",
  assay = "RNA"
)

update_data@active.ident <- factor(update_data@active.ident, levels = c("Cluster 1", "Cluster 2",
                                                                           "Cluster 3", "Cluster 4",
                                                                           "Cluster 5", "Cluster 6",
                                                                           "Cluster 8"))
names(update_data@active.ident) <- rownames(update_data@meta.data)

FeaturePlot(update_data, reduction="pca_new",
            features = c("Ihh", "Col10a1", "Mmp13", "Sp7", "Clec11a", "Col1a1", "Bglap"), 
            min.cutoff = "q10", max.cutoff = "q90", label = T,
            ncol = 2)

FeaturePlot(update_data, reduction="pca_new",
            features = c("Cre", "Sox9", "Acan", "Runx2","Spp1", "Alpl", "Ibsp"), 
            min.cutoff = "q10", max.cutoff = "q90", label = T,
            ncol = 2)

FeaturePlot(update_data, reduction="pca_new",
            features = c("Pdgfra", "Ly6a", "Lepr", "Cxcl12"), 
            min.cutoff = "q10", max.cutoff = "q90", label = T,
            ncol = 2)



pdf("PLOTS/example_realData_bone/trajs_cls_comparison_afterQC.pdf", onefile = TRUE, width = 20, height = 5)
par(mfrow=c(1, 4))
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(7)

obj <- PCA10_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[as.numeric(factor(update_data@active.ident))]
plot(dimred, col = plotcol, pch=16, main="PCA_10%genes_SS",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(x="topleft",legend = levels(update_data@active.ident),
       col = colors, pch=16)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')


##

load(file="RDATA/example_realData_bone/step0_m3_cls_afterQC.RData")
obj <- PCA10_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[as.numeric(factor(update_data@active.ident))]
plot(dimred, col = plotcol, pch=16, main="PCA_10%genes_M3",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(x="topleft",legend = levels(update_data@active.ident),
       col = colors, pch=16)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')


# load("RDATA/res_bone_SSPC_afterQC.RData")
# load("RDATA/example_realData_bone/step0_ss_orig_afterQC.RData")

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(7)

obj <- ss_orig_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[as.numeric(factor(update_data@active.ident))]
plot(dimred, col = plotcol, pch=16, main="from paper_SS",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(x="topright",legend = levels(update_data@active.ident),
       col = colors, pch=16)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')



##

# load("RDATA/res_bone_SSPC_afterQC.RData")
# load("RDATA/example_realData_bone/step0_m3_orig_afterQC.RData")

library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(7)

obj <- orig_eval
pse <- obj$pse
dimred <- obj$Embedding
fitLine <- obj$fitLine
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[as.numeric(factor(update_data@active.ident))]
plot(dimred, col = plotcol, pch=16, main="from paper",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
legend(x="topright",legend = levels(update_data@active.ident),
       col = colors, pch=16)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')


# dev.off()

library(Seurat)
pdf("PLOTS/example_realData_bone/markers_cls_comparison_afterQC.pdf", onefile = TRUE, width = 12, height = 12)
VlnPlot(object = update_data, features = c("Ihh", "Col10a1", "Mmp13", "Sp7", "Clec11a", "Col1a1", "Bglap"), cols = colors)
VlnPlot(object = update_data, features = c("Cre", "Sox9", "Acan", "Runx2","Spp1", "Alpl", "Ibsp"), cols = colors)
VlnPlot(object = update_data, features = c("Pdgfra", "Ly6a", "Lepr", "Cxcl12"), cols = colors)
dev.off()




# PT vs Expression

load("RDATA/example_realData_bone/step0_clean_afterQC.RData")
library(gridExtra)
library(ggplot2)
# select genes
topg <- "Ly6a"

obj <- orig_eval
pse <- obj$pse

gene <- topg[1]
Gene_PT <- data.frame(Time=obj$pse, Expresssion=norm_counts[which(rownames(norm_counts)==gene),], Group=rep("PT", length(obj$pse[,1])))
p1 <- ggplot(data=Gene_PT, aes(x=pse, y=Expresssion)) +
  geom_point(alpha=0.3, shape=16) +
  geom_smooth(method = "loess", se = F) +
  ggtitle(paste("obj_", gene, sep = "")) +
  theme_bw()


