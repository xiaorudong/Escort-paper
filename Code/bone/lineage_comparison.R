
rm(list = ls())


update_data <- subdata
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORResEval_v4.R")

# plot PCA10_SS, orig_ss

par(mfrow=c(1, 1))
load(file="RDATA/example_realData_bone/step0_ss_cls_afterQC.RData")
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)


pdf("PLOTS/example_realData_bone_forPaperFig6_B.pdf", onefile = TRUE, width = 4.5, height = 2.8)
par(mfrow=c(1, 2), mar=c(5,4,2,1))


pca_obj <- PCA10_eval
pca_pse <- pca_obj$pse
pca_dimred <- pca_obj$Embedding
pca_fitLine <- pca_obj$fitLine
pca_pse$Ave <- rowMeans(pca_pse, na.rm = T)
plotcol <- colors[cut(pca_pse$Ave, breaks=100)]
plot(pca_dimred, col = plotcol, pch=16, cex=.8,
     cex.lab=1, cex.axis=1)
segments(as.data.frame(pca_fitLine)$x0, as.data.frame(pca_fitLine)$y0, 
         as.data.frame(pca_fitLine)$x1, as.data.frame(pca_fitLine)$y1, 
         lwd=2, col='black')


load(file = "RDATA/example_realData_bone/step0_ss_orig_afterQC.RData")
orig_obj <- ss_orig_eval
orig_pse <- orig_obj$pse
orig_dimred <- orig_obj$Embedding
orig_fitLine <- orig_obj$fitLine
orig_pse$Ave <- rowMeans(orig_pse, na.rm = T)
plot(orig_dimred, col = plotcol, pch=16, cex=.8,
     cex.lab=1, cex.axis=1)
segments(as.data.frame(orig_fitLine)$x0, as.data.frame(orig_fitLine)$y0, 
         as.data.frame(orig_fitLine)$x1, as.data.frame(orig_fitLine)$y1, 
         lwd=2, col='black')

dev.off()

load("RDATA/res_bone_SSPC_afterQC.RData")
dim(subdata[[]])
dim(pca_pse)

pca_pse <- pca_pse[rownames(subdata[[]]),]

subdata$pt_lineage1 <- pca_pse$Lineage1
subdata$pt_lineage2 <- pca_pse$Lineage2

subdata$old_cls <- factor(subdata$old_cls, levels = paste0("Cluster ", c(1:6,8)))
par(mfrow=c(1,2))
boxplot(subdata$pt_lineage2 ~ subdata$old_cls)
boxplot(subdata$pt_lineage1 ~ subdata$old_cls)

head(subdata[[]])
subdata$subclust4 <- "Other"
subdata$subclust4[which(subdata$old_cls == "Cluster 4" & is.na(subdata$pt_lineage1) & !is.na(subdata$pt_lineage2))] <- "OnlyLineage2"
subdata$subclust4[which(subdata$old_cls == "Cluster 4" & is.na(subdata$pt_lineage2) & !is.na(subdata$pt_lineage1))] <- "OnlyLineage1"
subdata$subclust4[which(subdata$old_cls == "Cluster 4" & !is.na(subdata$pt_lineage2) & !is.na(subdata$pt_lineage1))] <- "Both"

subdata$subclust4 <- "Other"
subdata$subclust4[which(subdata$old_cls == "Cluster 4" & subdata$pt_lineage2 > .4)] <- "OnlyLineage2"
subdata$subclust4[which(subdata$old_cls == "Cluster 4" & subdata$pt_lineage1 < .5)] <- "OnlyLineage1"

table(subdata$old_cls)
table(subdata$subclust4)

library(Seurat)
Idents(subdata) <- "subclust4"
subdata2 <- subset(subdata, subclust4!= "Other")
subdata2

myde <- FindMarkers(subdata2, ident.1 = "OnlyLineage2", ident.2 = "OnlyLineage1", verbose = FALSE)

head(myde, 11)

siglin2g <- (subset(myde, avg_log2FC > 1))
dim(siglin2g); head(siglin2g,10)
VlnPlot(subdata2, features = "Taf7")

# Go to the line that has: ########### Run everything below to make Ruby's plots work in the plotg function
## !!!
plotg("Taf7")

write.table(rownames(siglin2g), file="~/DE_clust4_lineage2.txt", quote = F, row.names = F, col.names = F)

lin1g <- (subset(myde, avg_log2FC < -1))
dim(lin1g); head(lin1g)
plotg("Mmp13")
plotg("Alpl")
plotg("Spp1")
Idents(subdata) <- "old_cls"
VlnPlot(subdata, features = "Mmp13")

VlnPlot(subdata2, features = "Mmp13")

write.table(rownames(lin1g), file="~/DE_clust4_lineage1.txt", quote = F, row.names = F, col.names = F)
### Gives a good idea but not really enough cells for robust analysis here.


###### try something simple:

Idents(subdata) <- "old_cls"
myde2 <- FindMarkers(subdata, ident.1 = "Cluster 3", ident.2 = "Cluster 4", verbose = FALSE)

head(myde2, 11)

siglin2g <- (subset(myde2, avg_log2FC > 1))
dim(siglin2g)
head(siglin2g)
plotg("Col11a1")
VlnPlot(subdata, features = "Col11a1")

write.table(rownames(siglin2g), file="~/DE_clust4_clust3.txt", quote = F, row.names = F, col.names = F)


lin1g <- (subset(myde2, avg_log2FC < -1))
head(lin1g)
dim(lin1g)
plotg("Hspa8")
plotg("Col10a1")
plotg("Adipoq")

Idents(subdata) <- "old_cls"

VlnPlot(subdata, features = "Tnfsf11")
VlnPlot(subdata, features = "Sox9")
VlnPlot(subdata, features = "Hspa8")

VlnPlot(subdata2, features = "Adipoq")
VlnPlot(subdata2, features = "Ibsp")

plotg("Apoc1")

write.table(rownames(lin1g), file="~/DE_clust4_clust4.txt", quote = F, row.names = F, col.names = F)

  
  
########### Run everything below to make Ruby's plots work in the plotg function 





plotg<-function(gene) {
  Gene_PT <- data.frame(pca_obj$pse, Expresssion=norm_counts[which(rownames(norm_counts)==gene),])
  Gene_PT_long_pca <- Gene_PT %>% pivot_longer(!Expresssion, names_to = "Lineage", values_to = "PT")
  Gene_PT_long_pca$From <- rep("PCA_SS", nrow(Gene_PT_long_pca))
  
  Gene_PT <- data.frame(orig_obj$pse[,1:2], Expresssion=norm_counts[which(rownames(norm_counts)==gene),])
  Gene_PT_long_orig <- Gene_PT %>% pivot_longer(!Expresssion, names_to = "Lineage", values_to = "PT")
  Gene_PT_long_orig$From <- rep("Orig_SS", nrow(Gene_PT_long_orig))
  
  Gene_PT_long <- as.data.frame(rbind(Gene_PT_long_pca, Gene_PT_long_orig))
  
  
print(ggplot(data=Gene_PT_long, aes(x=PT, y=Expresssion, color = Lineage)) +
    geom_point(alpha=0.3, shape=16) + 
    scale_color_manual(values=c("#000000", "#bc3a3a")) +
    geom_smooth(method = "gam", se = F) +
    facet_grid(cols = vars(From))+
    ggtitle(paste(gene, sep = "--")) +
    theme_bw())

}

head(pca_pse)
pca_pse <- pca_pse[,1:2]
save(norm_counts,pca_pse, file="RDATA/example_realData_bone/dataToRun_scLane.Rdata")


