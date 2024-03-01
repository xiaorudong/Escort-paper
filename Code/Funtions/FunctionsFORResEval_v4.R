
library(scran)
### DR methods evaluation by comparing estimted PT with true PT: 
# obj: a trajectory evaluation object from TEvalobj() function
# dimred: A data frame. a 2D embedding from a DR method
# pse: A data frame. estimated pseudotime from TI. Each column contains the PT for one lineage. 
# fitLine: the fitted curve line got from TI.
# TruePT: true PT
# output is a list containing: 
#           + AllMSE: MSE from all lineages
#           + AllCOR: COR from all lineages
#           + MeanMSE: Mean MSE from all lineages
#           + MeanCOR: Mean COR from all lineages

DRMethodsEval <- function(obj=NULL, dimred=NULL, pse=NULL, fitLine=NULL, PLOTname=NULL, TruePT, NeedTimeFlip=F) {
  
  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  if (!require("scales")) install.packages("scales")
  
  library(RColorBrewer)
  library(scales)
  
  if (!is.null(obj)) {
    dimred <- obj$Embedding
    pse <- obj$pse
    fitLine <- obj$fitLine
    PLOTname <- deparse(substitute(obj))
  }
  pse <- as.data.frame(pse)
  
  par(mfrow=c(2, 2))
  
  scdata_time <- TruePT
  
  a <- match(rownames(pse), names(scdata_time))
  scdata_time <- scdata_time[a]
  scdata_time <-  scdata_time/max(scdata_time, na.rm=TRUE)
  
  if (NeedTimeFlip) {
    pse <- 1-pse
  }
  
  nlineage <- ncol(pse)
  mse_vec <- c()
  cor_vec <- c()
  mae_vec <- c()
  for (i in 1:nlineage) {
    lineage_pse <- pse[!is.na(pse[,i]),i]
    names(lineage_pse) <- rownames(pse)[!is.na(pse[,i])]
    
    lineage_dimred <- dimred[rownames(dimred) %in% names(lineage_pse),]
    lineage_pse <- lineage_pse[match(rownames(lineage_dimred), names(lineage_pse))]
    lineage_time <- scdata_time[names(scdata_time) %in% names(lineage_pse)]
    
    # compare PTs
    mse_diff <- mse(lineage_pse, lineage_time)
    cor_diff <- cor(lineage_pse, lineage_time, method="kendall")
    mse_vec <- c(mse_vec, mse_diff)
    cor_vec <- c(cor_vec, cor_diff)
    
    # compare orders
    mae_diff <- mean(abs(rank(lineage_pse)-rank(lineage_time)))
    mae_vec <- c(mae_vec, mae_diff)
  }
  
  finalres_mse <- mean(mse_vec)
  finalres_cor <- mean(cor_vec)
  finalres_mae <- mean(mae_vec)
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  pse$Ave <- rowMeans(pse, na.rm = T)
  
  plotcol <- colors[cut(pse$Ave, breaks=100)]
  plot(dimred, col = alpha(plotcol, 0.7), pch=16, main=paste(PLOTname, "Estimated PT", sep = ": "), 
       xlab=paste0("MSE: ", round(finalres_mse,3)))
  segments(x0 = fitLine$x0,                   
           y0 = fitLine$y0,
           x1 = fitLine$x1,
           y1 = fitLine$y1, lwd = 3)
  
  plotcol <- colors[cut(scdata_time, breaks=100)]
  plot(dimred, col = alpha(plotcol, 0.7), pch=16, main="True PT", 
       xlab=paste0("MAE: ", round(finalres_mae,3), "\n", 
                   "TauCOR: ", round(finalres_cor,3)))
  
  df <- data.frame(PT=pse$Ave, True=scdata_time)
  df <- df[order(scdata_time), ]
  
  plot(df$True, df$PT, ylab="Estimated PT", xlab="Scaled True Time", 
       col = alpha("grey", 0.7), pch=16, xlim=c(0, 1), ylim=c(0, 1))
  abline(coef = c(0,1), col="black", lty=2)
  
  return(list(AllMSE=mse_vec, AllCOR=cor_vec, AllMAE=mae_vec, 
              summary_df=data.frame(MeanMSE=finalres_mse, MeanCOR=finalres_cor, MeanMAE=finalres_mae)))
  
}




### MSE calculation
# pseudotime: estimated PT
# truetime: true PT
# output is MSE between estimated PT and true PT

mse <- function(pseudotime, truetime)
{
  if (cor(pseudotime, truetime, method="kendall") < 0) pseudotime <- rev(pseudotime)
  sum((pseudotime - truetime)^2) / length(truetime)
}
