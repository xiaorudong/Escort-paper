
library(tibble)

score_cal <- function(scoredf, SimiRetaincutoff=0.806, GOFcutoff=0.53, URATEcutoff=0.019) {
  if (all(is.na(scoredf))) {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("very bad(1)", nrow(final_scoredf))
    final_scoredf <- tibble::rownames_to_column(final_scoredf, "Row.names")
    return(final_scoredf)
  }
  if(sum(scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff)>0) {
    df <- scoredf[scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff,]
    if (nrow(df)==1) {
      df$score <- 999
      df$ranking <- 999
      df$decision <- "Only choice"
      df <- df[,5:7]
      final_scoredf <- merge(scoredf, df, by="row.names", all=T)
    } else {
      df_norm <- data.frame(Scaled_GOF = (GOFcutoff-df$GOF)/GOFcutoff, Scaled_USHAPE=(URATEcutoff-df$USHAPE)/URATEcutoff)
      df_norm$score <- rowSums(df_norm)*df$SimiRetain
      df_norm <- round(df_norm, 3)
      rownames(df_norm) <- rownames(df)
      final_scoredf <- merge(scoredf, df_norm, by="row.names", all=T)
      final_scoredf$ranking[!is.na(final_scoredf$score)] <- rank(-final_scoredf$score, na.last = NA)
      final_scoredf$decision <- ifelse(final_scoredf$score>0, "good", "bad")
      final_scoredf$decision[is.na(final_scoredf$decision)] <- "very bad"
    }
  }else {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("very bad(1)", nrow(final_scoredf))
    final_scoredf <- tibble::rownames_to_column(final_scoredf, "Row.names")
  }
  return(final_scoredf)
}
