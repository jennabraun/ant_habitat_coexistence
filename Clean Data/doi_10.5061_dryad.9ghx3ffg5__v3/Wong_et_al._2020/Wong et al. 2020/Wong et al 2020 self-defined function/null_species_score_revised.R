null_species_score <- function(df,data.clean,all=F) {
  
  if (all == F) {
  df <-df[df$value <1/3 | df$value > 3,]
  }
  
temp.sp.result <- list()
for (i in 1: max(df$j)) {
  temp.sp.df <- df[df$j == i, ]
  dummy <- data.frame(species_summary2(temp.sp.df, data.clean))
  dummy$sp <- rownames(dummy)
  temp.sp.result[[i]] <- dummy
}

final.sp.null <- do.call(rbind,temp.sp.result)
final.sp.null$diff <- final.sp.null$pos_str + final.sp.null$neg_str

average.diff <- c(by(final.sp.null$diff, final.sp.null$sp, mean))
sd.diff <- c(by(final.sp.null$diff, final.sp.null$sp, sd))
result <- data.frame(average.diff,sd.diff)

return(result)
}