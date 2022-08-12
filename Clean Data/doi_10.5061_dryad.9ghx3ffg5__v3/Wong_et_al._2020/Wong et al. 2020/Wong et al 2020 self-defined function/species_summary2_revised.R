species_summary2 <- function(df,data.clean) {
  require(plyr)
  
  name.list <- colnames(data.clean)
  
  pos.df <- df[df$log.OR > 0, ]
  neg.df <- df[df$log.OR < 0, ]
  
  link.list <- list()

  for (j in 1:length(name.list)){
    link.list[[j]] <- c(link_str =sum(df$log.OR[which((df$sp1 %in% name.list[[j]]))]),
                      link = length(df$log.OR[which((df$sp1 %in% name.list[[j]]))]),
                      pos_link = length(pos.df$Estimate[which((pos.df$sp1 %in% name.list[[j]]))]),
                      neg_link = length(neg.df$Estimate[which((neg.df$sp1 %in% name.list[[j]]))]),
                      pos_str = sum(pos.df$log.OR[which((pos.df$sp1 %in% name.list[[j]]))]),
                      neg_str = sum(neg.df$log.OR[which((neg.df$sp1 %in% name.list[[j]]))]),
                      LWD = mean(df$log.OR[which((df$sp1 %in% name.list[[j]]))]),
                      kurtosis = kurtosis(df$log.OR[which((df$sp1 %in% name.list[[j]]))]))
  }
  
 species_result <- do.call(rbind,link.list)
 rownames(species_result) <- colnames(data.clean)
  
  return(species_result)
}