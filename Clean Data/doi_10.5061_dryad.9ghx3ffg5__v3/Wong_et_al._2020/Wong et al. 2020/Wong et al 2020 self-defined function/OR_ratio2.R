or.asymmetric2 <- function (dataset, run.check = FALSE) 
{
  require(sppairs)
  if (run.check) 
    dataset <- or.check(dataset)
  for (i in 1:2) {
    dataset[, i] <- factor(dataset[, i], levels = c(0, 1), 
                           labels = c("0", "1"))
  }
  count.table <- as.matrix(table(dataset))
  both.pres <- count.table[2, 2] +0.5
  ApresBabs <- count.table[2, 1] +0.5
  Bpres <- sum(count.table[, 2]) +0.5
  Babs <- sum(count.table[, 1]) +0.5
  odds.ratio <- (both.pres/ApresBabs)/(Bpres/Babs)
  return(odds.ratio)
}
