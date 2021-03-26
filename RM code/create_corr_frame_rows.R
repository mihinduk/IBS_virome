#creates data frame with column pairs, correlation, and p-value
create_corr_frame_rows <- function(df1, df2, method="spearman"){
  #calculate number of rows we will need in our final data frame
  rows <- nrow(df1)*nrow(df2)
  #create factor vectors for x and y columns (vertices)
  #set levels for x and y ahead of time
  x <- factor(levels=unique(rownames(df1)))
  y <- factor(levels=unique(rownames(df2)))
  #create numeric vectors for correlation and p-value columns
  correlation <- numeric(rows)
  pval <- numeric(rows)
  index <- 1
  for(i in rownames(df1)){
    for(j in rownames(df2)){
      a <- as.numeric(df1[i,])
      b <- as.numeric(df2[j,])
      cor <- cor.test(a, b, method=method)
      x[index] <- i
      y[index] <- j
      correlation[index] <- cor$estimate
      pval[index] <- cor$p.value
      index <- index + 1
    }
    print(i)
  }
  dat <- data.frame(x=x, y=y, correlation=correlation, pval=pval)
  dat$qval <- p.adjust(dat$pval, method = "fdr")
  dat
}
