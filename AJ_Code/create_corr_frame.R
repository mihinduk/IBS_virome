#creates data frame with column pairs, correlation, and p-value
create_corr_frame <- function(df1, df2, method="spearman"){
  #calculate number of rows we will need in our final data frame
  rows <- ncol(df1)*ncol(df2)
  #create factor vectors for x and y columns (vertices)
  #set levels for x and y ahead of time
  x <- factor(levels=unique(colnames(df1)))
  y <- factor(levels=unique(colnames(df2)))
  #create numeric vectors for correlation and p-value columns
  correlation <- numeric(rows)
  pval <- numeric(rows)
  index <- 1
  for(i in colnames(df1)){
    for(j in colnames(df2)){
      a <- as.numeric(df1[,i])
      b <- as.numeric(df2[,j])
      cor <- cor.test(a, b, method=method, exact = F)
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