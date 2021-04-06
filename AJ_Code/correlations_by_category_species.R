# repeat for each IBS category

setwd("~/Dropbox/IBS/virome/R AJ/")
source("../../../../Documents/Projects/How-to/create_corr_frame.R")
source("assess_species_sums.R")

ave_items_2 <- read.delim("data/ave_items_2.txt", check.names = F)

run_cor_by_group <- function(ave_contigs, ave_items_2, group){
  subset <- meta[meta$cohort == group,]
  keep <- unique(as.character(subset$participantID))
  
  sub_contigs <- ave_contigs[rownames(ave_contigs)%in% keep,]
  sub_asa <- ave_items_2[rownames(ave_items_2) %in% keep,]
  
  # limit to presence in more than 3 people
  presence <- sub_contigs
  presence[presence > 0] <- 1
  presence <- presence[,colSums(presence) > 3]
  
  sub_contigs <- sub_contigs[,colnames(presence)]
  
  # limit to food groups in more than 3 people
  presence <- sub_asa
  presence[presence > 0] <- 1
  presence <- presence[, colSums(presence) > 3]
  
  sub_asa <- sub_asa[,colnames(presence)]
  
  corr <- create_corr_frame(sub_contigs, sub_asa)
  corr <- corr[corr$qval <= 0.2,]
  corr <- corr[order(corr$correlation, decreasing = T),]
  
  corr$cohort <- group
  
  return(corr)
}

corr_C <- run_cor_by_group(ave_contigs, ave_items_2, group = "C")
corr_D <- run_cor_by_group(ave_contigs, ave_items_2, group = "D")
corr_H <- run_cor_by_group(ave_contigs, ave_items_2, group = "H")

write.table(corr_C, file = "output/correlations_L2_species_C_sigs.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(corr_D, file = "output/correlations_L2_species_D_sigs.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(corr_H, file = "output/correlations_L2_species_H_sigs.txt", sep = "\t", quote = F, row.names = F, col.names = T)




