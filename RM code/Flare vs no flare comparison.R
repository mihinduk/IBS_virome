# Compare flare no flare with Mann Whitney

#requires running this to get the files loaded
#setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
#source('loading virome and multi-omics data.R')

require(ggplot2)
require(shapes)

virome_df <- df_contigs_ord

#subset only for the contigs used in the manuscript
groups_to_include <- c("Kingdom_undefined", "not classified", "Viruses")
contig_names_to_include <- as.character(ann_data$OTU[which(ann_data$Kingdom %in% groups_to_include)])

virome_df <- virome_df[contig_names_to_include,]

virome_meta_char <- as.data.frame(sapply(virome_meta, as.character, 2), stringsAsFactors = F)
meta <- virome_meta_char[order(virome_meta_char$SampleID),]

meta$Cohort <- as.character(meta$cohort)
meta$Cohort[meta$Cohort == "C"] <- "IBS-C"
meta$Cohort[meta$Cohort == "D"] <- "IBS-D"
meta$Cohort[meta$Cohort == "H"] <- "Healthy"

#split timepoint 1 and 2
colnames(virome_df)
virome_df_ord <- virome_df[, order(colnames(virome_df))]
virome_df_ord_T1 <- virome_df_ord[,seq(1, ncol(virome_df_ord),2)]
virome_df_ord_T2 <- virome_df_ord[,seq(2, ncol(virome_df_ord),2)]

#split metadata from the 2 timepoints
meta_T1 <- meta[seq(1, nrow(meta),2),]
meta_T2 <- meta[seq(2, nrow(meta),2),]
#meta_T1$SampleID
#meta_T2$SampleID

#subset rows that values in both (some are nonzero)
rows_to_remove <- unique(c(which(rowSums(virome_df_ord_T2 == 0) == ncol(virome_df_ord_T2)),
which(rowSums(virome_df_ord_T1 == 0) == ncol(virome_df_ord_T1))))
virome_df_ord_T1_filt <- virome_df_ord_T1[-rows_to_remove,]
virome_df_ord_T2_filt <- virome_df_ord_T2[-rows_to_remove,]
#dim(virome_df_ord_T1_filt)
#dim(virome_df_ord_T2_filt)

#----------------------------------------------------------------------
#get flare subjects and subset the 2 matrices; all Flares are in T2

flare_subjects <- meta_T2$subjectID[grep("Flare", meta_T2$SampleID)]

T1_flare_data <- virome_df_ord_T1_filt[,which(meta_T1$subjectID %in% flare_subjects)]
T2_flare_data <- virome_df_ord_T2_filt[,which(meta_T2$subjectID %in% flare_subjects)]

#combine the 2 and subset ones with more than 90% zeroes
groups <- c(rep("pre-flare", ncol(T1_flare_data)), rep("flare", ncol(T1_flare_data)))
flare_data_combined <- cbind(T1_flare_data, T2_flare_data)

#remove >75% zeroes
cutoff <- 0.75*ncol(flare_data_combined)
#e.g. species that are not present in 75% or more are removed
flare_data_combined_sub <- flare_data_combined[rowSums(flare_data_combined == 0) <= cutoff,]

category <- as.factor(groups)
x <- flare_data_combined_sub
ix1 <- category == levels(category)[1]
ix2 <- category == levels(category)[2]
pvals <- apply(x,1,function(taxon) wilcox.test(taxon[ix1], taxon[ix2],exact=FALSE)$p.value)
stats <- apply(x,1,function(taxon) wilcox.test(taxon[ix1], taxon[ix2],exact=FALSE)$statistic)
qvals <- p.adjust(pvals,'fdr')

min(pvals) 
min(qvals) #FDR minimum of 1 

hist(pvals, breaks=25)
abline(v=0.05)

which(pvals <= 0.05) #none below 0.01 
#contig_5 contig_2250 
#6         266 

#----------------------------------------------------------------------













