#compare with clinical symptoms from the longitudinal -omics paper

require(compositions)
require(vegan)
#----------------------------------------------------------------------

setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
source("create_corr_frame_rows.R") #correlation function


#----------------------------------------------------------------------
#quality of life
#general_meta$Effect_on_life_in_general_0_100

#number of stools
agg_total <- read.csv("../multi-omics data supplement/aggregate of nr stools per 7 days.csv", stringsAsFactors = F)

#Bristol stool scale per BM
agg_bss <- read.csv("../multi-omics data supplement/aggregate BSS.csv", stringsAsFactors = F)

#pain at passage does not make so much sense because it is linked to individual samples and we do not 
#have all samples in the virome


#----------------------------------------------------------------------
#----------------------------------------------------------------------
#quality of life

#subject average
general_meta_sub <- general_meta[which(!is.na(general_meta$Effect_on_life_in_general_0_100)),]
qol_collapsed <- aggregate(general_meta_sub$Effect_on_life_in_general_0_100, list(general_meta_sub$study_id), mean)
names(qol_collapsed) <- c("subjectID", "qol")
#hist(qol_collapsed$qol)
#add cohort data
qol_collapsed$cohort <- sapply(qol_collapsed$subjectID, function(x) unique(general_meta$cohort[which(general_meta$study_id == x)]))
#

#subset to only include the subjects that are in the virome study
qol_col_sub <- qol_collapsed[qol_collapsed$subjectID %in% virome_meta$subjectID,]
qol_col_sub <- qol_col_sub[-which(qol_col_sub$subjectID == "10007612"), ]
#set the Healthy subjects at 0 because no effect on QOL
#10007612 is listed as healthy with IBS symptoms; set to 0
healthy_ids <- unique(virome_meta_ord$subjectID[virome_meta_ord$Cohort == "Healthy"])
qol_col_sub_comb <- rbind(qol_col_sub,  as.data.frame(cbind("subjectID"=healthy_ids,
                               "qol"= rep(0, length(healthy_ids)),
                               "cohort"=rep("Healthy", length(healthy_ids)))))
#boxplot(as.numeric(qol_col_sub_comb$qol) ~ qol_col_sub_comb$cohort)


#subset collapsed contig data and reorder in the same order as qol_col_sub_comb
df_collapsed_sub_qol <- df_collapsed_sub[,colnames(df_collapsed_sub) %in% qol_col_sub_comb$subjectID]
df_collapsed_sub_qol <- df_collapsed_sub_qol[,order(colnames(df_collapsed_sub_qol))]

qol_col_sub_comb_ord <- qol_col_sub_comb[order(qol_col_sub_comb$subjectID),]

#correlation code 
y_data <- clr(df_collapsed_sub_qol) #clr transform as well

#hist(as.numeric(x_data))
x <- as.numeric(qol_col_sub_comb_ord$qol)
names(x) <- qol_col_sub_comb_ord$subjectID

#correlation with alpha diversity
alpha <- diversity(t(df_collapsed_sub_qol), "shannon")
cor.test(x, alpha, method = "spearman")
#not significant
#p-value = 0.8239 rho 0.03227196

#requires both input files to have rownames
#remove rows with >80% zeroes by cohort
y_data_sub <- y_data[rowSums(y_data == 0) < 0.8*ncol(y_data), ]

cor_list <- vector(mode = "list", length = nrow(y_data_sub))
for (i in 1:nrow(y_data_sub)) {
  y <- as.numeric(clr(y_data_sub[i,]))
  cor_list[[i]][1] <- cor(x, y, use="na.or.complete", method="spearman") #there might be  NAs in the data
  cor_list[[i]][2] <- cor.test(x, y, use="na.or.complete", method="spearman")$p.value
}
cor_df_qol <- as.data.frame(do.call(rbind, cor_list), stringsAsFactors=F)
names(cor_df_qol) <- c("cor", "p.val")
cor_df_qol$fdr <- p.adjust(cor_df_qol$p.val, "fdr")
rownames(cor_df_qol) <- rownames(y_data_sub)

#nothing obvious
par(mfrow=c(2,1))
plot(cor_df_qol$cor, -log10(cor_df_qol$fdr))
abline(h=-log10(0.05))
plot(cor_df_qol$cor, -log10(cor_df_qol$p.val))
abline(h=-log10(0.05))

cor_df_qol[which(cor_df_qol$cor > 0.4 & cor_df_qol$p.val < 0.01),]
#               cor       p.val       fdr
#contig_1322 0.4235791 0.002175855 0.3426674; dark_matter sp. cat1053
#contig_1587 0.4104699 0.003069035 0.3426674; not classified sp. cat1538
#contig_696  0.4046478 0.003560181 0.3426674; dark_matter sp. cat1339


#subset only phages and redo the fdr correction








#----------------------------------------------------------------------
#----------------------------------------------------------------------
#number of stools

#boxplot(agg_total$mean_per_week ~ agg_total$Cohort, las=1)

df_collapsed_sub_agg <- df_collapsed_sub[,colnames(df_collapsed_sub) %in% agg_total$subject_ID]
df_collapsed_sub_agg <- df_collapsed_sub_agg[,order(colnames(df_collapsed_sub_agg))]

agg_total_sub <- agg_total[agg_total$subject_ID  %in% colnames(df_collapsed_sub_agg), ]
agg_total_sub_ord <- agg_total_sub[order(agg_total_sub$subject_ID),]

#cbind(agg_total_sub_ord$subject_ID, colnames(df_collapsed_sub_agg))

#correlation code 
y_data <- clr(df_collapsed_sub_agg) #clr transform as well
#hist(as.numeric(x_data))
x <- as.numeric(agg_total_sub_ord$mean_per_week)
names(x) <- agg_total_sub_ord$subject_ID

#correlation with alpha diversity
alpha <- diversity(t(df_collapsed_sub_agg), "shannon")
cor.test(x, alpha, method = "spearman")
#not significant
#p-value = 0.6616  rho 0.06343947

#requires both input files to have rownames
#remove rows with >80% zeroes by cohort
y_data_sub <- y_data[rowSums(y_data == 0) < 0.8*ncol(y_data), ]

cor_list <- vector(mode = "list", length = nrow(y_data_sub))
for (i in 1:nrow(y_data_sub)) {
  y <- as.numeric(clr(y_data_sub[i,]))
  cor_list[[i]][1] <- cor(x, y, use="na.or.complete", method="spearman") #there might be  NAs in the data
  cor_list[[i]][2] <- cor.test(x, y, use="na.or.complete", method="spearman")$p.value
}
cor_df_nstool <- as.data.frame(do.call(rbind, cor_list), stringsAsFactors=F)
names(cor_df_nstool) <- c("cor", "p.val")
cor_df_nstool$fdr <- p.adjust(cor_df_nstool$p.val, "fdr")
rownames(cor_df_nstool) <- rownames(y_data_sub)

#nothing significant at fdr
par(mfrow=c(2,1))
plot(cor_df_nstool$cor, -log10(cor_df_nstool$fdr))
abline(h=-log10(0.05))
plot(cor_df_nstool$cor, -log10(cor_df_nstool$p.val)) 
abline(h=-log10(0.05))


cor_df_nstool[which(cor_df_nstool$cor > 0.4 & cor_df_nstool$p.val < 0.01),]
#                cor        p.val       fdr
#contig_1659 0.474466 0.0004993763 0.1922599; dark_matter sp. cat1128

min(cor_df_nstool$fdr)



#----------------------------------------------------------------------
#----------------------------------------------------------------------
#Bristol stool scale per BM

dim(agg_bss)
head(agg_bss)

#subject average BSS per week
agg_bss_collapsed <- aggregate(agg_bss$meanBSS_per_week, list(agg_bss$subject_ID), mean)
names(agg_bss_collapsed) <- c("subjectID", "agg_bss")
#hist(agg_bss_collapsed$agg_bss)
#add cohort data
agg_bss_collapsed$cohort <- sapply(agg_bss_collapsed$subjectID, function(x) unique(general_meta$cohort[which(general_meta$study_id == x)]))
#table(agg_bss_collapsed$cohort)
#boxplot(as.numeric(agg_bss_collapsed$agg_bss) ~ agg_bss_collapsed$cohort)

#subset collapsed contig data and reorder in the same order as df_collapsed_sub_aggbss
df_collapsed_sub_aggbss <- df_collapsed_sub[,colnames(df_collapsed_sub) %in% agg_bss_collapsed$subjectID]
df_collapsed_sub_aggbss <- df_collapsed_sub_aggbss[,order(colnames(df_collapsed_sub_aggbss))]

agg_bss_collapsed_sub <- agg_bss_collapsed[agg_bss_collapsed$subjectID  %in% colnames(df_collapsed_sub_aggbss), ]
agg_bss_collapsed_sub_ord <- agg_bss_collapsed_sub[order(agg_bss_collapsed_sub$subjectID),]

cbind(agg_bss_collapsed_sub_ord$subjectID, colnames(df_collapsed_sub_aggbss))


#correlation code 
y_data <- clr(df_collapsed_sub_aggbss) #clr transform as well
#hist(as.numeric(x_data))
x <- as.numeric(agg_bss_collapsed_sub_ord$agg_bss)
names(x) <- agg_bss_collapsed_sub_ord$subjectID

#correlation with alpha diversity
alpha <- diversity(t(df_collapsed_sub_aggbss), "shannon")
cor.test(x, alpha, method = "spearman")
#not significant
#p-value = 0.3714 rho 0.1288355


#requires both input files to have rownames
#remove rows with >80% zeroes by cohort
y_data_sub <- y_data[rowSums(y_data == 0) < 0.8*ncol(y_data), ]

cor_list <- vector(mode = "list", length = nrow(y_data_sub))
for (i in 1:nrow(y_data_sub)) {
  y <- as.numeric(clr(y_data_sub[i,]))
  cor_list[[i]][1] <- cor(x, y, use="na.or.complete", method="spearman") #there might be  NAs in the data
  cor_list[[i]][2] <- cor.test(x, y, use="na.or.complete", method="spearman")$p.value
}
cor_df_aggbss <- as.data.frame(do.call(rbind, cor_list), stringsAsFactors=F)
names(cor_df_aggbss) <- c("cor", "p.val")
cor_df_aggbss$fdr <- p.adjust(cor_df_qol$p.val, "fdr")
rownames(cor_df_aggbss) <- rownames(y_data_sub)

#nothing fdr corrected
par(mfrow=c(2,1))
plot(cor_df_aggbss$cor, -log10(cor_df_aggbss$fdr))
abline(h=-log10(0.05))
plot(cor_df_aggbss$cor, -log10(cor_df_aggbss$p.val))
abline(h=-log10(0.05))



cor_df_aggbss[which(cor_df_aggbss$cor > 0.4 & cor_df_aggbss$p.val < 0.01),]
#               cor        p.val       fdr
#contig_17   0.4352634 0.0015828019 0.9789799; crAss-like viruses sp. cat1631
#contig_2024 0.4717016 0.0005441772 1.0000000; not classified sp. cat1929

min(cor_df_aggbss$fdr)








