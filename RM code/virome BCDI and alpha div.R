#compare virome data with the BCDI analysis from the longitudinal omics paper

#----------------------------------------------------------------------
require(lmerTest)
require(psycho) 
require(vegan)
require(emmeans)
require(compositions)

setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
source("create_corr_frame_rows.R")

setwd("~/Dropbox/IBS/virome/R RM/")

#----------------------------------------------------------------------

#head(general_meta)
#general_meta$BCDI

x_data <- df_contigs_ord #need the non collapsed
#colnames(x_data)

#effectively removes all the healthy because they were not classified according to BCDI
general_meta_sub <- general_meta[general_meta$SampleID %in% colnames(x_data) & general_meta$BCDI %in% c("normal","irregular"),]
dim(general_meta_sub)

x_data_sub <- x_data[,colnames(x_data) %in% general_meta_sub$SampleID]  
x_data_order <- x_data_sub[,order(colnames(x_data_sub))]
general_meta_sub_order <- general_meta_sub[order(general_meta_sub$SampleID),]

#remove more than 80% zeroes
x_data_order_sub <- x_data_order[rowSums(x_data_order == 0) < 0.8*ncol(x_data_order), ]
dim(x_data_order_sub) #639  68

group <- general_meta_sub_order$BCDI
random_group <- general_meta_sub_order$study_id

pval_list <- list()
for (i in 1:nrow(x_data_order_sub)) {
  x <- as.numeric(x_data_order_sub[i,])
  fit <- lmer(x ~ group + (1|random_group), control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
  fit_emm <- emmeans(fit, "group")
  pval_list[[i]] <- as.numeric(pwpm(fit_emm)[1,2])
}
lmer_res_df <- as.data.frame(cbind(rownames(x_data_order_sub), as.numeric(unlist(pval_list))), stringsAsFactors = F)
names(lmer_res_df) <- c("contig_name", "lmer_pvalue")
lmer_res_df$lmer_pvalue <- as.numeric(lmer_res_df$lmer_pvalue)
lmer_res_df$fdr <- p.adjust(lmer_res_df$lmer_pvalue, method = "fdr") 

hist(lmer_res_df$lmer_pvalue, breaks=50)
#nothing associated to irregularity in the microbiome
min(lmer_res_df$lmer_pvalue)
min(lmer_res_df$fdr)



#----------------------------------------------------------------------
#microbial alpha diversity correlated with virome

#computing diversity and testing difference between the cohorts
x_data <- clr(df_contigs_collapsed) #collapsed data
y_genus_sub <- df_list_taxcollapsed$genus[,names(x_data)]

div_shannon <- diversity(y_genus_sub, index = "shannon", MARGIN = 2)
div_cohort <- as.character(sapply(names(div_shannon), function(x) {unique(general_meta$Cohort[which(general_meta$study_id == x)])}))
TukeyHSD(aov(div_shannon ~ div_cohort))
boxplot(div_shannon ~ div_cohort, las=1)

#D lower diversity in this cohort compared to C, not compared with Healthies

#----------------------------------------------------------------------
#three contigs are correlated with alpha diversity but plots don't look very convincing

cor_res_list_div <- create_corr_frame_rows(x_data, t(as.data.frame(div_shannon)), method="spearman")

plot(cor_res_list_div$correlation, -log10(cor_res_list_div$qval))
abline(h=-log10(0.25))
abline(v=c(-0.4, 0.4))

cor_res_list_div[which(abs(cor_res_list_div$correlation) >= 0.4 & cor_res_list_div$qval <= 0.25),]

par(mar=c(3,3,2,2))
plot(div_shannon, clr(x_data["contig_1502",]))
plot(div_shannon, clr(x_data["contig_1659",]))
plot(div_shannon, clr(x_data["contig_2009",]))

#----------------------------------------------------------------------
#now inspect average phage alpha diveristy with average bacterial alpha diversity

div_shannon_vir <- diversity(df_contigs_collapsed, index = "shannon", MARGIN = 2)
div_shannon_bac <- diversity(y_genus_sub, index = "shannon", MARGIN = 2)

#they are in the same order
cor.test(div_shannon_vir, div_shannon_bac, method = "spearman") #cor 0.185, p-value = 0.1978. Trend but not convincing

#----------------------------------------------------------------------






