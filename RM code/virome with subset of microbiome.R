#compare with transcriptome work from the longitudinal omics paper

require(compositions)
#----------------------------------------------------------------------

setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
source("create_corr_frame_rows.R") #correlation function

#loading data
source('loading virome and multi-omics data.R')


#----------------------------------------------------------------------
#load transcriptome data both timepoints

ids_used_by_Sambhwaw_tx <- read.table("../data_with_ids/virome_rnaseq_overlapping_samples.txt")
#ids_used_by_Sambhwaw_tx$V1

#get cohort n
table(sapply(ids_used_by_Sambhwaw_tx$V1, function(x) unique(general_meta$cohort[which(general_meta$study_id == x)])))
#Healthy   IBS-C   IBS-D 
#10      11       7 


#----------------------------------------------------------------------
#get the subset of microbiome data

#species
y_data_s <- clr(df_list_taxcollapsed$species[,colnames(df_list_taxcollapsed$species) %in% ids_used_by_Sambhwaw_tx$V1])
#colnames(y_data)

#genus to try
y_data_g <- clr(df_list_taxcollapsed$genus[,colnames(df_list_taxcollapsed$genus) %in% ids_used_by_Sambhwaw_tx$V1])

#Filter rare species, only keep those in at least 20% of the samples
y_data_s <- y_data_s[rowSums(y_data_s > 0) >= ncol(y_data_s)*0.2,]
y_data_g <- y_data_g[rowSums(y_data_g > 0) >= ncol(y_data_g)*0.2,]


#----------------------------------------------------------------------
#get the virome data
## Filter rare viral contigs -- only keep contigs found in at least 20% of samples

x_data <- df_contigs_collapsed[,colnames(df_contigs_collapsed) %in% ids_used_by_Sambhwaw_tx$V1]
unique(colnames(y_data_s) == colnames(x_data))

contig_names_to_include <- as.character(ann_data$OTU[which(ann_data$Kingdom %in% groups_to_include)])
x_data_sel <- x_data[rownames(x_data) %in% contig_names_to_include,]
dim(x_data_sel) #1211 28

#subset contigs of interest
select <- rowSums(x_data_sel > 0) >= ncol(x_data_sel)*0.2
x_data_sub <- x_data_sel[select,]
dim(x_data_sub) 
#406 by 28


#same order
unique(colnames(x_data_sub) == colnames(y_data_s))


#----------------------------------------------------------------------
#correlate microbiome with virome abundances; species

#adapt correlation code 
x_data <- clr(x_data_sub) #clr transform as well
#hist(as.numeric(x_data))
#hist(as.numeric(x_data))

#requires both input files to have rownames
cor_res <- create_corr_frame_rows(x_data, clr(y_data_s), method="spearman")
dim(cor_res)

#NAs
length(which(is.na(cor_res$pval))) / nrow(cor_res) *100 # ~3% is na
cor_res_sel <- cor_res[!is.na(cor_res$pval),]

min(cor_res_sel$qval)

plot(cor_res_sel$correlation, -log10(cor_res_sel$qval))
abline(h=-log10(0.25))

cor_res_sel[which(cor_res_sel$qval <= 0.25),]

write.csv(cor_res_sel, file = "../transcriptome/cor_all_to_check_species.csv", row.names = F)


#----------------------------------------------------------------------
#to read in and check later
cor_res_sel <- read.csv("../transcriptome/cor_all_to_check_species.csv", header = T)


#----------------------------------------------------------------------
#inspect the contigs coming out of Sambhawa's analysis with the transcriptome

vir_tx <- read.csv("../transcriptome/Sambhawa to use/significant correlations with annotation.csv")
dim(vir_tx)
table(vir_tx$taxa)
length(unique(vir_tx$taxa))

contigs_of_int <- unique(vir_tx$taxa) #21
cor_res_sel_sign_coi <- cor_res_sel[which(cor_res_sel$x %in% contigs_of_int), ]

#make and add the full viral taxonomy to gene - virome correlations
viral_tax_col <- apply(ann_data, 1, function(x) paste(x, collapse=" "))

vir_tx$full_tax <- sapply(vir_tx$taxa, function(x) viral_tax_col[names(viral_tax_col) == x])
write.csv(vir_tx, file = "../transcriptome/Sambhawa to use/significant correlations with annotation full tax.csv", row.names = F)


plot(cor_res_sel_sign_coi$correlation, -log10(cor_res_sel_sign_coi$pval))
abline(h=-log10(0.05))
abline(h=-log10(0.01))

#using 0.01 cutoff makes sense
cor_sel_for_plot <- cor_res_sel_sign_coi[which(cor_res_sel_sign_coi$pval <= 0.01), ]
dim(cor_sel_for_plot)
table(cor_sel_for_plot$x)

#add the full viral taxonomy
cor_sel_for_plot$full_tax <- sapply(cor_sel_for_plot$x, function(x) viral_tax_col[names(viral_tax_col) == x])

write.csv(cor_sel_for_plot, file = "../transcriptome/Sambhawa to use/cor_sel_for_plot_species.csv", row.names = F)

#----------------------------------------------------------------------
#inspect all contigs coming out of Sambhawa's analysis

vir_tx_all <- read.csv("../transcriptome/Sambhawa to use/IBS_all_samples_gene_virome_w_taxonomy_FDR_0.1.csv")

nrow(vir_tx_all) #548
length(unique(vir_tx_all$gene)) #548
length(unique(vir_tx_all$taxa)) #116
sort(table(vir_tx_all$taxa), decreasing = T)

#add the full taxonomy
vir_tx_all$full_tax <- sapply(vir_tx_all$taxa, function(x) viral_tax_col[names(viral_tax_col) == x])

write.csv(vir_tx_all, file = "../transcriptome/Sambhawa to use/IBS_all_samples_gene_virome_w_taxonomy_FDR_0.1.csv", row.names = F)





#----------------------------------------------------------------------
#----------------------------------------------------------------------
#same for genus

#requires both input files to have rownames
cor_res <- create_corr_frame_rows(x_data, clr(y_data_g), method="spearman")
dim(cor_res)

#NAs
length(which(is.na(cor_res$pval))) / nrow(cor_res) *100 # ~3% is na
cor_res_sel <- cor_res[!is.na(cor_res$pval),]

min(cor_res_sel$qval)

plot(cor_res_sel$correlation, -log10(cor_res_sel$qval))
abline(h=-log10(0.25))

cor_res_sel[which(cor_res_sel$qval <= 0.25),]

write.csv(cor_res_sel, file = "../transcriptome/cor_all_to_check_genus.csv", row.names = F)

#----------------------------------------------------------------------