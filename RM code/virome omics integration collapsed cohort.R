#Ruben Mars, PhD 20201125

#----------------------------------------------------------------------
require(compositions)
require(reshape2)
require(openxlsx)
library(dplyr)
library(data.table)

#----------------------------------------------------------------------
setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
source("create_corr_frame_rows.R")

#----------------------------------------------------------------------
#use the following datasets:

#collapsed contigs virome; df_contigs_collapsed_sel

#correlate with:

#collapsed microbiome species level; #dim(df_list_taxcollapsed$species)
#collapsed metabolomics; dim(scaled_mx_df_t)
#collapsed KEGG modules; dim(collapsed_kegg)


#----------------------------------------------------------------------
#subset everything by cohort and put in the same order

x_data <- clr(df_collapsed_sub) #clr transform as well
#hist(as.numeric(x_data))
y_data_list <- list(clr(df_list_taxcollapsed$species), scaled_mx_df_t, clr(collapsed_kegg))
names(y_data_list) <- c("metagenomics_species", "merged_metabolomics", "KEGG_terms")

general_meta$Cohort <- as.character(general_meta$Cohort)
general_meta$Cohort[general_meta$Cohort == "C"] <- "IBS-C"
general_meta$Cohort[general_meta$Cohort == "D"] <- "IBS-D"
general_meta$Cohort[general_meta$Cohort == "H"] <- "Healthy"

general_meta$study_id <- as.character(general_meta$study_id)

#get cohorts
virome_cohort <- as.character(sapply(colnames(x_data), function(x) {unique(general_meta$Cohort[which(general_meta$study_id == x)])}))
table(virome_cohort)
#Healthy   IBS-C   IBS-D 
#16      17      17 

#get indexes by cohort
ind_healthy <- virome_cohort == "Healthy"
ind_C <- virome_cohort == "IBS-C"
ind_D <- virome_cohort == "IBS-D"

cohort_ind_list <- list(ind_healthy, ind_C, ind_D)
names(cohort_ind_list) <- c("Healthy", "IBS-C", "IBS-D")

#index proper colnames as such
colnames(x_data)[ind_healthy]

#----------------------------------------------------------------------
#iterate over the input data; y_data_list
#iteratite over the cohort list; cohort_ind_list
#iterate over the rows; cross correlation code

#----------------------------------------------------------------------
#Below nested loop took 5 days to run on macbook pro in Sept 2020 so be careful

cor_res_list <- list()

for (i in 1:length(y_data_list)) {
  y_data_temp_full <- y_data_list[[i]]
  cor_res_list_cohort <- list()
  for (coh in 1:length(cohort_ind_list)) {
    print(coh)
    cohort_subjects_temp <- colnames(x_data)[cohort_ind_list[[coh]]]
    x_data_temp <- x_data[, colnames(x_data)[cohort_ind_list[[coh]]]]
    y_data_temp_sub <- y_data_temp_full[,cohort_subjects_temp]
    #remove rows with >80% zeroes by cohort
    x_data_sub <- x_data_temp[rowSums(x_data_temp == 0) < 0.8*ncol(x_data_temp), ]
    y_data_sub <- y_data_temp_sub[rowSums(y_data_temp_sub == 0) < 0.8*ncol(y_data_temp_sub), ]
    cor_res_list_cohort[[coh]] <- create_corr_frame_rows(x_data_sub, y_data_sub, method="spearman")
  }
  cor_res_list[[i]] <- cor_res_list_cohort
  print(i)
}
names(cor_res_list) <- names(y_data_list)


#could have done this naming in the loop but it takes very long to run
names(cor_res_list$metagenomics_species) <- names(cohort_ind_list)
names(cor_res_list$merged_metabolomics) <- names(cohort_ind_list)
names(cor_res_list$KEGG_terms) <- names(cohort_ind_list)


