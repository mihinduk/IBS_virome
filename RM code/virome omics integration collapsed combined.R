#Ruben Mars, PhD 20201202

#do not run seperately by cohort
#----------------------------------------------------------------------
require(compositions)
require(reshape2)
require(openxlsx)
library(dplyr)
library(data.table)

#----------------------------------------------------------------------
setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
source("create_corr_frame_rows.R") #correlation function

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

#----------------------------------------------------------------------
#iterate over the input data; y_data_list
#iterate over the cohort list; cohort_ind_list
#iterate over the rows; cross correlation code

#----------------------------------------------------------------------
# ! Can take couple hours / a day to run

cor_res_list_comb <- list()

for (i in 1:length(y_data_list)) {
  print(i)
  y_data_temp_full <- y_data_list[[i]]
  x_data_temp <- x_data[, colnames(x_data)]
  y_data_temp_sub <- y_data_temp_full[,colnames(x_data)]
  #remove rows with >80% zeroes by cohort
  x_data_sub <- x_data_temp[rowSums(x_data_temp == 0) < 0.8*ncol(x_data_temp), ]
  y_data_sub <- y_data_temp_sub[rowSums(y_data_temp_sub == 0) < 0.8*ncol(y_data_temp_sub), ]
  cor_res_list_comb[[i]] <- create_corr_frame_rows(x_data_sub, y_data_sub, method="spearman")
}
names(cor_res_list_comb) <- names(y_data_list)

