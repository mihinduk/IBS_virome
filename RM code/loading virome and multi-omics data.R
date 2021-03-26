#Ruben Mars, PhD 20201124

#----------------------------------------------------------------------
require(readxl)
#----------------------------------------------------------------------
#loading virome contig data and collapse
setwd("/Users/m210320/Dropbox/IBS/virome/data_with_ids/")

file_name <- "2020_09_10_IBS_phage_contigs_all_samples.txt"

df_contigs <- as.data.frame(read.csv(file_name, sep = "\t", row.names = 1), stringsAsFactors=F)
colnames(df_contigs) <- gsub("X", "", colnames(df_contigs))

df_contigs_ord <- df_contigs[,order(colnames(df_contigs))]

#----------------------------------------------------------------------
#metadata for the virome
virome_meta <- as.data.frame(read.csv("virome_meta.csv", sep = ","), stringsAsFactors=F)

virome_meta$Cohort <- as.character(virome_meta$cohort)
virome_meta$Cohort[virome_meta$Cohort == "C"] <- "IBS-C"
virome_meta$Cohort[virome_meta$Cohort == "D"] <- "IBS-D"
virome_meta$Cohort[virome_meta$Cohort == "H"] <- "Healthy"

#order in the same way
virome_meta_ord <- virome_meta[order(virome_meta$SampleID),]
#as.character(virome_meta_ord$SampleID)

#table(virome_meta$Cohort)

#----------------------------------------------------------------------
#collapse contig data based on subject
unique(colnames(df_contigs_ord) == as.character(virome_meta_ord$SampleID)) #these are in the same order

df_contigs_collapsed <- t(apply(df_contigs_ord, 1, function(xx) sapply(split(xx, virome_meta_ord$subjectID),mean)))
collapsed_cohort <- as.character(sapply(colnames(df_contigs_collapsed), function(x) unique(virome_meta_ord$Cohort[virome_meta_ord$subjectID == x])))

#write.csv(df_contigs_collapsed, file = "collapsed virome data.csv")

#----------------------------------------------------------------------
#subset contig data to only include the viruses and dark matter

#load contig annotation data
ann_data <- read.table("2020_08_19_Ruben_IBS_contig_taxonomy_table.txt", sep = "\t", header = T)
#dim(ann_data)
#table(ann_data$Kingdom) #Kingdom_undefined not classified Viruses
groups_to_include <- c("Kingdom_undefined", "not classified", "Viruses")
rownames(ann_data) <- ann_data$OTU

#some phage could be labeled under Kingdom Bacteria but most should be under viruses

#names
contig_names_to_include <- as.character(ann_data$OTU[which(ann_data$Kingdom %in% groups_to_include)])

df_contigs_collapsed_sel <- df_contigs_collapsed[rownames(df_contigs_collapsed) %in% contig_names_to_include,]
dim(df_contigs_collapsed_sel) #1211 50

#write.csv(df_contigs_collapsed_sel, file = "collapsed virome data sel.csv")

#test something that is expected to be different; Norwalk virus contig_1890 and contig_1168 
stripchart(df_contigs_collapsed_sel["contig_1168",] ~ collapsed_cohort)
stripchart(df_contigs_collapsed_sel["contig_1890",] ~ collapsed_cohort)

#table(df_contigs_collapsed_sel["contig_1168",]) #41 0s
#table(df_contigs_collapsed_sel["contig_1890",]) #48 zeroes

#----------------------------------------------------------------------
#number of unique species

length(unique(ann_data$Species)) 
#1552 species, consistent with Kathies number; collapsed by this before running statistics
#check how many after the subsetting

#----------------------------------------------------------------------
#remove contigs not present in at least 10% of people

cutoff <- 0.9*ncol(df_contigs_collapsed_sel)
#e.g. species that are not present in 90% or more are removed
df_collapsed_sub <- df_contigs_collapsed_sel[rowSums(df_contigs_collapsed_sel == 0) <= cutoff,]
dim(df_collapsed_sub)
#699  50

ann_data_sub <- ann_data[rownames(df_collapsed_sub),]
length(unique(ann_data_sub$Species))
#417 species out of 699 contigs; correlations should be done on these ideally..


#----------------------------------------------------------------------
#loading multi omics data; microbiome
#has multiple tabs so will do separate
setwd("../multi-omics data supplement")

#general metadata
general_meta <- as.data.frame(read.csv("Final_Cleaned_Condensed_Metadata.csv", sep = ","), stringsAsFactors=F)
general_meta$unique_id_nr <- paste(general_meta$study_id, general_meta$Timepoint_full, sep="_")

general_meta$cohort <- as.character(general_meta$Cohort)
general_meta$cohort[general_meta$cohort == "C"] <- "IBS-C"
general_meta$cohort[general_meta$cohort == "D"] <- "IBS-D"
general_meta$cohort[general_meta$cohort == "H"] <- "Healthy"


#microbiome all samples
file_name <- "IBS_microbiome_phylogenetic_levels_counts.xlsx"
excel_sheet_names <- excel_sheets(file_name) #we need 1 and 4 for stool taxatable and kegg orthology
df_list_tax <- lapply(excel_sheet_names, function(x) as.data.frame(read_excel(file_name, sheet=x), stringsAsFactors=F))
names(df_list_tax) <- excel_sheet_names

#set taxa column to rownames and remove, plus change the colnames #10007546_2 format #update to #"9.T.3" format

for (i in 1:length(df_list_tax)) {
  rownames(df_list_tax[[i]]) <- df_list_tax[[i]]$taxa
  df_list_tax[[i]] <- df_list_tax[[i]][,-1] #remove first column
  
  old_colnames <- colnames(df_list_tax[[i]])
  new_colnames <- sapply(old_colnames, function(x)  as.character(general_meta$SampleID[which(general_meta$unique_id_nr == x)]))
  colnames(df_list_tax[[i]]) <- new_colnames
}
#lapply(df_list_tax, colnames)


#microbiome collapsed samples
file_name <- "IBS_microbiome_phylogenetic_levels_counts_collapsed.xlsx"
excel_sheet_names <- excel_sheets(file_name) #we need 1 and 4 for stool taxatable and kegg orthology
df_list_taxcollapsed <- lapply(excel_sheet_names, function(x) as.data.frame(read_excel(file_name, sheet=x), stringsAsFactors=F))
names(df_list_taxcollapsed) <- excel_sheet_names

for (i in 1:length(df_list_taxcollapsed)) {
  rownames(df_list_taxcollapsed[[i]]) <- df_list_taxcollapsed[[i]]$taxa
  df_list_taxcollapsed[[i]] <- df_list_taxcollapsed[[i]][,2:ncol(df_list_taxcollapsed[[i]])]
}


#----------------------------------------------------------------------
#KEGG terms; collapsed
k_map <- read.table("ko-enzyme-annotations.txt", sep='\t', comment='')
collapsed_kegg <- t(as.data.frame(read.csv("keggc_collapsed_KOterms.tsv", sep="\t"), stringsAsFactors=F))
temp_colnames <- colnames(collapsed_kegg) #Subject 1 etc format
temp_colnames <- as.numeric(gsub("Subject_", "", temp_colnames))
new_colnames <- sapply(temp_colnames, function(x)  unique(as.character(general_meta$study_id[which(general_meta$ID_on_tube == x)])))
colnames(collapsed_kegg) <- new_colnames


#----------------------------------------------------------------------
#loading multi omics data; metabolomics collapsed

scaled_mx_data <- as.data.frame(read.csv("collapsed_mx_data.csv", stringsAsFactors = F))
rownames(scaled_mx_data) <- scaled_mx_data$subjectID
scaled_mx_df_t <- t(scaled_mx_data[2:ncol(scaled_mx_data)]) #25 metabolites, only NMR and BA methods

#----------------------------------------------------------------------
#checking Lactobacillus virus LBR48 and quality of life metrics

temp_data <- df_contigs_collapsed_sel[as.character(ann_data$OTU[which(ann_data$Species == "Lactobacillus virus LBR48")]),]

#temp_data[temp_data > 0]
#10007546     10007556     10007558     10007570     10007582     10007595     10007601 
#5.110382e-01 4.555824e+00 1.821976e+01 2.049087e+00 1.021840e+04 1.420399e+01 5.073330e-01 

#all IBS-C
#virome_meta_ord[virome_meta_ord$subjectID %in% names(temp_data[temp_data > 0]),]

#names(temp_data[temp_data > 0])


#----------------------------------------------------------------------
#inspect counts of 

#table(df_contigs_collapsed_sel["contig_1890",])

virome_cohort <- as.character(sapply(colnames(df_contigs_collapsed_sel), function(x) {unique(general_meta$Cohort[which(general_meta$study_id == x)])}))
temp_df <- as.data.frame(cbind(df_contigs_collapsed_sel["contig_1168",], 
                                df_contigs_collapsed_sel["contig_1168",] > 0, 
                                virome_cohort))
names(temp_df) <- c("norm_abundance", "presence_bool", "Cohort")
#table(temp_df$V2, temp_df$virome_cohort)
temp_df$subject_ID <- rownames(temp_df)

write.csv(temp_df, "Norovirus GII Cambria0299 abundance.csv", row.names = F)


