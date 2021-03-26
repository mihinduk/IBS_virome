#interpreting and plotting correlation results

require(stringr)
library(gplots)
library(dplyr)
library(dendextend)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
#BiocManager::install("keggorthology")
library(keggorthology)
require(circlize)


#----------------------------------------------------------------------

setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")

#load correlation results
#save(file="cor_res_list_comb", cor_res_list_comb)

#merged by cohort
load("cor_res_list_comb")
lapply(cor_res_list_comb, dim)

#$metagenomics_species
#[1] 318395      5
#$merged_metabolomics
#[1] 9625    5
#$KEGG_terms
#[1] 2099405       5

#separated by cohort
load("cor_res_list")
lapply(cor_res_list, length)
lapply(cor_res_list$metagenomics_species, dim)
lapply(cor_res_list$merged_metabolomics, dim)
lapply(cor_res_list$KEGG_terms, dim)

#----------------------------------------------------------------------
#load contig annotation data
ann_data <- read.table("../data_with_ids/2020_08_19_Ruben_IBS_contig_taxonomy_table.txt", sep = "\t", header = T)
rownames(ann_data) <- ann_data$OTU

#----------------------------------------------------------------------
#extract significant ones per cohort and make correlation networks.

#----------------------------------------------------------------------
#subset only the groups_to_include in the networks; 

#if wanting to subset the metabolomics

cor_res_list_comb_sub <- cor_res_list_comb
#metabolites_to_keep <- c("acetate", "propionate", "butyrate", "hypoxanthine", "cholic.acid", "chenodeoxycholic.acid")
#cor_res_list_comb_sub$merged_metabolomics <- cor_res_list_comb$merged_metabolomics[cor_res_list_comb$merged_metabolomics$y %in% metabolites_to_keep,]
#lapply(cor_res_list_comb_sub, dim) #9625 to 2310

#currently no subsetting
cor_res_list_subset <- cor_res_list


#----------------------------------------------------------------------
#iterate over the combined results
#get ones with FDR 0.25 (have explored 0.5 before)

q_cutoff <- 0.25
p_cutoff <- 0.001
cor_cutoff <- 0.0

#no difference when using 0 or 0.4 as a cutoff
cor_res_list_sign_comb <- cor_res_list_comb_sub
for (i in 1:length(cor_res_list_comb_sub)) {
  cor_res_list_comb_sub_temp <- cor_res_list_comb_sub[[i]]
  sign_rows <- which(cor_res_list_comb_sub_temp$qval <= q_cutoff & abs(cor_res_list_comb_sub_temp$correlation) >= cor_cutoff)
  #sign_rows <- which(cor_res_list_comb_sub_temp$pval <= p_cutoff & abs(cor_res_list_comb_sub_temp$correlation) >= cor_cutoff)
  cor_res_list_sign_comb[[i]] <-  cor_res_list_comb_sub_temp[sign_rows,]
}
lapply(cor_res_list_sign_comb, dim)

#$metagenomics_species
#[1] 62  5
#$merged_metabolomics
#[1] 1 5
#$KEGG_terms
#[1] 719   5


#optionally be more strict with the KEGG terms (only do for plotting)
cor_res_list_sign_comb$KEGG_terms <- cor_res_list_sign_comb$KEGG_terms[which(cor_res_list_sign_comb$KEGG_terms$qval < 0.1),]

#cor_res_list_sign_comb$metagenomics_species
#cor_res_list_sign_comb$merged_metabolomics
#cor_res_list_sign_comb$KEGG_terms


#filtering of contigs by phage or eukaryotic

#input only has 385 contigs
included_contigs <- unique(cor_res_list_comb$metagenomics_species$x)
#385
#subset ann data for actually included contigs
ann_data_sub <- ann_data[ann_data$OTU %in% included_contigs,]


#phage contigs
phage_phyla <- c("Hofneiviricota", "Phixviricota", "Uroviricota")
phage_contig_names <- as.character(ann_data_sub$OTU[which(ann_data_sub$Phylum %in% phage_phyla)])
length(phage_contig_names) #99

#eukaryotic viruses
#eukaryotic host: Kitrinoviricota and Pisuviricota
#plant: Virgaviridae, Potyviridae and Bromoviridae
euk_phyla <- c("Kitrinoviricota", "Pisuviricota", "Virgaviridae", "Potyviridae", "Bromoviridae")
euk_contig_names <- as.character(ann_data_sub$OTU[which(ann_data_sub$Phylum %in% euk_phyla)])
length(euk_contig_names)
#9

#dark matter
unknown_contig_names <- ann_data_sub$OTU[grep("dark_matter|not classified", ann_data_sub$Element)]
length(unknown_contig_names)
#277


#inspect eukaryotic viruses among overall sign correlations
cbind(unlist(lapply(cor_res_list_sign_comb, function(x) table(as.character(x$x) %in% euk_contig_names))))
#metagenomics_species.FALSE   61
#metagenomics_species.TRUE     1
#merged_metabolomics.FALSE     1
#KEGG_terms.FALSE            716
#KEGG_terms.TRUE               3

#1/62 species
#3/719 enzymes

#inspect classified phage
cbind(unlist(lapply(cor_res_list_sign_comb, function(x) table(as.character(x$x) %in% phage_contig_names))))
#metagenomics_species.FALSE   50
#metagenomics_species.TRUE    12 # n of correlations that are phage among all of them
#merged_metabolomics.FALSE     1
#KEGG_terms.FALSE            669
#KEGG_terms.TRUE              50

#if testing for significance take into account total number of input species. 
#Not significant as 99/9 = >11X more input phage species

#12/62 species
#50/719 enzymes


#inspect unclassified
cbind(unlist(lapply(cor_res_list_sign_comb, function(x) table(as.character(x$x) %in% unknown_contig_names))))
#metagenomics_species.FALSE   13
#metagenomics_species.TRUE    49
#merged_metabolomics.TRUE      1
#KEGG_terms.FALSE             53
#KEGG_terms.TRUE             666


#----------------------------------------------------------------------
#----------------------------------------------------------------------
#extracting significant ones for the cohort-specific correlations

q_cutoff <- 0.25
p_cutoff <- 0.001
cor_cutoff <- 0.0

cor_res_list_sign <- cor_res_list_subset
for (i in 1:length(cor_res_list_subset)) {
  cor_res_list_sub_temp <- cor_res_list_subset[[i]]
  for (j in 1:length(cor_res_list_sub_temp)) {
    cor_res_list_sub_temp_cohort_temp <- cor_res_list_sub_temp[[j]]
    sign_rows <- which(cor_res_list_sub_temp_cohort_temp$qval <= q_cutoff & abs(cor_res_list_sub_temp_cohort_temp$correlation) >= cor_cutoff)
    #sign_rows <- which(cor_res_list_sub_temp_cohort_temp$pval <= p_cutoff & abs(cor_res_list_sub_temp_cohort_temp$correlation) >= cor_cutoff)
    
    cor_res_list_sign[[i]][[j]] <-  cor_res_list_sub_temp_cohort_temp[sign_rows,] 
  }
}

#use 0.1 for KEGG terms to make consistent with the overall integration (only do for plotting)
cor_res_list_sign$KEGG_terms$Healthy <- cor_res_list_sign$KEGG_terms$Healthy[which(cor_res_list_sign$KEGG_terms$Healthy$qval < 0.1), ]
cor_res_list_sign$KEGG_terms$`IBS-C` <- cor_res_list_sign$KEGG_terms$`IBS-C`[which(cor_res_list_sign$KEGG_terms$`IBS-C`$qval < 0.1), ]
cor_res_list_sign$KEGG_terms$`IBS-D` <- cor_res_list_sign$KEGG_terms$`IBS-D`[which(cor_res_list_sign$KEGG_terms$`IBS-D`$qval < 0.1), ]

lapply(cor_res_list_sign, function(x) { lapply(x, dim)})

#$metagenomics_species
#$metagenomics_species$Healthy
#[1] 4 5
#$metagenomics_species$`IBS-C`
#[1] 6 5
#$metagenomics_species$`IBS-D`
#[1] 15  5

#$merged_metabolomics
#$merged_metabolomics$Healthy
#[1] 0 5
#$merged_metabolomics$`IBS-C`
#[1] 0 5
#$merged_metabolomics$`IBS-D`
#[1] 0 5

#$KEGG_terms
#$KEGG_terms$Healthy
#[1] 32  5        #39 at FDR 0.25
#$KEGG_terms$`IBS-C`
#[1] 7 5          #9 at FDR 0.25
#$KEGG_terms$`IBS-D`
#[1] 16  5        #18 at FDR 0.25


#No eukaryotic viruses in here either
lapply(cor_res_list_sign, function(x) lapply(x, function (xx) { unique(euk_contig_names %in% as.character(xx$x)) }))

lapply(cor_res_list_sign, function(x) lapply(x, function (xx) { unique(phage_contig_names %in% as.character(xx$x)) }))


#----------------------------------------------------------------------
#----------------------------------------------------------------------

#inspecting type of KEGG modules

#load from R RM folder
keggmapfp <- paste("ko-enzyme-annotations.txt")
k_map <- read.table(keggmapfp, sep='\t', comment='') #6451 keggs
k_map_unique <- k_map %>% distinct(V2, V3, V4, V5, .keep_all = TRUE)
dim(k_map)

sign_KEGG_terms <- as.character(k_map$V5[k_map$V1 %in% cor_res_list_sign_comb$KEGG_terms$y])
#sign_KEGG_terms[grep("purine", tolower(sign_KEGG_terms))]
#sign_KEGG_terms[grep("xantine", tolower(sign_KEGG_terms))]

#----------------------------------------------------------------------
#load higher level categories of the KEGG terms

KEGG_ann <- read.csv(file = "KEGG_annotation.csv", sep = ",")

#possibly explore some kind of pathway enrichment / inspection of the KEGG categories?


#----------------------------------------------------------------------
#expand names of the KO terms and contigs

#remove the EC number, str_split_fixed from stringr package to only split on the first space
sign_KEGG_terms_names <- as.character(sapply(sign_KEGG_terms, function(x) str_split_fixed(x, "  ", n=2)[2]))

#trim the longest name of the KEGG term
sign_KEGG_terms_names[3] <- "saccharopine dehydrogenase"
sign_KEGG_terms_names[13] <- "23S rRNA-methyltransferase"
sign_KEGG_terms_names[18] <- "16S rRNA-methyltransferase"
sign_KEGG_terms_names[20] <- "acetaldehyde dehydrogenase"
sign_KEGG_terms_names[25] <- "23S rRNA-methyltransferase"
sign_KEGG_terms_names[30] <- "beta-ketoacyl synthase II"


#overwrite y necessary for below
cor_res_list_sign_comb$KEGG_terms$y <- paste(cor_res_list_sign_comb$KEGG_terms$y, sign_KEGG_terms_names)


comb_contig_names <- lapply(cor_res_list_sign_comb, function(x) paste(x$x, ann_data[as.character(x$x),"Element"]))
cor_res_list_sign_comb[[1]]$comb_contig_names <- comb_contig_names[[1]]
cor_res_list_sign_comb[[2]]$comb_contig_names <- comb_contig_names[[2]]
cor_res_list_sign_comb[[3]]$comb_contig_names <- comb_contig_names[[3]]


#----------------------------------------------------------------------
#for the cohort-specific correlations

#getting the correlations to eventually merge with the circos file below; no metabolites significant

circos_data_comb_species <- do.call(rbind, lapply(cor_res_list_sign$metagenomics_species, function(x) x[,c("x", "y", "correlation")]))
#add column with cohort from rownames
circos_data_comb_species$Cohort <- as.character(sapply(rownames(circos_data_comb_species), function(x) strsplit(x, "\\.")[[1]][1]))

circos_data_comb_kegg <- do.call(rbind, lapply(cor_res_list_sign$KEGG_terms, function(x) x[,c("x", "y", "correlation")]))
#add column with cohort from rownames
circos_data_comb_kegg$Cohort <- as.character(sapply(rownames(circos_data_comb_kegg), function(x) strsplit(x, "\\.")[[1]][1]))


#----------------------------------------------------------------------
#update names KEGG ids
sign_KEGG_terms_coh <- as.character(k_map$V5[k_map$V1 %in% circos_data_comb_kegg$y])
#remove the EC number, str_split_fixed from stringr package to only split on the first space
sign_KEGG_terms_coh_names <- as.character(sapply(sign_KEGG_terms_coh, function(x) str_split_fixed(x, "  ", n=2)[2]))

#trim the longest name of the KEGG term
sign_KEGG_terms_coh_names[25] <- "site-specific DNA-methyltransferase"

#overwrite y necessary for below
circos_data_comb_kegg$y <- paste(circos_data_comb_kegg$y, sign_KEGG_terms_coh_names)

circos_data_coh <- circos_data_comb_species[,c("x", "y", "correlation", "Cohort")]
circos_data_coh <- rbind(circos_data_coh, circos_data_comb_kegg[,c("x", "y", "correlation", "Cohort")])

#update names of the contigs
circos_data_coh$comb_contig_names <- sapply(circos_data_coh$x, function(x) paste(x, ann_data[as.character(x),"Element"]))


#----------------------------------------------------------------------
#circos plot for the interrelations?; 
#combine them all into one circle; color by cohort?

circos.clear()

circos_data <- cor_res_list_sign_comb$metagenomics_species[,c("comb_contig_names", "y", "correlation")]
circos_data <- rbind(circos_data, cor_res_list_sign_comb$merged_metabolomics[,c("comb_contig_names", "y", "correlation")])
circos_data <- rbind(circos_data, cor_res_list_sign_comb$KEGG_terms[,c("comb_contig_names", "y", "correlation")])
circos_data$Cohort <- "Pooled"

#inspect shared interactions

circos_data_coh[which(paste(circos_data_coh$comb_contig_names, circos_data_coh$y) %in% paste(circos_data$comb_contig_names, circos_data$y)),]
circos_data[which(paste(circos_data$comb_contig_names, circos_data$y) %in% paste(circos_data_coh$comb_contig_names, circos_data_coh$y)),]

#the following are in both groups a cohort and a combined correlation; be mindful of this when making the plot

#x                            y correlation  Cohort                   comb_contig_names
#Healthy.128406 contig_1975    s__Flavonifractor_plautii  -0.8676471 Healthy contig_1975 dark_matter sp. cat1192
#IBS-C.12920    contig_1082 s__Subdoligranulum_variabile   0.8750000   IBS-C contig_1082 dark_matter sp. cat1013
#IBS-D.79269    contig_1577  s__Ruminococcus_bicirculans   0.8480392   IBS-D contig_1577 dark_matter sp. cat1106
#IBS-D.79273    contig_1577 s__Ruminococcus_flavefaciens   0.9164647   IBS-D contig_1577 dark_matter sp. cat1106

#comb_contig_names                            y correlation Cohort
#11541  contig_1082 dark_matter sp. cat1013 s__Subdoligranulum_variabile   0.7080862 Pooled
#70174  contig_1577 dark_matter sp. cat1106  s__Ruminococcus_bicirculans   0.7016760 Pooled
#70178  contig_1577 dark_matter sp. cat1106 s__Ruminococcus_flavefaciens   0.5561934 Pooled
#131057 contig_1975 dark_matter sp. cat1192    s__Flavonifractor_plautii  -0.6026355 Pooled

#----------------------------------------------------------------------
#-------------------- --------------------------------------------------
#prep separately for all cohorts separately and the combined circos plot


#full cohort
names(circos_data) <- c("from", "to", "value", "cohort")
circos_data$from <- as.character(circos_data$from)
circos_data$to <- as.character(circos_data$to)
circos_data$value <- as.numeric(circos_data$value)
circos_data$cohort <- as.character(circos_data$cohort)


type_list <- list(c("Metagenomic species"), 
                  c("Metabolites"), 
                  c("KEGG enzymes"))
names(type_list) <- c("Metagenomic species", "Metabolites", "KEGG enzymes")

input_cols <- c(brewer.pal(12, "Paired"), brewer.pal(6, "Set2"))
type_list_colors <- c(input_cols[10], input_cols[6], input_cols[4])
sel_cols_plot <- type_list_colors
sel_colors <- c("#666666", "#D95F02", "#0072B2","#000000")  #H, C, D, combined

#correlation dataset types; hardcoded
dataset_group <- c(rep("Metagenomic species", 62), c("Metabolites"), rep("KEGG enzymes", 69))

#names should be the same as the name!
names_grid_col <- c(circos_data$from, circos_data$to) #first nrow(circos_data) should be white (contigs)
grid.col_input <- c(rep("white", nrow(circos_data)), rep("black", nrow(circos_data)))
names(grid.col_input) <- names_grid_col

source("circos plots virome combined only.R")


#----------------------------------------------------------------------
#all cohorts separately

circos_data <- circos_data_coh
circos_data$x <- circos_data$comb_contig_names

names(circos_data) <- c("from", "to", "value", "cohort")
circos_data$from <- as.character(circos_data$from)
circos_data$to <- as.character(circos_data$to)
circos_data$value <- as.numeric(circos_data$value)
circos_data$cohort <- as.character(circos_data$cohort)


type_list <- list(c("Metagenomic species"),
                  c("KEGG enzymes"))
names(type_list) <- c("Metagenomic species", "KEGG enzymes")


#correlation dataset types; hardcoded
dataset_group <- c(rep("Metagenomic species", 25), rep("KEGG enzymes", 80-25))

#names should be the same as the name!
names_grid_col <- c(circos_data$from, circos_data$to) #first nrow(circos_data) should be white (contigs)
grid.col_input <- c(rep("white", nrow(circos_data)), sel_colors[factor(circos_data$cohort)])
names(grid.col_input) <- names_grid_col

source("circos plots virome cohort only.R")


#----------------------------------------------------------------------
#combined with cohort-specific correlations into one plot but this becomes too big for the manuscript

circos_data <- cor_res_list_sign_comb$metagenomics_species[,c("comb_contig_names", "y", "correlation")]
circos_data <- rbind(circos_data, cor_res_list_sign_comb$merged_metabolomics[,c("comb_contig_names", "y", "correlation")])
circos_data <- rbind(circos_data, cor_res_list_sign_comb$KEGG_terms[,c("comb_contig_names", "y", "correlation")])
circos_data$Cohort <- "Pooled"

circos_data <- rbind(circos_data, circos_data_coh[,c("comb_contig_names", "y", "correlation", "Cohort")])

names(circos_data) <- c("from", "to", "value", "cohort")
circos_data$from <- as.character(circos_data$from)
circos_data$to <- as.character(circos_data$to)
circos_data$value <- as.numeric(circos_data$value)
circos_data$cohort <- as.character(circos_data$cohort)


type_list <- list(c("Metagenomic species"), 
                  c("Metabolites"), 
                  c("KEGG enzymes"))
names(type_list) <- c("Metagenomic species", "Metabolites", "KEGG enzymes")

input_cols <- c(brewer.pal(12, "Paired"), brewer.pal(6, "Set2"))
type_list_colors <- c(input_cols[10], input_cols[6], input_cols[4])
sel_cols_plot <- type_list_colors
sel_colors <- c("#666666", "#D95F02", "#0072B2","#000000")  #H, C, D, combined


#correlation dataset types; hardcoded
dataset_group <- c(rep("Metagenomic species", 62), c("Metabolites"), rep("KEGG enzymes", 69), 
                   rep("Metagenomic species", 25), rep("KEGG enzymes", 55))


#names should be the same as the name!
names_grid_col <- c(circos_data$from, circos_data$to) #first 212 (contigs) should be white
grid.col_input <- c(rep("#FFFFFF", nrow(circos_data)), sel_colors[factor(circos_data$cohort)])
names(grid.col_input) <- names_grid_col

source("circos plots virome combined and cohort.R")


#----------------------------------------------------------------------
#only a single metabolite at FDR 0.25


#if the plot is needed input data has to be corrected



#cor_res_list_sign_comb$merged_metabolomics
#              x                     y correlation         pval      qval
#899 contig_1222 chenodeoxycholic.acid  -0.5716687 2.066864e-05 2.066864e-05

#rownames(y_data_list$merged_metabolomics) #metabolite 24

#y_temp <- y_data_list$merged_metabolomics[24,colnames(x_data)]
#x_temp <- x_data['contig_1222',colnames(x_data)]
#r2 <- -0.5716687
#pval <- 2.066864e-05
#fdr <- 2.066864e-05

#looks OK, not great
#pdf("CDCA contig 1222 cor.pdf", width=3.5, height=3)
#par(mar=c(3,4,2.5,1))
#plot(as.numeric(x_temp), as.numeric(y_temp), xlab="", ylab="", las=1, axes=F, pch=16, 
#     col=rgb(0.6,0.6,0.6,0.6))
#axis(1, cex.axis=0.7, line=0.2)
#axis(2, las=1, cex.axis=0.7, line=0.2)

#mtext("clr(contig abundance)", side = 1, line=2, cex=0.8)
#mtext("scaled metabolite abundance", side = 2.5, line=3.2, cex=0.8)
#title("CDCA vs. contig_1222\nundefined viral sequence", line=0.4, cex.main=0.7)

#stats_temp <- paste("R2=", round(r2, 3),"\nFDR=", round(fdr,5), "\np=", round(pval,5), sep="")
#legend("topleft", stats_temp, bty = "n", cex=0.5)
#dev.off()


#----------------------------------------------------------------------