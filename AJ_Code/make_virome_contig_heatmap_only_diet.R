# make a heatmap visualization of average virus content in diet subset
setwd("~/Dropbox/IBS/virome/R AJ/")

require(pheatmap)
require(RColorBrewer)
require(viridis)


# load the metadata
meta <- read.delim("data/meta.txt")
meta$participantID <- as.character(meta$participantID)

# load contigs
ave_contigs <- read.delim("data/ave_contigs.txt")

# load ffq to add HEI variable to this plot
ave_ffq <- read.delim("data/ave_ffq.txt")


# make annotation variables by participantID
ann_col <- meta[meta$participantID %in% rownames(ave_contigs),]
ann_col <- ann_col[,colnames(ann_col) %in% c("participantID", "disease", "disease2")]
ann_col <- unique(ann_col)

rownames(ann_col) <- ann_col$participantID
ann_col <- merge(ann_col, ave_ffq, by = 0 )

# just keep the HEI score
ann_col <- ann_col[colnames(ann_col) %in% c("participantID","disease","HEIScore")]
rownames(ann_col) <- ann_col$participantID
ann_col <- ann_col[,!colnames(ann_col) == "participantID"]
ann_col$disease <- as.character(ann_col$disease)
ann_col$disease[ann_col$disease == "constipation"] <- "IBS-C"
ann_col$disease[ann_col$disease == "diarrhea"] <- "IBS-D"
ann_col$disease[ann_col$disease == "healthy"] <- "Healthy"



# make just a presence absense plot with the sum of each person's two samples
presence <- t(ave_contigs)
presence[presence > 0] <- 1

ann_col <- ann_col[colnames(presence),]
# keep only differentially present contigs

# drop anything in everyone or in noone. Not doing this breaks Fisher's
#presence <- presence[!rowSums(presence) == 34,]
presence <- presence[!rowSums(presence) == 0,]

# check that all the Katie ID'd contigs are still there.
k_sigs <- read.delim("data/Kathie_sigs.txt", sep= "\t")
setdiff(k_sigs$Contig, rownames(presence)) # all there 

keeps <- as.character(k_sigs$Contig)

# subset to differentially abundant contigs
presence <- as.data.frame(presence[keeps,])

#reorder the annotation
ann_col <- ann_col[order(ann_col$disease, ann_col$HEIScore),]
presence <- presence[,rownames(ann_col)]

# create annotationrows to show kingdom origin
taxonomy <- read.delim("../data_with_ids/2020_08_19_Ruben_IBS_contig_taxonomy_table.txt", row.names = 1)

plants <- c("contig_1706", "contig_2338", "contig_660", "contig_1135","contig_1225", "contig_740", 
            "contig_1598", "contig_1229", "contig_1900", "contig_2352", "contig_2107","contig_2403")
taxonomy$`Plant origin` <- "No"
taxonomy[plants,"Plant origin"] <- "Yes"

taxonomy <- taxonomy[rownames(presence),]
taxonomy <- droplevels(taxonomy)

taxonomy$Kingdom <- as.character(taxonomy$Kingdom)
taxonomy[taxonomy$Kingdom == "Kingdom_undefined" |
                   taxonomy$Kingdom == "not classified",]$Kingdom <- "Dark matter"

taxonomy$Kingdom <- as.factor(taxonomy$Kingdom)

ann_rows <- taxonomy[,colnames(taxonomy) %in% c("Kingdom", "Plant origin", "Element"), drop = F]

#differentitally present contigs
write.table(ann_rows, "data/differntially_present_contigs.txt", sep = "\t", row.names = T, col.names = T)


source("Overall_contig_diet_correlations.R")
source("correlations_by_category.R")

ann_rows$`Diet Correlation All` <- "No"
#sig overall diet
sig_overall <- as.character(corr_g_L2_sigs$x)
ann_rows[sig_overall,"Diet Correlation All"] <- as.character(corr_g_L2_sigs$y)

# sig group diet
sig_C <- as.character(corr_C$x)
ann_rows$`Diet Correlation IBS-C` <- "No"
ann_rows[sig_C, "Diet Correlation IBS-C"] <- as.character(corr_C$y)

sig_D <- as.character(corr_D$x)
ann_rows$`Diet Correlation IBS-D` <- "No"
ann_rows[sig_D, "Diet Correlation IBS-D"] <- as.character(corr_D$y)


rownames(k_sigs) <- as.character(k_sigs$Contig)
k_sigs <- k_sigs[,"Category", drop= F]

ann_rows <- merge(ann_rows, k_sigs, by = 0, all.y = T)

colnames(ann_rows)[colnames(ann_rows) == "Category"] <- "Family"

ann_rows$`Family` <- as.character(ann_rows$`Family`)

#ann_rows[is.na(ann_rows$`Family`),]$`Family` <- "Not Significant"

ann_rows$`Family` <- factor(ann_rows$`Family`, 
                                                  levels = c("Podoviridae",
                                                             "Siphoviridae",
                                                             "Microviridae",
                                                             "Myoviridae",
                                                             "Potyviridae", 
                                                             "Virgaviridae",
                                                             "Dark matter"))
rownames(ann_rows) <- ann_rows$Row.names

ann_rows <- ann_rows[,c("Diet Correlation All", 
                        "Diet Correlation IBS-C",
                        "Diet Correlation IBS-D",
                        "Plant origin",
                        "Element",
                        "Family", 
                        "Kingdom")]


ann_rows <- ann_rows[ann_rows$`Diet Correlation All` != "No" | 
                       ann_rows$`Diet Correlation IBS-C` != "No" | 
                       ann_rows$`Diet Correlation IBS-D` != "No",]


presence_lim <- presence[rownames(ann_rows),]
rownames(presence_lim) <- paste0(rownames(presence_lim)," ", ann_rows$Element)

rownames(ann_rows) <- paste0(rownames(ann_rows)," ",ann_rows$Element)
ann_rows <- ann_rows[,c("Diet Correlation All", 
                        "Diet Correlation IBS-C",
                        "Diet Correlation IBS-D",
                        "Plant origin",
                        "Kingdom")]

ann_rows <- droplevels(ann_rows)

colnames(ann_col) <- c("Disease","HEI Score")

# Specify colors
ann_colors = list(
  `HEI Score` = c("White", "Dark Green"),
  Disease = c(Healthy = "#666666", `IBS-C` = "#D95F02",`IBS-D` = "#0072B2"),
  `Diet Correlation All` = c(`L1_Fruits;L2_Dried_fruits` = "Red", No = "White"),
  `Diet Correlation IBS-C` = c(#`L1_Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds;L2_Nuts_nut_butters_and_nut_mixtures` = "Brown",
                               `L1_Fats_Oils_and_Salad_Dressings;L2_Salad_dressings` = "Grey",
                               `L1_Fruits;L2_Other_fruits` = "Red",
                               `L1_Milk_and_Milk_Products;L2_Cheeses` = "Yellow",
                               `L1_Grain_Product;L2_Crackers_and_salty_snacks_from_grain` = "Orange",
                               `L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks` = "Blue",
                               #`L1_Vegetables;L2_Tomatoes_and_tomato_mixtures` = "Red",
                               No = "White"),
  `Diet Correlation IBS-D` = c(`L1_Vegetables;L2_Deepyellow_vegetables` = "Green", 
                                No = "White"),
  `Plant origin` = c(No = "White", Yes = "Black"),
  # `Family` = c(Potyviridae = "Black", 
  #              Microviridae = "Blue",
  #              Myoviridae = "Green",
  #              `Dark matter` = "Grey"),
  Kingdom = c(`Dark matter` = "Grey", Viruses = "Black")
)


myheatmap <- pheatmap(presence_lim, 
                      color = c("white", "#9370DB"), 
                      annotation_col = ann_col,
                      annotation_row = ann_rows,
                      annotation_names_row = T,
                      cluster_cols = F,
                      cluster_rows = T,
                      legend = F,
                      annotation_colors = ann_colors,
                      annotation_legend = F)



# code to save pheatmap output
save_pheatmap_png <- function(x, filename, width=3000, height=1200, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(myheatmap, "output/presenceheatmap_diet.png")

myheatmap <- pheatmap(presence_lim, 
                      color = c("white", "#9370DB"), 
                      annotation_col = ann_col,
                      annotation_row = ann_rows,
                      annotation_names_row = T,
                      cluster_cols = F,
                      cluster_rows = T,
                      legend = F,
                      annotation_colors = ann_colors,
                      annotation_legend = T)

# code to save pheatmap output
save_pheatmap_png <- function(x, filename, width=3000, height=3000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(myheatmap, "output/presenceheatmap_legend.png")



# plot the individual correlations

# Diarrhea
plot(ave_items_2$`L1_Vegetables;L2_Deepyellow_vegetables`, log10(ave_contigs$contig_2317))

# Constipation
plot(ave_items_2$`L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks`, log10(ave_contigs$contig_309))

corr_C$x
corr_C$y
