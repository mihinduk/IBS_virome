# visualize T1 vs T2 procrustes

#requires running this to get the files loaded
#setwd("/Users/m210320/Dropbox/IBS/virome/R RM/")
#source('loading virome and multi-omics data.R')


require(vegan)
require(ape)
require(ggplot2)
require(shapes)


virome_df <- df_contigs_ord

#subset only for the contigs used in the manuscript
groups_to_include <- c("Kingdom_undefined", "not classified", "Viruses")
contig_names_to_include <- as.character(ann_data$OTU[which(ann_data$Kingdom %in% groups_to_include)])

virome_df <- virome_df[contig_names_to_include,]
dim(virome_df)
#[1] 1211  100

virome_meta_char <- as.data.frame(sapply(virome_meta, as.character, 2), stringsAsFactors = F)
meta <- virome_meta_char[order(virome_meta_char$SampleID),]

meta$Cohort <- as.character(meta$cohort)
meta$Cohort[meta$Cohort == "C"] <- "IBS-C"
meta$Cohort[meta$Cohort == "D"] <- "IBS-D"
meta$Cohort[meta$Cohort == "H"] <- "Healthy"


#split timepoint 1 and 2
#colnames(virome_df)
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

#get the Bray Curtis distance separately by timepoint
beta_div_T1 <- vegdist(t(virome_df_ord_T1_filt), method = "bray") #as.matrix would make it accessible
beta_div_T2 <- vegdist(t(virome_df_ord_T2_filt), method = "bray") #as.matrix would make it accessible

# get pcoa vectors
pcoa_T1 <- as.data.frame(pcoa(beta_div_T1)$vectors)
pcoa_T2 <- as.data.frame(pcoa(beta_div_T2)$vectors)

# confirm the order is correct and both sets are the same
#rownames(pcoa_T1)
#rownames(pcoa_T2)

#-------------------------------------------------------------------------------
#try to get distances by group

#code from Robin Shields Cutler

pcoa_A <- pcoa_T1
pcoa_B <- pcoa_T2
metaA <- meta_T1
metaB <- meta_T2
thing <- 'subjectID'


set.seed(11)
crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)  # Run Procrustes
set.seed(11)
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif
A_crust <- data.frame(crusty$X)  # Recover the first group's coordinates
B_crust <- data.frame(crusty$Yrot)  # Recover the second group's coordinates

colnames(B_crust) <- colnames(A_crust)
ncoords = as.numeric(ncol(A_crust))
A_crust <- cbind(A_crust, metaA) #was previously merge(B_crust, metaB, by=0)
B_crust <- cbind(B_crust, metaB)
sample_ids <- A_crust[,"subjectID"]  # Get all the unifying participant IDs
real_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))  # Initialize the dataframe
colnames(real_dist) <- c('sampleID_A', 'sampleID_B', 'distance')

# Loop through each participant to get the multidimensional distance between their rotated points
for (i in 1:length(sample_ids)) {
  ix <- as.character(sample_ids[i])
  A_ix <- A_crust[A_crust[,thing] == ix, 1:ncoords]  # Keep all the PC axes
  B_ix <- B_crust[B_crust[,thing] == ix, 1:ncoords]
  AB_mat <- rbind(A_ix, B_ix)
  AB_dist <- matrix(dist(AB_mat, method = 'euclidean'))  # Calculates the distance
  real_dist[i,1] <- as.character(rownames(A_ix))  # Fill in the dataframe
  real_dist[i,2] <- as.character(rownames(B_ix))
  real_dist[i,3] <- as.numeric(AB_dist[1])
}

A_crust$data_type <- 'T1'; B_crust$data_type <- 'T2'
A_crust$realperm <- 'real'; B_crust$realperm <- 'real'

# The overall "Procrustes Distance"
real_pro_dist <- sqrt(sum(real_dist$distance^2))
cat(real_pro_dist)

#get the "Procrustes Distance" by group

#add cohort and flare destination to the df
real_dist$Cohort <- as.character(sapply(real_dist$sampleID_A, function(x) meta$Cohort[which(meta$SampleID == x)]))
flare_rows <- grep("Flare", real_dist$sampleID_B)

real_dist$Cohort_Flare <- real_dist$Cohort
real_dist$Cohort_Flare[flare_rows] <- paste(real_dist$Cohort_Flare[flare_rows], "Flare")


#make this into a proper boxplot
setwd("/Users/m210320/Dropbox/IBS/virome/R RM/figure_output/")

png("Euclidean Procrustes distance by group.png", width=4*300, height=3.5*300, res = 300)
par(mar=c(5,5,1,0))
b <- boxplot(as.numeric(real_dist$distance) ~ real_dist$Cohort_Flare, outline=F, xlab="", ylab="", 
             varwidth=T, axes=F)
stripchart(as.numeric(real_dist$distance) ~ real_dist$Cohort_Flare, vertical=T, 
           method = "jitter", jitter= 0.1, add = TRUE, pch=16, cex=1, col="black")
mtext("Euclidean Procrustes Distance", 2, line=3)
axis(2, las=1)
text(1:5, par('usr')[3], labels = b$names, srt = 45, adj = c(1.1,1.1), xpd = TRUE)
dev.off()

#TukeyHSD(aov(real_dist$distance ~ real_dist$Cohort_Flare)) #not significant as is visible


#-------------------------------------------------------------------------------
#IN PART REPEATED FROM ABOVE


# Set up the procrustes and plot

# procrustes
pro <- procrustes(pcoa_T1, pcoa_T2, symmetric = T)
pro_test <- protest(pcoa_T1, pcoa_T2, perm = 999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$ID <- rownames(beta_pro)
beta_pro$type <- "T1"
trans_pro$ID <- rownames(beta_pro) # hacky but should work, need to match the names
trans_pro$type <- "T2"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif)
plot <- rbind(beta_pro, trans_pro)

ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#fe9700", "#5f86b7")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=ID), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.05, y = -0.1, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

# save figure
setwd("/Users/m210320/Dropbox/IBS/virome/R RM/figure_output//")
ggsave("procrustes all samples.pdf", dev = "pdf", width = 4, height = 4 )



#facet plot
real_dist$Cohort_Flare[flare_rows] <- paste(real_dist$Cohort_Flare[flare_rows], "Flare")


real_dist$Cohort_Flare_pooled <- real_dist$Cohort_Flare
real_dist$Cohort_Flare_pooled[flare_rows] <- "Flare"
plot$groups <- real_dist$Cohort_Flare_poole

ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  facet_grid(. ~ groups) +
  scale_color_manual(values = c("#fe9700", "#5f86b7")) +
  theme_classic() +
  geom_line(aes(x= Axis.1, y=Axis.2, group=ID), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  #annotate("text", x = 0.00, y = -0.1, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

# save figure
ggsave("procrustes all samples facet.pdf", dev = "pdf", width = 5, height = 2.5)


