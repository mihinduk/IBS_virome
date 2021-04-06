# run correlations between diet and contigs
# run correlations on the sums per person data
setwd("~/Dropbox/IBS/virome/R AJ/")
source("../../../How-to/create_corr_frame.R")

ave_contigs <- read.delim("data/ave_contigs.txt", check.names = F)
ave_items_2 <- read.delim("data/ave_items_2.txt", check.names = F)
meta <- read.delim("data/meta.txt", check.names = F)

# subset contigs to just those differentially present
dif_present <- read.delim("data/differntially_present_contigs.txt")

ave_contigs <- ave_contigs[,rownames(dif_present)]

#limit to contigs in at least 3 people per group

get_contigs_per_group <- function(ave_contigs, ave_items_2, group){
  subset <- meta[meta$cohort == group,]
  keep <- unique(as.character(subset$participantID))
  
  sub_contigs <- ave_contigs[rownames(ave_contigs)%in% keep,]
  sub_asa <- ave_items_2[rownames(ave_items_2) %in% keep,]
  
  # limit to presence in more than 3 people
  n = 3
  presence <- sub_contigs
  presence[presence > 0] <- 1
  presence <- presence[,colSums(presence) >= n]
  
  sub_contigs <- sub_contigs[,colnames(presence)]
  
  contigs <- colnames(sub_contigs)
  
  # limit to food groups in more than 3 people
  n = 3
  presence <- sub_asa
  presence[presence > 0] <- 1
  presence <- presence[, colSums(presence) >= n]
  
  sub_asa <- sub_asa[,colnames(presence)]
  
  foodgroups <- colnames(sub_asa)
  
  names <- list(foodgroups = foodgroups, contigs = contigs)
  
  return(names) 
}

C_names <- get_contigs_per_group(ave_contigs, ave_items_2,group = "C")
D_names <- get_contigs_per_group(ave_contigs, ave_items_2,group = "D")
H_names <- get_contigs_per_group(ave_contigs, ave_items_2,group = "H")

contigs_keep <- c(C_names$contigs, D_names$contigs, H_names$contigs)
contigs_keep <- unique(contigs_keep) # 891

ave_contigs <- ave_contigs[,contigs_keep] 

items_keep <- c(C_names$foodgroups, D_names$foodgroups, H_names$foodgroups)
items_keep <- unique(items_keep) #32

ave_items_2 <- ave_items_2[,items_keep] #26

corr_g_L2 <- create_corr_frame(ave_contigs, ave_items_2, method = "spearman")
write.table(corr_g_L2, file = "output/correlations_L2_contigs.txt", sep = "\t", quote = F, row.names = F, col.names = T)


corr_g_L2_sigs <- corr_g_L2[corr_g_L2$qval <= 0.25,]
corr_g_L2_sigs <- corr_g_L2_sigs[order(corr_g_L2_sigs$correlation, decreasing = T),]

write.table(corr_g_L2_sigs, file = "output/correlations_L2_contigs_sigs.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# plot

require(ggplot2)

mycolors <- c(Healthy = "#666666", `IBS-C` = "#D95F02",`IBS-D` = "#0072B2")

# Overall
meta_plot <- meta[,colnames(meta) %in% c("participantID", "cohort")]
meta_plot <- unique(meta_plot)

# check
#meta_plot$participantID == rownames(ave_contigs)

plot <- data.frame(cohort = meta_plot$cohort,
                   Dried_fruits = ave_items_2$`L1_Fruits;L2_Dried_fruits`, 
                   contig_660 = ave_contigs$contig_660,
                   contig_1706 = ave_contigs$contig_1706)

plot$cohort <- as.factor(plot$cohort)
plot$cohort <- factor(plot$cohort, labels = c("IBS-C", "IBS-D", "Healthy"))

plot$cohort

ggplot(data = plot, aes(x = Dried_fruits, y = log10(contig_660),color = cohort)) +
  geom_point(size = 2) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  annotate("text", label = paste0("Spearman rank, q = ", round(corr_g_L2_sigs$qval[1], 3)), 
           x = 6, y = 3.1, size = 2) +
  xlab("Dried fruits, grams") +
  ylab("log10(contig_660)\nSquash vein yellowing virus")

ggsave(filename = "output/correlation_dried_fruit_contig_660.pdf", device = "pdf", height = 3, width = 4)


ggplot(data = plot, aes(x = Dried_fruits, y = log10(contig_1706), color = cohort)) +
  geom_point(size = 2) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  annotate("text", label = paste0("Spearman rank, q = ", round(corr_g_L2_sigs$qval[2], 3)), 
           x = 6, y = 3.6, size = 2) +
  xlab("Dried fruits, grams") +
  ylab("log10(contig_1706)\nPotyviridae")

ggsave(filename = "output/correlation_dried_fruit_contig_1706.pdf", device = "pdf", height = 3, width = 4)


# check without outliers
plot <- plot[-25,]

p1 <- cor.test(plot$Dried_fruits, plot$contig_660, method = "spearman", exact = F)$p.value
p2 <- cor.test(plot$Dried_fruits, plot$contig_1706, method = "spearman", exact = F)$p.value

ggplot(data = plot, aes(x = Dried_fruits, y = log10(contig_660),color = cohort)) +
  geom_point(size = 2) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  annotate("text", label = paste0("Spearman rank, p = ", round(p1, 4)), 
           x = 6, y = 0.5, size = 2) +
  xlab("Dried fruits, grams") +
  ylab("log10(contig_660)\nSquash vein yellowing virus") 

ggsave(filename = "output/correlation_dried_fruit_contig_660_rm_out.pdf", device = "pdf", height = 3, width = 4)


ggplot(data = plot, aes(x = Dried_fruits, y = log10(contig_1706), color = cohort)) +
  geom_point(size = 2) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  annotate("text", label = paste0("Spearman rank, p = ", round(p2, 3)), 
           x = 6, y = 0.5, size = 2) +
  xlab("Dried fruits, grams") +
  ylab("log10(contig_1706)\nPotyviridae")

ggsave(filename = "output/correlation_dried_fruit_contig_1706_rm_out.pdf", device = "pdf", height = 3, width = 4)

