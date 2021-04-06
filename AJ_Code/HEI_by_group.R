# is HEI different by group

setwd("~/Dropbox/IBS/virome/R AJ/")

require(ggplot2)
require(ggpubr)
require(ggbeeswarm)

# load the metadata
meta <- read.delim("data/meta.txt")
# load ffq to add HEI variable to this plot
ave_ffq <- read.delim("data/ave_ffq.txt")

map <- meta[c("participantID", "cohort")]
map <- map[!duplicated(map),]
rownames(map) <- map$participantID

ave_ffq <- merge(map, ave_ffq, by = 0)

colnames(ave_ffq)

ave_ffq$cohort <- as.factor(ave_ffq$cohort)
ave_ffq$cohort <- factor(ave_ffq$cohort, labels = c("IBS-C", "IBS-D", "Healthy"))
ave_ffq$cohort <- factor(ave_ffq$cohort, levels = c("Healthy", "IBS-C", "IBS-D"))

# difference in kcals?
kcal <- kruskal.test(ave_ffq$VioFFQ.Calories~ave_ffq$cohort)$p.value
aggregate(ave_ffq$VioFFQ.Calories, by = list(map$cohort), FUN = mean)


pvals_nutr <- apply(ave_ffq[,50:ncol(ave_ffq)], 2, function(x) {kruskal.test(x~ave_ffq$cohort)$p.value})
qvals <- p.adjust(pvals_nutr)
qvals <- qvals[qvals <0.25]
pvals <- sort(pvals_nutr[pvals_nutr <0.05])


# there are no differences in HEI and no differences in any other variables measured by ffq




# Calories from alcohol are not significant after adjustment for multiple comparisions.
# nominal pvalue 0.009 before adjustment.
ggplot(ave_ffq, aes(x = cohort, y = ave_ffq$A_CAL)) + geom_boxplot()+ geom_point() 
                 
my_comparisons <- list( c("IBS-C", "IBS-D"), c("IBS-C", "Healthy"), c("IBS-D", "Healthy") )
mycolors <- c(Healthy = "#666666", `IBS-C` = "#D95F02",`IBS-D` = "#0072B2")


ggplot(ave_ffq, aes(x = cohort, y = A_CAL, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_x_discrete(label = c("Healthy", "IBS-C", "IBS-D")) +
  scale_fill_manual(values = mycolors) +
  ggbeeswarm::geom_beeswarm() +
  theme_bw() +
  xlab("") +
  ylab("Alcohol (kilocalories)") +
  ylim(c(0,350)) +
  stat_compare_means(comparisons = my_comparisons, 
                     method.args = list(exact = F),
                     label.y = c(250, 275, 310), 
                     size = 3) +
  annotate("text", label = paste0("Kruskal-Wallis, p = ",round(pvals[1],3)), 
           x = 1.5, y = 350, size = 3) 

ggsave("output/alcohol.pdf", dev = "pdf", height = 3, width = 4)



# not different

# what was HEIscore for each group

aggregate(ave_ffq$HEIScore, by = list(map$cohort), FUN = mean)


# are there differences in ASA groups?
ave_asa <- read.delim("data/ave_items_2.txt")

ave_asa <- merge(map, ave_asa, by = 0)

ave_asa$cohort <- as.factor(ave_asa$cohort)
ave_asa$cohort <- factor(ave_asa$cohort, labels = c("IBS-C", "IBS-D", "Healthy"))
ave_asa$cohort <- factor(ave_asa$cohort, levels = c("Healthy", "IBS-C", "IBS-D"))


pvals_asa <- apply(ave_asa[,4:ncol(ave_asa)], 2, function(x) {kruskal.test(x~ave_asa$cohort)$p.value})
qvals <- p.adjust(pvals_asa)
qvals <- qvals[qvals <0.25]
pvals <- sort(pvals_asa[pvals_asa <0.05])

pairwise.wilcox.test(ave_asa$L1_Fruits.L2_Dried_fruits, ave_asa$cohort)
my_comparisons <- list(c("IBS-C", "Healthy"), c("IBS-D", "Healthy"))

ggplot(ave_asa, aes(x = cohort, y = L1_Fruits.L2_Dried_fruits, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(label = c("Healthy", "IBS-C", "IBS-D")) +
  scale_fill_manual(values = mycolors) +
  ggbeeswarm::geom_beeswarm() +
  theme_bw() +
  xlab("") +
  ylab("Dried Fruit") +
  ylim(c(0,39)) +
  stat_compare_means(size = 3) +
  stat_compare_means(comparisons = my_comparisons, 
                     method.args = list(exact = F),
                     label.y = c(32, 35), 
                     size = 3) 

ggsave("output/driedfruit.pdf", dev = "pdf", height = 3, width = 4)

pairwise.wilcox.test(ave_asa$L1_Milk_and_Milk_Products.L2_Creams_and_cream_substitutes, ave_asa$cohort, exact = F)
my_comparisons <- list(c("IBS-C", "Healthy"), c("IBS-D", "Healthy"), c("IBS-C", "IBS-D"))

ggplot(ave_asa, aes(x = cohort, y = L1_Milk_and_Milk_Products.L2_Creams_and_cream_substitutes, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(label = c("Healthy", "IBS-C", "IBS-D")) +
  scale_fill_manual(values = mycolors) +
  ggbeeswarm::geom_beeswarm() +
  theme_bw() +
  xlab("") +
  ylab("Creams and cream substitutes") +
  ylim(c(0,39)) +
  stat_compare_means(size = 3) +
  stat_compare_means(comparisons = my_comparisons, 
                     method.args = list(exact = F),
                     label.y = c(25, 30, 35), 
                     size = 3) 

ggsave("output/cream.pdf", dev = "pdf", height = 3, width = 4)

# are there differences in ASA groups?
ave_asa <- read.delim("data/ave_items_1.txt")

ave_asa <- merge(map, ave_asa, by = 0)

ave_asa$cohort <- as.factor(ave_asa$cohort)
ave_asa$cohort <- factor(ave_asa$cohort, labels = c("IBS-C", "IBS-D", "Healthy"))
ave_asa$cohort <- factor(ave_asa$cohort, levels = c("Healthy", "IBS-C", "IBS-D"))


pvals_asa <- apply(ave_asa[,4:ncol(ave_asa)], 2, function(x) {kruskal.test(x~ave_asa$cohort)$p.value})

ggplot(ave_asa, aes(x = cohort, y = L1_Milk_and_Milk_Products, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(label = c("Healthy", "IBS-C", "IBS-D")) +
  scale_fill_manual(values = mycolors) +
  ggbeeswarm::geom_beeswarm() +
  theme_bw() +
  xlab("") +
  ylab("Milk and milk products") +
  ylim(c(0, 330)) +
  stat_compare_means(size = 3) +
  stat_compare_means(comparisons = my_comparisons, 
                     method.args = list(exact = F),
                     label.y = c(240, 270 , 300), 
                     size = 3) 

ggsave("output/milk.pdf", dev = "pdf", height = 3, width = 4)


# compare ASA totals
ave_asa_totals <- read.delim("data/ave_asa_totals.txt")

map <- meta[c("participantID", "cohort")]
map <- map[!duplicated(map),]
rownames(map) <- map$participantID

ave_asa_totals <- merge(map, ave_asa_totals, by = 0)

colnames(ave_asa_totals)

ave_asa_totals$cohort <- as.factor(ave_asa_totals$cohort)
ave_asa_totals$cohort <- factor(ave_asa_totals$cohort, labels = c("IBS-C", "IBS-D", "Healthy"))
ave_asa_totals$cohort <- factor(ave_asa_totals$cohort, levels = c("Healthy", "IBS-C", "IBS-D"))

# difference in kcals?
kcal <- kruskal.test(ave_asa_totals$KCAL~ave_asa_totals$cohort)$p.value
aggregate(ave_asa_totals$KCAL, by = list(map$cohort), FUN = mean)

pvals_nutr <- apply(ave_asa_totals[,4:(ncol(ave_asa_totals)-1)], 2, function(x) {kruskal.test(x~ave_asa_totals$cohort)$p.value})
qvals <- p.adjust(pvals_nutr)
qvals <- qvals[qvals <0.25]
pvals <- sort(pvals_nutr[pvals_nutr <0.05])

