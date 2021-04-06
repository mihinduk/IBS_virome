# preprocess files so all come from people with both 2 ffqs, 2 dietary records, and 2 virome samples
setwd("~/Dropbox/IBS/virome/data_with_ids/")

# Load the IBS virome mapping data
meta <- read.delim("virome_meta.csv", sep = ",")

# load the diet map
map <- read.delim("../../../IBS_not_shared/Mapping/virome_diet_map.txt")

# rename map rownames
rownames(map) <- map$Timepoint.ID
map <- map[1]

# Load the IBS virome contig data
contigs <- read.delim("2020_09_10_IBS_phage_contigs_all_samples.txt", check.names = F)
contig_taxonomy <- read.delim("2020_08_19_Ruben_IBS_contig_taxonomy_table.txt", check.names = F)


# import the asa24 data (totals and items)
asa_totals_raw <- read.delim("../../../IBS_not_shared/Diet raw/Final ASA Totals.txt", check.names = F, row.names = 1)
asa_items <- read.delim("../../../IBS_not_shared/Food_tree_output/ibs.dhydrt.otu.txt", check.names = F, row.names = 1)
taxonomy <- asa_items$taxonomy
asa_items <- asa_items[,!colnames(asa_items) == "taxonomy"]

# start with ASA 24 totals
rownames(asa_totals_raw) <- gsub("\\.", "_", rownames(asa_totals_raw))
asa_totals_raw <- merge(map,asa_totals_raw, by = 0)
rownames(asa_totals_raw) <- asa_totals_raw$Virome_ID
asa_totals_raw <- asa_totals_raw[,-c(1,2)] # bad coding, but going quickly

# clean totals, drop NA rows
asa_totals_raw <- asa_totals_raw[!is.na(rowSums(asa_totals_raw)),] # this removes an additional sample from the set (n = 85)


# repeat with ASA 24 items (columns are samples)
colnames(asa_items) <- gsub("\\.", "_", colnames(asa_items))
asa_items <- t(asa_items)
asa_items <- merge(map, asa_items, by = 0, stringsAsFactors = FALSE)
rownames(asa_items) <- asa_items$Virome_ID
asa_items <- asa_items[,-c(1,2)]
asa_items <- t(asa_items)
asa_items <- as.data.frame(asa_items)
rownames(asa_items) <- taxonomy

# remove anything with no information
asa_items <- asa_items[,colSums(asa_items) != 0]

# subset asa and virome to the pairs that match
both <- intersect(colnames(asa_items), colnames(contigs))
asa_items <- asa_items[,both]
asa_totals_raw <- asa_totals_raw[both,]
contigs <- contigs[,both]

meta <- meta[meta$SampleID %in% both,]

# which participants also have ffqs
ffq_raw <- read.delim("../../../IBS_not_shared/Diet raw/Final FFQ.txt", check.names = F)

# people wiht ffqs 
ffq_ppl <- ffq_raw$`Study ID`

# asa24 people
asa_ppl <- meta[meta$SampleID %in% colnames(asa_items),]$subjectID

# virome people
vir_ppl <- meta[meta$SampleID %in% colnames(contigs),]$subjectID


# 100 people have virome data, 85 have asa 24, 136 have ffqs
# we need 2 ffqs and 2 asa24s for each person
ffq_ppl <- as.data.frame(table(ffq_ppl))
ffq_ppl_2 <- as.character(droplevels(ffq_ppl[ffq_ppl$Freq == 2,]$ffq_ppl)) #57

asa_ppl <- as.data.frame(table(asa_ppl))
asa_ppl_2 <- as.character(droplevels(asa_ppl)[asa_ppl$Freq == 2,]$asa_ppl) #37

vir_ppl <- as.data.frame(table(vir_ppl))
vir_ppl_2 <- as.character(droplevels(vir_ppl)[vir_ppl$Freq == 2,]$vir_ppl) #37


# find the intersection of all three strings
ids_keep <- Reduce(intersect,list(ffq_ppl_2,asa_ppl_2,vir_ppl_2)) #34

sample_ids <- as.character(droplevels(meta[meta$subjectID %in% ids_keep,]$SampleID)) #68

# reduce each set to these people

ffq <- ffq_raw[ffq_raw$`Study ID` %in% ids_keep,]

asa_items <- asa_items[,sample_ids]
asa_totals_raw <- asa_totals_raw[sample_ids,]
contigs <- contigs[,sample_ids]


# subset map to the samples in the virome dataset (86 samples)
map <- map[map$Virome_ID %in% sample_ids,,drop = FALSE]



# remove anything with no information
asa_items <- asa_items[,colSums(asa_items) != 0]

# Now need to get the items file at different levels of the food_tree taxonomy for easier comparison 
#summarize
split = strsplit(rownames(asa_items),";")
foodStrings = sapply(split,function(x) paste(x[1:3],collapse=";"))
asa_items_3 = rowsum(asa_items,foodStrings)

foodStrings = sapply(split,function(x) paste(x[1:2],collapse=";"))
asa_items_2 = rowsum(asa_items,foodStrings)

foodStrings = sapply(split,function(x) paste(x[1:1],collapse=";"))
asa_items_1 = rowsum(asa_items,foodStrings)


# find the overlap between the three sets

# now subset the metadata to just these individuals
meta <- meta[meta$SampleID %in% sample_ids,]
meta <- droplevels(meta)

# now subset the diet files
asa_items_1 <- asa_items_1[,colnames(asa_items_1) %in% meta$SampleID]
asa_items_2 <- asa_items_2[,colnames(asa_items_2) %in% meta$SampleID]
asa_items_3 <- asa_items_3[,colnames(asa_items_3) %in% meta$SampleID]

# subset the ffq files
ffq <- ffq[ffq$`Study ID` %in% ids_keep,]

# subset asa totals
asa_totals <- asa_totals_raw[rownames(asa_totals_raw) %in% meta$SampleID,]

# save these subsets
write.table(meta, file = "../R AJ/data/meta.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(contigs, file = "../R AJ/data/contigs.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(asa_items_1, file = "../R AJ/data/asa_items_1.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(asa_items_2, file = "../R AJ/data/asa_items_2.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(asa_items_3, file = "../R AJ/data/asa_items_3.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(ffq, file = "../R AJ/data/ffq.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(asa_totals, file = "../R AJ/data/asa_totals.txt", sep = "\t", quote = F, col.names = T, row.names = T )



