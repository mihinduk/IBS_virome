# average ASA24
# average datasets
setwd("~/Dropbox/IBS/virome/R AJ/")

# load the metadata
meta <- read.delim(file = "data/meta.txt", check.names = F, stringsAsFactors = T)

# function to average by person
# df is a dataframe with samples as columns
# id is the labels for the columns in the original df
# par_id is the grouping factor/variable to identify pariticpants
# method is the method of aggregation
averagebyperson <- function(df, id, par_id, method) {
  id <- droplevels(id)
  if(sum(colnames(df) == id) == length(colnames(df))) {
    
    par_id <- as.factor(par_id)
    
    ave <- apply(df, 1, function(x){
      y <- aggregate(x, by = list(par_id), FUN = method);
      return(y$x)
    })
    rownames(ave) <- unique(par_id)
    ave <- as.data.frame(ave)
    return(ave)
  }
}

# create average contigs
# load the contigs
asa_items_1 <- read.delim(file = "data/asa_items_1.txt", check.names = F)
asa_items_2 <- read.delim(file = "data/asa_items_2.txt", check.names = F)
asa_items_3 <- read.delim(file = "data/asa_items_3.txt", check.names = F)


ave_items_1 <- averagebyperson(df = asa_items_1, 
                               id = meta$SampleID, 
                               par_id = meta$participantID, 
                               method = mean)

ave_items_2 <- averagebyperson(df = asa_items_2, 
                               id = meta$SampleID, 
                               par_id = meta$participantID, 
                               method = mean)


ave_items_3 <- averagebyperson(df = asa_items_3, 
                               id = meta$SampleID, 
                               par_id = meta$participantID, 
                               method = mean)

# write files
write.table(ave_items_1, file = "../R AJ/data/ave_items_1.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(ave_items_2, file = "../R AJ/data/ave_items_2.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(ave_items_3, file = "../R AJ/data/ave_items_3.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# read in ASA totals
asa_totals <- read.delim(file = "data/asa_totals.txt", check.names = F)
# subset to the correct samples
asa_totals <- asa_totals[colnames(asa_items_1),]

# transpose
asa_totals <- as.data.frame(t(asa_totals))

ave_asa_totals <- averagebyperson(df = asa_totals,
                                  id = meta$SampleID,
                                  par_id = meta$participantID,
                                  method = mean)

write.table(ave_asa_totals, file = "../R AJ/data/ave_asa_totals.txt", sep = "\t", quote = F, col.names = T, row.names = T)

