# average datasets
setwd("~/Dropbox/IBS/virome/R AJ/")

# load the metadata
meta <- read.delim(file = "data/meta.txt", check.names = F)

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
contigs <- read.delim(file = "data/contigs.txt", check.names = F)

# drop the bacteria, eukaryota, and archea
taxonomy <- read.delim("../data_with_ids/2020_08_19_Ruben_IBS_contig_taxonomy_table.txt")
keep_contigs <- as.character(taxonomy[grep("Kingdom_undefined|not classified|Viruses", taxonomy$Kingdom),]$OTU)

contigs <- contigs[keep_contigs,]

ave_contigs <- averagebyperson(df = contigs, 
                               id = meta$SampleID, 
                               par_id = meta$participantID, 
                               method = mean)


# write file
write.table(ave_contigs, file = "../R AJ/data/ave_contigs.txt", sep = "\t", quote = F, col.names = T, row.names = T)

