# average datasets
setwd("~/Dropbox/IBS/virome/R AJ/")

# load the metadata
meta <- read.delim(file = "data/meta.txt", check.names = F)

ffq <- read.delim(file = "data/ffq.txt", check.names = F)

# reduce features to the ones we care about for these purposes
ffq_feat <- colnames(ffq)
ffq_feat <- ffq_feat[-c(1,3:14,16:27)]

ffq <- ffq[,ffq_feat]

# average by person
particpantID <- as.factor(ffq$`Study ID`)

# average these values within each person
ave_ffq <- apply(t(ffq), 1, 
                 function(x){
                   y <- aggregate(x, by = list(particpantID), FUN = mean); 
                   return(y$x)
                 }
)


ave_ffq <- as.data.frame(par_ffq)

# match to the study ids used with the other files
map <- meta[,c("participantID", "subjectID")]
map <- map[!duplicated(map),]

# check
map$subjectID == par_ffq$`Study ID`

# use the participant id from the contigs and other file types as rownames
rownames(ave_ffq) <- map$participantID

write.table(ave_ffq, file = "../R AJ/data/ave_ffq.txt", sep = "\t", quote = F, col.names = T, row.names = T)



