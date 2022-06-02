suppressWarnings(suppressPackageStartupMessages(library("TCGAbiolinks")))
suppressWarnings(suppressPackageStartupMessages(library(readr)))
suppressWarnings(suppressPackageStartupMessages(library(SummarizedExperiment)))

#input <- commandArgs(1)
input <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/COADREAD_genes_normalized__data.txt"

suppressWarnings(suppressPackageStartupMessages(library(purrr)))
cancer_type <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
cancer_type <-strsplit(cancer_type, "\\.")[[1]]
cancer_type <- unlist(cancer_type)
cancer_type <- cancer_type[1]
cancer_type <-strsplit(cancer_type, "\\_")[[1]]
cancer_type <- cancer_type[1]

message("Reading input...")
originalDF <- suppressMessages(read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE))
originalDF <- originalDF[-1,]
originalDF <- originalDF[-c(1:29),]
geneIDS <- as.data.frame(sapply(strsplit(originalDF$`Hybridization REF`,split= "\\|"),'[',2))
library(stringr)
message("Filtering out unsufficiently quantified genes...")
to_del_two_thrids=which(rowSums(originalDF == "0.0000") > (1/4* (length(originalDF) - 1 )))
to_keep_all=which(rowSums(originalDF == "0.0000") == 0)
unfilteredDF <- originalDF
filteredDF_two_thirds <- subset(originalDF,  ! rownames(originalDF) %in% to_del_two_thrids )
filteredDF_all <- subset(originalDF,   rownames(originalDF) %in% to_keep_all )
row.names(filteredDF_two_thirds) <- filteredDF_two_thirds$`Hybridization REF`
filteredDF_two_thirds$`Hybridization REF` <- NULL
row.names(filteredDF_all) <- filteredDF_all$`Hybridization REF`
filteredDF_all$`Hybridization REF` <- NULL
row.names(unfilteredDF) <- unfilteredDF$`Hybridization REF`
unfilteredDF$`Hybridization REF` <- NULL
message("Dividing cancer and healthy samples...")
dfCancer_two_thirds <- subset(filteredDF_two_thirds, select = grepl("^....-..-....-0", names(filteredDF_two_thirds)))
dfHealthy_two_thirds <- subset(filteredDF_two_thirds, select = grepl("^....-..-....-1", names(filteredDF_two_thirds)))
dfCancer_all <- subset(filteredDF_all, select = grepl("^....-..-....-0", names(filteredDF_all)))
dfHealthy_all <- subset(filteredDF_all, select = grepl("^....-..-....-1", names(filteredDF_all)))
dfCancer_unfil <- subset(unfilteredDF, select = grepl("^....-..-....-0", names(unfilteredDF)))
dfHealthy_unfil <- subset(unfilteredDF, select = grepl("^....-..-....-1", names(unfilteredDF)))
row.names(dfCancer_unfil) <- row.names(unfilteredDF)
row.names(dfHealthy_unfil) <- row.names(unfilteredDF)
row.names(dfCancer_two_thirds) <- row.names(filteredDF_two_thirds)
row.names(dfHealthy_two_thirds) <- row.names(filteredDF_two_thirds)
row.names(dfCancer_all) <- row.names(filteredDF_all)
row.names(dfHealthy_all) <- row.names(filteredDF_all)

df_two_thirds <- merge(dfHealthy_two_thirds, dfCancer_two_thirds, by= 0)
row.names(df_two_thirds) <- df_two_thirds$Row.names
df_two_thirds$Row.names <- NULL
df_all <- merge(dfHealthy_all, dfCancer_all, by= 0)
row.names(df_all) <- df_all$Row.names
df_all$Row.names <- NULL
df_unfil <- merge(dfHealthy_unfil, dfCancer_unfil, by= 0)
row.names(df_unfil) <- df_unfil$Row.names
df_unfil$Row.names <- NULL

#fdat_two_thirds<- t(as.data.frame(strsplit(row.names(df_two_thirds),split= "\\|")))
#fdat_all<- t(as.data.frame(strsplit(row.names(df_all),split= "\\|")))
fdat_two_thirds <- as.data.frame(sapply(strsplit(row.names(df_two_thirds),split= "\\|"),'[',2))
#fdat_two_thirds <- cbind(fdat_two_thirds, "bs")
#fdat_two_thirds$"bs" <- fdat_two_thirds$`sapply(strsplit(row.names(df_two_thirds), split = "\\|"), "[", 2)`
fdat_two_thirds$"bs"<- sapply(strsplit(row.names(df_two_thirds),split= "\\|"),'[',2)
fdat_all <- as.data.frame(sapply(strsplit(row.names(df_all),split= "\\|"),'[',2))
fdat_all$"bs" <- sapply(strsplit(row.names(df_all),split= "\\|"),'[',2)
fdat_unfil <- as.data.frame(sapply(strsplit(row.names(df_unfil),split= "\\|"),'[',2))
fdat_unfil$"bs" <- sapply(strsplit(row.names(df_unfil),split= "\\|"),'[',2)

pdat_two_thirds <- as.data.frame(colnames(df_two_thirds))
pdat_two_thirds$"tumor"[pdat_two_thirds$`colnames(df_two_thirds)` %in% colnames(dfHealthy_two_thirds)] <- 0 
pdat_two_thirds$"tumor"[pdat_two_thirds$`colnames(df_two_thirds)` %in% colnames(dfCancer_two_thirds)] <- 1 

pdat_all <- as.data.frame(colnames(df_all))
pdat_all$"tumor"[pdat_all$`colnames(df_all)` %in% colnames(dfHealthy_all)] <- 0 
pdat_all$"tumor"[pdat_all$`colnames(df_all)` %in% colnames(dfCancer_all)] <- 1 

pdat_unfil <- as.data.frame(colnames(df_unfil))
pdat_unfil$"tumor"[pdat_unfil$`colnames(df_unfil)` %in% colnames(dfHealthy_unfil)] <- 0 
pdat_unfil$"tumor"[pdat_unfil$`colnames(df_unfil)` %in% colnames(dfCancer_unfil)] <- 1 

write.table(df_two_thirds, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/exprs_three_quarters.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(pdat_two_thirds, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/pdat_three_quarters.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(fdat_two_thirds, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/fdat_three_quarters.txt", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(df_all, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/exprs_all.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(pdat_all, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/pdat_all.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(fdat_all, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/fdat_all.txt", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(df_unfil, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/exprs_unfil.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(pdat_unfil, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/pdat_unfil.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(fdat_unfil, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/fdat_unfil.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#write.table(geneIDS, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/gene_card_IDs.txt", row.names = F, col.names = F, quote = F, sep = "\t")

