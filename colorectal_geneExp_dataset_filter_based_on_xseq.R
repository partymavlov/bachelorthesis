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
row.names(originalDF) <- originalDF$`Hybridization REF`
originalDF$`Hybridization REF` <-NULL

notExpressed <- read.table("E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_not_expressed.txt", sep= "\t", header = F )
filteredDF <- subset(originalDF, !row.names(originalDF) %in% notExpressed$V1, drop = F)
row.names(filteredDF) <- subset(row.names(originalDF), !row.names(originalDF) %in% notExpressed$V1)

notExpCancer <- read.table("E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_expressed_only_healthy.txt", sep= "\t", header = F )
notExpHealthy <- read.table("E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_expressed_only_cancer.txt", sep= "\t", header = F )


dfCancer <- subset(filteredDF, select = grepl("^....-..-....-0", names(filteredDF)))
dfHealthy <- subset(filteredDF, select = grepl("^....-..-....-1", names(filteredDF)))
row.names(dfCancer) <- row.names(filteredDF)
row.names(dfHealthy) <- row.names(filteredDF)

dfCancer[row.names(dfCancer) %in% notExpCancer$V1, ] <- 0.0000
dfHealthy[row.names(dfHealthy) %in% notExpHealthy$V1, ] <- 0.0000

finalDF<- merge(dfHealthy, dfCancer, by=0)
row.names(finalDF) <- finalDF$Row.names
finalDF$Row.names <- NULL

write.table(format(finalDF, digits=4), file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_filtered_dataset.txt", sep = "\t", quote = F)
