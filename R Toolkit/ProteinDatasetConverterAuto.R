#commandArgs(1) <- "E:/Bachelor/Bachelor/new_data/RPPA/BLCA.rppa.txt"
#commandArgs(2) <- "E:/Bachelor/Bachelor/new_data/RPPA/BRCA.rppa.txt"
args<-commandArgs(trailingOnly = T)
#args<- "E:/Bachelor/Bachelor/new_data/RPPA/BLCA.rppa.txt"
for(i in 1:length(args))
{
  input <- args[i]
  cancer_type <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
  cancer_type <-strsplit(cancer_type, "\\.")[[1]]
  cancer_type <- unlist(cancer_type)
  cancer_type <- cancer_type[1]
  dataset<-read.delim(input, sep="\t", row.names=1,quote="",comment.char="")
  dataset<-as.matrix(dataset)

  colnames(dataset) <- paste(cancer_type ,1:length(colnames(dataset)), sep = "_")
  if(i==1)
  {
    total <- data.frame(Date=as.Date(character()), File=character(), User=character(), stringsAsFactors=FALSE)
    total <- dataset
  }
  else
  {
    total<- merge(total, dataset, by=0)
  }
  
}
total$Row.names <- sapply(strsplit(total$Row.names,split= "\\|"),'[',1)
total <- total[complete.cases(total), ]
total<- subset(total, !duplicated(total$Row.names))
suppressWarnings(library(biomaRt))
suppressWarnings(library(biomartr))
human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
uniprots = getBM(attributes = c("hgnc_symbol", "uniprot_gn"), filters = "hgnc_symbol", values = total$Row.names ,mart = human)
library(plyr)
uniprots <- ddply(uniprots, .(hgnc_symbol), summarize, uniprot_gn =  gsub(", ", ";", toString(uniprot_gn)))
colnames(uniprots)[1] <- "Row.names"
total <- merge(uniprots, total, by= "Row.names")
total<- total[!total$uniprot_gn=="",]
colnames(total)[1] <- "GeneCard_ID"
colnames(total)[2] <- "Uniprot_ID"
write.table(total, file = "E:/Bachelor/Bachelor/new_data/RPPA/merged.txt", sep = "\t", row.names = F, quote = F)