library(readr)
library(stringr)
input1 <- "E:/Bachelor/Bachelor/new_data/Gene DEA/BLCA_results.txt"
input2 <- "E:/Bachelor/Bachelor/new_data/Gene DEA/BRCA_results.txt"

#commandArgs(1)
suppressWarnings(suppressPackageStartupMessages(library(purrr)))
cancer_type1 <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input1))
cancer_type1 <-strsplit(cancer_type1, "\\.")[[1]]
cancer_type1 <- unlist(cancer_type1)
cancer_type1 <- cancer_type1[1]
cancer_type1 <-strsplit(cancer_type1, "\\_")[[1]]
cancer_type1 <- cancer_type1[1]

cancer_type2 <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input2))
cancer_type2 <-strsplit(cancer_type2, "\\.")[[1]]
cancer_type2 <- unlist(cancer_type2)
cancer_type2 <- cancer_type2[1]
cancer_type2 <-strsplit(cancer_type2, "\\_")[[1]]
cancer_type2 <- cancer_type2[1]
#originalDF <- read_delim("E:/Bachelor/Bachelor/new_data/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
#                                                                                                   "\t", escape_double = FALSE, trim_ws = TRUE)
message("Reading input...")
df1 <- suppressMessages(read_delim(input1, "\t", escape_double = FALSE, trim_ws = TRUE))
genes1 <- df1$mRNA
df2 <- suppressMessages(read_delim(input2, "\t", escape_double = FALSE, trim_ws = TRUE))
genes2 <- df2$mRNA
genesOverlap <- intersect(genes1, genes2)

exprs1<- suppressMessages(read_delim(paste("E:/Bachelor/Bachelor/new_data/", cancer_type1, ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE))
exprs2<- suppressMessages(read_delim(paste("E:/Bachelor/Bachelor/new_data/", cancer_type2, ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", sep = ""), "\t", escape_double = FALSE, trim_ws = TRUE))
geneNid1<- t(as.data.frame(strsplit(exprs1$`Hybridization REF`, "\\|")))
geneNid2<- t(as.data.frame(strsplit(exprs2$`Hybridization REF`, "\\|")))
exprs1$`Hybridization REF`<- sapply(strsplit(exprs1$`Hybridization REF`, "\\|"), `[`, 2)
exprs1<-subset(exprs1, exprs1$`Hybridization REF` %in% genesOverlap)
exprs2$`Hybridization REF`<- sapply(strsplit(exprs2$`Hybridization REF`, "\\|"), `[`, 2)
exprs2<-subset(exprs2, exprs2$`Hybridization REF` %in% genesOverlap)

pdat1 <- vector(mode="numeric", length=0)
for(i in 1: ncol(exprs1))
{
  pdat1<- c(pdat1, paste("cond0.rep", i, sep = ""))
}
for(i in 1: ncol(exprs2))
{
  pdat1<- c(pdat1, paste("cond1.rep", i, sep = ""))
}
pdat2 <- c(rep(0, ncol(exprs1)), rep(1, ncol(exprs2)))
library(reshape2)
pdat<- suppressMessages(melt(data.frame(pdat1,pdat2)))
pdat$variable <- NULL

exprs<- suppressMessages(plyr::join(exprs1, exprs2))
exprs$`Hybridization REF`<- NULL
fdat<- genesOverlap
exprs_path <- "E:/Bachelor/Bachelor/tmp/exprs.txt"
pdat_path <- "E:/Bachelor/Bachelor/tmp/pdat.txt"
fdat_path <- "E:/Bachelor/Bachelor/tmp/fdat.txt"
output_path <-  paste("E:/Bachelor/Bachelor/tmp/" , cancer_type1, "-", cancer_type2, "_DEAresults.txt" , sep = "")
degenes_path <- paste("E:/Bachelor/Bachelor/tmp/" , cancer_type1, "-", cancer_type2, "_DEgenes.txt" , sep = "")
message("Writing expression matrices for DEA...")
write.table(exprs,exprs_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(pdat,pdat_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(fdat,fdat_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
de.method <- "limma"
#out.file <- commandArgs()[10]
#library("limma")
message("Starting DEA...")
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser))))
message("Reading data ...")
eset <- read.eset(exprs, pdat, fdat, data.type= "ma")

#message("Removing genes with low read count ...")
#eset <- eset[rowSums(exprs(eset), na.rm=TRUE) > ncol(eset),]
#message("Normalizing values ...")
#eset <- normalize(eset, norm.method = "quantile", within = FALSE, data.type = "ma")

message("DE analysis ...")
eset <- de.ana(eset, de.method=de.method, padj.method="none")

de.tbl <- fData(eset)[,sapply(c("FC.COL","ADJP.COL"), config.ebrowser)]
de.tbl <- cbind(de.tbl, p.adjust(de.tbl[,2], method="BH"))
de.tbl <- cbind(featureNames(eset), de.tbl)
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL", "ADJ.PVAL")

write.table(de.tbl, file=output_path, row.names=FALSE, quote=FALSE, sep="\t")
message("DEA table written to ", output_path)
onlyde<- de.tbl[de.tbl$ADJ.PVAL < 0.05,]
write.table(onlyde, file=degenes_path, row.names=FALSE, quote=FALSE, sep="\t")
message("DE genes written to ", degenes_path)
#invisible(file.remove(exprs_path_two_thirds))
#invisible(file.remove(pdat_path_two_thirds))
#invisible(file.remove(fdat_path_two_thirds))
