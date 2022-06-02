
if(length(commandArgs()) != 10) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <results.file> <differentiallyExpressed.file>")
stopifnot(length(commandArgs()) == 11)

message("Loading EnrichmentBrowser")
print(.libPaths())
suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser)))

exprs.file <- commandArgs()[6]
pdat.file <- commandArgs()[7]
fdat.file <- commandArgs()[8]
de.method <- commandArgs()[9]
out.file <- commandArgs()[10]
out.file1 <- commandArgs()[11]

#exprs.file <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/exprs.txt"
#pdat.file <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/pdat.txt"
#fdat.file <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/fdat.txt"
#de.method <- "limma"
#out.file <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/out_differentially expressed.txt"

message("Reading data ...")
eset <- read.eset(exprs.file, pdat.file, fdat.file, data.type= "ma")

#message("Removing genes with low read count ...")
#eset <- eset[rowSums(assay(eset), na.rm=TRUE) > ncol(eset),]
#message("Normalizing values ...")
#eset <- normalize(eset, norm.method = "quantile", within = FALSE, data.type = "ma")

message("DE analysis ...")
eset <- de.ana(eset, de.method=de.method, padj.method="none")

#de.tbl <- fData(eset)[,sapply(c("FC.COL","ADJP.COL"), config.ebrowser)]
de.tbl <- rowData(eset)
de.tbl <- cbind(de.tbl, p.adjust(de.tbl[,4], method="BH"))
de.tbl$ENTREZID <- NULL
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL", "limma.STAT","ADJ.PVAL")
dea.tbl <- de.tbl[abs(de.tbl$log2FC) > 1 & de.tbl$ADJ.PVAL <= 0.05,]
write.table(de.tbl, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
write.table(dea.tbl, file=out.file1, row.names=FALSE, quote=FALSE, sep="\t")
message("DE tables written to ", out.file)
