input <- "E:/Bachelor/Bachelor/new_data/RPPA/Broad/BRCA.rppa.txt"

df <- suppressMessages(read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE))
library(data.table)
df<- subset(df, unique(df$Composite.Element.REF))
#df$Composite.Element.REF<- sapply(strsplit(df$Composite.Element.REF, "\\|"), `[`, 1)
fdat<- t(as.data.frame(strsplit(df$Composite.Element.REF, "\\|")))
dfHealthy <- subset(df, select = grepl("^....-..-....-1", names(df)))

#patient<-substr(colnames(dfHealthy), 0, 12)
dfCancer <- subset(df, select = grepl("^....-..-....-0", names(df)))
dfCancer_paired_names <- vector(mode="numeric", length=0)
for(i in 1: ncol(dfHealthy))
{
  patient<-substr(colnames(dfHealthy)[i], 0, 12)
  canc_barcode <- names(subset(dfCancer, select = grepl(paste(patient, "-01", sep = ""), names(dfCancer))))
  #message(canc_barcode)
  if(!is.null(canc_barcode))
  {
    dfCancer_paired_names<- c(dfCancer_paired_names, canc_barcode)
  }
}
dfCancer<- subset(dfCancer, select = dfCancer_paired_names)

exprs<-merge(dfHealthy, dfCancer, by = 'row.names')
exprs$Row.names<-NULL

pdat1 <- vector(mode="numeric", length=0)
for(i in 1: ncol(dfHealthy))
{
  pdat1<- c(pdat1, paste("cond0.rep", i, sep = ""))
}
for(i in 1: ncol(dfCancer))
{
  pdat1<- c(pdat1, paste("cond1.rep", i, sep = ""))
}
pdat2 <- c(rep(0, ncol(dfHealthy)), rep(1, ncol(dfCancer)))
library(reshape2)
pdat<- suppressMessages(melt(data.frame(pdat1,pdat2)))
pdat$variable <- NULL

exprs_path <- "E:/Bachelor/Bachelor/tmp/exprs.txt"
pdat_path <- "E:/Bachelor/Bachelor/tmp/pdat.txt"
fdat_path <- "E:/Bachelor/Bachelor/tmp/fdat.txt"
output_path <-  paste("E:/Bachelor/Bachelor/tmp/" , "BRCA.rppa", "_DEAresults.txt" , sep = "")
degenes_path <- paste("E:/Bachelor/Bachelor/tmp/" , "BRCA.rppa", "_DEproteins.txt" , sep = "")
message("Writing expression matrices for DEA...")
write.table(exprs,exprs_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(pdat,pdat_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(fdat,fdat_path,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
de.method <- "limma"

message("Starting DEA...")
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser))))
message("Reading data ...")
eset <- read.eset(exprs_path, pdat_path, fdat_path, data.type= "ma")

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
