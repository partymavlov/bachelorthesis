source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
biocLite("affy")  # Methods for Affymetrix Oligonucleotide Arrays
biocLite("hgu133a.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
biocLite("hgu133acdf")
library(GEOquery)
setwd("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/11 Cancer Cell Lines/expression")
#getGEOSuppFiles("GSM180457")


cels = list.files("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/11 Cancer Cell Lines/expression/a549", pattern = "cel", ignore.case = TRUE)
cels

library(affy)
library(hgu133a.db)
library(hgu133acdf)

setwd("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/11 Cancer Cell Lines/expression/a549")
raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
data.rma.norm = rma(raw.data)
rma = exprs(data.rma.norm)
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")
tt = cbind(row.names(rma), rma)
colnames(tt) = c("ProbID", sub(".cel", "", colnames(rma), ignore.case = TRUE))
rownames(tt) = NULL
tt
annot <- data.frame(ENSEMBL=sapply(contents(hgu133aENSEMBL), paste, collapse=", "), ENTREZID=sapply(contents(hgu133aENTREZID), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
annot<- cbind(rownames(annot), annot)
rownames(annot) <- NULL
colnames(annot) <- c("Prob" , "ENSEMBL", "ENTREZID" , "DESC")
head(annot)
comb = merge(annot, tt, by.x="Prob", by.y= "ProbID")
head(comb)
write.table(comb, file = "exprs.txt", quote = FALSE, sep = "\t", row.names = FALSE)







# If multiple probe sets corresponded to the same gene, then the expression
# values of these probe sets were averaged.
comb2 <- subset(comb, select = -c(ProbeSetID))
comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)
out <- aggregate(. ~ EntrezGene, data = comb2, mean)

# Format values to 5 decimal places
out = format(out, digits = 5)
out[1:5, 1:5]
write.table(out, file = "GSE1297.RMA.txt", quote = FALSE, sep = "\t", row.names = FALSE)