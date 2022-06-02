library(limma)
library(edgeR)
library(Biobase)
pdatFile <- read.delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/pdat.txt", sep = "\t", header = F, col.names = c("label", "condition"))
pdat1<- as.data.frame(pdatFile$label)
fdatFile <- read.delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/fdat.txt", sep = "\t", header = F, col.names = c("geneSymbol", "geneSymbol2"))
fdatFile<- as.data.frame(fdatFile$geneSymbol)
exprs <- read.delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/exprs.txt", sep = "\t", header = F)
colnames(exprs) <- pdat$label
row.names(exprs) <- fdat$`fdat$geneSymbol`
pdatFile[pdatFile == "0"] <- "NT"
pdatFile[pdatFile == "1"] <- "T"
counts <- read.delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/colorectal_90tumorVS30healthy.txt", sep = "\t", row.names = 1)
dge <- DGEList(counts=counts)
u <- unique(pdatFile$condition)
f <- factor(pdatFile$condition, levels=u)
design <- model.matrix(~0+f)
colnames(design) <- u
#keep <- filterByExpr(dge, design)
#dge <- dge[keep,,keep.lib.sizes=FALSE]
#logCPM <- cpm(dge, log=T)
pdat<- AnnotatedDataFrame(pdat1)
fdat<- AnnotatedDataFrame(fdatFile)
protodat <- AnnotatedDataFrame(pdat1)
sampleNames(protodat) <- pdatFile$label
featureNames(fdat) <- fdatFile$`fdatFile$geneSymbol`
sampleNames(pdat) <- pdatFile$label
exprs <- as.matrix(exprs)
dge<- ExpressionSet(assayData = exprs, phenoData = pdat, featureData = fdat, protocolData = protodat)
dge
a<- fData(dge)
de.tbl <- fData(dge)[,sapply(c("FC.COL","ADJP.COL"), config.ebrowser)]
exprs(dge)
fit <- lmFit(dge, design)
fit <- eBayes(fit, trend=TRUE)
out<- topTable(fit, coef=ncol(design), sort="none",n=Inf, lfc = 1.0, p.value = 0.05)
write.table(out, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/DEA/results4_filtered.txt", sep = "\t")
