#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
library(readr)
input <- "E:/Bachelor/Bachelor/new_data/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
                                                                                                            
             #commandArgs(1)
suppressWarnings(suppressPackageStartupMessages(library(purrr)))
cancer_type <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
cancer_type <-strsplit(cancer_type, "\\.")[[1]]
cancer_type <- unlist(cancer_type)
cancer_type <- cancer_type[1]
#originalDF <- read_delim("E:/Bachelor/Bachelor/new_data/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
#                                                                                                   "\t", escape_double = FALSE, trim_ws = TRUE)
message("Reading input...")
originalDF <- suppressMessages(read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE))
originalDF = originalDF[-1,]
library(stringr)
message("Filtering out unsufficiently quantified genes...")
to_del_two_thrids=which(rowSums(originalDF == "0.0000") > (2/3* (length(originalDF) - 1 )))
to_keep_all=which(rowSums(originalDF == "0.0000") == 0)
message("Dividing cancer and healthy probes...")
filteredDF_two_thirds <- subset(originalDF,  ! rownames(originalDF) %in% to_del_two_thrids )
#genes_two_thirds <- sapply(strsplit(filteredDF_two_thirds$`Hybridization REF`,split= "\\|"),'[',2)
genes_two_thirds <- as.data.frame(sapply(strsplit(filteredDF_two_thirds$`Hybridization REF`,split= "\\|"),'[',2))

filteredDF_all <- subset(originalDF,   rownames(originalDF) %in% to_keep_all )
genes_all <- sapply(strsplit(filteredDF_all$`Hybridization REF`,split= "\\|"),'[',2)
dfCancer_two_thirds <- plyr::join(genes_two_thirds,subset(filteredDF_two_thirds, select = grepl("^....-..-....-01", names(filteredDF_two_thirds))))
dfHealthy_two_thirds <- subset(filteredDF_two_thirds, select =  grepl("^....-..-....-11", names(filteredDF_two_thirds)))
dfCancer_all <- subset(filteredDF_all, select = grepl("^....-..-....-01", names(filteredDF_all)))
dfHealthy_all <- subset(filteredDF_all, select = grepl("^....-..-....-11", names(filteredDF_all)))
suppressWarnings(suppressPackageStartupMessages(library(plyr)))
message("Preparing expression matrices for DEA...")
exprs_two_thirds <- suppressMessages(plyr::join(dfCancer_two_thirds, dfHealthy_two_thirds))
pdat_two_thirds_1 <- vector(mode="numeric", length=0)
for(i in 1: ncol(dfCancer_two_thirds))
{
  pdat_two_thirds_1<- c(pdat_two_thirds_1, paste("cond0.rep", i, sep = ""))
}
for(i in 1: ncol(dfHealthy_two_thirds))
{
  pdat_two_thirds_1<- c(pdat_two_thirds_1, paste("cond1.rep", i, sep = ""))
}
pdat_two_thirds_2 <- c(rep(0, ncol(dfCancer_two_thirds)), rep(1, ncol(dfHealthy_two_thirds)))
library(reshape2)
pdat_two_thirds <- suppressMessages(melt(data.frame(pdat_two_thirds_1,pdat_two_thirds_2)))
pdat_two_thirds$variable <- NULL
fdat_two_thirds <- genes_two_thirds
exprs_path_two_thirds <- "E:/Bachelor/Bachelor/tmp/exprs.txt"
pdat_path_two_thirds <- "E:/Bachelor/Bachelor/tmp/pdat.txt"
fdat_path_two_thirds <- "E:/Bachelor/Bachelor/tmp/fdat.txt"
output_path <-  paste("E:/Bachelor/Bachelor/tmp/" , cancer_type, "_DEAresults_third.txt" , sep = "")
degenes_path <- paste("E:/Bachelor/Bachelor/tmp/" , cancer_type, "_DEgenes_third.txt" , sep = "")
message("Writing expression matrices for DEA...")
write.table(exprs_two_thirds,exprs_path_two_thirds,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(pdat_two_thirds,pdat_path_two_thirds,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(fdat_two_thirds,fdat_path_two_thirds,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
de.method <- "limma"
#out.file <- commandArgs()[10]
#library("limma")
message("Starting first DEA...")
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser))))
message("Reading data ...")
eset <- read.eset(exprs_path_two_thirds, pdat_path_two_thirds, fdat_path_two_thirds, data.type= "ma")

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

exprs_all <- suppressMessages(plyr::join(dfCancer_all, dfHealthy_all))
pdat_all_1 <- vector(mode="numeric", length=0)
for(i in 1: ncol(dfCancer_all))
{
  pdat_all_1<- c(pdat_all_1, paste("cond0.rep", i, sep = ""))
}
for(i in 1: ncol(dfHealthy_all))
{
  pdat_all_1<- c(pdat_all_1, paste("cond1.rep", i, sep = ""))
}
pdat_all_2 <- c(rep(0, ncol(dfCancer_all)), rep(1, ncol(dfHealthy_all)))
pdat_all <- suppressMessages(melt(data.frame(pdat_all_1,pdat_all_2)))
pdat_all$variable <- NULL
fdat_all <- genes_all
exprs_path_all <- "E:/Bachelor/Bachelor/tmp/exprs.txt"
pdat_path_all <- "E:/Bachelor/Bachelor/tmp/pdat.txt"
fdat_path_all <- "E:/Bachelor/Bachelor/tmp/fdat.txt"
output_path <-  paste("E:/Bachelor/Bachelor/tmp/" , cancer_type, "_DEAresults_all.txt" , sep = "")
degenes_path <- paste("E:/Bachelor/Bachelor/tmp/" , cancer_type, "_DEgenes_all.txt" , sep = "")
write.table(exprs_all,exprs_path_all,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(pdat_all,pdat_path_all,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
write.table(fdat_all,fdat_path_all,sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE)
de.method <- "limma"
#out.file <- commandArgs()[10]
#library("limma")
message("Starting second DEA...")
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser))))
message("Reading data ...")
eset <- suppressWarnings(read.eset(exprs_path_all, pdat_path_all, fdat_path_all, data.type= "ma"))

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
invisible(file.remove(exprs_path_all))
invisible(file.remove(pdat_path_all))
invisible(file.remove(fdat_path_all))


#dev.off()
#End of working part

# 
# 
# 
# 
# 
# exprs<-
# genes<-read.table("E:/Bachelor/Bachelor/11 Cancer Cell Lines/differential expression analysis/A549-GAMG/fdat.txt")
# pdat<-read.table("E:/Bachelor/Bachelor/11 Cancer Cell Lines/differential expression analysis/A549-GAMG/pdat.txt")
# design <- cbind(WT=c(1,1,1,0,0,0),MU=c(0,0,0,1,1,1))
# #data<-cbind(genes, exprs)
# #data
# #dim(data)
# #dim(design)
# fit<-lmFit(exprs, design)
# cont.matrix <- makeContrasts(MUvsWT=MU-WT, levels=design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# fit2 <- eBayes(fit2)
# topTable(fit2, adjust="BH")
# 
# 
# design <- cbind(WT=1,MUvsWT=c(0,0,0,1,1,1))
# fit <- lmFit(exprs, design)
# fit <- eBayes(fit)
# topTable(fit, coef="MUvsWT", adjust="BH")
# 
# pval_two_thrids = matrix(data=NA,nrow=nrow(dfCancer_two_thirds),ncol=11)
# pval_all = matrix(data=NA,nrow=length(datasetNormalized[,1]),ncol=11)
# fc = matrix(data=NA,nrow=length(datasetNormalized[,1]),ncol=11)
# for(i in 1:11)
# {
#   # preparing the design matrix by considering replicate samples as well as samples that were filtered as being outliers.
#   Treatment = rep(0,33)
#   # labeling which samples will be tested against.
#   ch=((i*3):(i*3+2))-2
#   Treatment[ch] = 1
#   # removing outlier samples from the design matrix.
#   design = model.matrix(~ -1+factor(Treatment[c(-3,-21,-30)]))
#   # renaming the design matrix.
#   colnames(design) = c("control","treatment")
#   fit = lmFit(datasetNormalized, design=design) # Fit the original matrix to the above design.
#   contrastsMatrix = makeContrasts("control-treatment",levels = design)
#   fit2 = contrasts.fit(fit, contrasts = contrastsMatrix) # Making the comparison.
#   fit2 = eBayes(fit2) # Moderating the t-tests by empirical Bayes smoothing method.
#   # extracting the results:
#   a = topTable(fit2, coef = "control-treatment",number=length(datasetNormalized[,1]),adjust.method="fdr")
#   pval[,i] = a[rownames(datasetNormalized),5]
#   fc[,i] = a[rownames(datasetNormalized),2]  
# }
# # renaming the row and column names for the limma results
# ensemblid=comps[,2]
# names(ensemblid)=comps[,3]
# celllinespval = pval
# rownames(celllinespval) = ensemblid[rownames(datasetNormalized)]
# colnames(celllinespval) = c("A549","GAMG",	"HEK293",	"HeLa",	"HepG2",	"Jurkat",	"K562",	"LnCap",	"MCF7",	"RKO",	"U2OS")
# celllinesfc = fc
# rownames(celllinesfc) = ensemblid[rownames(datasetNormalized)]
# colnames(celllinesfc) = c("A549","GAMG",	"HEK293",	"HeLa",	"HepG2",	"Jurkat",	"K562",	"LnCap",	"MCF7",	"RKO",	"U2OS")

