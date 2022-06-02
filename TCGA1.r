# You can define a list of samples to query and download providing relative TCGA barcodes.
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#source("https://bioconductor.org/biocLite.R")
#biocLite("TCGAbiolinks")
suppressWarnings(suppressPackageStartupMessages(library("TCGAbiolinks")))
suppressWarnings(suppressPackageStartupMessages(library(readr)))
suppressWarnings(suppressPackageStartupMessages(library(SummarizedExperiment)))
#install.packages("devtools")
#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")

#input <- commandArgs(1)
input <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/COADREAD_genes_normalized__data.txt"
#commandArgs(1)
suppressWarnings(suppressPackageStartupMessages(library(purrr)))
cancer_type <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input))
cancer_type <-strsplit(cancer_type, "\\.")[[1]]
cancer_type <- unlist(cancer_type)
cancer_type <- cancer_type[1]
cancer_type <-strsplit(cancer_type, "\\_")[[1]]
cancer_type <- cancer_type[1]
#originalDF <- read_delim("E:/Bachelor/Bachelor/new_data/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
#                                                                                                   "\t", escape_double = FALSE, trim_ws = TRUE)
message("Reading input...")
originalDF <- suppressMessages(read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE))
originalDF <- originalDF[-1,]
originalDF <- originalDF[-c(1:29),]
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
dfCancer_all <- subset(filteredDF_all, select = grepl("^....-..-....-0", names(filteredDF_all)))
dfHealthy_all <- subset(filteredDF_all, select = grepl("^....-..-....-1", names(filteredDF_all)))
dfCancer_paired_names <- vector(mode="numeric", length=0)
for(i in 1: ncol(dfHealthy_all))
{
  patient<-substr(colnames(dfHealthy_all)[i], 0, 12)
  canc_barcode <- names(subset(dfCancer_all, select = grepl(patient, names(dfCancer_all))))
  if(!is.null(canc_barcode))
  {
    dfCancer_paired_names<- c(dfCancer_paired_names, canc_barcode)
  }
}
dfHealthy_paired_names <- colnames(dfHealthy_all)
for(i in dfCancer_paired_names){
  patient<-substr(i, 0, 12)
}
message("Generating query...")
query <- GDCquery(project = paste("TCGA-", cancer_type, sep = ""), 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "normalized_results",
                  barcode = c(dfCancer_paired_names, dfHealthy_paired_names), 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
message("Downloading query data...")
GDCdownload(query, method = "api")

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE)
#BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)
#BRCARnaseq_CorOutliers
#TCGAvisualize_BarPlot(BRCARnaseq_CorOutliers, Subtype = cancer_type, filename = paste("E:/Bachelor/Bachelor/tmp/", canc_barcode, "_sample_cluster_plot.png", sep = ""))
message("Normalizing and filtering expression data...")
#dataNorm <- TCGAanalyze_Normalization(tabDF = BRCAMatrix, geneInfo =  geneInfo)

# quantile filter of genes
#dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                  method = "quantile", 
#                                  qnt.cut =  0.25)
dataFilt <- TCGAanalyze_Filtering(tabDF = BRCAMatrix, method = "quantile", qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))
message("Performing DEA...")
# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])
message("Writing results...")
write.table(dataDEGsFiltLevel, file = paste("E:/Bachelor/Bachelor/CPTAC/Colorectal/Gene expression/", cancer_type, "_results.txt", sep = ""),sep = "\t")
      