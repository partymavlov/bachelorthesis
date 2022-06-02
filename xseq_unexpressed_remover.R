library(xseq)
library(readr)
#data(mut, expr, cna.call, cna.logr, net)

input <- "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/COADREAD_genes_normalized__data.txt"
originalDF <- suppressMessages(read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE))
originalDF <- originalDF[-1,]
originalDF <- originalDF[-c(1:29),]
exprs <- originalDF
row.names(exprs) <- exprs$`Hybridization REF`
exprs$`Hybridization REF` <-NULL

dfCancer_unfil <- subset(exprs, select = grepl("^....-..-....-0", names(exprs)))
dfHealthy_unfil <- subset(exprs, select = grepl("^....-..-....-1", names(exprs)))
row.names(dfCancer_unfil) <- row.names(exprs)
row.names(dfHealthy_unfil) <- row.names(exprs)

exprsCancer <- t(dfCancer_unfil)
exprsHealthy <- t(dfHealthy_unfil)

data(exprsHealthy)

ExprH <- log2(as.matrix.data.frame(exprsHealthy))
ExprC <- as.data.frame(log2(as.numeric(exprsCancer[,1:20502])))
weight = EstimateExpression(ExprC)
EstimateExpression(ExprC, show.plot = FALSE, loglik = TRUE,
                   xlab = "Expression", ylab = "Density")


id = weight[weight >= 0.8] 
id = id & !is.na(id)

##

setwd("/home/users/e.kataka/CancerPairedData")
patient_sigs = list.files(path = "../CancerPairedData/", pattern = "*siggenes",full.names = TRUE)
list_of_profiles <- lapply(patient_sigs, read.table, header = F)
merged_patient_sigs = Reduce(function(...) merge(..., all=T), list_of_profiles) # 7121



LIVERPROFILES=list.files(path = "~/CancerPairedData/Allother11cancersfilt", pattern = "*-LIHC$",full.names = TRUE,recursive = T)

lihccancer <- exprs[-grep("-11A-", exprs, fixed=T)]
lihchealthy <- exprs[-grep("-01A-", exprs, fixed=T)]



dataMergeC <- data.frame("isoform_id"= character(0))
for(f in lihccancer){ 
  ReadInMergeC <- read.table(file=f, header=T, na.strings="NULL")
  colnames(ReadInMergeC)[2] <- substr(f, 61, 76)
  dataMergeC <- merge(dataMergeC, ReadInMergeC, by="isoform_id",all=T)
}

dataMergeC 
rownames(dataMergeC) <- dataMergeC[,1]
dataMergeC[,1] <- NULL




library(reshape2)
lihccancerset <- melt(dataMergeC)
lihcexpr <- dcast(lihccancerset, variable~isoform_id , value.var='value', fill=0)
rownames(lihcexpr) <- lihcexpr[,1]
lihcexpr[,1] <- NULL
lihcexpr[1:5,1:5]


exprsHealthy <- exprsHealthy[apply(exprsHealthy, 1, function(x) !all(x==0)),]

rm(lihcexpr2)


library(xseq)
#exprsHealthy[1:51, 1:20502] <- sapply(exprsHealthy[1:51, 1:20502], as.numeric)
exprsHealthy <- apply(exprsHealthy, 2, as.numeric)
exprsCancer <- apply(exprsCancer, 2, as.numeric)

exprsHealthy <- exprsHealthy[1:51,1:20502] + 1
exprsHealthy[1:51,1:20502]  <- lapply(exprsHealthy[1:51,1:20502], log2)

exprsCancer <- exprsCancer[1:382,1:20502] + 1
exprsCancer[1:382,1:20502]  <- lapply(exprsCancer[1:382,1:20502], log2)

exprsHealthy <- apply(exprsHealthy, 2, log2)

exprsCancer <- apply(exprsCancer, 2, log2)


weightlihc = EstimateExpression(exprsCancer)

weightlihc <- EstimateExpression(exprsCancer, show.plot = T, loglik = TRUE,xlab = "Expression", ylab = "Density")

write.table(weightlihc, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/xseq_results_cancer.txt", sep = "\t", quote = F)

weightsCancer <- as.data.frame(weightlihc)
weightsHealthy <- read.table(file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_results_healthy.txt", sep = "\t")
weightsCancer <- subset(weightsCancer, weightsCancer$weightlihc >= 0.8)
weightsHealthy <- subset(weightsHealthy, weightsHealthy$x >= 0.8)
weightsOnlyC <- subset(weightsCancer, !(row.names(weightsCancer) %in% row.names(weightsHealthy)))
weightsOnlyH <- subset(weightsHealthy, !(row.names(weightsHealthy) %in% row.names(weightsCancer)))
weightsInNone <- as.data.frame(subset(row.names(exprs), !(row.names(exprs) %in% row.names(weightsHealthy))))
weightsInNone <- subset(weightsInNone, !(weightsInNone$`subset(row.names(exprs), !(row.names(exprs) %in% row.names(weightsHealthy)))` %in% row.names(weightsCancer)))
write.table(weightsInNone, file = "E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq_not_expressed.txt", sep = "\t", quote = F, row.names = F)

dataMergeH <- data.frame("isoform_id"= character(0))
for(f in lihchealthy){ 
  ReadInMergeH <- read.csv(file=f, header=T, na.strings="NULL")
  dataMergeH <- merge(dataMergeH, ReadInMergeH, 
                      by=isoform_id, all=T)
}


substr(f, 61, 76)


########

lihc_h <- lihc_h[-1,]
lihc_h$Hybridization.REF <- sub("^(\\?\\|).*", "", lihc_h$Hybridization.REF)
lihc_h <- lihc_h[-c(1:29), ]
colnames(lihc_h)[1] <- "gene_id"
lihccancerset <- melt(lihc_h, id.vars = "gene_id")
rm(lihccancerset)
#rm(lihcsets)
lihcexprG <- dcast(lihccancerset, variable~gene_id, value.var='value', fill=0)
rownames(lihcexprG) <- lihcexprG[,1]
lihcexprG[,1] <- NULL
lihcexprG[1:5,1:5]


#lihcexpr2 <- lihcexpr[apply(lihcexpr, 1, function(x) !all(x==0)),]

rm(lihcexpr2)


library(xseq)


#lihcexprG <-  as.numeric(gsub(",", ".", lihcexprG[,c(1:20502)]))
#USE READ.CSV to avoid errors on as.num
lihcexprG[,c(1:20531)] <- sapply(lihcexprG[,c(1:20531)], as.numeric)

lihcexprG[1:5,1:5]

lihcexprG <- lihcexprG + 1

lihcexprG[,c(1:20531)]  <- lapply(lihcexprG[,c(1:20531)], log2)

weightlihc = EstimateExpression(lihcexpr)
EstimateExpression(lihcexprG, show.plot = T, loglik = TRUE,xlab = "Expression", ylab = "Density")



weightlihc = EstimateExpression(expr, show.plot = T, loglik = TRUE,
                                xlab = "Expression", ylab = "Density")

##
lihc_C[1:5,1:5]

lihc_C <- lihc_C[-1,]
colnames(lihc_C)[1] <- "gene_id"
lihccancerset <- melt(lihc_C, id.vars = "gene_id")
rm(lihccancerset)
#rm(lihcsets)
lihcexprG <- dcast(lihccancerset, variable~gene_id, value.var='value', fill=0)
rownames(lihcexprG) <- lihcexprG[,1]
lihcexprG[,1] <- NULL
lihcexprG[1:5,1:5]


#lihcexpr2 <- lihcexpr[apply(lihcexpr, 1, function(x) !all(x==0)),]

rm(lihcexpr2)


library(xseq)


#lihcexprG <-  as.numeric(gsub(",", ".", lihcexprG[,c(1:20502)]))
#USE READ.CSV to avoid errors on as.num
lihcexprG[,c(1:20531)] <- sapply(lihcexprG[,c(1:20531)], as.numeric)

lihcexprG[1:5,1:5]

lihcexprG <- lihcexprG + 1

lihcexprG[,c(1:20531)]  <- lapply(lihcexprG[,c(1:20531)], log2)

weightlihc = EstimateExpression(lihcexpr)
EstimateExpression(lihcexprG, show.plot = T, loglik = TRUE,xlab = "Expression", ylab = "Density")