library(biomaRt)
library(biomartr)

human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes<- read.csv("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/Reprogramming Dataset/Gene_names_unique.txt",header=F)$V1
genes

#orth =  getBM(c("mgi_id","human_ensembl_gene","human_orthology_type","human_percent_identity"),
#        filters="mgi_id",values =
#        genes,
#        mart = mouse)

#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes ,mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
genes<- data.frame(genes[,1], genes[,2])
genes
write.table(genes, file="C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/Reprogramming Dataset/Gene_names_ensmbl+cards.txt", sep="\t",row.names=FALSE,col.names=FALSE, quote = FALSE)
