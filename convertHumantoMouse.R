#sbiocLite('biomartr')
library(biomaRt)
library(biomartr)

human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes<- read.csv("/home/p/pavlovma/Bachelor/Bachelor/Ensembl_Genes_new.txt",header=F)$V1
genes

#orth =  getBM(c("mgi_id","human_ensembl_gene","human_orthology_type","human_percent_identity"),
#        filters="mgi_id",values =
#        genes,
#        mart = mouse)

#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes ,mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
genes<- data.frame(genes[,1], genes[,2])
genes
write.table(genes, file="/home/p/pavlovma/Bachelor/Bachelor/Reprogramming Dataset/Gene_names_ensmbl+cards_v2.txt", sep="\t",row.names=FALSE,col.names=FALSE, quote = FALSE)
