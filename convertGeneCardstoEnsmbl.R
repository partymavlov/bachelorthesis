library(biomaRt)
library(biomartr)

human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genes<- read.csv("E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/gene_card_IDs.txt",header=F)
#genes

#orth =  getBM(c("mgi_id","human_ensembl_gene","human_orthology_type","human_percent_identity"),
#        filters="mgi_id",values =
#        genes,
#        mart = mouse)

#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes = getBM(attributes = c("entrezgene", "ensembl_gene_id"), filters = "entrezgene", values = genes ,mart = human)
genes
genes<- data.frame(genes[,1], genes[,2])
genes
write.table(genes, file="E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/geneIDs_ensmbl+card.txt", sep="\t",row.names=FALSE,col.names=FALSE, quote = FALSE)

