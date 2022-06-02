library(readr)
dep <- read_delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/ProtExp/DEA/differentially_expressed.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
deg <- read_delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/GeneExp/DEA/xseq/xseq_differentially_expressed.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)
complexes <- read_delim("E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/input-files/Protein-Complexes-UniprotIds.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
ensToCard <-read.table("E:/Bachelor/Bachelor/CPTAC/Colorectal/geneIDs_ensmbl+card.txt", 
                       "\t", header = F, sep = "\t")
degs<- deg$GENE.ID
diffExed <- dep[, -c(3:4)]
diffExed <- cbind(diffExed, 'Regulation')
downreg <- vector(mode="numeric", length=0)
upreg <- vector(mode="numeric", length=0)

for(i in 1:nrow(diffExed)){
  if(diffExed$log2FC[i] > 0){
    diffExed$"Regulation"[i] <-  "up"
  }
  else{
    diffExed$"Regulation"[i] <-  "down"
  }  
  
}
diffExed <- diffExed[, -c(2,4)]
diffExed<- cbind(diffExed, 'Transcriptome level inclusion')
for(i in 1:nrow(diffExed)){
  if(diffExed$GENE.ID[i] %in% degs){
    diffExed$"Transcriptome level inclusion"[i] <-  '+'
  }
  else{
    diffExed$"Transcriptome level inclusion"[i] <-  '-'
  }  
  
}
diffExed <- diffExed[, -4]
occurances <- as.data.frame(table(complexes$EnsemblGeneID))

diffExed <- cbind(diffExed, "Part of # Predefined Complexes")

for(i in 1:nrow(diffExed)){
  if(diffExed$GENE.ID[i] %in% occurances$Var1){
    diffExed$"Part of # Predefined Complexes"[i] <-  occurances$Freq[occurances$Var1 == diffExed$GENE.ID[i]]
  }
  else{
    diffExed$"Part of # Predefined Complexes"[i] <- 0
  }  
  
}
diffExed <- diffExed[, -5]


colnames(ensToCard)[2] <- 'GENE.ID'

diffExed <- merge(diffExed, ensToCard, by = 'GENE.ID')
diffExed <- cbind(diffExed$V1, diffExed[, c(1:5)])
colnames(diffExed)[1] <- 'GeneSymbol'
colnames(diffExed)[2] <- 'EnsemblID'
colnames(diffExed)[3] <- 'LIMMA p-value'
diffExed <- cbind(diffExed[,1:2], diffExed[,4], diffExed[,3], diffExed[,5:6])
colnames(diffExed)[3] <- 'Regulation'
colnames(diffExed)[4] <- 'LIMMA p-value'

write.table(diffExed, file="E:/Bachelor/Bachelor/CPTAC/Colorectal/COADREAD_results.txt", row.names=FALSE, quote=FALSE, sep="\t")

for(i in 1:nrow(xseq_differentially_expressed)){
  if(xseq_differentially_expressed$log2FC[i] > 0){
    xseq_differentially_expressed$"Regulation"[i] <-  "up"
  }
  else{
    xseq_differentially_expressed$"Regulation"[i] <-  "down"
  }  
  
}
