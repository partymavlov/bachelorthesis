library(readr)
library(stringi)


colorectal_protein_db <- read_delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/ProtExp/DEA/DEA_results.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE)

colorectal_DE_protein_db <- read_delim("E:/Bachelor/Bachelor/CPTAC/Colorectal/ProtExp/DEA/differentially_expressed.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)

breast_protein_db <- read_delim("E:/Bachelor/Bachelor/TCGA/RPPA/Broad/BRCA.rppa.txt", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)

breast_DE_protein_db <- read_delim("E:/Bachelor/Bachelor/TCGA/BRCA_results.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)

colorectal_IDs <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/colorectal_IDs.txt", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)

breast_IDs <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/breast_IDs.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)

colorectal_IDs<- colorectal_IDs[!duplicated(colorectal_IDs$EnsemblID),]
breast_IDs<- breast_IDs[!duplicated(breast_IDs$GeneCardID),]

colorectal_protein_db <- colorectal_protein_db$GENE.ID
for(i in 1:length(colorectal_protein_db)){
  if(colorectal_protein_db[i] %in% colorectal_IDs$EnsemblID){
    colorectal_protein_db[i] <- colorectal_IDs$UniProtID[colorectal_protein_db[i] == colorectal_IDs$EnsemblID]
  }
  }
colorectal_DE_protein_db <- data.frame(colorectal_DE_protein_db$GENE.ID, colorectal_DE_protein_db$log2FC, stringsAsFactors = F)
for(i in 1:nrow(colorectal_DE_protein_db)){
  if(colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID[i] %in% colorectal_IDs$EnsemblID){
    colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID[i] <- colorectal_IDs$UniProtID[colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID[i] == colorectal_IDs$EnsemblID]
  }
  if(colorectal_DE_protein_db$colorectal_DE_protein_db.log2FC[i] >0){
    colorectal_DE_protein_db$colorectal_DE_protein_db.log2FC[i] <- "up"
  }
  else{
    colorectal_DE_protein_db$colorectal_DE_protein_db.log2FC[i] <- "down"
  }
}
breast_protein_db<- breast_protein_db$Composite.Element.REF 
breast_protein_db <- sapply(strsplit(breast_protein_db, "\\|"), `[`, 1)
for(i in 1:length(breast_protein_db)){
  if(breast_protein_db[i] %in% breast_IDs$GeneCardID){
    breast_protein_db[i] <- breast_IDs$UniProtID[breast_protein_db[i] == breast_IDs$GeneCardID]
  }
}
breast_DE_protein_db<-data.frame(breast_DE_protein_db$Gene_Symbol, breast_DE_protein_db$Regulation, stringsAsFactors = F)
for(i in 1:length(breast_DE_protein_db)){
  if(breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol[i] %in% breast_IDs$GeneCardID){
    breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol[i] <- breast_IDs$UniProtID[breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol[i] == breast_IDs$GeneCardID]
  }
}

write.table(breast_protein_db, file = "E:/Bachelor/Bachelor/CORUM complexes/breast.txt", sep = "\t", quote = F, row.names = F)
write.table(breast_DE_protein_db, file = "E:/Bachelor/Bachelor/CORUM complexes/breastDE.txt", sep = "\t", quote = F, row.names = F)
write.table(colorectal_protein_db, file = "E:/Bachelor/Bachelor/CORUM complexes/colorectal.txt", sep = "\t", quote = F, row.names = F)
write.table(colorectal_DE_protein_db, file = "E:/Bachelor/Bachelor/CORUM complexes/colorectalDE.txt", sep = "\t", quote = F, row.names = F)

prot_complexes <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/nonRedundantComplexes.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)

#prot_comps <- prot_comps[prot_comps$Organism == "Human",]
#prot_comps_sufficient_members <- prot_comps[lengths(strsplit(prot_comps$`subunits(UniProt IDs)`, ";")) > 4,]
#write.table(prot_comps, file = "E:/Bachelor/Bachelor/CORUM complexes/humanComplexes.txt", sep = "\t", quote = F, row.names = F)

complex_scores <- data.frame(ComplexID = prot_complexes$ComplexID, ComplexName = prot_complexes$ComplexName, Members = stri_join_list(strsplit(prot_complexes$`subunits(UniProt IDs)`, ";"), sep = " "), QuantifiedInColorectal = character(nrow(prot_complexes)), PercentageQuantifiedInColorectal = integer(nrow(prot_complexes)), QuantifiedInBreast = character(nrow(prot_complexes)), PercentageQuantifiedInBreast = integer(nrow(prot_complexes)), stringsAsFactors = FALSE)

for(i in 1:nrow(complex_scores)){
  members <- strsplit(prot_complexes$`subunits(UniProt IDs)`[i], ";")[[1]]
  size <- length(members)
  quantifiedColorectal <- ""
  quantifiedsizeColorectal <- character(0)
  quantifiedBreast <- ""
  quantifiedsizeBreast <- character(0)
  for(j in 1:length(members)){
  if(members[j] %in% colorectal_protein_db){
    quantifiedColorectal <- paste(quantifiedColorectal, members[j], sep = " ")
    quantifiedsizeColorectal <- c(quantifiedsizeColorectal, members[j])
    }
  if(members[j] %in% breast_protein_db){
    quantifiedBreast <- paste(quantifiedBreast, members[j], sep = " ")
    quantifiedsizeBreast <- c(quantifiedsizeBreast, members[j])
    }
  }
  quantifiedColorectal<- substring(quantifiedColorectal, 2)
  complex_scores$QuantifiedInColorectal[i] <- quantifiedColorectal
  complex_scores$PercentageQuantifiedInColorectal[i] <- (length(quantifiedsizeColorectal)/size)
  quantifiedBreast<- substring(quantifiedBreast, 2)
  complex_scores$QuantifiedInBreast[i] <- quantifiedBreast
  complex_scores$PercentageQuantifiedInBreast[i] <- (length(quantifiedsizeBreast)/size)
}

write.table(complex_scores, file = "E:/Bachelor/Bachelor/CORUM complexes/complexQuantification.txt", sep = "\t", quote = F, row.names = F)

quantifiedColorectal025 <- complex_scores[complex_scores$PercentageQuantifiedInColorectal > 0.25,]
quantifiedColorectal0375 <- complex_scores[complex_scores$PercentageQuantifiedInColorectal > 0.375,]
quantifiedColorectal05 <- complex_scores[complex_scores$PercentageQuantifiedInColorectal > 0.5,]
quantifiedBreast025 <- complex_scores[complex_scores$PercentageQuantifiedInBreast > 0.25,]
quantifiedBreast0375 <- complex_scores[complex_scores$PercentageQuantifiedInBreast > 0.375,]
quantifiedBreast05 <- complex_scores[complex_scores$PercentageQuantifiedInBreast > 0.5,]

write.table(quantifiedColorectal025, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedColorectal025.txt", sep = "\t", quote = F, row.names = F)
write.table(quantifiedColorectal0375, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedColorectal0375.txt", sep = "\t", quote = F, row.names = F)
write.table(quantifiedColorectal05, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedColorectal05.txt", sep = "\t", quote = F, row.names = F)
write.table(quantifiedBreast025, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedBreast025.txt", sep = "\t", quote = F, row.names = F)
write.table(quantifiedBreast0375, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedBreast0375.txt", sep = "\t", quote = F, row.names = F)
write.table(quantifiedBreast05, file = "E:/Bachelor/Bachelor/CORUM complexes/quantifiedBreast05.txt", sep = "\t", quote = F, row.names = F)

