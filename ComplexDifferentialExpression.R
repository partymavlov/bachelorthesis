quantifiedColorectal <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/quantifiedColorectal05.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)

quantifiedBreast <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/quantifiedBreast05.txt", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)


differentiallyExpressedColorectal <- data.frame(ComplexID = quantifiedColorectal$ComplexID, ComplexName = quantifiedColorectal$ComplexName, Members = quantifiedColorectal$Members,
                             QuantifiedInColorectal = quantifiedColorectal$QuantifiedInColorectal, DifferentiallyExpressedInColorectal = character(nrow(quantifiedColorectal)), 
                             DifferentiallyExpressedUpRegulated = character(nrow(quantifiedColorectal)), DifferentiallyExpressedDownRegulated = character(nrow(quantifiedColorectal)),
                             PercentageDifferentiallyExpressedInColorectal = integer(nrow(quantifiedColorectal)), stringsAsFactors = FALSE)
differentiallyExpressedBreast <- data.frame(ComplexID = quantifiedBreast$ComplexID, ComplexName = quantifiedBreast$ComplexName, Members = quantifiedBreast$Members,
                                                QuantifiedInBreast = quantifiedBreast$QuantifiedInBreast, DifferentiallyExpressedInBreast = character(nrow(quantifiedBreast)),
                                                DifferentiallyExpressedUpRegulated = character(nrow(quantifiedBreast)), DifferentiallyExpressedDownRegulated = character(nrow(quantifiedBreast)),
                                                PercentageDifferentiallyExpressedInBreast = integer(nrow(quantifiedBreast)), stringsAsFactors = FALSE)

for(i in 1:nrow(quantifiedColorectal)){
  members <- strsplit(quantifiedColorectal$QuantifiedInColorectal[i], " ")[[1]]
  size <- length(members)
  differentiallyExpressed <- ""
  differentiallyExpressedSize <- character(0)
  upRegulated <- 0
  downRegulated <- 0

  for(j in 1:length(members)){
    if(members[j] %in% colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID){
      differentiallyExpressed <- paste(differentiallyExpressed, members[j], sep = " ")
      differentiallyExpressedSize <- c(differentiallyExpressedSize, members[j])
      if(colorectal_DE_protein_db$colorectal_DE_protein_db.log2FC[colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID == members[j]] == 'up'){
        upRegulated <- upRegulated + 1
      }
      if(colorectal_DE_protein_db$colorectal_DE_protein_db.log2FC[colorectal_DE_protein_db$colorectal_DE_protein_db.GENE.ID == members[j]] == 'down'){
        downRegulated <- downRegulated + 1
    }
    }
    }
  differentiallyExpressed<- substring(differentiallyExpressed, 2)
  differentiallyExpressedColorectal$DifferentiallyExpressedInColorectal[i] <- differentiallyExpressed
  differentiallyExpressedColorectal$PercentageDifferentiallyExpressedInColorectal[i] <- (length(differentiallyExpressedSize)/size)
  differentiallyExpressedColorectal$DifferentiallyExpressedUpRegulated[i] <- upRegulated
  differentiallyExpressedColorectal$DifferentiallyExpressedDownRegulated[i] <- downRegulated
}

for(i in 1:nrow(quantifiedBreast)){
  members <- strsplit(quantifiedBreast$QuantifiedInBreast[i], " ")[[1]]
  size <- length(members)
  differentiallyExpressed <- ""
  differentiallyExpressedSize <- character(0)
  upRegulated <- 0
  downRegulated <- 0
  
  for(j in 1:length(members)){
    if(members[j] %in% breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol){
      differentiallyExpressed <- paste(differentiallyExpressed, members[j], sep = " ")
      differentiallyExpressedSize <- c(differentiallyExpressedSize, members[j])
      if(breast_DE_protein_db$breast_DE_protein_db.Regulation[breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol== members[j]] == 'up'){
        upRegulated <- upRegulated + 1
      }
      if(breast_DE_protein_db$breast_DE_protein_db.Regulation[breast_DE_protein_db$breast_DE_protein_db.Gene_Symbol == members[j]] == 'down'){
        downRegulated <- downRegulated + 1
      }
    }
  }
  differentiallyExpressed<- substring(differentiallyExpressed, 2)
  differentiallyExpressedBreast$DifferentiallyExpressedInBreast[i] <- differentiallyExpressed
  differentiallyExpressedBreast$PercentageDifferentiallyExpressedInBreast[i] <- (length(differentiallyExpressedSize)/size)
  differentiallyExpressedBreast$DifferentiallyExpressedUpRegulated[i] <- upRegulated
  differentiallyExpressedBreast$DifferentiallyExpressedDownRegulated[i] <- downRegulated
}

write.table(differentiallyExpressedColorectal, file = "E:/Bachelor/Bachelor/CORUM complexes/DEA results/differentialExpressionOfComplexesColorectal05_withReg.txt", sep = "\t", quote = F, row.names = F)
write.table(differentiallyExpressedBreast, file = "E:/Bachelor/Bachelor/CORUM complexes/DEA results/differentialExpressionOfComplexesBreast05_withReg.txt", sep = "\t", quote = F, row.names = F)

significantlyExpressedColorectal <- differentiallyExpressedColorectal[differentiallyExpressedColorectal$PercentageDifferentiallyExpressedInColorectal > 0.2,]
significantlyExpressedColorectal <- significantlyExpressedColorectal[order(-significantlyExpressedColorectal$PercentageDifferentiallyExpressedInColorectal),]

significantlyExpressedBreast<- differentiallyExpressedBreast[differentiallyExpressedBreast$PercentageDifferentiallyExpressedInBreast > 0.2,]
significantlyExpressedBreast<- significantlyExpressedBreast[order(-significantlyExpressedBreast$PercentageDifferentiallyExpressedInBreast),]

write.table(significantlyExpressedColorectal, file = "E:/Bachelor/Bachelor/CORUM complexes/DEA results/differentiallyExpressedInColorectal05_withReg.txt", sep = "\t", quote = F, row.names = F)
write.table(significantlyExpressedBreast, file = "E:/Bachelor/Bachelor/CORUM complexes/DEA results/differentiallyExpressedInBreast05_withReg.txt", sep = "\t", quote = F, row.names = F)
