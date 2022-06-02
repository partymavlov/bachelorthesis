library(readr)
results <- read_delim("E:/Bachelor/Bachelor/new_data/RPPA/Broad/converted/BRCA_results_full.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
deg <- read_delim("E:/Bachelor/Bachelor/new_data/Gene DEA/BRCA_results.txt", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)
complexes <- read_delim("E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/input-files/Protein-Complexes-UniprotIds.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
degs<- deg$mRNA
diffExed <- subset(results, results$`LIMMA p-value T1-C1` <= 0.05)
diffExed <- diffExed[, -c(5:8)]
downreg <- vector(mode="numeric", length=0)
upreg <- vector(mode="numeric", length=0)

for(i in 1:nrow(BRCA_results_full)){
    if(BRCA_results_full$`Fold-change T1-C1`[i] > 0){
    upreg<- cbind(BRCA_results_full$Analyte_ID[i], upreg)
  }
  else{
    downreg<- cbind(BRCA_results_full$Analyte_ID[i], downreg)
    }  
  
}
occurances <- as.data.frame(table(complexes$EnsemblGeneID))

diffExed <- cbind(diffExed, "Part of # Predefined Complexes")


test <- diffExed
test<- data.frame(stringsAsFactors=FALSE)
write.table(diffExed, file="E:/Bachelor/Bachelor/new_data/BRCA_results.txt", row.names=FALSE, quote=FALSE, sep="\t")

counts <- as.data.frame(table(diffExed$`Part of # Predefined Complexes`))
counts <- c(904, 172, 28, 8, 1, 0, 3, 3)
barplot(counts, names.arg=c("in 0 protein complexes (904)", "in 1 protein complexes (172)", "in 2 protein complexes (28)", "in 3 protein complexes (8)", "in 4 protein complexes (1)", "in 5 protein complexes (0)", "in 6 protein complexes (3)", "in 7 protein complexes (3)"), col = rainbow(8), las= 2, ylab = " # of proteins" , ylim = c(0, 1119))
par(mar=c(12,6,1,2))
