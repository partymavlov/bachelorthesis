#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
library(readr)
originalDF <- read_delim("E:/Bachelor/Bachelor/new_data/BLCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", 
                                                                                                   "\t", escape_double = FALSE, trim_ws = TRUE)
library(stringr)
dfCancer <- subset(originalDF, select = grepl("^....-..-....-01", names(originalDF)))
dfHealthy <- subset(originalDF, select = grepl("^....-..-....-11", names(originalDF)))
for(i in 1:length(originalDF))
{
  if(str_count(originalDF[i], "0.0000") > 2/3* length(originalDF[i]))
    originalDF <- originalDF[-originalDF[i],]
}


