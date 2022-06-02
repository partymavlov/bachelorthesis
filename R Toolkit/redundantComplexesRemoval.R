library(readr)
library(stringi)

quantified_complexes <- read_delim("E:/Bachelor/Bachelor/CORUM complexes/humanBigComplexes.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T)
sizes <- numeric(604)
for(i in 1:nrow(quantified_complexes)){
  members <- strsplit(quantified_complexes$`subunits(UniProt IDs)`[i], ";")[[1]]
  size <- length(members)
sizes[i] <- size
}
quantified_complexes <- cbind(quantified_complexes, sizes)
quantified_complexes <- quantified_complexes[order(-quantified_complexes$sizes),]

for(i in 1:(nrow(quantified_complexes)-1)){
  membersFirst <- strsplit(quantified_complexes$`subunits(UniProt IDs)`[i], ";")[[1]]
  for(j in (i+1):nrow(quantified_complexes)){
    membersSecond <- strsplit(quantified_complexes$`subunits(UniProt IDs)`[j], ";")[[1]]
   if(length(intersect(membersFirst, membersSecond)) > (length(membersFirst)/2)){
     quantified_complexes <- quantified_complexes[-j,]
  }
 }
}
quantified_complexes$sizes <- NULL
write.table(quantified_complexes, file = "E:/Bachelor/Bachelor/CORUM complexes/nonRedundantComplexes.txt", sep = "\t", quote = F, row.names = F)
