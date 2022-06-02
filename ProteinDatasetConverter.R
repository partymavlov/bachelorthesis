##  This source code file is part of the analysis of variable protein complexes (VariableComplexes).
##  Copyright (C) 2016  Murat Iskar, Alessandro Ori
## 
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
## 
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##  Please check the publication "Spatiotemporal variation of mammalian protein complex stoichiometries", XXX
##  and the supplementary website: www.bork.embl.de/variable_complexes
##
## preprocessing_id_matching.R
## v1
## Murat Iskar <muratiskar at gmail.com>
## Alessandro Ori <alessandro.ori at embl.de> 
##  20.1.2016
## 
## 
## input files: 
## Geiger_etal_Mol_Cell_Proteomics_2012_11_cell_lines_proteomic_dataset.csv
## Protein-Complexes-UniprotIds.txt
##
## output files: 
## 11-cell-line-dataset-only-Protein-complexes.txt
## matched-protein-complexes-11-cell-line-dataset.txt
## Outlier-detection-correlation-of-samples.pdf
##

options(max.print=200)
options(stringsAsFactors=FALSE)
#install.packages("gplots")
# to use heatmap.2 function
library("gplots")


########################################################
########################################################
# INPUTS
########################################################
#grep("^TCGA\\...\\.....\\.1", colnames(dataset1))
# proteomic dataset, log2 normalized, rownames = protein ids ";" separated



input1 <- "E:/Bachelor/Bachelor/new_data/RPPA/BLCA.rppa.txt"
cancer_type1 <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input1))
cancer_type1 <-strsplit(cancer_type1, "\\.")[[1]]
cancer_type1 <- unlist(cancer_type1)
cancer_type1 <- cancer_type1[1]
dataset1<-read.delim(input1, sep="\t", row.names=1,quote="",comment.char="")
dataset1<-as.matrix(dataset1)
colnames(dataset1) <- paste(cancer_type1 ,1:length(colnames(dataset1)), sep = "_")
input2<- "E:/Bachelor/Bachelor/new_data/RPPA/BRCA.rppa.txt"
cancer_type2 <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input2))
cancer_type2 <-strsplit(cancer_type2, "\\.")[[1]]
cancer_type2 <- unlist(cancer_type2)
cancer_type2 <- cancer_type2[1]
dataset2<-read.delim(input2, sep="\t", row.names=1,quote="",comment.char="")
dataset2<-as.matrix(dataset2)
colnames(dataset2) <- paste(cancer_type2 ,1:length(colnames(dataset2)), sep = "_")
total<- merge(dataset1, dataset2, by=0)
total$Row.names <- sapply(strsplit(total$Row.names,split= "\\|"),'[',1)
total <- total[complete.cases(total), ]
total<- subset(total, !duplicated(total$Row.names))
#genecards <- as.data.frame(sapply(strsplit(total$Row.names,split= "\\|"),'[',1))
library(biomaRt)
library(biomartr)
human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
uniprots = getBM(attributes = c("hgnc_symbol", "uniprot_gn"), filters = "hgnc_symbol", values = total$Row.names ,mart = human)
uniprots<- subset(uniprots, !duplicated(uniprots$hgnc_symbol))
colnames(uniprots)[1] <- "Row.names"
total <- merge(uniprots, total, by= "Row.names")
total<- total[!total$uniprot_gn=="",]
write.table(total, file = "E:/Bachelor/Bachelor/new_data/RPPA/merged.txt", sep = "\t")

# reading in complex definitions
comps<-read.table("E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/input-files/Protein-Complexes-UniprotIds.txt",sep="\t",quote="",comment.char="",header=TRUE)


##########
# missing value (NA) removal function
##########
#ns: number of samples
#nrep: number of replicates per sample
#maxna: maximum number of missing values allowed per sample

na.removal<-function(data,ns,nrep,maxna){
  NAcounts<-NULL
  p<-seq(1,(nrep*ns),by=nrep)
  for(i in 1:ns){
    tmp=apply(is.na(data[,p[i]:(p[i]+nrep-1)]),1,sum)
    NAcounts<-cbind(NAcounts,tmp)
  }
  datafiltered = data[apply(NAcounts,1,function(x) all(x <=maxna)),]
  datafiltered
}

# Filtering quantified proteins if there are multiple missing values in any cell line.
#defining number of samples
ns = 11
#number of replicates per sample
nrep = 3
#max number of missing value allowed per sample
maxna = 1
# filtering of proteins with missing values in any of the samples
dataset = na.removal(dataset,ns,nrep,maxna)

correlation = cor(dataset, method="pearson",use="pairwise.complete.obs")
pdf("E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/repro/Outlier-detection-correlation-of-samples.pdf", width=11, height=11)
heatmap.2(correlation,trace = "none",density.info="none", cexRow=0.9, cexCol=0.9, 
          colsep=c(seq(1,(nrep*ns),1)),
          rowsep=c(seq(1,(nrep*ns),1)),
          sepcolor='white',sepwidth=c(0.0125,0.02))
dev.off()

# samples A549.3 (3rd position), RKO.3 (30) and K562.3 (21) were identified as outliers and were discarded from further analysis.
#dataset = dataset[,c(-3,-21,-30)]

#preprocessing of the proteomics profiles
# profiles were centered around their trimmed mean.
medpoint = apply(dataset,2,function(x) mean(x,trim=0.5,na.rm=TRUE))
dataset = t(t(dataset)-medpoint)

##########################################################
comp.sel<-as.character(unique(comps$ProteinComplex))
comps2<-c()
data2<-c()
for(h in 1:length(comp.sel))
{
  print(h)
  subcomp = as.matrix(comps[comps$ProteinComplex==comp.sel[h],])
  subunitRedundantNames = subcomp[,"UniProtIDS"]
  for(m in 1:length(subunitRedundantNames))
  {
    subunitname = unique(unlist(strsplit(as.character(subunitRedundantNames[m]),":")))
    match = c()
    idmatch = c()
    for(s in 1:length(dataset[,1]))
    {
      subjectids = unique(unlist(strsplit(rownames(dataset)[s],";")))
      if(sum(subunitname%in%subjectids))
      {
        if(length(match)==0)
        {
          match = subunitname[subunitname%in%subjectids][1];
          idmatch = s;
        }
        else
        {
          print(c("reporting redundantcase in name matching:",h,m,s,subunitname))
        }
      }
    }
    if(length(match)>0)
    {
      subcomp[m,"UniProtIDS"] = match;
      rownames(dataset)[idmatch] = match;
      comps2 = rbind(comps2,subcomp[m,]);
    }
  }
}

# we only retain complexes having at least 5 members quantified.
comp.count = table(comps2[,"ProteinComplex"])
compsmin5members = comps2[comps2[,"ProteinComplex"]%in%names(comp.count)[comp.count>=5],]
# Finally, we also get the subset of the proteomics dataset that matches with the protein complexes.
ProteinCompDataset = dataset[rownames(dataset)%in%compsmin5members[,3],]

# writing output files that will be used in further analysis.
write.table(compsmin5members,"E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/repro/matched-protein-complexes-11-cell-line-dataset.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(ProteinCompDataset,"E:/Bachelor/Bachelor/source/Analysis_of_variable_protein_complexes/2_proteomics_datasets/11_cell_line_dataset/repro/11-cell-line-dataset-only-Protein-complexes.txt",sep="\t",quote=FALSE)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
