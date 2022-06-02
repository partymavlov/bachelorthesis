#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
# author: Mario Picciani
# date: 2017-13-04 03:46:16
# 
# descr: create statistics and plots for organism
# 
############################################################
#
# call: Rscript ratioComp.R <workingDirectory> <counttype>
#

if(length(commandArgs()) != 7) message("usage: ratioComp.R <workingDirectory> <mapper> ");
stopifnot(length(commandArgs()) == 7);
#
#setwd(commandArgs()[6]);
#data<-read.csv(paste(commandArgs()[7],"_counttypes_ratio.tsv",sep=""),sep="\t",header=TRUE,row.names=1)
#data<-data[order(rowSums(data),decreasing=TRUE),]
#png("ratioComp.png")
#par(mar=c(10.1,4.1,4.1,2.1))
library(RColorBrewer)
data <- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Genes/Homo_sapiens_top_genes_reps.txt"
                   , sep="\t", header=TRUE)
#data <- data[order(data$gene),]
data <- data[order(data$tissue),]
data <- reshape(data,idvar="gene",timevar="tissue",direction="wide")
rownames(data) <- data$gene
data$gene <- NULL
data

colnames(data) <- sapply(colnames(data), function(x) sub("psi.", "", x))
data.matrix <- data.matrix(data)
data.matrix <- t(data.matrix)
par(mar=c(9,4.1,4.1,4.1))
matplot(data.matrix,type="l",xaxt = "n", yaxt= "n",ylab="PSI-value",  ylim=c(0,150),col=palette(rainbow(5)), lty= 1,main="Top 5 genes in Homo sapiens over all replicates")
axis(1, at=1:ncol(data), labels=colnames(data),las=2)
axis(2, at=seq(0, 150, by=10))
legend("topright",rownames(data),lty=1,col=1:5)


dev.off()
