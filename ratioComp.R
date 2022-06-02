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
setwd(commandArgs()[6]);
data<-read.csv(paste(commandArgs()[7],"_counttypes_ratio.tsv",sep=""),sep="\t",header=TRUE,row.names=1)
data<-data[order(rowSums(data),decreasing=TRUE),]
png("ratioComp.png")
par(mar=c(10.1,4.1,4.1,2.1))
matplot(data,type="l",xaxt = "n",ylab="ratio",ylim=c(0,1),col=1:8,main="comparison between ratios over organisms")
axis(1, at=1:nrow(data), labels=row.names(data),las=2)
legend("topright",names(data),lty=1,col=1:8)
dev.off()
