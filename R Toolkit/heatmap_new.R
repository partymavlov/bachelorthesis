#!/home/proj/biosoft/software/R/R-3.3.0/bin/Rscript

#call: Rscript evalueaDE.R <condition1Counts.file> <condition2Counts.file> <limmaTbl.file> <edgeRTbl.file> <DESeqTbl.file> <output>

args <- commandArgs(trailingOnly=TRUE)
#edit
if(length(args) != 1){
	stop("usage: Rscript task1.R <condition1Counts.file> <condition2Counts.file> <limmaTbl.file> <edgeRTbl.file> <DESeqTbl.file> <output>")
}
install.packages("plotly")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("reshape2")
library(plotly)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

#correlation.file <- args[1]
data <- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Not_pooled/correlation_not_pooled_homo_sap_brain.txt"
                   , sep="\t", header=TRUE)
data <- data[order(data$Var1),]
data <- data[order(data$Var2),]
data <- reshape(data,idvar="Var1",timevar="Var2",direction="wide")
row.names(data) <- data$Var1
data$Var1 <- NULL
colnames(data) <- row.names(data)



data.matrix <- data.matrix(data)

melted_data <- melt(data.matrix)
head(melted_data)
melted_data
#stop()
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+  geom_tile()+  scale_fill_gradient2(low = "blue", high = "red", mid = "white",  limit = c(0,1), space = "Lab",  name="Correlation") +  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(angle = 45, vjust = 1, size = 10, hjust = 1), plot.title=element_text(hjust=0.5))+  coord_fixed()+ ggtitle("Gene set size used for correlation heatmap for mutually exclusive exons events")
#ggsave(filename="plot1.pdf", plot=p)

p <- plot_ly(     x = rownames(data.matrix), y = colnames(data.matrix), z = data.matrix, type = "heatmap" )

