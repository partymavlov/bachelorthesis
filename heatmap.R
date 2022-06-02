#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript

#install.packages("gplot")
#install.packages('RColorBrewer')
#install.packages("plotly")


x <- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Pooled_new/correlation_pooled_alternative_splicing_site_full.txt", sep="\t", header=TRUE)
#x <- x[order(x$orgX),]
#x <- x[order(x$orgY),]
library(gplots)
library(plotly)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
y <- reshape2(x,idvar="orgX",timevar="orgY",direction="wide")
row.names(y) <- y$orgX
#colnames(y) <- y$orgX
#colnames(y) <- c(colnames(y)[1], gsub("^[[:alpha:]]", "", colnames(y)[-1]))
y$orgX <- NULL
#rownames(y) <- y$orgX
#colnames(y) <- y$orgY
y <- data.matrix(y)
y
#org_heatmap <- heatmap(y, Rowv=NA, main="Correlation heatmap for alternative splicing site events", Colv="Rowv", col = heat.colors(256), scale="column", margins=c(10,10) )
#hm <- heatmap.2(y, col=heat.colors(256), Rowv=NA, main="Correlation heatmap for alternative splicing site events", Colv="Rowv",density.info="none", trace="none",  dendrogram="none",
#                symm=TRUE, scale="column", key.title = "Correlation", key.xlab = NA,key.ylab = NA, keysize=1, densadj = 0.2, margins=c(12,12), lhei=c(0.2,0.2,2), lmat=rbind(c(0,0,0),c(0,0,0),c(0,0,0)), 
#                lwid=c(0.3,4,0.3)) 

melted_data <- melt(y)
head(melted_data)
#stop()
ggplot(data = melted_data, aes(Var2, Var1, fill = value))+  geom_tile()+  scale_fill_gradient2(low = "red", high = "white", mid = "yellow",    midpoint = 0.5, limit = c(0,1), space = "Lab",    name="Correlation") +  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(angle = 45, vjust = 1, size = 1, hjust = 1))+  coord_fixed()

#par(mar=c(10.1,5.1,15.1,10.1)
#ggplot(y, aes(variable, Name)) 
#geom_tile(aes(fill = rescale), colour = "white")
#scale_fill_gradient(low = "white", high = "steelblue")

#res <- data.matrix(subset(y, select=-orgX))
#rownames(res) <- y$orgX
#colnames(res) <- y$orgY
#res <- res[order(res["orgX"]),]
#res <- res[order(res["orgY"]),]



dev.off()



