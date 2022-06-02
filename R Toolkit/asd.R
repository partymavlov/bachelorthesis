#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript
library(reshape2)
library("ComplexHeatmap")
library(RColorBrewer)
library(colorspace)
pal <-choose_palette()
data <- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Pooled_new/correlation_pooled_mutually_exclusive_exons_full_quant.txt"
                   , sep="\t", header=TRUE)
data <- data[order(data$Var1),]
data <- data[order(data$Var2),]
data

y <- reshape(data,idvar="Var1",timevar="Var2",direction="wide")
row.names(y) <- y$Var1
#colnames(y) <- y$orgX
#colnames(y) <- c(colnames(y)[1], gsub("^[[:alpha:]]", "", colnames(y)[-1]))
y$Var1 <- NULL
#rownames(y) <- y$orgX
#colnames(y) <- y$orgY
y <- data.matrix(y)
y
orgs<-rownames(y)
tissue = sub(".*_.*_","",orgs)
org =sub("_[^_]*$","",orgs)
legend<-data.frame(ORG<-org,T<-tissue)
write.table(legend,file="/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Pooled_new/legend",row.names=FALSE,quote=FALSE)
names<-unique(legend)
org_colors <- rainbow(8)
org_colors
#org_colors<- rainbow_hcl(8)
#names(org_colors)<-unique(legend[,1])
tissue_colors<-rainbow(length(unique(legend[,2])))
names(tissue_colors)<-unique(legend[,2])
#names(tissue_colors)<- unique(tissue)
ha2<-HeatmapAnnotation(df=legend,col=list(ORG=org_colors, T=tissue_colors),which="row")
ha1<-HeatmapAnnotation(df=legend,col=list(ORG=org_colors,T=tissue_colors),show_legend=F)
colnames(y)<-NULL
rownames(y)<-NULL
ht<-Heatmap(y,name="Num of genes used",top_annotation=ha1)
ha2 + ht

