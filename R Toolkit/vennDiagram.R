#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript

install.packages('VennDiagram')
library(VennDiagram)

#x<- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Mapper\ Comparison/mapper_comp_names.txt", sep="\t", header=TRUE)



par(mar=c(10,30,0,0))
a <- draw.pairwise.venn(area1 = 183, area2 = 42, cross.area = 35, 
                     category = c("down-regulated", "differentially expressed"), lty = "dashed", euler.d= T, scaled = T, cex = c(3, 3, 3), cat.cex = c(2, 2), ext.text = T, cat.pos= c(190,180), cat.dist = c(0.05, 0.05),
                 fill = c("blue", "red"), main="Protein regulation and differential expression intersection", main.pos= c(0, 0), height= 3000, width= 3000, resolution = 500);
require(gridExtra)
grid.arrange(gTree(children=a), top="Protein regulation and differential expression intersection")

jpeg('E:/Bachelor/Bachelor/rplot.jpg')
plot(a)
dev.off()
b<- draw.triple.venn(area1 = 43, area2 = 42, area3 = 183, n12 = 7, n23 = 35, n13 = 0, n123 = 0, category = c("up-regulated", "differentially expressed", "down-regulated"), lty = "dashed", euler.d= T, scaled = T, cex = c(3, 3, 3, 3, 3, 3, 3), cat.cex = c(2, 2, 2), ext.text = T, cat.pos= c(190,160, 160), cat.dist = c(0.05, 0.05, 0.05),
                     fill = c("blue", "red", "yellow"), main="Protein regulation and differential expression intersection", main.pos= c(0, 0), height= 3000, width= 3000, resolution = 500);
                   


a <- draw.quad.venn(area1 = 16, area2 = 10, area3 = 22, area4 = 4, n12 = 0, n13 = 14, n14 = 1, n23 = 7, n24 = 3, n34 = 0, n123 = 0, n124 = 0, n134 = 0, n234 = 0, n1234 = 0,
                        category = c("down-regulated", "up-regulated", "gene diff. exp.", "gene not diff. exp."), lty = "dashed", euler.d= T, scaled = T, cat.cex = c(2, 2, 2, 2), ext.text = T, cat.pos= c(190,180, 0, 0), cat.dist = c(0.05, 0.05, 0.05, 0.05),
                        fill = c("blue", "red" , "green", "orange"), main="Protein regulation and differential expression intersection", main.pos= c(0, 0), height= 3000, width= 3000, resolution = 500)
