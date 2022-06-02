#!/home/proj/biosoft/software/R/R-3.2.2/bin/Rscript

x <- read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/distance_matrix_testis.txt"
                , sep="\t", header=TRUE)
x <- reshape(x,idvar="orgX",timevar="orgY",direction="wide")
row.names(x) <- x$orgX
#rownames(x) <- x$orgX
x$orgX <- NULL
#x <- data.matrix(x)
x
#methods for dist: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
y<-dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#methods for hclust: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
y<-hclust(y, method = "average", members = NULL )
par(mar=c(2, 4, 4, 2))
plot(y, labels = NULL, hang = 0.2, check = FALSE, yaxt= "n", axes = TRUE, frame.plot = FALSE, main = "Phylogenetic tree (PSI-value method)",
     ann = TRUE, xlab = NA, sub = NA , ylab = "Distance (testis tissue correlation dependant)")
axis(2, at=seq(0.2, 1.2, by=0.2))
#main("Phylogenetic tree (PSI-value method)")

