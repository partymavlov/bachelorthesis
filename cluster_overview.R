
library("ComplexHeatmap")
data<-read.table("hisat_nfrags_filtered_ratios.matrix",sep="\t",header=TRUE)
names<-read.csv("hisat_nfrags_20_names.tsv",sep="\t",header=TRUE)
data$X<-NULL
my_matrix<-as.matrix(data)
colnames(my_matrix)<-NULL
rownames(my_matrix)<-NULL
#rownames(my_matrix)=colnames(my_matrix)
#rownames(my_matrix)<-colnames(my_matrix)
#ORG<-as.vector(c("Rattus","Bos","Homo","Macaca"))
#T1<-as.vector(c("liver","kidney","thyroid","brain"))
#T2<-as.vector(c("liver","brain","cerebellum","lung"))
#names<-data.frame(ORG,T1,T2)
human<-which(names[,1]=="Homo_sapiens")
names<-names[human,]
unique_orgs<-unique(names[,1])
unique_orgs
length(unique_orgs)
org_cols<-rainbow(length(unique_orgs))
org_cols
names(org_cols)<-unique_orgs

unique_t1<-as.vector(unique(names[,2]))
unique_t2<-as.vector(unique(names[,3]))
unique_ts<-unique(unlist(list(unique_t1,unique_t2)))
ts<-rainbow(length(unique_ts))
names(ts)<-unique_ts
ts
t1_col<-ts[unique_t1]
names(t1_col)<-unique_t1
t2_col<-ts[unique_t2]
names(t2_col)<-unique_t2
t1_col
t2_col
ha1<-HeatmapAnnotation(df=names,col=list(ORG=org_cols,T1=t1_col,T2=t2_col),show_legend = F)
ha2<-HeatmapAnnotation(df=names,col=list(ORG=org_cols,T1=t1_col,T2=t2_col),which = "row")
ht1 <-Heatmap(my_matrix[human,human],top_annotation = ha1,name="ratio")
ha2 + ht1


