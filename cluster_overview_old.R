

setwd("C:/Users/Mario/Desktop/GOBI_tmp")
data<-read.table("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Pooled_new/matrix",sep="\t",header=TRUE)
names<-read.csv("/home/p/pavlovma/Java2/Java/Eclipse/Gobi/output/assignmentBlock/Correlations/Pooled_new/legend",sep="\t",header=TRUE)
data$X<-NULL
my_matrix<-as.matrix(data)
my_matrix
#rownames(my_matrix)=colnames(my_matrix)
colnames(my_matrix)<-NULL
rownames(my_matrix)<-NULL
#ORG<-as.vector(c("Rattus","Bos","Homo","Macaca"))
#T1<-as.vector(c("liver","kidney","thyroid","brain"))
#T2<-as.vector(c("liver","brain","cerebellum","lung"))
#names<-data.frame(ORG,T1,T2)
names
unique_orgs<-unique(names[,1])
unique_orgs
length(unique_orgs)
org_cols<-rainbow(length(unique_orgs))
org_cols
names(org_cols)<-unique_orgs

unique_t1<-as.vector(unique(names[,2]))
unique_ts<-rainbow(length(unique_t1))
names(ts)<-unique_t
ts
t_col<-ts[unique_t]
names(t_col)<-unique_t
t_col

ha1<-HeatmapAnnotation(df=names,col=list(ORG=org_cols,T=t_col),show_legend = F)
ha2<-HeatmapAnnotation(df=names,col=list(ORG=org_cols,T=t_col),which = "row")
ht1 <-Heatmap(my_matrix,top_annotation = ha1,name="ratio")
ha2 + ht1


