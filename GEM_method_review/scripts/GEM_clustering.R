##This process is stochastic(contains randomized elements), so this enables us to get the same output each time we run our code.
set.seed(100)

require(vegan)
require("ggplot2")

par(mfrow=c(1,1))

#chose flux_maxes_single_cell.csv created by GEM_output_example.py
ranges<-read.csv(file.choose(),row.names=1)
#remove first column, because it contains row names which are already stored in data frame
ranges<-ranges[,-1]

#swap rows/columns of data frame so samples are rows and variables are columns
ranges<-as.data.frame(t(as.matrix(ranges)))

#replaces NA values with 0
ranges[is.na(ranges)]<-0

wss <- (nrow(ranges)-1)*sum(apply(ranges,2,var))
for (i in 2:10) wss[i] <- sum(kmeans(ranges,
                                     centers=i)$withinss)

###Plot helps choose number of clusters, K, but it is hard to interpret. Will likely change methodology later
plot(1:10, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",main=paste("Single Cell GEM Fit")) 



set.seed(100)
##Next 20 lines of code creates labels dataframe which associates each sample with its molecular subtype of breast cancer
set_cluster<-kmeans(ranges,5)
set_cluster$cluster <- as.factor(set_cluster$cluster)

u<-as.data.frame(set_cluster$cluster)

groups<-c("BC01","BC02","BC03","BC03LN","BC04","BC05","BC06","BC07","BC07LN","BC08","BC09","BC10","BC11")
types<-c("ER+","ER+","ER+ HER2+","Lymph metastatis BC03","HER2+","HER2+","HER2+","TNBC","Lymph metastasis BC07","TNBC","TNBC","TNBC","TNBC")
associations<-as.data.frame(rbind(groups,types))

samples<-row.names(ranges)


labels<-c("id","type")
for (i in samples){
  group<-unlist(strsplit(i,"_"))[1]
  type<-as.character(associations[,associations[1,]==group])[2]
  labels<-rbind(labels,c(i,type))
}

labels<-as.data.frame(labels)
colnames(labels)<-c("id","subtype")
labels<-labels[-c(1),]

all_types<-cbind(labels,u)
row.names(all_types)<-row.names(u)

#performs bray curtis clustering
mds <- cmdscale(vegdist(ranges, method='bray'))
new_mds <- cbind(data.frame(mds),all_types$subtype,all_types$`set_cluster$cluster`)
colnames(new_mds)<-c("x1","x2","Molecular Subtype","Cluster")

#plots bray curtis clustering in GGplot2, assining colors by cluster membership and shape by breast cancer subtype
#GGplot is extremely hard to read but produces high-quality graphs
p<-ggplot(new_mds, aes(new_mds$x1,new_mds$x2,group=new_mds$Cluster)) +
  geom_point(aes(color=new_mds$Cluster,shape=new_mds$`Molecular Subtype`),size=2) +
  scale_shape_manual(values = c(seq(0,13,1))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),plot.title = element_text(hjust = 0.5))
print(p+labs(title="Single Cell GEM Metabolic Potential Clustering",x="",y="",
             color="Cluster",shape="Molecular Subtype"))