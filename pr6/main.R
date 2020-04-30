mydata<-read.csv("ex4.3.csv")
rownames(mydata)<-mydata[,1]
mydata<-mydata[,-1]
library(cluster)
# Compute with agnes (make sure you have the package cluster)
hc1 <- agnes(mydata, method = "ward")
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram of agnes")
rect.hclust(hc1, k = 3)
#先求样本之间两两相似性
result <- dist(mydata, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
#进行初步展示
fviz_dend(result_hc, cex = 0.6)
res2<-fviz_dend(result_hc, k = 5, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07","#2ca25f","#dd1c77"),
          color_labels_by_k = TRUE, 
          rect = TRUE          
)
summary(result_hc)

heatmap(as.matrix(result), labRow = F, labCol = F)

re <- cutree(result_hc, k = 5)

mds = cmdscale(result, k = 2, eig = T)

X <- mds$points[, 1]
Y <- mds$points[, 2]

p = ggplot(data.frame(X, Y ), aes(X, Y ))
library(ggplot2)
p + geom_point(size = 3, alpha = 0.8, aes(colour =factor(re)))

                                         
library(factoextra) 
b<-scale(mydata) 
set.seed(123)  
#确定最佳聚类个数，使用组内平方误差和法  

fviz_nbclust(b,kmeans,method="wss")+geom_vline(xintercept=5,linetype=2) 

library(cluster)
res<-kmeans(b,5)
res1<-cbind(mydata,res$cluster)
fviz_cluster(res,data=mydata[,1:ncol(mydata)-1])

head(mydata)


data<-read.csv("ex6.7.csv")

rownames(data)<-data[,1]
data<-data[,-1]
head(data)
pr<-princomp(data,cor=TRUE)
screeplot(pr,type="lines")
summary(pr,loadings=TRUE) 

library("FactoMineR")
library("factoextra")
data<-read.csv("ex6.7.csv")
rownames(data)<-data[,1]
data<-data[,-1]
res.pca <- PCA(data, graph = FALSE)
print(res.pca)
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(res.pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
fviz_pca_var(res.pca, col.var = "black")
library("corrplot")
corrplot(var$cos2, is.corr=FALSE)
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)


