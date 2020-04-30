data<-read.csv("ex7.11.csv")
pr<-princomp(data,cor=TRUE)
screeplot(pr,type="lines")
summary(pr,loadings=TRUE) 
res<-predict(pr)[,1:3]
rank(predict(pr)[,1])
mycol<-c("#FF0000","#FFFF00","#008B8B","#7FFFD4","#FF00FF","#0000FF",
         "#8A2BE2","#A52A2A","#000000","#7FFF00","#80000040","#FF7F50","#6495ED")
plot(res,main= 'PCA plot',xlab= "PC1",ylab= "PC2",
     type= 'p',col=mycol,cex=2,lwd=2)  #画点�?

legend("topright",title="Sample", legend = rownames(res),
       pch=1,col = mycol,cex=0.5,ncol=2)
library(cluster)
# Compute with agnes (make sure you have the package cluster)
hc1 <- agnes(res, method = "ward")
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram of agnes")
rect.hclust(hc1, k = 3)

data<-read.csv("ex8-6.txt")
data<-read.table("ex8-6.txt")
data2<-data[,6:10]
data<-data[,-c(6:10)]
dat<-cbind(t(data),t(data2))
data<-as.data.frame(t(dat))
library(psych)
principal(data,nfactors = 4,rotate = "varimax",covar = F,scores = TRUE)

library("corrplot")
mr <- fa(data,4,rotate="varimax",fm="pa")  # principal axis factor analysis
load<-mr[["loadings"]][1:5,]
corrplot(load, is.corr=FALSE) 
factanal(data, 3, scores="none", rotation="varimax")
Bartl
