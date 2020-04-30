ex1<-read.csv("ex6.7.csv")
rownames(ex1)=ex1[,1]
ex1<-ex1[,-1]
Bartl<-factanal(ex1, 3, scores="Bartlett", rotation="none")
head(Bartl[["scores"]])
Regre<-factanal(ex1, 3, scores="regression", rotation="none")
head(Regre[["scores"]])


rawdata<-read.csv("data.csv",stringsAsFactors = FALSE)
rawdata$Weight<-as.numeric(rawdata$Weight)
rawdata$Height<-as.numeric(rawdata$Height)
rawdata<-na.omit(rawdata) 
data<-rawdata[,-c(1:3)]
str(data)
library(psych)
fa.parallel(data)
pc<-principal(data,nfactors=7) 
pc
pc_score<-pc[["scores"]]

library(factoextra) 
b<-scale(pc_score)
set.seed(123)  
fviz_nbclust(b,kmeans,method="wss")+geom_vline(xintercept=3,linetype=2)

res<-kmeans(b,3)
res1<-cbind(rawdata[,1:3],res$cluster)
fviz_cluster(res,data=data)

Bartl<-factanal(data[,1:14], 3, scores="none", rotation="none")
pr<-principal(data,nfactors = 7,rotate = "varimax",covar = F,scores = TRUE)
load<-as.data.frame(pr[["loadings"]])
library("corrplot")
mr <- fa(data,7,rotate="varimax",fm="pa")  # principal axis factor analysis
load<-mr[["loadings"]][1:15,]
corrplot(load, is.corr=FALSE) 
mr_score<-mr[["scores"]]
regredata<- as.data.frame(cbind(mr_score,rawdata$Age)) 
regre <-lm(V8~PA3+PA4+PA2+PA1+PA6+PA5+PA7, data=regredata) 
summary(regre)


factor.analy1<-function(S, m){
  p<-nrow(S); diag_S<-diag(S); sum_rank<-sum(diag_S)
  rowname<-paste("X", 1:p, sep="")
  colname<-paste("Factor", 1:m, sep="")
  A<-matrix(0, nrow=p, ncol=m,
            dimnames=list(rowname, colname))
  eig<-eigen(S)
  for (i in 1:m)
    A[,i]<-sqrt(eig$values[i])*eig$vectors[,i]
  h<-diag(A%*%t(A))
  rowname<-c("SS loadings","Proportion Var","Cumulative Var")
  B<-matrix(0, nrow=3, ncol=m,
            dimnames=list(rowname, colname))
  for (i in 1:m){
    B[1,i]<-sum(A[,i]^2)
    B[2,i]<-B[1,i]/sum_rank
    B[3,i]<-sum(B[1,1:i])/sum_rank
  }
  method<-c("Principal Component Method")
  list(mehod=method, loadings=A,
       var=cbind(common=h, spcific=diag_S-h), B=B)
} 
S<-cor(ex1)
factor.analy1(S,3)
