rawdata<-read.csv("data.csv")
rawdata[rawdata$group==1,][1,]
rawdata[rawdata$group==2,][1:5,]
rawdata[rawdata$group==3,][1,]
rawdata[rawdata$group==4,][1,]
rawdata[rawdata$group==5,][1,]
data<-rawdata[,-1]
data<-data[,1:4]
head(data)

library("FactoMineR")
library("factoextra")

res.pca <- PCA(data, graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.character(rawdata$group), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#2E9FDF","#dd1c77"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             mean.point = TRUE
)

fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             group.ind = as.character(rawdata$group), # color by groups
             legend.title = "Groups",
             addEllipses = TRUE,
             mean.point = TRUE)

varcomp <- function(covmat,n) {
  if (is.list(covmat)) {
    if (length(covmat) < 2)
      stop("covmat must be a list with at least 2 elements")
    ps <- as.vector(sapply(covmat,dim))
    if (sum(ps[1] == ps) != length(ps))
      stop("all covariance matrices must have the same dimension")
    p <- ps[1] # 样本维度p
    m <- length(covmat) # 样本组数m
    if (length(n) == 1)
      nl <- rep(n,m)
    else if (length(n) == m)
      nl <- n # n存储每组样本的样本数量
    else
      stop("n must be equal length(covmat) or 1")
    
    DNAME <- deparse(substitute(covmat)) # data.name
  }
  
  else
    stop("covmat must be a list")
  Sl <- covmat # S_l是第l组样本协方差矩阵
  # (n_l-1)*S_l
  S1 <- lapply(1:length(covmat),function(i,mat,n) { (n[i]-1) * mat[[i]] },mat=covmat,n=nl)
  S <- matrix(colSums(matrix(unlist(S1),ncol=p^2,byrow=T)),ncol=p)
  S_pooled <- 1/sum(nl-1)*S # S_pooled是合并的样本协方差矩阵
  
  
  detSl <- sapply(Sl,det)
  det_S_pooled <- det(S_pooled)
  # lambda4 <- prod((detSl/det_S_pooled)^((nl-1)/2)) 
  # 这个很容易Nan，直接取对数进行运算
  ln_lambda4 <- sum((nl-1)/2*log(detSl/det_S_pooled))
  b <- (sum(1/(nl-1)) - 1/sum(nl-1))*(2*p^2+3*p-1)/(6*(p+1)*(m-1))
  
  STATISTIC <- -2*(1-b)*ln_lambda4
  mu <- p*(p+1)*(m-1)/2
  PVAL <- 1 - pchisq(STATISTIC,mu)

  names(STATISTIC) <- "corrected lambda*"
  names(mu) <- "df"
  RVAL <- structure(list(statistic = STATISTIC, parameter = mu,p.value = PVAL, data.name = DNAME, method = "Equality of Covariances Matrices Test"),class="htest")
  return(RVAL)
}
covmat <- vector("list",5) 
for (i in 1:5){
  covmat[[i]]<-cov(data[rawdata$group==i,])
}
n=c()
for(i in 1:5){
  n[i]=length(rawdata$group==i)
}
varcomp(covmat,n)

x1<-data[rawdata$group==1,];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-data[rawdata$group==2,];x2m<-apply(x2,2,mean);S2<-var(x2)
n=nrow(x1);m=nrow(x2);A1<-(n-1)*S1;A2<-(m-1)*S2;
Dsq=(n+m-2)%*%t(x1m-x2m)%*%solve(A1+A2)%*%(x1m-x2m)
p=2
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]

gnum<-length(unique(rawdata$group))
covmat <- vector("list",gnum) 
for (i in 1:gnum){
  covmat[[i]]<-cov(data[rawdata$group==i,])
}
meanmat <- vector("list",gnum) 
for (i in 1:gnum){
  meanmat[[i]]<-apply(data[rawdata$group==i,],2,mean)
}

group<-rawdata$group
x=data2
r=2
discri <- function(x,gnum,meanmat,covmat){
  
  re<-c()
  for(r in 1:nrow(x)){
    dist<-c()
    for(i in 1:gnum){
      dist[i]<-as.matrix(x[r,]-meanmat[[i]])%*%solve(covmat[[i]])%*%t(x[r,]-meanmat[[i]]) 
    }
    re[r]<-which(dist==min(dist))
  }
  return(re)
}
re<-discri(data,gnum,meanmat,covmat)
table(re==rawdata$group)/nrow(rawdata)

data2<-read.csv("data2.csv")
# 提取相应变量数据
data2<-cbind(data2$总胆红素.T.BIL.,
             data2$白蛋白.Alb.,data2$碱性磷酸酶.ALP.,data2$谷丙转氨酶.ALT.)
colnames(data2)<-c("BIL","ALB","ALP","ALT")
data2<-as.data.frame(na.omit(data2))
redata2<-discri(data2,gnum,meanmat,covmat)
result<-as.data.frame(cbind(data2,redata2))
table(redata2)
write.csv(result,"result.csv")
