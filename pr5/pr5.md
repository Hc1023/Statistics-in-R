# 判别分析

之前主要做了两类的判别，下面进行$k$类的判别分析。
不同之处在于组间均值显著性差异与组间方差显著性差异的不同。

利用判别分析的统计方法建立肝胆疾病的判别函数
（a）检验组间均值显著性差异
（b）分析组间方差显著性差异
（c）建立判别函数 
（d）分析判别效果
并应用于“体检数据”中，根据体检资料分析是否有得肝胆疾病的可能性。

数据的临床诊断指标为总胆红素（T-BIL），白蛋白（Alb），碱性磷酸酶（ALP），谷丙转氨酶（ALT），样本分为五组，1代表慢肝，2比较杂糅，3代表慢胆，4代表急胆，5代表正常。
```r
rawdata<-read.csv("data.csv")
```

- 观察数据

```r
> rawdata[rawdata$group==1,][1,]
  临床诊断  BIL  Alb ALP ALT group
1     肝炎 97.5 40.1 189 446     1
> rawdata[rawdata$group==2,][1:5,]
    临床诊断   BIL  Alb ALP  ALT group
109     肝炎 140.8 39.0 241  761     2
110     肝炎 127.2 37.0 126 3443     2
111     急肝  42.8 30.3 187  623     2
112     急肝 112.7 33.7 132 1232     2
113     肝炎 292.9 29.4 206 1364     2
> rawdata[rawdata$group==3,][1,]
    临床诊断  BIL  Alb ALP ALT group
124     慢胆 18.3 45.5  71 108     3
> rawdata[rawdata$group==4,][1,]
    临床诊断  BIL  Alb ALP ALT group
183     急胆 20.9 46.5  81  34     4
> rawdata[rawdata$group==5,][1,]
    临床诊断  BIL  Alb ALP ALT group
217     正常 19.8 51.9 118  22     5
```

- 处理data

```r
data<-rawdata[,-1]
data<-data[,1:4]
```

## Step1 $k$个多元正态总体的均值方差分析

我们要对$k$类总体进行判别分析，首先希望这$k$类总体有均值显著性差异，而若它有方差显著性差异那么一般的线性判别法将忽略掉方差信息，也会影响它的准确性。

首先可以提取它们前两个主成分直观地感受不同肝胆疾病的数据差异。

```r
library("FactoMineR")
library("factoextra")
res.pca <- PCA(data, graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.character(rawdata$group), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07","#2E9FDF","#dd1c77"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
```
![2.png](https://upload-images.jianshu.io/upload_images/17587926-1c16877fddbe9258.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

每个椭圆中心都有点代表它们的二维均值，直观上看它们是有差异的，且它们的方差也有差异。

### 多总体协方差阵的检验

检验
$$H_0: \Sigma_{1}=\Sigma_{2}=\ldots=\Sigma_{m}=\Sigma$$

构造检验统计量
$$\lambda_{4}=\prod_{l=1}^{m}\left(\left|\mathbf{S}_{l}\right| /\left|\mathbf{S}_{p o o l e d}\right|\right)^{\left(n_{l}-1\right) / 2}$$

其中，$n_l$是第$l$组的样本量，$S_l$是第$l$组
样本协方差矩阵，$S_{pooled}$是合并的样本协方差矩阵，满足

$$\mathbf{S}_{p o o l e d}=\frac{1}{\sum_{l=1}^{m}\left(n_{l}-1\right)} \sum_{l=1}^{m}\left(n_{l}-1\right) \mathbf{S}_{l}$$

当样本容量$n=\Sigma_{l=1}^m n_l$很大时，在$H_0$为真时有以下近似分布

$$
-2(1-b) \ln \lambda_{4} \sim_{a p p r o x} \chi^{2}(\nu), \nu=p(p+1)(m-1) / 2
$$

其中 $b=\left(\sum_{l=1}^{m}\left(n_{l}-1\right)^{-1}-\left(\sum_{l=1}^{m}\left(n_{l}-1\right)\right)^{-1}\right)\left(\frac{2 p^{2}+3 p-1}{6(p+1)(m-1)}\right)$

由此我们可以自己写出进行**多元变量多总体协方差的检验**的函数，以下代码的字符与运算公式基本与前面的数学公式保持一致。
一个需要注意的事情是$\lambda_4$是幂次乘积直接对它去做运算极有可能不是0就是Inf，这是十分不明智的，所以应该直接取对数运算。总之在写代码时还是需要一些小skills，请体谅一下计算机的运算精度。

```r
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
```

- 函数应用

首先计算好每个组别的协方差阵将其存为list，每个组别的样本量存为一个向量。

```r
covmat <- vector("list",5) 
for (i in 1:5){
  covmat[[i]]<-cov(data[rawdata$group==i,])
}
n=c()
for(i in 1:5){
  n[i]=length(rawdata$group==i)
}
```

Perfectly done!

```
> varcomp(covmat,n=n)

	Equality of Covariances Matrices Test

data:  covmat
corrected lambda* = 9661.6, df = 40, p-value < 2.2e-16
```

其$p$值远小于$0.01$，说明方差有显著性差异。

### 多总体均值的检验

**多元方差分析**

这里我们需要在方差不相等情况下进行均值的检验。

又没有现成的函数可以用，我又要写理论写函数了，天越来越冷了，我越来越懒了，我们直接抓两个组别出来比一下。

```r
x1<-data[rawdata$group==1,];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-data[rawdata$group==2,];x2m<-apply(x2,2,mean);S2<-var(x2)
n=nrow(x1);m=nrow(x2);A1<-(n-1)*S1;A2<-(m-1)*S2;
Dsq=(n+m-2)%*%t(x1m-x2m)%*%solve(A1+A2)%*%(x1m-x2m)
p=2
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]
```

$p=1.575946e-13$，故在显著水平0.01时两总体的均值有显著差异。

是否有现成的函数或包比如说ANOVA做**多元变量多总体而且方差不齐性的均值检验**有待探索。
好了，上面的这些工作其实都是为了保证后面的判别有意义。

## Step2 多总体判别

### 建立判别函数

利用**马氏距离**进行多总体距离判别
$$d^2(X,G)=(X-\mu)'\Sigma^{-1}(X-\mu)$$

- 多总体距离判别

把$X$判给距离最近的那一个，设$i=l$时，若
$$d_l^2(X)=min_{i=1,...,k}\{d_i^2(X)\}$$则$X\in G_l$

```r
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
discri <- function(x,gnum,meanmat,covmat){
  re<-c()
  r=nrow(x)
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
```

判断正确的概率为56.1%，由于是五类判别，样本数据集也不多，还过得去吧。

```r
> table(re==rawdata$group)/nrow(rawdata)
    FALSE      TRUE 
0.4389535 0.5610465 
```

可以稍微看看它的结果，我发现它很容易把第五类错判为第三类，再反观之前的PCA图发现两者确实重合度比较高。所以说如果被判为慢胆，很容易是误诊说不定你就是正常人？？？
好吧其实还是可以考虑换其它判别方法再做尝试。或者有其它诊断指标，进一步鉴定……

## Step3 应用于待测数据

然后我们将其应用于“体检数据”，进行医生肝胆病乱诊断系列

```r
data2<-read.csv("data2.csv")
# 提取相应变量数据
data2<-cbind(data2$总胆红素.T.BIL.,
             data2$白蛋白.Alb.,data2$碱性磷酸酶.ALP.,data2$谷丙转氨酶.ALT.)
colnames(data2)<-c("BIL","ALB","ALP","ALT")
data2<-as.data.frame(na.omit(data2))
redata2<-discri(data2,gnum,meanmat,covmat)
result<-as.data.frame(cbind(data2,redata2))
write.csv(result,"result.csv")
```

```r
> table(redata2)
redata2
  1   3   4   5 
  4 159  10  49 
```

哈哈哈根据前面的分析判断为第三类（慢胆）的人极有可能是正常人，得出结论：看病不要找统计学家……
尴尬而不失礼貌地微笑，结束啦