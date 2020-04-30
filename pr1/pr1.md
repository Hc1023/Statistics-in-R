# R基础和统计图示

## CLT模拟

（1）利用R软件生成随机数，分析如下几个参数对中心极限定理(CLT)的影响：
（a）样本容量 n：是否n的取值越大，CLT的拟合效果越好？
（b）模拟次数 m: 分析不同模拟次数对于CLT拟合效果的影响；
（c）随机数的分布：在相同样本容量和模拟次数下，取自不同分布的随机数对CLT有什么影响？

中心极限定理（CLT）$$\{X_i\}\ i.i.d \sim\ (\mu,\sigma^2)$$
记$S_n=\Sigma_{i=1}^{n}X_i$，则满足 $$W_n=\frac{S_n-E(S_n)}{\sqrt{var(S_n)}}\xrightarrow[n\rightarrow\infty]{d}N(0,1)$$

```R
require(graphics)
rm(list=ls())

fclt<-function(n,m,fnum){
  w<-vector(mode="numeric",length=m)
  fn=switch(fnum,rt,runif)
  for (i in 1:m){
    # wnorm=rnorm(m)
    if (fnum == 1){# t分布
      df=5;x=rt(1:n, df=df)
      w[i]=sum(x)/(sqrt(n*df/(df-2)))
    }else if (fnum == 2){# [0,1]均匀分布
      x=runif(n, min = 0, max = 1)
      w[i]=(sum(x)-n*0.5)/(sqrt(n/12))
    }else if(fnum == 3){# 卡方分布
      df=5;x=rchisq(n,df=df,ncp=0)
      w[i]=(sum(x)-n*df)/(sqrt(n*2*df))
    }else if(fnum == 4){# F分布
      df1=3;df2=5;x=rf(n,df1=df1,df2=df2,ncp=0)
      v=2*df2^2*(df1+df2-2)/(df1*(df2-2)^2*(df2-4))
      w[i]=(sum(x)-n*(df2/(df2-2)))/(sqrt(n*v))
    }
  }
  ## Have a look at the densities
  plot(density(w))
  ## hist
  hist(w)
  ## Perform the test
  # print(unlist(shapiro.test(wnorm))[1:2])
  print(unlist(shapiro.test(w))[1:2])
  
}
# fclt(样本容量n, 模拟次数m, 随机变量fnum)
# fnum 1: t分布 2: [0,1]均匀分布 3: 卡方分布 4: F分布

fclt(30000,100,4)
```
(a) 取F分布为例，取m=100，n较小是，p值小，否定原来的“为正态分布”的原假设，确实是n越大，CLT效果越好，这也正印证了中心极限定理
```R
> fclt(10,100,4)
           statistic.W                p.value 
   "0.885283766689354" "3.05790289028539e-07" 
> fclt(100,100,4)
           statistic.W                p.value 
   "0.835925783099004" "3.70441187997276e-09" 
> fclt(500,100,4)
           statistic.W                p.value 
   "0.946722423424597" "0.000507830934013186" 
> fclt(1000,100,4)
        statistic.W             p.value 
"0.979001382417503" "0.111096689923914" 
> fclt(5000,100,4)
        statistic.W             p.value 
"0.989184759111393" "0.599160300921035" 
```

- n=10
![n=10](https://upload-images.jianshu.io/upload_images/17587926-0ecb0be06d6bb8d9.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![n=10](https://upload-images.jianshu.io/upload_images/17587926-e52123de07eaf4ed.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- n=5000
![n=5000](https://upload-images.jianshu.io/upload_images/17587926-35aa64455c301ce0.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![n=5000](https://upload-images.jianshu.io/upload_images/17587926-03e948b4bae6487f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

(b)从F分布的随机变量上看，m变大似乎CLT效果变差了，然而其它分布随机变量似乎没有很大的影响
```
> fclt(1000,50,1)
        statistic.W             p.value 
"0.986368604263566" "0.828481828547294" 
> fclt(1000,5000,1)
        statistic.W             p.value 
"0.999637322739418" "0.518279189942571" 
> fclt(1000,50,2)
        statistic.W             p.value 
"0.979515531909071" "0.531714677872245" 
> fclt(1000,5000,2)
        statistic.W             p.value 
"0.999702799571169"  "0.71507223205894" 
> fclt(1000,50,3)
        statistic.W             p.value 
"0.973350997552889" "0.315087007830422" 
> fclt(1000,5000,3)
        statistic.W             p.value 
"0.999570017722336" "0.344273646107321" 
> fclt(1000,50,4)
        statistic.W             p.value 
"0.981613517259601" "0.621724898190725" 
> fclt(1000,5000,4)
           statistic.W                p.value 
   "0.964348978725153" "1.94456354298482e-33" 
```
(c)选取了1: t分布 2: [0,1]均匀分布 3: 卡方分布 4: F分布 进行模拟，都分别印证了CLT。m合适时，当n变大都趋于正态分布。
可以直接从上面的代码结果中看出，n=1000时四种随机变量都趋于正态。而随机变量本身若接近正态分布，则拟合程度更好，在更小的n时就趋于正态，对于这四个随机变量可以发现F分布需要较大的n才能通过假设检验，而前三种在较小的n就能通过假设检验。


## 数据分析图示
![](https://upload-images.jianshu.io/upload_images/17587926-77fcdb8e18775109.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

（1）	建立 R数据集
```R
# write.csv(NA, file = "MyData.csv")
mydata = read.csv("MyData.csv",row.names = 1)
```
（2）	单变量数据包括集中趋势的特征值：均值、众数、中位数等；离散趋势的特征值： 标准差、方差，(0.05, 0.95, 0.025, 0.975,0.1, 0.9)分位数。
（这里结果比较琐碎不放上）
```R
# mean
apply(mydata,2,mean)
# mode
# v=c(1,3,2,2,3,5)
# get one mode
Mode1 <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# what to do when the number of modes is more than one
Mode2 <- function(v) {
  a <- table(v)
  as.numeric(names(a)[a == max(a)])
}
# however sometimes the data is just scattered
Mode3 <- function(v) {
  a <- table(v)
  b <- as.numeric(names(a)[a == max(a)])
  if(length(b)>2){
    print("Oops, the data is scattered.")
  }else{
    return(b)
  }
}
apply(mydata,2,Mode3)# Hah, this data has no mode!
apply(iris[,1:4],2,Mode3)# iris is always a good one to analyze

apply(mydata,2,median)# median
apply(mydata,2,sd)# sd
apply(mydata,2,var)# var

myquantile<-function(x){
  return(quantile(x, probs = c(0.05, 0.95, 0.025, 0.975,0.1, 0.9)))
}
apply(mydata,2,myquantile)# quantile
```
- 轮廓图
```R
# contour
library(lattice)

parallel(~mydata, mydata, groups = rownames(mydata),
         horizontal.axis = FALSE, scales = list(x = list(rot = 90)))
# some other contour plot
v <- ggplot(faithfuld, aes(waiting, eruptions, z = density))
v + geom_contour()
# look better maybe
v + geom_raster(aes(fill = density)) +
  geom_contour(colour = "white")
```
![轮廓图](https://upload-images.jianshu.io/upload_images/17587926-daa211ba1a994407.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```R
# Libraries
library(GGally)

# Data set is provided by R natively
data <- data.frame(mydata,rownames(mydata))
colnames(data)[7]<-"city"

# Plot
ggparcoord(data,
           columns = 1:6, 
           groupColumn = 7
)
```
![1005.png](https://upload-images.jianshu.io/upload_images/17587926-000332bb645ffa5d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 雷达图
```R
# radar plot
# devtools::install_github("ricardo-bion/ggradar", 
#                         dependencies = TRUE)
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
mydata_radar <- mydata %>%
  as_tibble(rownames = "group") %>%
  mutate_at(vars(-group),rescale)
ggradar(mydata_radar)
```
![雷达图](https://upload-images.jianshu.io/upload_images/17587926-317082644964be1f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 调和曲线图

假设 $X_r$ 是$p$维数据的第$r$个观测值，即
$$X_r=(x_{r1},x_{r2},...,x_{rp})'$$

则调和曲线为
$$f_r(t)=\frac{x_{r1}}{\sqrt{2}}+x_{r2}sin t+x_{r3}cos t+x_{r4}sin 2t+x_{r5}cos 2t+...,\ -\pi\leq t\leq \pi$$
```R
# harmonic graph
x = as.matrix(mydata)
t = seq(-pi, pi, pi/30)
m = nrow(x)
n = ncol(x)
f = matrix(0, m, length(t))
for (i in 1:m) {
  f[i, ] = x[i, 1]/sqrt(2)
  for (j in 2:n) {
    if (j%%2 == 0)
      f[i, ] = f[i, ] + x[i, j] * sin(j/2 * t)
    else f[i, ] = f[i, ] + x[i, j] * cos(j%/%2 * t)
  }
}
plot(c(-pi, pi), c(min(f), max(f)), type = "n", main = "The harmonic graph of data",
     xlab = "t", ylab = "f(t)")
for (i in 1:m) lines(t, f[i, ],col = i)
legend(x = -3, y = 400, rownames(mydata),
       lty = 1, col = c(1:m),cex=0.8,
       y.intersp=0.15,x.intersp = 0.05,
       bty="n", ncol=2,lwd=2)
```
![调和曲线图](https://upload-images.jianshu.io/upload_images/17587926-dfc9556da45497dc.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```R
library(andrews)
andrews(data,clr=7,ymax=4.5)
legend( "topright",legend = data$city,
        col=unique(data[,7]),lty=1,
        cex=0.8,ncol=2,bty="n")
```
![1005-2.png](https://upload-images.jianshu.io/upload_images/17587926-90fbb4b3a84f71e9.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 散布图矩阵
```R
# scatterplot-matrix
pairs(mydata,main="scatterplot-matrix",
      pch = 21, bg=c(1:m))
```
![散布图](https://upload-images.jianshu.io/upload_images/17587926-3661b4784aee21bd.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

分析：从各种图中能更直观地发现数据之间的比较和关系，显然的是杭州基本从各方面碾压了其它城市，但这也是自然的。