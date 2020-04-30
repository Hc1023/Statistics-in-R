# 正态检验和置信区域

## 实验一
- 数据预处理
使用boxplot剔除了极端情况下的异常值并删除了有缺损的数据
```R
{
  sample<-read.csv("sample.csv",na.strings=c(""),stringsAsFactors=FALSE)
  #str(sample)
  sample<-na.omit(sample[-1,])# 对na数据的处理
  sample[,c(3:7,10:12)]<-as.numeric(unlist(sample[,c(3:7,10:12)]))
  #str(sample)
  bmi<-sample$weight/(0.01*sample$height)^2
  sample<-data.frame(sample,bmi)
  x<-data.frame(sample$bmi,sample$FPG,sample$sbp,
                sample$dbp,sample$TG,sample$HDL.C)
  # 剔除异常值
  sapply(x,function(X){
    boxplot.stats(sample$bmi)$out
  })
  # 观察发现可能只有bmi需要处理
  outlier_location <- which(sample$bmi%in%boxplot.stats(sample$bmi)$out)
  todel <- (sort(unique(unlist(outlier_location))))
  sample<-sample[-todel,]
  x<-x[-todel,]
}
```

判断是否患有代谢综合征
```R
{
  b1<-as.numeric(sample$bmi>=25)
  b2<-as.numeric(sample$FPG>=6.1)
  b3<-as.numeric(sample$sbp>=140|sample$dbp>=90)
  b4_1<-sample$TG>=1.7
  b4_2<-sample$gender == "F" & sample$HDL.C<1
  b4_3<-sample$gender == "M" & sample$HDL.C<0.9
  b4<-as.numeric(b4_1 | b4_2 | b4_3)
  b<-(b1+b2+b3+b4>=3)
}
```

检验相关数据的相关性
```R
{
  xmean<-apply(x, 2, mean) # 样本均值
  xdispersion<-(nrow(x)-1)*var(x) # 样本离差阵
  xvar<-var(x) # 样本协方差
  xcor<-cor(x) # 样本相关阵
}
```
检验数据的正态性
- 单变量正态性检验
对每个变量进行Q-Q Plot正态性检验
```R
# one dimension Normal Q-Q Plot
{
  par(mfrow=c(3,2))
  qqnorm(x$sample.bmi, main = "bmi",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
  qqnorm(x$sample.FPG, main = "FPG",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
  qqnorm(x$sample.sbp, main = "sbp",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
  qqnorm(x$sample.dbp, main = "dbp",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
  qqnorm(x$sample.TG, main = "TG",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
  qqnorm(x$sample.HDL.C, main = "HDL-C",
         xlab = "Theoretical Normal Quantiles",
         ylab = "Sample Normal Quantiles",pch=19)
}
```
![1027.png](C:\Users\Lenovo\Desktop\学习\多元统计\pr\pr3\1027.png)
各变量的QQ图基本落于一条直线上，我们可以认为他们符合正态性分布

- p维变量正态性检验
```R
# p-dimension Normal Q-Q plot
p_qqplot<-function(x){
  S<-var(x)
  xmean<-as.matrix(apply(x, 2, mean))
  D.2<-apply(x,1,function(X){
    t(X-xmean)%*%solve(S)%*%(X-xmean)
  })
  D.2<-sort(D.2) # 从小到大进行排序，马氏距离
  t<-1:nrow(x)
  p<-(t-0.5)/nrow(x)
  chi.2<-qchisq(p, df=ncol(x), lower.tail = TRUE) # TRUE代表左侧面积为p
  data<-data.frame(D.2,chi.2)
  # 剔除异常值
  {
    outlier_location <- sapply(data,function(X){
      which(X%in%boxplot.stats(X)$out)
    })
    todel <- (sort(unique(unlist(outlier_location))))
    data<-data[-todel,]
  }
  colnames(data)<-c("Mahalanobis_distance","chi_square_quantile")
  ggplot(data, aes(x=Mahalanobis_distance, y=chi_square_quantile))+ 
    geom_point() + stat_smooth(method=lm, level=0.95)+ # 95%的置信度
    labs(title="p-dimension Q-Q plot")+
    geom_abline(intercept=0,slope=1)
}
p_qqplot(x)
```
![1027-2.png](C:\Users\Lenovo\Desktop\学习\多元统计\pr\pr3\1027-2.png)

```R
# p-dimension Normal Q-Q plot
p_qqplot<-function(x){
  S<-var(x)
  xmean<-as.matrix(apply(x, 2, mean))
  D.2<-apply(x,1,function(X){
    t(X-xmean)%*%solve(S)%*%(X-xmean)
  })
  D.2<-sort(D.2) # 从小到大进行排序，马氏距离
  t<-1:nrow(x)
  p<-(t-0.5)/nrow(x)
  chi.2<-qchisq(p, df=ncol(x), lower.tail = TRUE) # TRUE代表左侧面积为p
  data<-data.frame(D.2,chi.2)
  # 剔除异常值
  {
    outlier_location <- sapply(data,function(X){
      which(X%in%boxplot.stats(X)$out)
    })
    todel <- (sort(unique(unlist(outlier_location))))
    data<-data[-todel,]
  }
  colnames(data)<-c("Mahalanobis_distance","chi_square_quantile")
  ggplot(data, aes(x=Mahalanobis_distance, y=chi_square_quantile))+ 
    geom_point() + stat_smooth(method=lm, level=0.95)+ # 95%的置信度
    labs(title="p-dimension Q-Q plot")+
    geom_abline(intercept=0,slope=1)
}
p_qqplot(x)
```

人群患代谢综合症的比例
```R
rate<-table(b)[2]/length(b)
```
人群患代谢综合征的比例为13.83%

 计算患代谢综合症的群体与没有患代谢综合症群体各类指标
将数据分为患代谢综合症的群体与没有患代谢综合症群体，其中x1为没患代谢综合征，x2为患代谢综合征群体
分析他们的体重指数、血压、血指、血糖等等指标的均值和置信区间分析差异
```R
xx<-split(x=x, f=b)
x1<-xx[[1]];x2<-xx[[2]]# 1没患代谢综合征 2患代谢综合征
# t检验
for (j in 1:ncol(x)){
  print(t.test(as.matrix(x1[,j]), as.matrix(x2[,j])))
}
```

- 体重指数
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = -8.5064, df = 44.759, p-value = 6.671e-11
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -5.033321 -3.105867
sample estimates:
mean of x mean of y 
 23.16380  27.23339 
```
他们的体重指数有显著性差异，患有代谢综合征的体重指数均值为27.23，要高于没有患代谢综合征的群体的均值23.16。
t统计量的置信区间
95 percent confidence interval:
 -5.033321 -3.105867

- 血糖
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = -5.1512, df = 27.678, p-value = 1.894e-05
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -2.0401547 -0.8788101
sample estimates:
mean of x mean of y 
 5.393210  6.852692 
```
他们的血糖有显著性差异，患有代谢综合征的血糖均值为6.85，要高于没有患代谢综合征的群体的均值5.39。
t统计量的置信区间
95 percent confidence interval:
 -2.0401547 -0.8788101

- 收缩压
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = -5.8399, df = 38.67, p-value = 8.895e-07
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -22.04958 -10.70256
sample estimates:
mean of x mean of y 
 115.2778  131.6538 
```
他们的收缩压有显著性差异，患有代谢综合征的收缩压均值为131.65，要高于没有患代谢综合征的群体的均值115.28。
t统计量的置信区间
95 percent confidence interval:
 -22.04958 -10.70256

- 舒张压
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = -7.0634, df = 38.079, p-value = 1.996e-08
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -17.054270  -9.456651
sample estimates:
mean of x mean of y 
 73.97531  87.23077 
```
他们的舒张压有显著性差异，患有代谢综合征的舒张压均值为87.23，要高于没有患代谢综合征的群体的均值73.98。
t统计量的置信区间
95 percent confidence interval:
 -17.054270  -9.456651

- 甘油三酯
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = -4.0374, df = 27.604, p-value = 0.0003881
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -2.7678213 -0.9037836
sample estimates:
mean of x mean of y 
 1.739198  3.575000 
```
他们的甘油三酯有显著性差异，患有代谢综合征的甘油三酯均值为3.58，要高于没有患代谢综合征的群体的均值1.74。
t统计量的置信区间
95 percent confidence interval:
 -2.7678213 -0.9037836

- 胆固醇
```R
	Welch Two Sample t-test

data:  as.matrix(x1[, j]) and as.matrix(x2[, j])
t = 3.4114, df = 36.825, p-value = 0.001583
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.08593825 0.33745207
sample estimates:
mean of x mean of y 
 1.625926  1.414231 
```
他们的胆固醇有显著性差异，患有代谢综合征的胆固醇均值为1.41，要低于没有患代谢综合征的群体的均值1.63。
t统计量的置信区间
95 percent confidence interval:
 0.08593825 0.33745207

（年龄）对患代谢综合症的影响分析，可以探讨下面的问题
（a）不同年龄组的人患代谢综合症的比例，检验显著性差异
（b）不同年龄组是否患代谢综合症群体各类指标的均值估计和置信区间或区域
- 将其以50岁为界限分为高年龄组和低年龄组
```R
a=matrix(data=NA,nrow=nrow(sample),ncol=1)
for (i in 1:nrow(sample)){
  #if(sample$age[i]<=30){a[i]="<=30"}
  if(sample$age[i]<=50){a[i]="l"}
  if(sample$age[i]>50){a[i]="h"}
  # if(sample$age[i]>70){a[i]=">70"}
}
t.test(b~a)
```
```R
	Welch Two Sample t-test

data:  b by a
t = 2.0781, df = 133.56, p-value = 0.03961
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.005342706 0.216234347
sample estimates:
mean in group h mean in group l 
     0.20253165      0.09174312 

```
高年龄段患代谢综合征的比例为20.25%，低年龄段患代谢综合征的比例为9.17%；高年龄段患代谢综合征比例高于低年龄段，并显著。
- 分为四个年龄组，<=30,30-50,50-70,>70
需要进行多元方差分析
```R
a=matrix(data=NA,nrow=nrow(sample),ncol=1)
for (i in 1:nrow(sample)){
  if(sample$age[i]<=30){a[i]="<=30"}
  if(sample$age[i]>30 & sample$age[i]<=50){a[i]="30-50"}
  if(sample$age[i]>50 & sample$age[i]<=70){a[i]="50-70"}
  if(sample$age[i]>70){a[i]=">70"}
}
library(agricolae)
data<-data.frame(as.numeric(b),a)
colnames(data)<-c("x","a")
model <- aov(x ~ a,data) # 先进行方差分析
print(summary(model))
out <- LSD.test(model, "a", p.adj="none" ) # 进行多重比较，不矫正P值
print(out$groups)# 查看每个组的label
plot(out) # 可视化展示
```
```R
> print(summary(model))
             Df Sum Sq Mean Sq F value Pr(>F)
a             3   0.65  0.2165   1.831  0.143
Residuals   184  21.75  0.1182               
> print(out$groups)# 查看每个组的label
               x groups
50-70 0.20779221      a
30-50 0.09375000      b
<=30  0.07692308      b
>70   0.00000000      b
```
虽然50-70年龄段患代谢综合征比例高于其它年龄段，但差异不显著

- 选取50-70年龄段的人群进行分析
```R
samp<-sample[a=="50-70",]
b<-b[a=="50-70"]
xsamp<-x[a=="50-70",]
xx<-split(x=xsamp, f=b)
x1<-xx[[1]];x2<-xx[[2]]# 1没患代谢综合征 2患代谢综合征
# t检验
for (j in 1:ncol(x)){
  print(t.test(as.matrix(x1[,j]), as.matrix(x2[,j])))
}
```
这与之前是否患代谢综合征各数据指标的比较是类似的

- 分析患代谢综合症的年龄差异
```R
> # t检验
> t.test(as.matrix(sample1$age), as.matrix(sample2$age))

	Welch Two Sample t-test

data:  as.matrix(sample1$age) and as.matrix(sample2$age)
t = -2.5479, df = 40.769, p-value = 0.0147
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -8.651375 -1.000097
sample estimates:
mean of x mean of y 
 46.63580  51.46154 

> # hotellingT检验
> HotellingsT2(as.matrix(sample1$age), as.matrix(sample2$age))

	Hotelling's two sample T2-test

data:  as.matrix(sample1$age) and as.matrix(sample2$age)
T.2 = 4.3001, df1 = 1, df2 = 186, p-value = 0.03949
alternative hypothesis: true location difference is not equal to c(0)
```
患代谢综合的平均年龄为51岁略高于未患代谢综合征群体的47岁，并有显著性差异

## 实验二
(II)数据ex2.1：给出了27,名糖尿病人血清总胆固醇(x1), 甘油（x2）,空腹胰岛素(x3),糖化血红蛋白(x4),空腹血糖(y)的测量值。
```R
library(readxl)
diabetes<-read_excel("ex2.1.xls")[,-1]
```
（1）试建立血糖(y)与其他指标的线性回归方程，并进行分析；
先观察一下数据
```
> head(diabetes)
# A tibble: 6 x 5
     x1    x2    x3    x4     y
  <dbl> <dbl> <dbl> <dbl> <dbl>
1  5.68  1.9   4.53   8.2  11.2
2  3.79  1.64  7.32   6.9   8.8
3  6.02  3.56  6.95  10.8  12.3
4  4.85  1.07  5.88   8.3  11.6
5  4.6   2.32  4.05   7.5  13.4
6  6.05  0.64  1.42  13.6  18.3
```
检查一下各变量之间的相关性
```R
> cor(diabetes)
           x1          x2          x3         x4          y
x1  1.0000000  0.63150583 -0.35479471  0.4152708  0.5585251
x2  0.6315058  1.00000000 -0.03863221  0.2189743  0.4585096
x3 -0.3547947 -0.03863221  1.00000000 -0.3297787 -0.5101213
x4  0.4152708  0.21897432 -0.32977870  1.0000000  0.6096420
y   0.5585251  0.45850963 -0.51012130  0.6096420  1.0000000
```
可以看到总胆固醇(x1), 甘油（x2），糖化血红蛋白(x4)与空腹血糖有明显的正相关性，而空腹胰岛素(x3)与空腹血糖有负相关性
```R
library(car)
scatterplotMatrix(diabetes,main='Scatter Plot Matrix')
```
![1027-3.png](C:\Users\Lenovo\Desktop\学习\多元统计\pr\pr3\1027-3.png)

从图中可以看出空腹血糖含量随着总胆固醇(x1), 甘油（x2），糖化血红蛋白(x4)的增加而增加，随空腹胰岛素(x3)增加而减少，其余各变量间似乎也有一定的相互依赖关系。

- 建立多元线性回归模型
自变量为血清总胆固醇(x1), 甘油（x2）,空腹胰岛素(x3),糖化血红蛋白(x4), 因变量为空腹血糖(y)
```R
fit<-lm(y~x1+x2+x3+x4,data=diabetes)
summary(fit)
```
```R
Call:
lm(formula = y ~ x1 + x2 + x3 + x4, data = diabetes)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.6268 -1.2004 -0.2276  1.5389  4.4467 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)  
(Intercept)   5.9433     2.8286   2.101   0.0473 *
x1            0.1424     0.3657   0.390   0.7006  
x2            0.3515     0.2042   1.721   0.0993 .
x3           -0.2706     0.1214  -2.229   0.0363 *
x4            0.6382     0.2433   2.623   0.0155 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.01 on 22 degrees of freedom
Multiple R-squared:  0.6008,	Adjusted R-squared:  0.5282 
F-statistic: 8.278 on 4 and 22 DF,  p-value: 0.0003121
```
回归系数的含义: 当其它自变量保持不变时，某一子变量增加时，因变量所增加的值。
对于p<0.05的显著性水平，空腹胰岛素(x3),糖化血红蛋白(x4)与y有显著的线性相关性，而另两个变量并没有显著的线性相关性。由Multiple R-squared: 0.6008可知，自变量解释了60.08%的方差。由Residual standard error: 2.01知估计标准误差为2.01，说明用以上四个预测变量来估计空腹血糖时，平均的估计误差为2.01。

这时我们将不显著相关的变量去掉重新进行线性回归的拟合
```R
fit<-lm(y~x3+x4,data=diabetes)
summary(fit)
```
```R
Call:
lm(formula = y ~ x3 + x4, data = diabetes)

Residuals:
   Min     1Q Median     3Q    Max 
-3.753 -1.675 -0.150  1.446  5.609 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)   
(Intercept)   6.3784     2.6705   2.388   0.0251 * 
x3           -0.2764     0.1244  -2.222   0.0360 * 
x4            0.7947     0.2505   3.173   0.0041 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.198 on 24 degrees of freedom
Multiple R-squared:  0.4788,	Adjusted R-squared:  0.4354 
F-statistic: 11.03 on 2 and 24 DF,  p-value: 0.0004014
```
```R
qqPlot(fit,id.method='identify',simulate = TRUE,
       labels=row.names(diabetes),main='Q-Q plot')
```
所有的点都在直线附近，并都落在置信区间内，这表明正态性假设符合得很完美。

（2）（x1, x2, x3, x4）是否服从多元正态？
使用之前已经写好的函数
```R
p_qqplot(diabetes[,1:4])
```
![1027-4.png](https://upload-images.jianshu.io/upload_images/17587926-5b99519163fc8361.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

（x1, x2, x3, x4）的QQ图近似落于一条过原点的直线，服从多元正态分布。

（x1,x2）与（x3,x4）是否相互独立？
```
n<-nrow(diabetes)
A<-(n-1)*var(diabetes[,1:4])
A11<-(n-1)*var(diabetes[,1:2])
A22<-(n-1)*var(diabetes[,3:4])
p=4;p1=2;p2=2;
b=n-3/2-(p^3-p1-p2)/(3*(p^2-p1-p2))
f=1/2*(p*(p+1)-p1*(p1+1)-p2*(p2+1))
e=-b*log(det(A)/(det(A11)*det(A22)))
chi.2<-qchisq(1-0.05, df=f, lower.tail = TRUE)
```
 e=7.710149，chi.2=9.487729，e<chi.2说明（x1,x2）与（x3,x4）不相互独立。