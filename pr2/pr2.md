# 参数估计和假设检验

实验内容 ：附表中的数据sample.xls进行分析。记X1=BMI, X2=FPG, X3=SBP, X4=DBP, X5=TG, X6=HDL-C，并构成一个向量。X=(X1, X2, X3, X4, X5, X6)
- 数据预处理
```R
rm(list = ls())
sample<-read.csv("sample.csv",na.strings=c(""),stringsAsFactors=FALSE)
str(sample)
sample<-na.omit(sample[-1,])# 对na数据的处理
sample[,c(3:7,10:12)]<-as.numeric(unlist(sample[,c(3:7,10:12)]))
str(sample)
bmi<-sample$weight/(0.01*sample$height)^2
sample<-data.frame(sample,bmi)
x<-data.frame(sample$bmi,sample$FPG,sample$sbp,
              sample$dbp,sample$TG,sample$HDL.C)
```
- 如何处理缺失数据
```R
sample<-na.omit(sample[-1,])# 对na数据的处理
```
删除了有数据缺失的病人样本
- X=(X1, X2, X3, X4, X5, X6)
X1=BMI, X2=FPG, X3=SBP, X4=DBP, X5=TG, X6=HDL-C

![1007.png](https://upload-images.jianshu.io/upload_images/17587926-0fa655e3f362765c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

(1)	分析X各变量之间的相关性？
```R
cormatrix<-cor(sample[,c(3:7,10:12)], method ="pearson" )
```
![1007-2.png](https://upload-images.jianshu.io/upload_images/17587926-b974105604de3d64.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

可以发现一些变量如sbp与dbp之间有明显的正相关性

(2)	分析患代谢综合症的比例有没有性别差异，与吸烟或喝酒是否有关？
女性F，男性M；
戒烟以及无吸烟记为0，吸烟记为1；
不喝酒记为0，喝酒记为1
- 判断是否患代谢综合征
 中华医学会糖尿病学分会（CDS）建议代谢综合征的诊断标准: 具备以下4项中的3项及以上即为代谢综合症:
	(1)	 超重:  　BMI>= 25.0  Kg/M^2(体重/身高平方) ;
	(2)	高血糖:FPG>= 6.1mmol/L(110mg/dl)或 2hPG>=7.8  mmol/L(140mg/dl), 或已确诊糖尿病并治疗者;
	(3)	高血压:　收缩压 SBP >=e 140 mmHg  或 舒张压DBP>= 90mmHg,  或已确诊高血压并治疗者;
	(4)	空腹血:  甘油三脂 TG>=1.7 mmol/L(110mg/dl)  或   HDL-C  <0.9  mmol/L( 35 mg/dl)（男）,  <1.0 mmol/L( 39  mg/dl) （女）.

```R
{
  bmi<-sample$weight/(0.01*sample$height)^2
  b1<-as.numeric(bmi>=25)
  b2<-as.numeric(sample$FPG>=6.1)
  b3<-as.numeric(sample$sbp>=140|sample$dbp>=90)
  b4_1<-sample$TG>=1.7
  b4_2<-sample$gender == "F" & sample$HDL.C<1
  b4_3<-sample$gender == "M" & sample$HDL.C<0.9
  b4<-as.numeric(b4_1 | b4_2 | b4_3)
  b<-(b1+b2+b3+b4>=3)
}
```
- 是否有性别差异
```R
> t.test(b~sample$gender)

	Welch Two Sample t-test

data:  b by sample$gender
t = -4.8393, df = 153.34, p-value = 3.148e-06
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.2705478 -0.1136899
sample estimates:
mean in group F mean in group M 
     0.01449275      0.20661157 
```
通过t检验发现有显著的性别差异，男性患有代谢综合征的概率高于女性

- 与吸烟是否有关
```R
> t.test(b~as.numeric(unlist(sample$smoke)))

	Welch Two Sample t-test

data:  b by as.numeric(unlist(sample$smoke))
t = -2.3243, df = 85.197, p-value = 0.02249
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.26166093 -0.02039035
sample estimates:
mean in group 0 mean in group 1 
     0.09230769      0.23333333 
```
通过t检验发现有较显著的差异，抽烟群体患代谢综合征的概率大于非抽烟群体

```R
> t.test(b~as.numeric(unlist(sample$drunk)))

	Welch Two Sample t-test

data:  b by as.numeric(unlist(sample$drunk))
t = -3.0949, df = 108.29, p-value = 0.002506
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.27958855 -0.06128102
sample estimates:
mean in group 0 mean in group 1 
     0.06956522      0.24000000 
```
通过t检验发现有较显著的差异，喝酒群体患代谢综合征的概率大于不喝酒群体

(3)	分年龄(小于等于30， 30~50, 50~70, 70以上)，分析X中的各个指标是否有年龄上的差异？
- 年龄分组
```R
a=matrix(data=NA,nrow=nrow(sample),ncol=1)

for (i in 1:nrow(sample)){
  if(sample$age[i]<=30){a[i]="<=30"}
  if(sample$age[i]>30 & sample$age[i]<=50){a[i]="30-50"}
  if(sample$age[i]>50 & sample$age[i]<=70){a[i]="50-70"}
  if(sample$age[i]>70){a[i]=">70"}
}
```
- 分析差异
```R
library(agricolae)
for(i in 1:6){
  data<-data.frame(x[,i],a)
  colnames(data)<-c("x","a")
  model <- aov(x ~ a,data) # 先进行方差分析
  print(summary(model))
  out <- LSD.test(model, "a", p.adj="none" ) # 进行多重比较，不矫正P值
  print(out$groups)# 查看每个组的label
  plot(out) # 可视化展示
}
```
- X1=BMI
```R
             Df Sum Sq Mean Sq F value Pr(>F)
a             3   2957   985.6   0.595  0.619
Residuals   186 308277  1657.4               
             x groups
50-70 31.43127      a
30-50 23.53183      a
<=30  22.85015      a
>70   21.98399      a
```
无显著差异

- X2=FPG
```R
             Df Sum Sq Mean Sq F value  Pr(>F)   
a             3  12.69   4.229   4.155 0.00706 **
Residuals   186 189.31   1.018                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
             x groups
50-70 5.884615      a
30-50 5.440103      b
>70   5.285000      b
<=30  5.066923      b
```
有较显著的差异，50-70年龄段会略高一些

- X3=SBP
```R
             Df Sum Sq Mean Sq F value   Pr(>F)    
a             3   5957  1985.8   8.162 3.92e-05 ***
Residuals   186  45254   243.3                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
             x groups
>70   157.0000      a
50-70 121.3846      b
30-50 114.7423      c
<=30  109.0769      c
```
年龄组越高，SBP越高

- X4=DBP
```R
             Df Sum Sq Mean Sq F value  Pr(>F)   
a             3   1728   576.1   4.875 0.00275 **
Residuals   186  21979   118.2                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
             x groups
>70   80.50000      a
50-70 78.80769      a
30-50 74.24742     ab
<=30  68.23077     ab
```
较SBP,DBP对年龄组的差异并不明显，但也有差异显著性，年龄高DBP高

- X5=TG
```R
             Df Sum Sq Mean Sq F value Pr(>F)
a             3   12.8   4.252   1.737  0.161
Residuals   186  455.3   2.448               
             x groups
30-50 2.155979      a
50-70 1.957308     ab
<=30  1.238462      b
>70   0.835000      b
```
无显著差异

- X6=HDL-C
```R
             Df Sum Sq Mean Sq F value Pr(>F)
a             3  0.208 0.06941   0.608  0.611
Residuals   186 21.232 0.11415               
             x groups
>70   1.875000      a
<=30  1.624615      a
50-70 1.606410      a
30-50 1.577423      a
```
无显著差异

```R
xmean<-apply(x, 2, mean) # 样本均值
xdispersion<-(ncol(x)-1)*var(x) # 样本离差阵
xvar<-var(x) # 样本协方差
xcor<-cor(x) # 样本相关阵
```

(4)	计算X样本均值、样本离差阵、样本协方差和样本相关阵

- 样本均值向量$\bar X$:
$$\bar X= \frac 1 n\Sigma_{i=1}^n X_{(i)}=(\bar{x}_1,\dots,\bar{x}_p)'=\frac{1}{n}X'\bf 1_n$$
其中
$$\bar x_i=\frac 1 n\Sigma_{\alpha=1}^nx_{\alpha i}\ (i=1,2,\dots,p).$$

```R
> xmean
  sample.bmi   sample.FPG   sample.sbp   sample.dbp    sample.TG 
   26.711822     5.595421   117.526316    75.773684     1.997737 
sample.HDL.C 
    1.595684
```

- 样本离差阵（又称交叉乘积阵）$A$:
$$A=\Sigma_{\alpha=1}^n(X_{(\alpha)}-\bar X)(X_{(\alpha)}-\bar X)'=X'X-n\bar X\bar X'$$
$$=X'[I_n-\frac 1 n \bf1_n 1_n']\it X\rm\overset{def}=(a_{ij})_{p}$$
其中
$$ a_{ij}=\Sigma_{\alpha=1}^n(x_{\alpha i}-\bar x_i)(x_{\alpha j}-\bar x_j)\ (i,j=1,2,\dots,p)$$

![1017.png](https://upload-images.jianshu.io/upload_images/17587926-40f9088022213d1d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 样本协方差$S$:
$$S=\frac{1}{n-1}A=(s_{ij})_{p\times p}\ \ (\ or\ S^*=\frac1n A )$$
其中
$$s_{ii}=\frac{1}{n-1}\Sigma_{a=1}^n(x_{ai}-\bar x_i)^2\ (i=1,2,\dots,p)$$
称为变量$X_i$的样本方差；样本方差的平方根$\sqrt{s_{ii}}$称为变量$X_i$的样本标准差。

![1007-3.png](https://upload-images.jianshu.io/upload_images/17587926-99c415aeda38c7b2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 样本相关阵$R$:
$$R=(r_{ij})_{p\times p}$$
其中
$$r_{ij}=\frac{s_{ij}}{\sqrt{s_{ii}}\sqrt{s_{jj}}}\overset{or}=\frac{a_{ij}}{\sqrt{a_{ii}}\sqrt{a_{jj}}}\ (i,j=1,2,\dots,p)$$

![1007-2.png](https://upload-images.jianshu.io/upload_images/17587926-1014da80770e1b6f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

(5)	分析X2, X3是否服从正态分布?
```R
> print(unlist(shapiro.test(unlist(x[2])))[1:2])
           statistic.W                p.value 
    "0.76565518467017" "4.00113250383294e-16" 
> print(unlist(shapiro.test(unlist(x[3])))[1:2])
        statistic.W             p.value 
 "0.98209695047987" "0.015692440430348" 
```
X2服从正态分布，并具有较高显著性。
X3服从正态分布的显著性不高。