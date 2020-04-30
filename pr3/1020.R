rm(list = ls())
library(ggplot2)
library(ICSNP)
# (I)
# 数据预处理
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

# metabolism syndrome
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

sam<-split(x=sample, f=b)
sample1<-sam[[1]];sample2<-sam[[2]]# 1患代谢综合征 2没患代谢综合征
# 1.	分析患代谢综合症的年龄差异
# t检验
t.test(as.matrix(sample1$age), as.matrix(sample2$age))
# hotellingT检验
HotellingsT2(as.matrix(sample1$age), as.matrix(sample2$age))

# 检验相关数据的正态性和相关性
{
  xmean<-apply(x, 2, mean) # 样本均值
  xdispersion<-(nrow(x)-1)*var(x) # 样本离差阵
  xvar<-var(x) # 样本协方差
  xcor<-cor(x) # 样本相关阵
}

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

# 人群患代谢综合症的比例
rate<-table(b)[2]/length(b)
# (d)	计算患代谢综合症的群体与没有患代谢综合症群体各类指标
#（体重指数、血压、血指、血糖等等指标的均值和置信区间分析差异
xx<-split(x=x, f=b)
x1<-xx[[1]];x2<-xx[[2]]# 1 not患代谢综合征 2患代谢综合征
# t检验
for (j in 1:ncol(x)){
  print(t.test(as.matrix(x1[,j]), as.matrix(x2[,j])))
}

# （年龄）对患代谢综合症的影响分析，可以探讨下面的问题
# （a）不同年龄组的人患代谢综合症的比例，检验显著性差异
# （b）不同年龄组是否患代谢综合症群体各类指标的均值估计和置信区间或区域
a=matrix(data=NA,nrow=nrow(sample),ncol=1)
for (i in 1:nrow(sample)){
  #if(sample$age[i]<=30){a[i]="<=30"}
  if(sample$age[i]<=50){a[i]="l"}
  if(sample$age[i]>50){a[i]="h"}
  # if(sample$age[i]>70){a[i]=">70"}
}
t.test(b~a)
# 高年龄段患代谢综合征比例高于低年龄段，并显著

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
out <- LSD.test(model, "a", p.adj="none" ) # 进行多重比较，不矫正P值
print(out$groups)# 查看每个组的label
plot(out) # 可视化展示
# 虽然50-70年龄段患代谢综合征比例高于其它年龄段，但不显著

# 选取50-70年龄段的人群进行分析
samp<-sample[a=="50-70",]
b<-b[a=="50-70"]
xsamp<-x[a=="50-70",]
xx<-split(x=xsamp, f=b)
x1<-xx[[1]];x2<-xx[[2]]# 1没患代谢综合征 2患代谢综合征
# t检验
for (j in 1:ncol(x)){
  print(t.test(as.matrix(x1[,j]), as.matrix(x2[,j])))
}

# (II)数据ex2.1：给出了27,名糖尿病人血清总胆固醇(x1), 
# 甘油（x2）,空腹胰岛素(x3),糖化血红蛋白(x4),空腹血糖(y)的测量值。
library(readxl)
diabetes<-read_excel("ex2.1.xls")[,-1]
# （1）试建立血糖(y)与其他指标的线性回归方程，并进行分析；
head(diabetes)
# 检查一下各变量之间的相关性
cor(diabetes)
# 可以看到总胆固醇(x1), 甘油（x2），糖化血红蛋白(x4)与空腹血糖有明显的正相关性
# 而空腹胰岛素(x3)与空腹血糖有负相关性
library(car)
scatterplotMatrix(diabetes,main='Scatter Plot Matrix')
# 从图中可以看出空腹血糖含量随着总胆固醇(x1), 甘油（x2），糖化血红蛋白(x4)的增加而增加
# 随空腹胰岛素(x3)增加而减少
# 其余各变量间似乎也有一定的相互依赖关系
# 建立多元线性回归模型
# 自变量为血清总胆固醇(x1), 甘油（x2）,空腹胰岛素(x3),糖化血红蛋白(x4),
# 因变量为空腹血糖(y)
fit<-lm(y~x1+x2+x3+x4,data=diabetes)
summary(fit)
# 回归系数的含义: 当其它自变量保持不变时，某一子变量增加时，因变量所增加的值。
# 对于p<0.05的显著性水平，空腹胰岛素(x3),糖化血红蛋白(x4)与y有显著的线性相关性，
# 而另两个变量并没有显著的线性相关性
# 由Multiple R-squared: 0.6008可知，自变量解释了60.08%的销售价格的方差。
# 由Residual standard error: 2.01知估计标准误差为2.01，说明用以上四个预测变量来估计空腹血糖时，平均的估计误差为2.01。
fit<-lm(y~x3+x4,data=diabetes)
summary(fit)
# 这时我们将不显著相关的变量去掉重新进行线性回归的拟合
qqPlot(fit,id.method='identify',simulate = TRUE,
       labels=row.names(diabetes),main='Q-Q plot')
# 所有的点都在直线附近，并都落在置信区间内，这表明正态性假设符合得很完美


# （2）（x1, x2, x3, x4）是否服从多元正态？
p_qqplot(diabetes[,1:4])
# （x1, x2, x3, x4）的QQ图近似落于一条过原点的直线，服从多元正态


# （x1,x2）与（x3,x4）是否相互独立？
n<-nrow(diabetes)
A<-(n-1)*var(diabetes[,1:4])
A11<-(n-1)*var(diabetes[,1:2])
A22<-(n-1)*var(diabetes[,3:4])
p=4;p1=2;p2=2;
b=n-3/2-(p^3-p1-p2)/(3*(p^2-p1-p2))
f=1/2*(p*(p+1)-p1*(p1+1)-p2*(p2+1))
e=-b*log(det(A)/(det(A11)*det(A22)))
chi.2<-qchisq(1-0.05, df=f, lower.tail = TRUE)
# > e
# [1] 7.710149
# > chi.2
# [1] 9.487729
# e<chi.2,说明（x1,x2）与（x3,x4）不相互独立













