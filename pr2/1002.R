rm(list = ls())
{
# library(readxl)
# sample<-read_excel("sample.xls")
sample<-read.csv("sample.csv",na.strings=c(""),stringsAsFactors=FALSE)
str(sample)
sample<-na.omit(sample[-1,])# 对na数据的处理
sample[,c(3:7,10:12)]<-as.numeric(unlist(sample[,c(3:7,10:12)]))
str(sample)
bmi<-sample$weight/(0.01*sample$height)^2
sample<-data.frame(sample,bmi)
x<-data.frame(sample$bmi,sample$FPG,sample$sbp,
              sample$dbp,sample$TG,sample$HDL.C)
}
# 1
cormatrix<-cor(x, method ="pearson" )
# 2
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

t.test(b~sample$gender)
t.test(b~as.numeric(unlist(sample$smoke)))
t.test(b~as.numeric(unlist(sample$drunk)))
# 3
a=matrix(data=NA,nrow=nrow(sample),ncol=1)

for (i in 1:nrow(sample)){
  if(sample$age[i]<=30){a[i]="<=30"}
  if(sample$age[i]>30 & sample$age[i]<=50){a[i]="30-50"}
  if(sample$age[i]>50 & sample$age[i]<=70){a[i]="50-70"}
  if(sample$age[i]>70){a[i]=">70"}
}

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
# 4
xmean<-apply(x, 2, mean) # 样本均值
xdispersion<-(ncol(x)-1)*var(x) # 样本离差阵
xvar<-var(x) # 样本协方差
xcor<-cor(x) # 样本相关阵
# 5
print(unlist(shapiro.test(unlist(x[2])))[1:2])
hist(unlist(x[2]))
print(unlist(shapiro.test(unlist(x[3])))[1:2])
hist(unlist(x[3]))

