library(MASS)
attach(UScereal)
shelf=factor(shelf) #转化为因子变量
y=cbind(calories,fat,sugars) #将因变量合并成一个矩阵 
aggregate(y,by=list(shelf),FUN=mean)  #求各类均值
fit=manova(y~shelf) 
summary(fit)  #多元方差分析
summary.aov(fit)  #对每个变量做单因素方差分析

rm(list = ls())
data<-read.csv("ex9.csv")
X=data[1:7,]
Y=data[8:14,]
Xmean<-apply(X,2,mean)
Ymean<-apply(Y,2,mean)
n=nrow(X)
m=nrow(Y)
A1<-(n-1)*var(X) 
A2<-(m-1)*var(Y)
Dsq=(n+m-2)%*%t(Xmean-Ymean)%*%solve(A1+A2)%*%(Xmean-Ymean)
p=3
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]


S=1/(n+m-2)*(A1+A2)
C1=solve(S)%*%Xmean
C10=log(0.5)-0.5*t(Xmean)%*%solve(S)%*%(Xmean)
C2=solve(S)%*%Ymean
C20=log(0.5)-0.5*t(Ymean)%*%solve(S)%*%(Ymean)
Y1=apply(data,1,function(x){x%*%C1+C10})
Y2=apply(data,1,function(x){x%*%C2+C20})
Y1>Y2

x=c(2.95,2.15,1.54)
x%*%C1+C10>x%*%C2+C20
