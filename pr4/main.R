rm(list = ls())
data1<-read.csv("ex5.2.csv")
head(data1)
gp<-data1$G # 1下雨 2不下雨
tr<-data1[,2:3]

x1<-tr[gp=="1",];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-tr[gp=="2",];x2m<-apply(x2,2,mean);S2<-var(x2)
# 是否有均值差异
n=nrow(x1);m=nrow(x2);A1<-(n-1)*S1;A2<-(m-1)*S2;
Dsq=(n+m-2)%*%t(x1m-x2m)%*%solve(A1+A2)%*%(x1m-x2m)
p=2
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]
# 故在显著水平0.01时两总体的均值有显著差异，
# 即讨论两个总体的判别是有意义的

dsq1<-function(x){
  return(t(x-x1m)%*%solve(S1)%*%(x-x1m))
}
dsq2<-function(x){
  return(t(x-x2m)%*%solve(S2)%*%(x-x2m))
}
d1<-apply(tr,1,dsq1)
d2<-apply(tr,1,dsq2)
re<-cbind(d1,d2)
regp<-apply(re,1,function(a){which(a==min(a),arr.ind=TRUE)})
redata<-cbind(data1,regp)
head(redata)
table(redata$G==redata$regp)
dsq1(c(8.1,2))<dsq2(c(8.1,2))
dsq1(c(7.5,3.5))<dsq2(c(7.5,3.5))

data2<-read.csv("ex5.4.csv")
gp<-data2$GROUP # 1类表示已液化，2类表示未液化
tr<-data2[,-1]

x1<-tr[gp=="1",];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-tr[gp=="2",];x2m<-apply(x2,2,mean);S2<-var(x2)
# 是否有均值差异
n=nrow(x1);m=nrow(x2);A1<-(n-1)*S1;A2<-(m-1)*S2;
Dsq=(n+m-2)%*%t(x1m-x2m)%*%solve(A1+A2)%*%(x1m-x2m)
p=ncol(tr)
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]
# 故在显著水平0.01时两总体的均值有显著差异，
# 即讨论两个总体的判别是有意义的

d1<-apply(tr,1,dsq1)
d2<-apply(tr,1,dsq2)
re<-cbind(d1,d2)
regp<-apply(re,1,function(a){which(a==min(a),arr.ind=TRUE)})
redata<-cbind(data2,regp)
f<-table(redata$GROUP==redata$regp)[[1]]
f/nrow(redata)







