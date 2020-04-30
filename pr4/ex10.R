rm(list = ls())
library(MASS)
data<-read.csv("ex10.csv")
tr<-data[1:17,]
gp=tr$group
tr=tr[,1:4]
x1<-tr[gp=="1",];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-tr[gp=="2",];x2m<-apply(x2,2,mean);S2<-var(x2)
x3<-tr[gp=="3",];x3m<-apply(x3,2,mean);S3<-var(x3) 
dsq1<-function(x){
  return(t(x-x1m)%*%solve(S1)%*%(x-x1m))
}
dsq2<-function(x){
  return(t(x-x2m)%*%ginv(S2)%*%(x-x2m))# S2矩阵奇异
}
dsq3<-function(x){
  return(t(x-x3m)%*%solve(S3)%*%(x-x3m))
}
tt=data[,1:4]
d1<-apply(tt,1,dsq1)
d2<-apply(tt,1,dsq2)
d3<-apply(tt,1,dsq3)
re<-cbind(d1,d2,d3)
regp<-apply(re,1,function(a){which(a==min(a),arr.ind=TRUE)})
redata<-cbind(data,regp)

Dsq1<-function(x){
  return(t(x-x1m)%*%solve(S1)%*%(x-x1m)+log(abs(det(S1))))
}
Dsq2<-function(x){
  return(t(x-x2m)%*%ginv(S2)%*%(x-x2m)+log(abs(det(S2))))
}
Dsq3<-function(x){
  return(t(x-x3m)%*%solve(S3)%*%(x-x3m)+log(abs(det(S3))))
}
d1<-apply(tt,1,Dsq1)
d2<-apply(tt,1,Dsq2)
d3<-apply(tt,1,Dsq3)
regp<-apply(re,1,function(a){which(a==min(a),arr.ind=TRUE)})
redata<-cbind(data,regp)



