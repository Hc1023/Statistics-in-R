rm(list=ls())
# write.csv(NA, file = "MyData.csv")
mydata = read.csv("MyData.csv",row.names = 1)
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

library(ggplot2)
# contour
# library(lattice)
# parallelplot(~mydata, mydata, groups = rownames(mydata),
#          horizontal.axis = FALSE, scales = list(x = list(rot = 90)))

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
# some other contour plot
v <- ggplot(faithfuld, aes(waiting, eruptions, z = density))
v + geom_contour()
# look better maybe
v + geom_raster(aes(fill = density)) +
  geom_contour(colour = "white")

# radar plot
# Example from https://github.com/ricardo-bion/ggradar
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
#
library(andrews)
andrews(data,clr=7,ymax=4.5)
legend( "topright",legend = data$city,
        col=unique(data[,7]),lty=1,
        cex=0.8,ncol=2,bty="n")
# scatterplot-matrix
pairs(mydata,main="scatterplot-matrix",
      pch = 21, bg=c(1:m))
