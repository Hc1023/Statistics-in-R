# 判别分析

## 理论
### 两总体协方差阵相等（但未知）时均值向量的检验
设$X_{\alpha}(\alpha=1,...,n)$为来自总体$X \sim N_p(\mu^{(1)},\Sigma)$的随机样本；$Y_{\alpha}(\alpha=1,...,m$为来自总体$X \sim N_p(\mu^{(2)},\Sigma)$的随机样本，且相互独立，$\Sigma$未知。
检验
$$H_0:\mu^{(1)}=\mu^{(2)},\ H_1:\mu^{(1)}\neq\mu^{(2)}$$取检验统计量
$$F=\frac{n+m-p-1}{(n+m-2)p}T^2$$

可计算样本均值$\bar X,\bar Y$，样本离差阵$A_1,A_2$，进一步计算得
$$D^2=(n+m-2)(\bar X-\bar Y)'(A_1+A_2)^{-1}(\bar X-\bar Y)$$
$$T^2=\frac{nm}{n+m}D^2$$ 
$$F_0=\frac{(n+m-2)-p+1}{(n+m-2)p}T^2\sim F(p,n+m-p-1)$$
给定显著水平$\alpha$，只需计算$p$值，即$p=P\{F\geq F_0\}$，若$p<\alpha$，即可否定$H_0$，即两组样本间存在显著性差异。在这种情况下，可能犯第一类错误，犯第一类错误的概率为$\alpha$.

### 马氏距离判别
设总体$G$为$m$元总体（考察m个指标），均值向量为$\mu=(\mu_1,\mu_2,...,\mu_m)'$，协方差阵为$\Sigma=(\sigma_{ij})_{m\times m}$，则样品$X=(x_1,x_2,...,x_m)'$与总体$G$的**马氏距离**定义为
$$d^2(X,G)=(X-\mu)'\Sigma^{-1}(X-\mu)$$
用估计量替代上式的距离和方差，均值用样本均值代，方差用组内协方差阵代替
$$S_i=\frac{1}{n_i-1}A_i$$其中
$$A_i=\Sigma_{t=1}^{n_i}(X_{(t)}^{(i)}-\bar{X^{i}})(X_{(t)}^{(i)}-\bar{X^{i}})'$$称为组间离差阵
当假定$\Sigma_1=\Sigma_2\overset{def}=\Sigma$时，反映分散性的协方差$\Sigma$的估计为
$$S=\frac{1}{n-k}\Sigma_{i=1}^{k}A_i$$
- $\Sigma_1=\Sigma_2$

$$\frac 1 2 (d_2^2(X)-d_1^2(X))=(X-\frac{1}{2}(\bar X^{(1)}+\bar X^{(2)}))'S^{-1}(\bar X^{(1)}-\bar X^{(2)})\overset{def}=W(X)$$
$W(X)>0$，判$X\in G_1$，否则判$X\in G_2$
所以计算马氏距离可改为对线性函数的计算
$Y_{i}(X)=\left(S^{-1} \bar{X}^{(i)}\right)^{\prime} X-\frac{1}{2}\left(\bar{X}^{(i)}\right)^{\prime} S^{-1} \bar{X}^{(i)}$

- $\Sigma_1\neq\Sigma_2$二次判别

$\begin{aligned} W(X)=& d^{2}\left(X, G_{2}\right)-d^{2}\left(X, G_{1}\right) \\=&\left(X-\mu^{(2)}\right)^{\prime} S_{2}^{-1}\left(X-\mu^{(2)}\right)-\left(X-\mu^{(1)}\right)^{\prime} S_{1}^{-1}\left(X-\mu^{(1)}\right) \\=& X^{\prime}\left(S_{2}^{-1}-S_{1}^{-1}\right) X-2\left(\mu^{(2)}-\mu^{(1)}\right)^{\prime} X \\ &+\left(\mu^{(2)^{\prime}} S_{2}^{-1} \mu^{(2)}-\mu^{(1)^{\prime}} S_{1}^{-1} \mu^{(1)}\right) \\ \cong & Z(X)+Z_{0} \end{aligned}$
其中$Z(X)$为$X$的二次函数，$Z_0$是一常数

注意：经检验不能拒绝协方差相等时，一般可采用线性判别，但仍可采用二次判别，二次判别的正确率高于线性判别；如果数据量很大，二次判别需要的时间远大于线性判别。

- 多总体距离判别

把$X$判给距离最近的那一个，设$i=l$时，若
$$d_l^2(X)=min_{i=1,...,k}\{d_i^2(X)\}$$则$X\in G_l$

### 贝叶斯判别及广义平方距离判别

当先验概率$\{q_i\}$总体分布的密度函数$\{f_i(x)\}$及损失函数$\{L(j|i)\}$给定时，贝叶斯判别的解满足：
$$G_l=\{y|h_l(y)<h_j(y),j\neq l,j=1,2,...,k\}$$其中
$h_l(y)=\Sigma_{i\neq l}q_if_i(y)L(l|i)$，进一步，如果$\Sigma_{i=1}^{k}\Sigma_{i\neq l}\int_{h_i(y)=h_l(y)}f_i(x)dx=0$，则$R_0\backslash\cup G_i$的概率测度为0。

在正态总体下，$G_i\sim N_m(\mu_i,\Sigma_i)$的贝叶斯判别法$L(j|i)=1-\delta_{ij}$，则贝叶斯判别的解为
$$G_l=\{y|q_lf_l(y)>q_jf_j(y),j\neq l\}$$

可以求得它的线性判别函数
$$ln\ q_if_i(X)=C_0+C_{i0}+X'C_i\overset{def}=C_0+Y_i(X)$$
其中$C_0$为与$G_i$无关的依赖于$X$的常数。
$Y_j(X)=C_{j0}+C_j'(X)$，并且称$Y_j(X)$为**线性判别函数**，称$C_j=S^{-1}\bar X^{(j)}$为**判别系数**，$C_{j0}=ln\ q_j-\frac 1 2[\bar X^{(j)}]'S^{-1}\bar X^{(j)}$为常数项。贝叶斯判别的解为最大的$Y_j$对应的下标。

这与用广义平方距离判别是一致的，定义样品$X$到总体$G_t(t=1,...,k)的$**广义平方距离**$D_t^2(X)$为
$$D_t^2(X)=d_t^2(X)+ln\ |S_t|-2ln\ |q_t|$$
将样本判给广义平方距离最小的对应总体。

我们之前所熟知的**后验概率**最大法判别
$$P(t|X)=P(X\in G_t|X已知)=\frac{q_tf_t(x)}{\Sigma_{i=1}^{k}q_if_i(x)}$$
其实就是贝叶斯判别法的一种特殊情况。

### Fisher判别
直观来讲，基本思想是**投影**将来自不同总体的样本间的距离拉大。
设从$m$元总体$G_t(t=1,2,...,k)$分别抽取样本如下

 $X_{(i)}^{(t)}=\left(x_{i 1}^{(t)}, \cdots, x_{i m}^{(t)}\right)\ \left(t=1, \cdots, k ; i=1, \cdots, n_{t}\right)$
投影后为
$G_{k}: \quad a^{\prime} X_{(1)}^{(k)}, \cdots, a^{\prime} X_{\left(n_{k}\right)}^{(k)}, \quad \bar{X}^{(k)}=\frac{1}{n_{k}} \sum_{j=1}^{n_{k}} X_{(j)}^{(k)}$

投影后为一元数据，其组间平方和为
$\begin{aligned} B_{0} &=\sum_{t=1}^{k} n_{t}\left(a^{\prime} \bar{X}^{(t)}-a^{\prime} \bar{X}\right)^{2} \\ &=a^{\prime}\left[\sum_{t=1}^{k} n_{t}\left(\bar{X}^{(t)}-\bar{X}\right)\left(\bar{X}^{(t)}-\bar{X}\right)^{\prime}\right] a \\ &=a^{\prime} B a \end{aligned}$
组内离差阵为
$\begin{aligned} A_{0} &=\sum_{t=1}^{k} 
 \sum_{j=1}^{n_{t}}\left(a^{\prime} \bar{X}_{(j)}^{(t)}-a^{\prime} \bar{X}^{(t)}\right)^{2} \\ &=a^{\prime} A a \end{aligned}$
使得投影后组间差异尽可能大，即使得$\frac{a^{\prime} B a}{a^{\prime} A a} \stackrel{d e f}{=} \Delta(a)$尽可能大。

不过按照老师的讲法，更常见的分母是各个总体的$\Sigma$直接求和。

## 实例及代码实现

（1）在天气预报中，常根据当天天气的湿温差(x1)和气温差（x2），来预测第二天是否下雨。试利用观测到的天气数据ex5.2（见附件），判断当今天测得(x1, x2)=(8.1,2.0) 或 (7.5, 3.5)时,明天的天气应判断为下雨还是不下雨？
- 读入数据，计算均值和协方差

```r
data1<-read.csv("ex5.2.csv")

gp<-data1$G # 1下雨 2不下雨
tr<-data1[,2:3]

x1<-tr[gp=="1",];x1m<-apply(x1,2,mean);S1<-var(x1) 
x2<-tr[gp=="2",];x2m<-apply(x2,2,mean);S2<-var(x2)
```
```r
> head(data1)
  G   x1   x2
1 1 -1.9  3.2
2 1 -6.9  0.4
3 1  5.2  2.0
4 1  5.0  2.5
5 1  7.3  0.0
6 1  6.8 12.7
```

- step1 是否有均值差异

```r
n=nrow(x1);m=nrow(x2);A1<-(n-1)*S1;A2<-(m-1)*S2;
Dsq=(n+m-2)%*%t(x1m-x2m)%*%solve(A1+A2)%*%(x1m-x2m)
p=2
Tsq=n*m/(n+m)*Dsq
FF=(n+m-p-1)/((n+m-2)*p)*Tsq
pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]
```
$p=0.00826$, 故在显著水平0.01时两总体的均值有显著差异， 即讨论两个总体的判别是有意义的。

- step2 马氏距离判别

```r
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
```
```r
> head(redata)
  G   x1   x2 regp
1 1 -1.9  3.2    1
2 1 -6.9  0.4    1
3 1  5.2  2.0    1
4 1  5.0  2.5    1
5 1  7.3  0.0    1
6 1  6.8 12.7    1
> table(redata$G==redata$regp)
FALSE  TRUE 
    2    18 
> table(redata$G==redata$regp)/nrow(redata)
FALSE  TRUE 
  0.1   0.9 
```
第一列为原来的类别，最后一列为数据回代后的判别结果，20个数据错判了2个，错判比例为0.1。

```r
> dsq1(c(8.1,2))<dsq2(c(8.1,2))
     [,1]
[1,] TRUE
> dsq1(c(7.5,3.5))<dsq2(c(7.5,3.5))
     [,1]
[1,] TRUE
```
两组数据均被判断为下雨。

（2）在研究沙基液化问题中，选取7个因子。现从已液化和未液化的地层中分别抽了12个和23个样本，其中1类表示已液化，2类表示未液化。试用距离判别法对原来的35个样本进行回代分类并分析误判情况（也就是对观测到的35个样本逐个进行判断，得到的判断结果与原先的分类是否一致，错误判断了多少个？错误判断所占的比例为多少？数据见附件ex5.4.）

这题与（1）同样的分析步骤，略写，直接上代码吧
```r
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
```
- 结果

```r
> pf(FF, df1=p, df2=n+m-p-1, lower.tail = F)# P[X > x]
             [,1]
[1,] 0.0001083051
> table(redata$GROUP==redata$regp)
FALSE  TRUE 
    1    34 
> f/nrow(redata)
[1] 0.02857143
```
35个数据中错判了1个，错判概率为0.0286。
