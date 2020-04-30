# 聚类分析和主成分分析

## 聚类分析
（1）现有16种饮料的热量、咖啡因含量、钠含量和价格的数据（见ex4.2),根据这4个变量对16种饮料进行聚类

```r
> head(mydata)
  热量x1 咖啡因含量x2 钠含量x3 价格x4
1  207.2          3.3     15.5    2.8
2   36.8          5.9     12.9    3.3
3   72.2          7.3      8.2    2.4
4   36.7          0.4     10.5    4.0
5  121.7          4.1      9.2    3.5
6   89.1          4.0     10.2    3.3
```
### 系统聚类

这里展示的是离差平方和法（WARD）进行系统聚类。它基于方差分析的思想，同类样品之间的离差平方和应当较小，不同类之间的离差平方和应当较大
$G_t$中样品的离差平方和为
$$W_t=\Sigma_{i=1}^{n_t}(X_{(i)}^{(t)}-\bar{X}^{(t)})' (X_{(i)}^{(t)}-\bar{X}^{(t)})$$
$k$个类的总组内离差平方和为
$$W=\Sigma_{t=1}^k W_t$$
当$k$固定式要选择使W达到极小的分类
代码一（比较丑）：

```r
library(cluster)
# Compute with agnes (make sure you have the package cluster)
hc1 <- agnes(mydata, method = "ward")
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram of agnes")
rect.hclust(hc1, k = 3)
```

![1213-2.png](https://upload-images.jianshu.io/upload_images/17587926-bd59b53da1d6e806.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

代码二：

```r
#用欧式距离代表样本之间两两相似性
result <- dist(mydata, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")
fviz_dend(result_hc, k = 3, 
          cex = 0.5, 
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect = TRUE          
)
```

![1213-4.png](https://upload-images.jianshu.io/upload_images/17587926-ffec691389ab3e6d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

与之前代码的结果是一致的，但是图比较好看嘿嘿嘿

经WARD系统聚类法，我们得到了聚类结果的谱系图，从中我们可以看到饮料分为三类的结果是{12,14,2,4,13,3,8,11},{1,10},{7,6,9,16,5,15}.

### 动态聚类

这里采用K-means聚类方法展示动态聚类法。
Step1 规定样品间的距离人为定出三个数：$K$（分类数），$C$（类间距离最小值）和$R$（类内距离最大值）；取前K个样品为**凝聚点**
Step2 计算这$K$个凝聚点两两之间的距离，小于$C$则合并这两个凝聚点，用这两个点的重心作为新的凝聚点重复Step2至所有凝聚点之间的距离$\geq C$.
Step3 剩下$n-K$个样品逐个归类，若与凝聚点距离$>R$则作为新的凝聚点否则归入距离最近的凝聚点所在的类。重新计算每类重心作为新凝聚点。重复Step2。
Step4 所有样品重复Step3且分类不一致是要重新计算重心。

```r
library(factoextra) 
b<-scale(mydata) 
set.seed(123)  
fviz_nbclust(b,kmeans,method="wss")+geom_vline(xintercept=3,linetype=2) 
```
使用组内平方误差和确定最佳聚类个数，预判分为三类或四类，这里分为三类效果较好（可看后面的聚类图）。

![1213-1.png](https://upload-images.jianshu.io/upload_images/17587926-3e8e0381ab2dfc2f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

```r
res<-kmeans(b,3)
res1<-cbind(mydata,res$cluster)
fviz_cluster(res,data=mydata[,1:ncol(mydata)-1])
```
![1213.png](https://upload-images.jianshu.io/upload_images/17587926-1d32fa34ac8b255c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

使用K-means动态聚类方法与Ward系统聚类法的聚类结果并不相同，K-means动态聚类方法将饮料分为三类的结果是{1,8,12,14},{4,9,10,16},{3,5,6,7,11,12,13,15}

（2）中国31个城市2011年的空气质量数据（见ex4.3）,根据这个数据对31个城市进行聚类分析。
```r
> head(mydata)
         可吸入颗粒物.PM10. 二氧化硫.SO2. 二氧化氮.NO2.
北京             0.11329589    0.02810959    0.05582740
天津             0.09287671    0.04203562    0.03830685
石家庄           0.09887397    0.05225205    0.04118356
太原             0.08433151    0.06391233    0.02293151
呼和浩特         0.07576164    0.05413973    0.03944110
沈阳             0.09648219    0.05859452    0.03272329
         空气质量达到及好于二级的天数.天.
北京                                  286
天津                                  320
石家庄                                320
太原                                  308
呼和浩特                              347
沈阳                                  332
         空气质量达到二级以上天数占全年比重...
北京                                  78.35616
天津                                  87.67123
石家庄                                87.67123
太原                                  84.38356
呼和浩特                              95.06849
沈阳                                  90.95890
```
这次先进行K-means动态聚类法。
首先我根据动态聚类法判定类别数，$k=5$时出现拐点，分成五类合适。

![1213_2.png](https://upload-images.jianshu.io/upload_images/17587926-e0790f72272c7363.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

绘制动态聚类结果

![1213_2-2.png](https://upload-images.jianshu.io/upload_images/17587926-8047b7b4f808fc74.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

根据可吸入颗粒物(PM10)、二氧化硫(SO2)、二氧化氮(NO2)、空气质量达到及好于二级的天数(天)、空气质量达到二级以上天数占全年比重(%)这几个变量，用K-means动态聚类法将城市分成五类，结果如下：

```r
> rownames(res1)[res1$`res$cluster`=="1"]
 [1] "天津"     "石家庄"   "太原"     "呼和浩特" "沈阳"     "南昌"    
 [7] "济南"     "重庆"     "贵阳"     "西宁"     "银川"
> rownames(res1)[res1$`res$cluster`=="2"]
[1] "长春" "上海" "杭州" "长沙" "广州" "南宁" "昆明"
> rownames(res1)[res1$`res$cluster`=="3"]
[1] "福州" "海口" "拉萨"
> rownames(res1)[res1$`res$cluster`=="4"]
[1] "兰州"     "乌鲁木齐"
> rownames(res1)[res1$`res$cluster`=="5"]
[1] "北京"   "哈尔滨" "南京"   "合肥"   "郑州"   "武汉"   "成都"  
[8] "西安"  
```

- WARD系统聚类法

![1213_2-3.png](https://upload-images.jianshu.io/upload_images/17587926-5920a6891e11f93d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

> 我精挑细选的配色，绘图的关键在配色（√）

使用WARD系统聚类法的聚类结果为
第一类：{"福州" "广州" "拉萨" "海口" "昆明"}
第二类：{"南宁" "贵阳" "长春" "呼和浩特" "南昌" "沈阳" "杭州"  "上海"  "长沙"}
第三类：{"太原" "合肥" "武汉" "西安" "西宁" "郑州" "哈尔滨" "南京" "天津" "石家庄" "济南" "重庆" "成都"}
第四类：{"兰州"}
第五类：{"北京" "乌鲁木齐"}

下面我们通过**热图**的方法发现类并确定类的个数。
用result表示两两样本间的欧式距离矩阵，距离越近相关性越高。

```r
result <- dist(mydata, method = "euclidean")
```

### 通过热图发现类

```r
heatmap(as.matrix(result), labRow = F, labCol = F)
```
在图中可以看出颜色越深，距离越近。通过观察感觉分成五类或六类OK，分成两类也似有道理但总觉得两类太少了点。

![heatmap.png](https://upload-images.jianshu.io/upload_images/17587926-fc97218f285fcdc3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

确定类的个数为五类

```r
result_hc <- hclust(d = result, method = "ward.D2")
re <- cutree(result_hc, k = 5)
```
采用多维标定，用第一和第二主成分，表示原本分类。

```r
mds = cmdscale(result, k = 2, eig = T)
X <- mds$points[, 1]
Y <- mds$points[, 2]
p = ggplot(data.frame(X, Y ), aes(X, Y ))
library(ggplot2)
p + geom_point(size = 3, alpha = 0.8, aes(colour =factor(re)))
```

![cluster.png](https://upload-images.jianshu.io/upload_images/17587926-ea360101d371601e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)


### 系统聚类法和动态聚类法的优缺点

系统聚类法：
优点：可解释性好，可以考虑先取K比较大的K-means后，合并阶段用系统聚类也许能产生更高质量的聚类。
缺点：时间复杂度高。

动态聚类法：
优点：适用于大样本的Q型聚类分析。大型数据一般较集中，异常值影响较弱。且算法简单高效，时间复杂度低。
缺点：依赖初始$K$个点的选取或说依赖于初始划分，容易落入局部最优。对噪声和离群值敏感。

动态聚类法的改进方法：为了检验聚类的稳定性，可用一个新的初始分类重新检验整个聚类算法。如最终分类与原来一样，则不必再进行计算；否则，须另行考虑聚类算法。

## 主成分分析

- 由于变量个数太多，且彼此有相关性，从而数据信息重叠。
- 当变量较多，在高维空间研究样本分布规律较复杂

于是我们希望，用较少的综合变量代替原来较多变量，又能尽可能多地反映原来数据的信息，并且彼此之间互不相关。

叮！这就孕育了**主成分分析**！

设$X=(X_1,X_2,...,X_p)'$为$p$维随机向量。称$Z_i=a_i'X$为$X$的第$i$主成分$(i=1,2,...,p)$，如果：
- $a_i'a_i=1\ (i=1,2,...,p)$;
- when $i>1$, $a_i'\Sigma a_j=0\ (j=1,....,i-1)$;
- $Var(Z_i)=max_{a_i'a_i=1,a_i'\Sigma a_j=0\ (j=1,....,i-1)}Var(\alpha' X)$

易见，$Var(Z_i)$越大表明$Z_i$蕴含信息越多，此时可见$a_i$需要有限制否则方差可以趋于无穷，常用限制为$a_i'a_i=1$. 而我们不想信息重复，于是又有$Cov(Z_i,Z_j)=a_i'\Sigma a_j=0\ (j\neq i)$. 这时我们不难理解上述对主成分的定义了。

一般地，求$X$的第$i$主成分可通过求$\Sigma$的第$i$大特征值所对应的特征向量得到。

*Thm* 设$X=(X_1,X_2,...,X_p)'$是$p$维随机向量，且$D(X)=\Sigma$，$\Sigma$的特征值为$\lambda_1\geq\lambda_2\geq...\geq\lambda_p\geq 0$，$a_1,a_2,...,a_p$为相应的单位正交特征向量，则$X$的第$i$主成分为
$$Z_i=a_i'X\ (i=1,2,...,p)$$

下面这张图就形象地展现了如何利用主成分分析将二维降至一维。
![pca1.png](https://upload-images.jianshu.io/upload_images/17587926-723fc28c4272912c.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

注意，当数据集中的变量高度相关时，PCA方法特别有用。 相关性表明数据中存在冗余。 由于这种冗余，PCA可用于将原始变量减少为较少数量的新变量（=主成分），从而解释了原始变量中的大多数方差。

![pca_red.png](https://upload-images.jianshu.io/upload_images/17587926-9f5488203c4897a4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

以下解释贡献率的概念，贡献率的大小又能说是解释了总体方差的多少。
我们称$\lambda_k / \Sigma_{i=1}^{p}\lambda_i$为主成分$Z_k$的**贡献率**；又称$\Sigma_{k=1}^{m}\lambda_k / \Sigma_{i=1}^{p}\lambda_i$为主成分$Z_1,...,Z_m(m<p)$的**累计贡献率**。
（3'）我国2010你那各地区城镇居民家庭平均每人全年消费数据如ex6.7所示，这些数据指标分别从食品（x1）,衣着，居住，医疗，交通，通信，教育，家政和耐用消费品来描述消费。试对该数据进行主成分分析。

```r
> head(data)
          食品    衣着    居住    医疗 交通通讯   教育 家庭服务 耐用消费品
北京   5561.54 1571.74 1286.32 1563.10  2293.23 809.25    84.71     548.55
天津   5005.09 1153.66 1528.28 1220.92  1567.87 715.24    45.50     467.75
河北   3155.40 1137.22 1097.41  808.88  1062.31 386.60    28.84     305.70
山西   2974.76 1137.71 1250.87  769.79   931.33 570.79    35.38     259.05
内蒙古 3553.48 1616.56 1028.19  869.71  1191.70 568.35    30.49     307.92
辽宁   4378.14 1187.41 1270.95  913.13  1295.70 670.13    30.40     235.46
```

```r
pr<-princomp(data,cor=TRUE)
screeplot(pr,type="lines")
```

![1213_3-1.png](https://upload-images.jianshu.io/upload_images/17587926-c8508c45ae27a328.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

根据**碎石图**转折点在2个主成分或3个主成分

主成分分析结果：
```r
> summary(pr,loadings=TRUE)
Importance of components:
                          Comp.1    Comp.2     Comp.3    Comp.4     Comp.5
Standard deviation     2.3213318 1.1100881 0.72943408 0.5464987 0.47643779
Proportion of Variance 0.6735727 0.1540369 0.06650926 0.0373326 0.02837412
Cumulative Proportion  0.6735727 0.8276096 0.89411886 0.9314515 0.95982558
                          Comp.6     Comp.7      Comp.8
Standard deviation     0.4351019 0.28334611 0.227588880
Proportion of Variance 0.0236642 0.01003563 0.006474587
Cumulative Proportion  0.9834898 0.99352541 1.000000000

Loadings:
           Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7 Comp.8
食品        0.358  0.396  0.158  0.288  0.503         0.282  0.522
衣着        0.257 -0.536  0.703        -0.130 -0.336         0.135
居住        0.374        -0.412 -0.570 -0.112 -0.512  0.224  0.198
医疗        0.275 -0.599 -0.336         0.600  0.148 -0.248       
交通通讯    0.393  0.292  0.137  0.120  0.166 -0.233  0.114 -0.795
教育        0.386         0.195 -0.466 -0.178  0.729  0.168       
家庭服务    0.396  0.264               -0.211        -0.837  0.152
耐用消费品  0.361 -0.205 -0.373  0.599 -0.503  0.114  0.251  
```

可以发现当选用前三个主成分时，已经解释了接近90%的方差，所以我认为选取三个主成分进行降维是合理的。

### FactorMineR

下面我们用FactoMineR包来对同样的数据作PCA，它的展示效果更好。

- 计算变量的主成分

```r
res.pca <- PCA(data, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
```

- 观察特征值，它代表了每个主成分的方差

```
> eig.val
      eigenvalue variance.percent cumulative.variance.percent
Dim.1 5.38858121       67.3572651                    67.35727
Dim.2 1.23229559       15.4036948                    82.76096
Dim.3 0.53207408        6.6509260                    89.41189
Dim.4 0.29866081        3.7332601                    93.14515
Dim.5 0.22699297        2.8374121                    95.98256
Dim.6 0.18931363        2.3664203                    98.34898
Dim.7 0.08028502        1.0035627                    99.35254
Dim.8 0.05179670        0.6474587                   100.00000
```

- 观察碎石图

```r
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
```

![scree.png](https://upload-images.jianshu.io/upload_images/17587926-216a269504218614.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

从PCA输出中提取变量结果的一种简单方法是使用函数get_pca_var（）[factoextra package]。 该函数提供了一个矩阵列表，其中包含活动变量的所有结果（坐标，变量与轴之间的相关性，余弦平方和贡献）

```r
> var <- get_pca_var(res.pca)
> var
Principal Component Analysis Results for variables
 ===================================================
  Name       Description                                    
1 "$coord"   "Coordinates for the variables"                
2 "$cor"     "Correlations between variables and dimensions"
3 "$cos2"    "Cos2 for the variables"                       
4 "$contrib" "contributions of the variables"
```
```r
> # Coordinates
> head(var$coord)
             Dim.1        Dim.2      Dim.3        Dim.4       Dim.5
食品     0.8309549 -0.439440979  0.1151560  0.157267242  0.23942030
衣着     0.5968496  0.594904015  0.5130757  0.013283797 -0.06211260
居住     0.8674876 -0.038553622 -0.3007570 -0.311292972 -0.05351647
医疗     0.6377140  0.664657743 -0.2453518  0.002703152  0.28597319
交通通讯 0.9113413 -0.323711002  0.1000397  0.065838532  0.07908669
教育     0.8970832 -0.007529241  0.1420381 -0.254423709 -0.08501927
> # Cos2: quality on the factore map
> head(var$cos2)
             Dim.1        Dim.2      Dim.3        Dim.4       Dim.5
食品     0.6904861 1.931084e-01 0.01326091 2.473299e-02 0.057322080
衣着     0.3562295 3.539108e-01 0.26324670 1.764593e-04 0.003857975
居住     0.7525348 1.486382e-03 0.09045477 9.690331e-02 0.002864013
医疗     0.4066792 4.417699e-01 0.06019753 7.307029e-06 0.081780663
交通通讯 0.8305430 1.047888e-01 0.01000794 4.334712e-03 0.006254705
教育     0.8047584 5.668948e-05 0.02017483 6.473142e-02 0.007228276
> # Contributions to the principal components
> head(var$contrib)
             Dim.1        Dim.2     Dim.3        Dim.4     Dim.5
食品     12.813875 15.670621251  2.492305  8.281295869 25.252800
衣着      6.610821 28.719634446 49.475573  0.059083501  1.699601
居住     13.965361  0.120618931 17.000408 32.445942285  1.261719
医疗      7.547054 35.849346573 11.313750  0.002446598 36.027840
交通通讯 15.413018  8.503545275  1.880930  1.451383010  2.755462
教育     14.934513  0.004600315  3.791734 21.673892623  3.184361
```
可以看出在第一主成分中，食品、居住、交通通讯和教育贡献率较大，而在第二主成分中医疗贡献率较大，在第三主成分中衣着贡献率较大。

- Correlation circle

The correlation between a variable and a principal component (PC) is used as the coordinates of the variable on the PC i.e. variables are represented by their correlations.

![varpca.png](https://upload-images.jianshu.io/upload_images/17587926-28f7ea4d67a7e894.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- Quality of representation

The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates) . 

其实吧……我觉得它的第一主成分就包含了足够多的信息哈哈，把它直接降至一维就很不错。

![cor.png](https://upload-images.jianshu.io/upload_images/17587926-6fa146a126f20cb3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

It’s also possible to create a bar plot of variables cos2 using the function fviz_cos2()[in factoextra]

```r
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:2)
```

![bar.png](https://upload-images.jianshu.io/upload_images/17587926-3ac6214879efbc55.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- A high cos2 indicates a good representation of the variable on the principal component. In this case the variable is positioned close to the circumference of the correlation circle.

- A low cos2 indicates that the variable is not perfectly represented by the PCs. In this case the variable is close to the center of the circle.

而我们发现大部分的变量的cos2均较高，这与这些变量在之前的相关圆中接近圆周是一致的。这也表明用两个主成分能很好地反应这些变量的信息。

下用颜色代表不同的gradient重新绘制 color variables by their cos2 values

![varpca_grad.png](https://upload-images.jianshu.io/upload_images/17587926-3ba8555726d8a078.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

后续关于PCA分析绘图可以看ref链接

reference: [http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials)