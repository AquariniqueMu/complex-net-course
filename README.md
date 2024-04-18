# 复杂网络的技术实践及分析

## 一、复杂网络构建与计算工具：Networkx

中文文档：[复杂网络软件 — NetworkX 2.8 文档 (osgeo.cn)](https://www.osgeo.cn/networkx/)

英文文档：[Tutorial — NetworkX 3.3 documentation](https://networkx.org/documentation/stable/tutorial.html)

NetworkX是一个用于创建、操作和研究复杂网络的结构、动力学和功能的Python包。它提供：

- 研究社会、生物和基础设施网络结构和动态的工具；
- 一种适用于多种应用的标准编程接口和图形实现；
- 为协作性、多学科项目提供快速发展环境；
- 与现有的数值算法和C、C++和FORTRAN代码的接口；
- 能够轻松处理大型非标准数据集。

> Networkx对图网络的节点、连边、属性等特征做了完善的处理，可以非常方便地执行图中常用的操作。

### Networkx的安装（已配置好Python环境和Jupyter）

```bash
pip install networkx
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple networkx
```

## 二、可视化工具：Gephi

**Gephi**是基于**Java**的制作关系网络图的工具，常用于网络分析的可视化环节。Networkx可以保存gexf文件，与gephi具有良好的衔接性，是绘制网络拓扑图的必备软件。





Gephi下载链接：[Download (gephi.org)](https://gephi.org/users/download/)

JAVA下载链接：[[Java Downloads | Oracle 中国](https://www.oracle.com/cn/java/technologies/downloads/#jdk22-windows)](https://www.java.com/zh-CN/)

文件格式：.gexf



## 三、节点中心性算法

### 1、度中心性

概念上最简单是度中心性，它定义为一个节点上事件的链接数量（即一个节点拥有的关系数量）。度可以解释为节点捕获的任何流经网络的东西（例如病毒或某些信息）的直接风险。在有向网络的情况下（关系有方向） ，我们通常定义两个独立的度中心性的度量，即 **入度 Indegree**和**出度 Outdegree**。因此，入度是指向该节点的关系数，出度是该节点指向其他节点的关系数。当关系与一些积极的方面如友谊或合作有关时，入度通常被解释为一种受欢迎的形式，而出度则被解释为一种合群的形式。
$$
degree(u)=\frac{N_u}{N-1}
$$
度中心性的数值一般为邻居节点的数量除以网络中的节点总数-1。

### 2、接近中心性

节点的接近中心性是所有 n-1 个可到达节点到节点的平均最短路径距离的倒数。反映在网络中某一节点与其他节点之间的接近程度。如果一个节点离其他的节点都很近，那么传递信息的时候就不需要依赖其他的节点，说明这个节点很重要。
$$
C(u) = \frac{n - 1}{\sum_{v=1}^{n-1} d(v, u)},
$$

> Linton C. Freeman: Centrality in networks: I. Conceptual clarification. Social Networks 1:215-239, 1979. https://doi.org/10.1016/0378-8733(78)90021-7

### 3、介数中心性

**介数中心性**（betweenness centrality，又译作**中间中心性**）是基于最短路径的衡量标准。对全连接网络图，其中任意两个节点均至少存在一个最短路径，在无权重网络图中该最短路径是路径包含边的数量求和，加权网络图中该最短路径则是路径包含边的权重求和。每个节点的介数中心性即为这些最短路径穿过该节点的次数。
$$
c_B(v) =\sum_{s,t \in V} \frac{\sigma(s, t|v)}{\sigma(s, t)}
$$

> Ulrik Brandes: A Faster Algorithm for Betweenness Centrality. Journal of Mathematical Sociology 25(2):163-177, 2001. https://doi.org/10.1080/0022250X.2001.9990249

### 4、PageRank中心性

PageRank 根据传入链接的结构计算图 G 中节点的排名。它最初是作为网页排名的算法。

**核心思想**是基于有向图上的随机游走模型，是一个一阶马尔可夫链。描述了一个随机游走者如何沿着图的边随机移动，从一个节点访问到另一个节点。在满足某些条件的前提下，这个随机游走过程最终会收敛到一个平稳分布。在这个平稳分布中，每个节点被访问的概率即为其 PageRank 值，这个值可以被解释为节点的重要性。

- PageRank 是递归定义的，这意味着一个页面的 PageRank 值部分地取决于链接到它的其他页面的 PageRank 值。因此，计算 PageRank 值通常需要迭代方法。

> Page, Lawrence; Brin, Sergey; Motwani, Rajeev and Winograd, Terry,  The PageRank citation ranking: Bringing order to the Web. 1999. http://dbpubs.stanford.edu:8090/pub/showDoc.Fulltext?lang=en&doc=1999-66&format=pdf

### 5、特征向量中心性

特征向量中心性是根据节点邻域的中心性来计算节点的中心性。节点 i 的特征向量中心度为特征向量的第 i 个元素
$$
Ax = \lambda x
$$
其中，A 是图`G`的邻接矩阵，其特征值为 λ 。根据Perron–Frobenius定理，如果 λ 是邻接矩阵的最大特征值，则存在唯一的解 x，其所有项均为正值。

### 6、K-shell 中心性

k-壳分解法(k-shell)确定网络中节点的位置, 将外围的节点层层剥去, 处于内层的节点拥有较高的影响力。可看成是一种 基于节点度的粗粒化排序方法。

网络中如果存在度为 1 的节点, 从度中心性的角度看它们就是最不重要的节点. 如果把这些度为 1 的节点及其所连接的边都去掉, 剩下的网络中会新出现一些度为 1 的节点, 再将这些度为 1 的节点去掉, 循环操作, 直到所剩的网络中没有度为 1 的节点为止. 此时, 所有被去掉的节点组成一个层, 称为 1-壳(记 为 ks=1)。





## 四、SIR传播模型



SIR模型最早由 *W. O. Kermack* 与 *McKendrick, A. G. McKendrick* 于1927年发表，逐渐发展为最成功、最著名的传染病传播模型之一。在 SIR 模型中，全体人口被划分成三类人群：

S: 易感个体的数量。当易感个体和传染个体发生 "传染接触 "时，易感个体就会感染疾病并转入传染区。
I：感染者的数量。这些个体已被感染，并能感染易感个体。
R：移除（免疫）或死亡个体的数量。这些个体受到感染后，要么已经痊愈并进入移除区，要么已经死亡。假定死亡人数与总人口相比可以忽略不计。这部分人也可称为 "康复者 "或 "抵抗者"。

将SIR传播模型和复杂网络相结合，传播过程可以基于网络结构进行模拟和分析，研究传染病在网络中的传播机制、传播速度、传播规模等问题，为疾病控制和预防提供科学依据。类似地，信息、舆论的传播也可以用这样的模型进行模拟。





## 五、案例实践

### 1、复杂网络构建与处理

- 添加节点和连边
- 观察网络信息
- 画出网络拓扑图
- 读取和保存各种类型的网络文件
- 提取网络的子网
- 提取网络的最大连通组件



### 2、节点中心性算法与传播模型模拟

- 度中心性算法
- 接近中心性算法
- 介数中心性算法
- Pagerank算法
- K-shell算法
- 特征向量中心性算法
- 对比各种算法的结果

### 3、Networkx + Gephi可视化

- 根据上一步的结果进行可视化，观察不同算法评价节点的差异