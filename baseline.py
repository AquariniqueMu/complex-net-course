'''
Description: 一些经典baseline算法的代码合集
Author: Junwen Yang
Date: 2022-10-13 09:42:40
LastEditTime: 2023-11-27 23:12:49
LastEditors: Junwen Yang
'''
import numpy as np
import networkx as nx
import math
import pandas as pd
import kshell as ks
import graph_create as gc
import tools_func as tf
from alive_progress import alive_bar
import warnings
warnings.filterwarnings('ignore')
import time

def GM(G, DIST_MAT=None):
    """
    经典重力模型中心性 :
    以度值为质量, 平均路径长度的一半为截断半径
    """
    # 获取节点的度值,key为节点,value为度值
    MASS = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        MASS[node] = G.degree(node)

    # 计算网络的平均路径长度, 取一半并四舍五入, 为截断半径
    AVER_DIST = tf.AVERAGE_DISTANCE(G)
    INT_AVER_DIST = int(AVER_DIST / 2 + 0.5)
    # print(INT_AVER_DIST)
    # 搜索所有节点的 [AVER_DIST四舍五入] 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, INT_AVER_DIST)

    # 根据计算的截断半径循环计算节点的中心性分数
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        # print(DICT_NEIGHBOR[node])
        for i in range(INT_AVER_DIST):
            # print(list(MASS.keys()))
            CENTRALITY[node] += sum(
                [MASS[target]*MASS[node] / (i+1)**2 for target in DICT_NEIGHBOR[node][i]])

    return CENTRALITY




def KGM(G, DIST_MAT=None):
    """
    K-shell重力模型中心性:
    截断半径为三阶邻居, 并最后做一次一阶增广
    """
    # 获取节点的K-shell值
    MASS = ks.kshell(G)

    # 搜索所有节点一二三阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, 3)

    # 在节点的三阶范围内计算重力中心性
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        CENTRALITY[node] = sum(
            [MASS[target]*MASS[node] for target in DICT_NEIGHBOR[node][0]]) + sum(
                [MASS[target]*MASS[node] / 4 for target in DICT_NEIGHBOR[node][1]]) + sum(
                    [MASS[target]*MASS[node] / 9 for target in DICT_NEIGHBOR[node][2]])

    # 在节点的一阶范围内计算增广重力中心性, 即所有一节邻居的重力中心性之和
    CENTRALITY_EX = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        CENTRALITY_EX[node] = sum(
            [CENTRALITY[target] for target in DICT_NEIGHBOR[node][0]])

    return CENTRALITY_EX


def MultiGM(G, DIST_MAT = None):
    """
    多属性重力模型中心性:
    """
    MASS_DEGREE = G.degree()
    MASS_KSHELL = ks.kshell(G)


def GMM(G, DIST_MAT = None):
    """
    广义重力模型中心性: 特征向量EC中心性的分数为系数, 度值为质量
    """
    # 计算节点度值和EC中心性
    MASS = G.degree()
    COEFFICIENT_EC = nx.eigenvector_centrality(G, max_iter=5000)

    # 计算网络的平均路径长度, 取一半并四舍五入, 为截断半径
    AVER_DIST = tf.AVERAGE_DISTANCE(G)
    INT_AVER_DIST = int(AVER_DIST / 2 + 0.5)

    # 搜索所有节点的 [AVER_DIST四舍五入] 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, INT_AVER_DIST)

    # 根据计算的截断半径循环计算节点的中心性分数
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        for i in range(INT_AVER_DIST):
            CENTRALITY[node] += sum(
                [COEFFICIENT_EC[node] * MASS[target] * MASS[node] / (i+1)**2 for target in DICT_NEIGHBOR[node][i]])

    return CENTRALITY


def GGC(G, DIST_MAT = None):
    """
    广义重力中心性: 质量为局部聚类系数的指数形式乘以度值, 聚类系数前有α作为调节参数
    """
    DEGREE = [a[1] for a in G.degree()]
    CLUSTER = nx.clustering(G)
    alpha = 2
    MASS = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        MASS[node] = math.exp((-1) * alpha * CLUSTER[node]) * G.out_degree(node)

    # 计算网络的平均路径长度, 取一半并四舍五入, 为截断半径
    AVER_DIST = tf.AVERAGE_DISTANCE(G)
    INT_AVER_DIST = int(AVER_DIST / 2 + 0.5)

    # 搜索所有节点的 [AVER_DIST四舍五入] 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, INT_AVER_DIST)

    # 根据计算的截断半径循环计算节点的中心性分数
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        for i in range(INT_AVER_DIST):
            CENTRALITY[node] += sum(
                [MASS[target] * MASS[node] / (i+1)**2 for target in DICT_NEIGHBOR[node][i]])

    return CENTRALITY


def KSGC(G):
    """
    计算k-shell based重力中心性: k-shell归一化指数为系数, 度值为质量
    """
    # 获取节点的度值
    MASS = G.degree()
    KSHELL = ks.kshell(G)
    KS_POLAR_DIFF = max(KSHELL.values()) - min(KSHELL.values())

    # 计算网络的平均路径长度, 取一半并四舍五入, 为截断半径
    AVER_DIST = tf.AVERAGE_DISTANCE(G)
    INT_AVER_DIST = int(AVER_DIST / 2 + 0.5)

    # 搜索所有节点的 [AVER_DIST四舍五入] 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, INT_AVER_DIST)

    # 计算中心性分数
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        for i in range(INT_AVER_DIST):
            CENTRALITY[node] += sum(
                [MASS[target] * MASS[node] * math.exp((KSHELL[node] - KSHELL[target]) / (KS_POLAR_DIFF)) / (i+1)**2 for target in DICT_NEIGHBOR[node][i]])

    return CENTRALITY


def EffG(G):
    """
    计算基于有效距离的重力中心性分数
    """
    MASS = G.degree()
    ADJ_MAT = nx.adjacency_matrix(G)
    ADJ_MAT[ADJ_MAT != 0] = 1

    # 搜索所有节点的 1 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, 1)

    G_copy = nx.DiGraph()
    G_copy.add_nodes_from(list(G.nodes()))
    # print(G_copy)
    for node in G.nodes():
        for neighbor in DICT_NEIGHBOR[node][0]:
            # ADJ_MAT[node, neighbor] = 1 - math.log2(1 / G.out_degree(node))
            G_copy.add_edge(node,neighbor,length=1 - math.log2(1 / G.out_degree(node)))
            # print( G_copy.get_edge_data(node,neighbor))


    SHORTEST_LENGTH_MAT = dict()
    with alive_bar(nx.number_of_nodes(G_copy)**1) as bar:
        for node in G_copy.nodes():
            nodeslist = nx.single_source_shortest_path_length(G_copy,node)
            del nodeslist[node]
            NODE_SHORTEST_LENGTH = dict()
            for target in nodeslist.keys():

                NODE_SHORTEST_LENGTH[target] = nx.shortest_path_length(G_copy,source=node,target=target,weight='length')
            SHORTEST_LENGTH_MAT[node] = NODE_SHORTEST_LENGTH
            bar()
    # print(SHORTEST_LENGTH_MAT[1])


    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        # print(SHORTEST_LENGTH_MAT[node].keys())
        for target in SHORTEST_LENGTH_MAT[node].keys():
            # print(target)
            CENTRALITY[node] += MASS[node] * MASS[target] / SHORTEST_LENGTH_MAT[node][target]**2

    # print(CENTRALITY)




    return CENTRALITY
    
    
def compute_probability(Source_G):
    """compute the infection probability
    # Arguments
        Source_G: a graph as networkx Graph
    Returns
        the infection probability computed by  formula: <k> / (<k^2> - <k>)
    """
    G = nx.DiGraph()
    G = Source_G
    degree_dict = dict(G.degree())
    k = 0.0
    k_pow = 0.0
    for i in degree_dict:
        k = k + degree_dict[i]
        k_pow = k_pow + degree_dict[i] * degree_dict[i]

    k = k / G.number_of_nodes()
    k_pow = k_pow / G.number_of_nodes()
    pro_cover = k / (k_pow - k)
    pro_sir=k/k_pow
    return pro_cover,pro_sir


def improved_wvoterank(G):
    # 初始化节点分数与投票能力
    voting_score_Sv = dict.fromkeys(G.nodes(), 0)
    voting_ability_VAv = dict.fromkeys(G.nodes(), 1)
    # 初始化一跳投票分数与二跳投票分数
    one_hop_score = dict.fromkeys(G.nodes(), 0)
    two_hop_score = dict.fromkeys(G.nodes(), 0)
    # 初始化一跳投票权重与二跳投票权重的衰减因子
    d = dict(nx.degree(G,weight=None))
    one_hop_sigma = 1/(sum(d.values())/len(G.nodes()) * 1)
    two_hop_sigma = 1/(sum(d.values())/len(G.nodes()) * 2)
    
    # 初始化每次投票的高分节点存储列表
    sort_list = []
    
    # 投票过程，每次投票选出最高分节点，将其投票能力置为0，将其一跳节点与二跳节点的投票能力减去相应的权重
    while len(sort_list) < len(G.nodes()):
        # 进行投票
        for node in G.nodes():
            for neighbor in G.neighbors(node):
                # 计算一阶邻居的投票分数
                one_hop_score[node] += voting_ability_VAv[neighbor] * G[node][neighbor]['weight']
                # voting_score_Sv[neighbor] += voting_ability_VAv[node]
                for neighbor2 in G.neighbors(neighbor):
                    # 计算二阶邻居的投票分数
                    two_hop_score[node] += voting_ability_VAv[neighbor2]* G[neighbor][neighbor2]['weight']
            # 计算每个节点本轮投票的得分
            voting_score_Sv[node] = math.sqrt(len(list(G.neighbors(node))) * one_hop_score[node] + two_hop_score[node])   
             
        # 选出本轮投票的最高分节点
        sort_list.append(max(voting_score_Sv, key=voting_score_Sv.get))
        
        # 对选出的高分节点及其一阶与二阶邻居节点的投票能力进行更新
        for node in sort_list:
            # 被选节点投票能力置为0
            voting_ability_VAv[node] = 0
            for neighbor in G.neighbors(node):
                # 一阶邻居投票能力减去一阶衰减因子
                voting_ability_VAv[neighbor] -= one_hop_sigma
                # 投票能力最低为0
                if voting_ability_VAv[neighbor] < 0:
                    voting_ability_VAv[neighbor] = 0
                for neighbor2 in G.neighbors(neighbor):
                    # 二阶邻居投票能力减去二阶衰减因子
                    voting_ability_VAv[neighbor2] -= two_hop_sigma
                    # 投票能力最低为0
                    if voting_ability_VAv[neighbor2] < 0:
                        voting_ability_VAv[neighbor2] = 0
    # 返回节点的投票分数
    return voting_score_Sv



def weighted_PageRank(G, alpha=0.85,theta=0.8, personalization=None, max_iter=100, tol=1.0e-6, nstart=None, weight='weight', dangling=None):
    if len(G) == 0: 
        return {}
    # 将图G转换为有向图 
    if not G.is_directed(): 
        D = G.to_directed() 
    else: 
        D = G 

    # 创建一个字典来存储每个节点的PageRank值 
    # 使用一个字典来存储每个节点的出度 
    W = nx.stochastic_graph(D, weight=weight) 
    N = W.number_of_nodes() 

    # 如果没有指定初始值，则使用均匀分布 
    if nstart is None: 
        x = dict.fromkeys(W, 1.0 / N) 
    else: 
        # 将初始值转换为字典 
        s = float(sum(nstart.values())) 
        x = dict((k, v / s) for k, v in nstart.items()) 

    if personalization is None: 

        # 使用均匀分布作为个性化向量 
        p = dict.fromkeys(W, 1.0 / N) 
    else: 
        missing = set(G) - set(personalization) 
        if missing: 
            raise NetworkXError('Personalization dictionary '
                                'must have a value for every node. '
                                'Missing nodes %s' % missing) 
        s = float(sum(personalization.values())) 
        p = dict((k, v / s) for k, v in personalization.items()) 

    if dangling is None: 

        # 如果没有指定悬挂节点，则使用均匀分布 
        dangling_weights = p 
    else: 
        missing = set(G) - set(dangling) 
        if missing: 
            raise NetworkXError('Dangling node dictionary '
                                'must have a value for every node. '
                                'Missing nodes %s' % missing) 
        s = float(sum(dangling.values())) 
        dangling_weights = dict((k, v/s) for k, v in dangling.items()) 
    dangling_nodes = [n for n in W if W.out_degree(n, weight=weight) == 0.0] 

    # 初始化迭代器 
    for _ in range(max_iter): 
        xlast = x 
        x = dict.fromkeys(xlast.keys(), 0) 
        danglesum = dict.fromkeys(xlast.keys(),0)
        for n in x:
            danglesum[n] = alpha * G.out_degree(weight=weight)[n] / sum([G.out_degree(weight=weight)[n] for n in G.nodes()]) 
        for n in x: 
            if n in dangling_nodes: 
                x[n] = danglesum[n] * (1.0 - alpha)
            else:
            # 个性化向量 
                for nbr in W[n]: 
                    x[nbr] += alpha * xlast[n] * ((theta * G[n][nbr]['weight'] / G.out_degree(weight=weight)[n]) + ((1-theta)/ G.out_degree(weight=None)[n]))
                x[n] += danglesum[n] * (1.0 - alpha)

        # 检查收敛性 
        err = sum([abs(x[n] - xlast[n]) for n in x]) 
        if err < N*tol: 
            return x 
    raise NetworkXError('pagerank: power iteration failed to converge '
                        'in %d iterations.' % max_iter)


def LGC(G:nx.Graph):
    """
    计算Laplacian重力中心性: 以度值组合取代之前的Laplace质量
    
    """
    # 获取节点的度值
    DEGREE = G.degree()
    # KSHELL = ks.kshell(G)
    # KS_POLAR_DIFF = max(KSHELL.values()) - min(KSHELL.values())
    MASS_LC = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        MASS_LC[node] = DEGREE[node] + DEGREE[node] ** 2
        for neighbor in G.neighbors(node):
            MASS_LC[node] += DEGREE[neighbor] * 2
    
    # 计算网络的平均路径长度, 取一半并四舍五入, 为截断半径
    AVER_DIST = tf.AVERAGE_DISTANCE(G)
    INT_AVER_DIST = int(AVER_DIST / 2 + 0.5)

    # 搜索所有节点的 [AVER_DIST四舍五入] 阶邻居
    DICT_NEIGHBOR = dict()
    for node in G.nodes():
        DICT_NEIGHBOR[node] = tf.neighbor_search(G, node, INT_AVER_DIST)

    # 计算中心性分数
    CENTRALITY = dict.fromkeys(G.nodes(), 0)
    for node in G.nodes():
        for i in range(INT_AVER_DIST):
            CENTRALITY[node] += sum(
                [MASS_LC[target] * MASS_LC[node] / (i+1)**2 for target in DICT_NEIGHBOR[node][i]])

    return CENTRALITY







def test():
    G = nx.Graph()
    G = nx.read_weighted_edgelist('./data/stormofswords.edgelist',create_using=nx.DiGraph(),nodetype=int)
    # print(EffG(G,'stormofswords.edgelist'))
    res = LGC(G)
    
    res_sort = sorted(res.items(), key=lambda x: x[1], reverse=True)
    print(res_sort)
# test()

