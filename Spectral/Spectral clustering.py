import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.algorithms import community

def Girvan_Newman_algorithm(G, weight):
    """ G는 원래 네트워크 g는 Edge를 한개씩 끊어나갈 네트워크 """
    g = G.copy()

    """ initial """
    step = 0                # step
    log_step = []           # step 기록
    log_modularity = []     # modularity 기록
    old_max_m = 0           # 이전 최대 modularity 기억
    max_g = g.copy()        # modularity가 최대일 때의 네트워크 여기서는 초기화 작업
    k = sorted(nx.connected_components(G), key=len, reverse=True)   # k 는 모두 연결되어있는 Community를 노드로 나타낸 값
    k_list = []
    for j in range(len(k)):
        k_list = k_list + [list(k[j])]
    max_k = k_list  # max_k 는 modularity가 최대일 때의 k 값 저장용
    m = community.modularity(G, communities=k, weight=weight)   # modularity
    max_m = m       # max_m은 modularity가 최대일 때 값 기록용
    max_step = 0    # max_step은 modularity가 최대일 때 step값 기록용

    """ Girvan-Newman algorithm """
    while len(g.edges()) > 0:
        k = sorted(nx.connected_components(g), key=len, reverse=True)  # 커뮤니티 추출
        m = community.modularity(G, communities=k, weight=weight)   # 추출된 커뮤니티의 modularity 계산
        if m > old_max_m:   # 이전 최대 modularity보다 현재 modularity가 높을 경우 기록
            max_g = g.copy()
            max_m = m
            k_list = []
            for j in range(len(k)):
                k_list = k_list + [list(k[j])]
            max_k = k_list
            max_step = step
            old_max_m = m
        log_step = log_step + [step]    # 로깅용
        log_modularity = log_modularity + [m]   # 로깅용
        print("step: ", step, "  modularity: ", m)

        """ remove edge """
        step = step + 1
        betweenness = nx.edge_betweenness_centrality(g, weight=weight)  # betweennes centrality 계산
        max_edge = max(betweenness, key=betweenness.get)    # betweeness centrality가 가장 큰 Edge 선택
        g.remove_edge(max_edge[0], max_edge[1])     # Edge 추출

    return log_step, log_modularity, max_g, max_m, max_k, max_step


def Add_Inner_edges(range, num):
    inner_edges = []
    while len(inner_edges) < num:
        new_edge = tuple(np.sort(np.random.choice(range, size=2, replace=None)))
        if new_edge not in inner_edges:
            inner_edges = inner_edges + [new_edge]
    return inner_edges


def Add_Outer_edges(Community_all, num):
    # 두 커뮤니티 선택
    outter_edges = []
    while len(outter_edges) < num:
        #group_choiced = np.random.choice(range(len(Community_all)), size=2, replace=None) # 범용적으로 커뮤니티 선택할 시
        if len(outter_edges) < 3:
            group_choiced = np.random.choice([0, 1], size=2, replace=None)
        elif len(outter_edges) < 6:
            group_choiced = np.random.choice([0, 2], size=2, replace=None)
        elif len(outter_edges) < 10:
            group_choiced = np.random.choice([1, 2], size=2, replace=None)
        new_edge = tuple(np.sort([np.random.choice(Community_all[group_choiced[0]], replace=None),
                    np.random.choice(Community_all[group_choiced[1]], replace=None)]))
        if new_edge not in outter_edges:
            outter_edges = outter_edges + [new_edge]
    return outter_edges


""" 본문 """

""" Network 생성 """
# 8명, 10명, 12명 의 3팀 생성, 총 30명

G = nx.Graph()
G.add_nodes_from(range(30))
# Community 설정
Community1 = range(0, 8)
Community2 = range(8, 18)
Community3 = range(18, 30)
Community_all = Community1, Community2, Community3

# 내부 연결 추가
Inner_edges_1 = Add_Inner_edges(range=Community1, num=15)
Inner_edges_2 = Add_Inner_edges(range=Community2, num=20)
Inner_edges_3 = Add_Inner_edges(range=Community3, num=25)
G.add_edges_from(Inner_edges_1)
G.add_edges_from(Inner_edges_2)
G.add_edges_from(Inner_edges_3)

# 외부 연결 추가
Outter_edges = Add_Outer_edges(Community_all, 10)
G.add_edges_from(Outter_edges)

""" 수동 Modularity 계산 """
E = np.array([[30, 3, 3],[3, 40, 4],[3, 4, 50]]) / 140
Modularity = np.trace(E) - np.sum(np.dot(E, E))


""" Girvan Newman algorithm"""
log_step, log_modularity, max_g, max_m, max_k, max_step = Girvan_Newman_algorithm(G, weight=None)   # weight='weight' : weighted network
# Algorithm 실행결과 플롯
fig = plt.figure()
plt.subplots_adjust(hspace=0.5, wspace=0.3)
plt.plot(log_step, log_modularity)
plt.xlabel('step')
plt.ylabel('modularity')
plt.title("unweighted")
plt.show()


""" Graph G Plot """
pos = nx.spring_layout(G)
fig = plt.figure(figsize=(7, 6))
node_size = 100
node_color_list = []
# max_k, 즉 Modularity 가 가장 높은 지점의 Community k 별로 색깔 설정
for i in range(len(G.nodes)):
    if i in max_k[0]:
        node_color_list = node_color_list + ['red']
    elif i in max_k[1]:
        node_color_list = node_color_list + ['yellow']
    elif i in max_k[2]:
        node_color_list = node_color_list + ['blue']
im = nx.draw_networkx_nodes(G, pos, node_color=node_color_list, node_size=node_size)
nx.draw_networkx_edges(G, pos)
nx.draw_networkx_labels(G, pos, font_size=10, font_color="black")
plt.xticks([])
plt.yticks([])
plt.show()