import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.algorithms import community

def Girvan_Newman_algorithm(G, weight):

    g = G.copy()

    step = 0
    log_step = []
    log_modularity = []
    old_max_m = 0
    max_g = g.copy()
    k = sorted(nx.connected_components(G), key=len, reverse=True)
    k_list = []
    for j in range(len(k)):
        k_list = k_list + [list(k[j])]
    max_k = k_list
    m = community.modularity(G, communities=k, weight=weight)
    max_m = m
    max_step = 0

    while len(g.edges()) > 0:
        k = sorted(nx.connected_components(g), key=len, reverse=True)
        m = community.modularity(G, communities=k, weight=weight)
        if m > old_max_m:

            max_g = g.copy()
            max_m = m
            k_list = []
            for j in range(len(k)):
                k_list = k_list + [list(k[j])]
            max_k = k_list
            max_step = step
            old_max_m = m
        log_step = log_step + [step]
        log_modularity = log_modularity + [m]
        print("step: ", step, "  modularity: ", m)

        step = step + 1
        betweenness = nx.edge_betweenness_centrality(g, weight=weight)
        max_edge = max(betweenness, key=betweenness.get)
        g.remove_edge(max_edge[0], max_edge[1])

    return log_step, log_modularity, max_g, max_m, max_k, max_step,k_list


src = []
dst = []
f = open('414.csv', 'r')
while True:
	line = f.readline()
	if not line: break
	term_list = line.split(',')
	src.append(int(term_list[0]))
	dst.append(int(term_list[1]))

num_node = max(max(src),max(dst)) + 1
num_edge = len(src)


adj = [[] for i in range(num_node)]
G = nx.Graph()
ego_num = num_node + 1
for i in range(num_edge):
	u = src[i]
	v = dst[i]
	adj[u].append(v)
	adj[v].append(u)
	G.add_edge(u,v)
	G.add_edge(ego_num,u)
	G.add_edge(ego_num,v)



log_step, log_modularity, max_g, max_m, max_k, max_step,k_list = Girvan_Newman_algorithm(G, weight=None)
fig = plt.figure()
plt.subplots_adjust(hspace=0.5, wspace=0.3)
plt.plot(log_step, log_modularity)
plt.xlabel('step')
plt.ylabel('modularity')
plt.title("unweighted")
plt.show()
print(k_list)

