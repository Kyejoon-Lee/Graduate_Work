import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
from sklearn.metrics.pairwise import kernel_metrics
import networkx as nx
import random
import pandas as pd
import openpyxl
import csv
from networkx.algorithms.community.centrality import girvan_newman
import itertools

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

print('#node: ' + str(num_node-1))
print('#edge: ' + str(num_edge))

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

k = 15
comp = girvan_newman(G)
limited = itertools.takewhile(lambda c: len(c) <= k, comp)
for communities in limited:
    print(tuple(sorted(c) for c in communities))
