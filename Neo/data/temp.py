'''
f = open('414.circles','r')
f_o = open('index.txt','w')
num_node = 0
node_index = {}
node_class = {}
cls = 0
dup_node = 0
while True:
	line = f.readline()
	if not line: break
	line_list = line.split('	')
	for n in line_list:
		node = int(n)
		try:
			exist = node_index[node]
		except:
			num_node = num_node + 1
			node_index[node] = num_node
			f_o.write(str(node) + ' : ' + str(num_node) + '\n')
	cls = cls + 1
f.close()
f_o.close()
print('#node: ' + str(len(node_index)))
for key,value in node_index.items():
	print(key)

'''
node_index = {}
num_node = 0

f = open('414.edges','r')
f_o = open('index.txt','w')
while True:
	line = f.readline()
	if not line: break
	line_list = line.split()
	u = int(line_list[0])
	v = int(line_list[1])
	try:
		exist = node_index[u]
	except:
		num_node = num_node + 1
		node_index[u] = num_node
		f_o.write(str(u) + '\n')
	try:
		exist = node_index[v]
	except:
		num_node = num_node + 1
		node_index[v] = num_node
		f_o.write(str(v) + '\n')
f.close()
f_o.close()
ego_node = num_node + 1
node_adj = [[] for _ in range(num_node + 2)]
f = open('414.edges','r')
while True:
	line = f.readline()
	if not line: break
	line_list = line.split()
	u = int(line_list[0])
	v = int(line_list[1])
	u_index = node_index[u]
	v_index = node_index[v]
	node_adj[u_index].append(v_index)
	node_adj[v_index].append(u_index)
f.close()

for i in range(1,num_node+1):
	node_adj[ego_node].append(i)
	node_adj[i].append(ego_node)

num_edge = 0
for i in range(1,num_node+2):
	node_adj[i] = list(set(node_adj[i]))
	num_edge = num_edge + len(node_adj[i])
num_edge = num_edge / 2

print('#edge: ' + str(int(num_edge)))

f = open('facebook-graph','w')
f.write(str(num_node) + ' ' + str(int(num_edge)) + '\n')
for i in range(1,num_node + 2):
	for j in node_adj[i]:
		f.write(str(j) + ' ')
	f.write('\n')

