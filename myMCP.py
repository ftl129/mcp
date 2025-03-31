'''
============================================================================
 Name        : A fast maximum clique algorithm based on network decomposition
 Author      : Tianlong Fan (tianlong.fan@unifr.ch)
 Description : Python implementation of the CUBIS-based MCP algorithm

 Copyright (C) 2025 by Tianlong Fan. All rights reserved.

 Please cite the following paper if used:
   Fan T, Jiang W, Zhang Y-C, Lü L. A fast maximum clique algorithm based on
   network decomposition for large sparse networks. Chaos, Solitons and Fractals,
   2025, 194:116239. https://doi.org/10.1016/j.chaos.2025.116239
============================================================================
'''

import time
import copy
import igraph as ig

# === Configuration ===


NetworkAddress = r'G:/Your/Path/To/Network.txt'  # Replace with your actual path
networkName = 'networkname'

LayerNum = 4  # Number of top core layers used in the first decomposition

# === Initialize structures ===
print("Network name:", networkName)
print("Top core layers:", LayerNum)

Coreness = dict()
NodesinSaShell = {}
MaxCliques = 0

IdCovNam = list()
igGraph = ig.Graph()

# === Load network ===
NodeNum = 0
edgeSet = set()
nodeSet = set()
with open(NetworkAddress) as file:
    while True:
        lines = file.readlines(100000)
        if not lines:
            break
        for line in lines:
            # input format of the network

            # example_network
            line = line[:-1]
            edge = list(line.split(','))
            ai = int(edge[0]) - 1
            aj = int(edge[1]) - 1
            unedge = [ai, aj]
            unedge.sort()
            edgeSet.add(tuple(unedge))
            nodeSet.add(ai)
            nodeSet.add(aj)

# === Reindex node IDs and construct graph ===
newID = {old: new for new, old in enumerate(nodeSet)}
igGraph.add_vertices(newID.values())
edges_mapped = [(newID[e[0]], newID[e[1]]) for e in edgeSet]
igGraph.add_edges(edges_mapped)
igGraph.to_undirected()

for i in igGraph.vs.indices:
    igGraph.vs[i]["name"] = i

NodeNum = igGraph.vcount()
EdgeNum = igGraph.ecount()
print("Total nodes:", NodeNum)
print("Total edges:", EdgeNum)

# === Heuristic clique search for pruning ===
def pre_find_cliques(G):
    global MaxCliques
    Q = [None]
    subg = set(G.vs.indices)
    cand = set(G.vs.indices)
    u = max(subg, key=lambda u: len(cand & set(G.neighbors(u))))
    ext_u = cand - set(G.neighbors(u))
    stack = []
    try:
        while True:
            if ext_u:
                q = ext_u.pop()
                cand.remove(q)
                if G.degree(q) + 1 <= MaxCliques:
                    continue
                Q[-1] = q
                adj_q = set(G.neighbors(q))
                subg_q = subg & adj_q
                if len(subg_q) + len(Q) <= MaxCliques:
                    continue
                if not subg_q:
                    MaxCliques = max(MaxCliques, len(Q))
                else:
                    cand_q = cand & adj_q
                    if cand_q:
                        stack.append((subg, cand, ext_u))
                        Q.append(None)
                        subg = subg_q
                        cand = cand_q
                        u = max(subg, key=lambda u: len(cand & set(G.neighbors(u))))
                        ext_u = cand - set(G.neighbors(u))
            else:
                Q.pop()
                subg, cand, ext_u = stack.pop()
    except IndexError:
        pass

# === Exact maximum clique search ===
def find_cliques(G):
    global MaxCliques
    if G.vcount() == 0:
        return
    adj = {u: {v for v in G.neighbors(u) if v != u} for u in G.vs.indices}
    Q = [None]
    subg = set(G.vs.indices)
    cand = set(G.vs.indices)
    u = max(subg, key=lambda u: len(cand & adj[u]))
    ext_u = cand - adj[u]
    stack = []
    try:
        while True:
            if ext_u:
                q = ext_u.pop()
                cand.remove(q)
                if Coreness[G.vs[q]["name"]] + 1 <= MaxCliques:
                    continue
                Q[-1] = q
                adj_q = adj[q]
                subg_q = subg & adj_q
                if len(subg_q) + len(Q) <= MaxCliques:
                    continue
                if not subg_q:
                    MaxCliques = max(MaxCliques, len(Q))
                else:
                    cand_q = cand & adj_q
                    if cand_q:
                        stack.append((subg, cand, ext_u))
                        Q.append(None)
                        subg = subg_q
                        cand = cand_q
                        u = max(subg, key=lambda u: len(cand & adj[u]))
                        ext_u = cand - adj[u]
            else:
                Q.pop()
                subg, cand, ext_u = stack.pop()
    except IndexError:
        pass

# === Step 1: Preprocessing ===
time000 = time.time()
PreDegree = igGraph.degree()
PreMaxDegNode = PreDegree.index(max(PreDegree))
subSet = set([PreMaxDegNode] + igGraph.neighbors(PreMaxDegNode))
randomselectsubNet = igGraph.subgraph(subSet)
pre_find_cliques(randomselectsubNet)

removeNodes = [i for i in igGraph.vs.indices if igGraph.degree(i) < MaxCliques]
igGraph.delete_vertices(removeNodes)
print("Nodes removed by pruning:", len(removeNodes))

# === Step 2: Compute coreness and shell structure ===
timeCore = time.time()
igCoreness = igGraph.shell_index(mode='ALL')

IdCovNam = [0] * len(igCoreness)
for i in igGraph.vs.indices:
    IdCovNam[i] = igGraph.vs[i]["name"]
    Coreness[igGraph.vs[i]["name"]] = igCoreness[i]

for nod in range(len(igCoreness)):
    NodesinSaShell.setdefault(igCoreness[nod], set()).add(nod)
NumShell = len(NodesinSaShell)

# === Step 3: Construct and search first shell ===
def parallelFormShellGraph_1(coreList):
    CoreNode = set().union(*(NodesinSaShell[i] for i in coreList))
    ShellGraph = igGraph.subgraph(CoreNode)
    find_cliques(ShellGraph)
    return ShellGraph.vcount(), ShellGraph.ecount()

# === Step 4: Construct and search second shell ===
def parallelFormShellGraph(coreList):
    if not coreList:
        return 0, 0
    maxCo = max(coreList)
    minCo = min(coreList)
    CoreNode = set().union(*(NodesinSaShell[i] for i in coreList))

    realCoreNode = set()
    for noCur in CoreNode:
        neiL = [i for i in igGraph.neighbors(noCur) if igCoreness[i] >= minCo]
        neiL.sort(key=lambda d: igGraph.degree(d))
        addFlag = True
        Cot = 0
        satisfied = 0
        for n in neiL:
            if len(set(neiL) & set(igGraph.neighbors(n))) > MaxCliques - 2:
                satisfied += 1
                if satisfied > MaxCliques - 1:
                    addFlag = True
                    break
            else:
                Cot += 1
                if Cot > len(neiL) - (MaxCliques - 1):
                    addFlag = False
                    break
        if addFlag:
            realCoreNode.add(noCur)
            realCoreNode.update([nei for nei in igGraph.neighbors(noCur) if igCoreness[nei] > maxCo])

    ShellGraph = igGraph.subgraph(realCoreNode)
    find_cliques(ShellGraph)
    return ShellGraph.vcount(), ShellGraph.ecount()

# === Step 5: Execute algorithm ===
coreLi = sorted(NodesinSaShell.keys(), reverse=True)
PcoreList = set()
time_start = time.time()

if len(coreLi) <= LayerNum:
    PcoreList = set(coreLi)
    coreLi.clear()
else:
    for _ in range(LayerNum):
        PcoreList.add(coreLi.pop(0))

nM1, eM1 = parallelFormShellGraph_1(PcoreList)
time_end1 = time.time()
MaxCliques1 = MaxCliques

PcoreList.clear()
minCore = 0
for i in coreLi:
    if i + 1 > MaxCliques:
        PcoreList.add(i)
    else:
        minCore = i
        break

if not PcoreList:
    print("Total core layers:", NumShell)
    print("Search ended early after first shell.")
    print("Maximum clique size:", MaxCliques)
    print("First shell node count = %d, ratio = %.6f" % (nM1, nM1 / float(NodeNum)))
    print("First shell edge count = %d, ratio = %.6f" % (eM1, eM1 / float(EdgeNum)))
    print("Time used for first shell:", round(time_end1 - time_start, 3), "seconds")
    exit("Search completed.")

nM2, eM2 = parallelFormShellGraph(PcoreList)
time_end2 = time.time()

# === Final Output ===
print("First shell nodes:", nM1)
print("First shell edges:", eM1)
print("Time for first search:", round(time_end1 - time_start, 3))

print("Second shell nodes:", nM2)
print("Second shell edges:", eM2)
print("Time for second search:", round(time_end2 - time_end1, 3))

print("Maximum clique size:", MaxCliques)
print("Total runtime:", round(time_end2 - time_start, 3))