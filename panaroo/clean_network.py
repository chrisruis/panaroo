import networkx as nx
from panaroo.cdhit import *
from panaroo.merge_nodes import merge_nodes
from panaroo.isvalid import del_dups
from collections import defaultdict, deque
from panaroo.cdhit import is_valid
from itertools import chain, combinations
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
from tqdm import tqdm
# import graph_tool.all as gt

# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=3, max_recursive=2):

    # fix trailing
    for i in range(max_recursive):
        bad_nodes = []
        removed = False
        for (node, val) in G.degree():
            if val <= 1:  # trailing node
                if G.nodes[node]['size'] < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)
            removed = True

        if not removed: break

    return G

def mod_bfs_edges(G, source, depth_limit=None):
    """Iterate over edges in a breadth-first search.
    Modified version of 'generic_bfs_edges' from networkx
    """
    neighbors = G.neighbors

    visited = {source}
    if depth_limit is None:
        depth_limit = len(G)
    queue = deque([(source, depth_limit, neighbors(source))])
    while queue:
        parent, depth_now, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent, child, depth_now
                visited.add(child)
                if depth_now > 1:
                    queue.append((child, depth_now - 1, neighbors(child)))
        except StopIteration:
            queue.popleft()

def max_clique(G):
    max_clique = []
    max_len = 0
    for clique in nx.find_cliques(G):
        if len(clique)>max_len:
            max_len = len(clique)
            max_clique = clique
    return max_clique

def collapse_families(G,
                      outdir,
                      family_threshold=0.7,
                      dna_error_threshold=0.99,
                      correct_mistranslations=False,
                      n_cpu=1,
                      quiet=False,
                      distances_bwtn_centroids=None, 
                      centroid_to_index=None):

    node_count = max(list(G.nodes())) + 10

    
    if correct_mistranslations:
        depths = [1, 2, 3]
        threshold = [0.99, 0.98, 0.95, 0.9]
    else:
        depths = [1, 2, 3]
        threshold = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]

    # precluster for speed
    if correct_mistranslations:
        cdhit_clusters = iterative_cdhit(
            G,
            outdir,
            thresholds=threshold,
            n_cpu=n_cpu,
            quiet=True,
            dna=True,
            word_length=7,
            accurate=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, dna_error_threshold, dna=True, n_cpu=n_cpu)

        # keep track of centroids for each sequence. Need this to resolve clashes
        seqid_to_index = {}
        for node in G.nodes():
            for sid in G.nodes[node]['seqIDs']:
                seqid_to_index[sid] = centroid_to_index[G.nodes[node]["longCentroidID"][1]]

    elif distances_bwtn_centroids is None:
        cdhit_clusters = iterative_cdhit(G,
                                         outdir,
                                         thresholds=threshold,
                                         n_cpu=n_cpu,
                                         quiet=True,
                                         dna=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, family_threshold, dna=False, n_cpu=n_cpu)
    for depth in depths:
        search_space = set(G.nodes())
        while len(search_space) > 0:
            # look for nodes to merge
            temp_node_list = list(search_space)
            removed_nodes = set()
            for node in temp_node_list:
                if node in removed_nodes: continue

                if G.degree[node] <= 2:
                    search_space.remove(node)
                    removed_nodes.add(node)
                    continue

                # find neighbouring nodes and cluster their centroid with cdhit
                neighbours = [
                    v for u, v in nx.bfs_edges(G, source=node, depth_limit=depth)
                ]
                if correct_mistranslations:
                    neighbours += [node]

                # find clusters
                index = []
                neigh_array = []
                for neigh in neighbours:
                    for sid in G.nodes[neigh]['centroid']:
                        if sid not in centroid_to_index:
                            print(centroid_to_index)
                            print(sid)
                        index.append(centroid_to_index[sid])
                        neigh_array.append(neigh)
                index = np.array(index, dtype=int)
                neigh_array = np.array(neigh_array)

                n_components, labels = connected_components(
                    csgraph=distances_bwtn_centroids[index][:, index],
                    directed=False,
                    return_labels=True)
                # labels = labels[index]
                for neigh in neighbours:
                    l = list(set(labels[neigh_array == neigh]))
                    if len(l)>1:
                        for i in l[1:]:
                            labels[labels==i] = l[0]

                clusters = [
                    del_dups(list(neigh_array[labels == i])) for i in np.unique(labels)
                ]

                for cluster in clusters:

                    # check if there are any to collapse
                    if len(cluster) <= 1: continue

                    # check for conflicts
                    members = []
                    for n in cluster:
                        for m in G.nodes[n]['members']:
                            members.append(m)

                    if (len(members) == len(set(members))):
                        # no conflicts so merge
                        node_count += 1
                        for neig in cluster:
                            removed_nodes.add(neig)
                            if neig in search_space: search_space.remove(neig)
                        temp_c = cluster.copy()
                        G = merge_nodes(
                            G,
                            temp_c.pop(),
                            temp_c.pop(),
                            node_count,
                            multi_centroid=(not correct_mistranslations))
                        while (len(temp_c) > 0):
                            G = merge_nodes(
                                G,
                                node_count,
                                temp_c.pop(),
                                node_count + 1,
                                multi_centroid=(not correct_mistranslations))
                            node_count += 1
                        search_space.add(node_count)
                    else:
                        # if correct_mistranslations:

                        # merge if the centroids don't conflict and the nodes are adjacent in the conflicting genome
                        # this corresponds to a mistranslation/frame shift/premature stop where one gene has been split 
                        # into two in a subset of genomes
                        
                        # build a mini graph of allowed pairwise merges
                        tempG = nx.Graph()
                        for nA, nB in itertools.combinations(cluster, 2):
                            mem_inter = sorted(G.nodes[nA]['members'].intersection(G.nodes[nB]['members']))
                            if len(mem_inter) > 0:
                                if distances_bwtn_centroids[centroid_to_index[G.nodes[nA]["longCentroidID"][1]], 
                                    centroid_to_index[G.nodes[nB]["longCentroidID"][1]]]==0:
                                    tempG.add_edge(nA, nB)
                                else:
                                    for imem in mem_inter:
                                        tempids = []
                                        for sid in sorted(G.nodes[nA]['seqIDs'] | G.nodes[nB]['seqIDs']):
                                            if int(sid.split("_")[0])==imem:
                                                tempids.append(sid)
                                        shouldmerge = True
                                        for sidA, sidB in itertools.combinations(tempids, 2):
                                            if abs(int(sidA.split("_")[2])-int(sidB.split("_")[2])) >= len(tempids):
                                                shouldmerge = False
                                            if distances_bwtn_centroids[seqid_to_index[sidA],seqid_to_index[sidB]]==1: 
                                                shouldmerge=False
                                        if shouldmerge:
                                            tempG.add_edge(nA, nB)
                            else:
                                tempG.add_edge(nA, nB)

                        # merge from largest clique to smallest
                        clique = max_clique(tempG)
                        while len(clique)>1:
                            node_count += 1
                            for neig in clique:
                                removed_nodes.add(neig)
                                if neig in search_space:
                                    search_space.remove(neig)
                            
                            temp_c = clique.copy()
                            G = merge_nodes(G,
                                            temp_c.pop(),
                                            temp_c.pop(),
                                            node_count,
                                            multi_centroid=(not correct_mistranslations),
                                            check_merge_mems=False)
                            while (len(temp_c) > 0):
                                G = merge_nodes(G,
                                                node_count,
                                                temp_c.pop(),
                                                node_count + 1,
                                                multi_centroid=(not correct_mistranslations),
                                                check_merge_mems=False)
                                node_count += 1
                            search_space.add(node_count)
                            tempG.remove_nodes_from(clique)
                            clique = max_clique(tempG)
                
                if node in search_space:
                    search_space.remove(node)

    return G, distances_bwtn_centroids, centroid_to_index

def collapse_paralogs(G, centroid_contexts, max_context=5, quiet=False): 
    
    # contexts [centroid] = [[node, member, contig, context], ...]
    node_count = max(list(G.nodes())) + 10

    # first sort by context length, context dist to ensure ties
    #  are broken the same way
    for centroid in centroid_contexts:
        centroid_contexts[centroid] = sorted(centroid_contexts[centroid])

    # set up for context search
    centroid_to_index = {}
    ncentroids=-1
    for node in G.nodes():
        centroid = G.nodes[node]['centroid'][0]
        if centroid not in centroid_to_index:
            ncentroids += 1
            centroid_to_index[centroid] = ncentroids
            centroid_to_index[G.nodes[node]['centroid'][0]] = ncentroids
        else:
            centroid_to_index[G.nodes[node]['centroid'][0]] = centroid_to_index[centroid]
    ncentroids += 1

    for centroid in tqdm(centroid_contexts):
        # calculate distance
        # d = 1 - 1/(abs(contextA-contextB))
        member_paralogs = defaultdict(list)
        for para in centroid_contexts[centroid]:
            member_paralogs[para[1]].append(para)

        ref_paralogs = max(member_paralogs.items(), key=lambda x: len(x[1]))[1]
        # for each paralog find its closest reference paralog
        cluster_dict = defaultdict(set)
        cluster_mems = defaultdict(set)
        for c, ref in enumerate(ref_paralogs):
            cluster_dict[c].add(ref[0])
            cluster_mems[c].add(ref[1])

        # print("here1")
        # row,col,data = zip(*((u,v,1) for u,v in G.edges()))
        # spG = csc_matrix((data + data, (row + col, col + row)), 
        #     shape=(node_count+2, node_count+2))
        # ids  = [ref[0] for ref in ref_paralogs] + [para[0] for para in centroid_contexts[centroid]]
        # spath = shortest_path(csgraph=spG, directed=False, indices=np.array(ids))
        # Gt = gt.Graph()
        # Gt.add_edge_list(G.edges())
        # print("here2")

        for para in centroid_contexts[centroid]:
            # print("here3")
            d_max = np.inf
            s_max = -np.inf
            best_cluster = None

            if para[1]==ref_paralogs[0][1]:
                # this is the reference so skip
                continue

            # first attempt by shortest path            
            for c, ref in enumerate(ref_paralogs):
                if para[1] in cluster_mems[c]:
                    #dont match paralogs of the same isolate
                    continue
                # d = spath[para[0], ref[0]]
                # d = gt.shortest_distance(Gt, para[0], ref[0])
                try:
                    d = nx.shortest_path_length(G, ref[0], para[0])
                except nx.NetworkXNoPath:
                    continue
                if d<d_max:
                    d_max = d
                    best_cluster = c

            # if this fails use context
            if d_max==np.inf:
                best_cluster = 0
                s_max = -np.inf
                para_context = np.zeros(ncentroids)
                for u, node, depth in mod_bfs_edges(G, para[0], max_context):
                    para_context[centroid_to_index[G.nodes[node]['centroid'][0]]] = depth
                for c, ref in enumerate(ref_paralogs):
                    if para[1] in cluster_mems[c]:
                        #dont match paralogs of the same isolate
                        continue
                    ref_context = np.zeros(ncentroids)
                    for u, node, depth in mod_bfs_edges(G, ref[0], max_context):
                        ref_context[centroid_to_index[G.nodes[node]['centroid'][0]]] = depth
                    s = np.sum(1/(1+np.abs((para_context - ref_context)[(para_context*ref_context)!=0])))
                    if s>s_max:
                        s_max=s
                        best_cluster = c
            
            cluster_dict[best_cluster].add(para[0])
            cluster_mems[best_cluster].add(para[1])
        
        # merge
        for cluster in cluster_dict:
            if len(cluster_dict[cluster])<2: continue
            temp_c = list(cluster_dict[cluster].copy())
            node_count += 1
            G = merge_nodes(G, temp_c.pop(), temp_c.pop(), node_count)            
            while (len(temp_c) > 0):
                G = merge_nodes(G, node_count, temp_c.pop(), node_count + 1)    
                node_count += 1
                

    return(G)

def merge_paralogs(G):

    node_count = max(list(G.nodes())) + 10

    # group paralog nodes by centroid
    paralog_centroid_dict = defaultdict(list)
    for node in G.nodes():
        if G.nodes[node]['paralog']:
            paralog_centroid_dict[G.nodes[node]['centroid'][0]].append(node)

    # merge paralog nodes that share the same centroid
    for centroid in paralog_centroid_dict:
        node_count += 1
        temp_c = paralog_centroid_dict[centroid]
        G = merge_nodes(G, temp_c.pop(), temp_c.pop(), node_count)
        while (len(temp_c) > 0):
            G = merge_nodes(G, node_count, temp_c.pop(), node_count + 1)
            node_count += 1

    return (G)


def clean_misassembly_edges(G, edge_support_threshold):

    bad_edges = set()
    max_weight = 0

    # remove edges with low support near contig ends
    for node in G.nodes():
        max_weight = max(max_weight, G.nodes[node]['size'])
        for neigh in G.neighbors(node):
            if G.nodes[neigh]['hasEnd']:
                if G[node][neigh]['weight'] < edge_support_threshold:
                    bad_edges.add((node, neigh))

    # remove edges that have much lower support than the nodes they connect
    for edge in G.edges():
        if float(G.edges[edge]['weight']) < (0.05 * min(
                int(G.nodes[edge[0]]['size']), int(G.nodes[edge[1]]['size']))):
            if float(G.edges[edge]['weight']) < edge_support_threshold:
                bad_edges.add(edge)

    for edge in bad_edges:
        if G.has_edge(edge[0], edge[1]):
            G.remove_edge(edge[0], edge[1])

    return (G)


def identify_possible_highly_variable(G,
                                      cycle_threshold_max=20,
                                      cycle_threshold_min=5,
                                      size_diff_threshold=0.5):

    # add family paralog attribute to nodes
    for node in G.nodes():
        G.nodes[node]['highVar'] = 0

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, list(sub_G.nodes())[0])
        complete_basis += [
            set(b) for b in basis if len(b) <= cycle_threshold_max
        ]

    # remove cycles that are too short
    complete_basis = [b for b in complete_basis if len(b) >= 3]

    # merge cycles with more than one node in common (nested)
    if len(complete_basis) < 1:
        return G

    merged_basis = [[1, set(complete_basis[0])]]
    for b in complete_basis[1:]:
        b = set(b)
        merged = False
        for i, mb in enumerate(merged_basis):
            if len(mb[1].intersection(b)) > 1:
                merged = True
                merged_basis[i][0] += 1
                merged_basis[i][1] |= b
        if not merged:
            merged_basis.append([1, b])

    for b in merged_basis:
        if b[0] < cycle_threshold_min: continue
        max_size = max([G.nodes[node]['size'] for node in b[1]])
        for node in b[1]:
            if G.nodes[node]['size'] < (size_diff_threshold * max_size):
                G.nodes[node]['highVar'] = 1

    return G
