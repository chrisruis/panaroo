import networkx as nx
from panaroo.cdhit import *
from panaroo.merge_nodes import merge_node_cluster 
from panaroo.isvalid import del_dups, max_clique
from collections import defaultdict, deque
from panaroo.isvalid import single_source_shortest_path_length_mod
from itertools import chain, combinations
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
from tqdm import tqdm
from time import time

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


def single_linkage(G, distances_bwtn_centroids, centroid_to_index, neighbours):
    index = []
    neigh_array = []
    for neigh in neighbours:
        for sid in G.nodes[neigh]['centroid']:
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
        if len(l) > 1:
            for i in l[1:]:
                labels[labels == i] = l[0]

    clusters = [
        del_dups(list(neigh_array[labels == i])) for i in np.unique(labels)
    ]

    return (clusters)


def collapse_families(G,
                      seqid_to_centroid,
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
        cdhit_clusters = iterative_cdhit(G,
                                         outdir,
                                         thresholds=threshold,
                                         n_cpu=n_cpu,
                                         quiet=True,
                                         dna=True,
                                         word_length=7,
                                         accurate=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, dna_error_threshold, dna=True, n_cpu=n_cpu)
    elif distances_bwtn_centroids is None:
        cdhit_clusters = iterative_cdhit(G,
                                         outdir,
                                         thresholds=threshold,
                                         n_cpu=n_cpu,
                                         quiet=True,
                                         dna=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, family_threshold, dna=False, n_cpu=n_cpu)

    print("calculated pairwise distances..")

    # keep track of centroids for each sequence. Need this to resolve clashes
    seqid_to_index = {}
    for node in G.nodes():
        for sid in G.nodes[node]['seqIDs']:
            if "refound" in sid:
                seqid_to_index[sid] = centroid_to_index[G.nodes[node]
                                                        ["longCentroidID"][1]]
            else:
                seqid_to_index[sid] = centroid_to_index[seqid_to_centroid[sid]]

    nonzero_dist = distances_bwtn_centroids.nonzero()
    nonzero_dist = set([(i, j)
                        for i, j in zip(nonzero_dist[0], nonzero_dist[1])])

    # create a dictionary to of node->member->centroids to improve performance
    # TODO: incoporate this into new node structure in major refactor
    node_mem_index_dict = defaultdict(lambda: defaultdict(set))
    for n in G.nodes():
        for sid in G.nodes[n]['seqIDs']:
            node_mem_index_dict[n][sid.split("_")[0]].add(seqid_to_index[sid])


    print("starting to collapse")
    tim_sl = 0
    tim_set_up_minig = 0
    tim_clique_w_merge=0

    for depth in depths:
        print("processing depth:", depth)
        search_space = set(G.nodes())
        while len(search_space) > 0:
            # look for nodes to merge
            temp_node_list = list(search_space)
            removed_nodes = set()

            if tim_sl>0: print("tim_sinlge link:", tim_sl)
            if tim_set_up_minig>0: print("tim_set_up_minig:", tim_set_up_minig)
            if tim_clique_w_merge>0: print("tim_clique_w_merge:", tim_clique_w_merge)
            tim_sl = 0
            tim_set_up_minig = 0
            tim_clique_w_merge=0

            for node in tqdm(temp_node_list):
                if node in removed_nodes: continue

                if G.degree[node] <= 2:
                    search_space.remove(node)
                    removed_nodes.add(node)
                    continue

                # find neighbouring nodes and cluster their centroid with cdhit
                neighbours = [
                    v
                    for u, v in nx.bfs_edges(G, source=node, depth_limit=depth)
                ] + [node]

                # find clusters
                tim=time()
                clusters = single_linkage(G, distances_bwtn_centroids,
                                          centroid_to_index, neighbours)
                tim_sl += time()-tim

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
                        G = merge_node_cluster(G, cluster, node_count,
                            multi_centroid=(not correct_mistranslations))
                        for n in cluster:
                            for mem in node_mem_index_dict[n]:
                                node_mem_index_dict[node_count][mem] |= node_mem_index_dict[n][mem]
                        search_space.add(node_count)
                    else:
                        # merge if the centroids don't conflict and the nodes are adjacent in the conflicting genome
                        # this corresponds to a mistranslation/frame shift/premature stop where one gene has been split
                        # into two in a subset of genomes

                        # build a mini graph of allowed pairwise merges
                        tim=time()

                        tempG = nx.Graph()
                        for nA, nB in itertools.combinations(cluster, 2):
                            mem_inter = G.nodes[nA]['members'].intersection(
                                G.nodes[nB]['members'])
                            if len(mem_inter) > 0:

                                shouldmerge = True
                                if len(
                                        set(G.nodes[nA]['centroid']).
                                        intersection(
                                            set(G.nodes[nB]['centroid']))) > 0:
                                    shouldmerge = False

                                    if shouldmerge:
                                        for imem in mem_inter:
                                            for sidA in node_mem_index_dict[nA][imem]:
                                                for sidB in node_mem_index_dict[nB][imem]:
                                                    if (
                                                        (sidA,
                                                        sidB) in nonzero_dist
                                                    ) or ((sidB,
                                                        sidA) in nonzero_dist):
                                                        shouldmerge = False
                                                        break
                                                if not shouldmerge: break
                                            if not shouldmerge: break

                                if shouldmerge:
                                    tempG.add_edge(nA, nB)
                            else:
                                tempG.add_edge(nA, nB)

                                # if shouldmerge:
                                #     idsA = defaultdict(list)
                                #     for sid in G.nodes[nA]['seqIDs']:
                                #         ssid = sid.split("_")
                                #         if ssid[0] in mem_inter:
                                #             idsA[ssid[0]].append(sid)

                                #     idsB = defaultdict(list)
                                #     for sid in G.nodes[nB]['seqIDs']:
                                #         ssid = sid.split("_")
                                #         if ssid[0] in mem_inter:
                                #             idsB[ssid[0]].append(sid)

                                #     for imem in mem_inter:
                                #         for sidA in set([
                                #                 seqid_to_index[sid]
                                #                 for sid in idsA[imem]
                                #         ]):
                                #             for sidB in set([
                                #                     seqid_to_index[sid]
                                #                     for sid in idsB[imem]
                                #             ]):
                                #                 if (
                                #                     (sidA,
                                #                      sidB) in nonzero_dist
                                #                 ) or ((sidB,
                                #                        sidA) in nonzero_dist):
                                #                     shouldmerge = False
                                #                     break
                                #             if not shouldmerge: break
                                #         if not shouldmerge: break

                        tim_set_up_minig += time()-tim
                        tim=time()

                        # merge from largest clique to smallest
                        sys.setrecursionlimit(max(len(tempG.nodes), 10000))
                        clique = max_clique(tempG)
                        while len(clique) > 1:
                            clique_clusters = single_linkage(
                                G, distances_bwtn_centroids, centroid_to_index,
                                clique)
                            for clust in clique_clusters:
                                if len(clust) <= 1: continue
                                node_count += 1
                                for neig in clust:
                                    removed_nodes.add(neig)
                                    if neig in search_space:
                                        search_space.remove(neig)

                                G = merge_node_cluster(G, clust, node_count,
                                        multi_centroid=(not correct_mistranslations),
                                        check_merge_mems=False)
                                for n in clust:
                                    for mem in node_mem_index_dict[n]:
                                        node_mem_index_dict[node_count][mem] |= node_mem_index_dict[n][mem]
                                search_space.add(node_count)
                            tempG.remove_nodes_from(clique)
                            clique = max_clique(tempG)
                        
                        tim_clique_w_merge += time()-tim

                if node in search_space:
                    search_space.remove(node)

    return G, distances_bwtn_centroids, centroid_to_index


def collapse_paralogs(G, centroid_contexts, max_context=2, quiet=False):

    node_count = max(list(G.nodes())) + 10

    # first sort by context length, context dist to ensure ties
    #  are broken the same way
    for centroid in centroid_contexts:
        centroid_contexts[centroid] = sorted(centroid_contexts[centroid])

    # set up for context search
    centroid_to_index = {}
    ncentroids = -1
    for node in G.nodes():
        centroid = G.nodes[node]['centroid'][0]
        if centroid not in centroid_to_index:
            ncentroids += 1
            centroid_to_index[centroid] = ncentroids
            centroid_to_index[G.nodes[node]['centroid'][0]] = ncentroids
        else:
            centroid_to_index[G.nodes[node]['centroid']
                              [0]] = centroid_to_index[centroid]
    ncentroids += 1

    for centroid in tqdm(centroid_contexts, disable=quiet):
        print("preparing..")
        print("Number of nodes in G:", len(G.nodes()))
        tim=time()
        # calculate distance
        member_paralogs = defaultdict(list)
        para_nodes = set()
        for para in centroid_contexts[centroid]:
            member_paralogs[para[1]].append(para)
            para_nodes.add(para[0])
        para_nodes = sorted(list(para_nodes), key=lambda x: G.nodes[x]['size'], reverse=True)

        ref_paralogs = max(member_paralogs.items(), key=lambda x: len(x[1]))[1]
        ref_paralogs = [ref[0] for ref in ref_paralogs]
        ref_paralogs = sorted(ref_paralogs, key=lambda x: G.nodes[x]['size'], reverse=True)

        ref_path_lengths = []
        for ref in ref_paralogs:
            ref_path_lengths.append(single_source_shortest_path_length_mod(G, ref, para_nodes))

        # for each paralog find its closest reference paralog
        cluster_dict = defaultdict(set)
        cluster_mems = defaultdict(set)
        for c, ref in enumerate(ref_paralogs):
            cluster_dict[c].add(ref)
            cluster_mems[c] = G.nodes[ref]['members']
        print("prepared in ", time()-tim)

        print("finding location..")
        tim_sp = 0
        n_sp_searches = 0
        tim_cont = 0
        for ip, para in enumerate(para_nodes):
            if para in ref_paralogs: continue # as this is the reference
            distances = []
            
            # calculate path distances
            tim=time()
            for c, ref in enumerate(ref_paralogs):
                n_sp_searches+=1
                if para in ref_path_lengths[c]:
                    distances.append((ref_path_lengths[c][para], c))
                else:
                    distances.append((np.inf, c))
                # try:
                #     distances.append((nx.shortest_path_length(G, ref, para), c))
                # except nx.NetworkXNoPath:
                #     distances.append((np.inf, c))
            
            distances = sorted(distances)
            tim_sp += time()-tim

            # check if a merge is possible
            best_cluster = None
            for dist, c in distances:
                if dist==np.inf: break
                if len(cluster_mems[c].intersection(G.nodes[para]['members']))==0:
                    best_cluster = c
                    cluster_mems[c] = cluster_mems[c] | G.nodes[para]['members']
                    break
            
            # if merge based on shortest path failed use context
            tim=time()
            if best_cluster is None:
                distances = []

                s_max = -np.inf
                para_context = np.zeros(ncentroids)
                for u, node, depth in mod_bfs_edges(G, para, max_context):
                    para_context[centroid_to_index[G.nodes[node]['centroid']
                                                   [0]]] = depth
                
                for c, ref in enumerate(ref_paralogs):
                    if para in G.nodes[ref_paralogs[c]]['members']:
                        #dont match paralogs of the same isolate
                        continue
                    ref_context = np.zeros(ncentroids)
                    for u, node, depth in mod_bfs_edges(
                            G, ref, max_context):
                        ref_context[centroid_to_index[G.nodes[node]['centroid']
                                                      [0]]] = depth
                    d = 1 - np.sum(1 / (1 + np.abs((para_context - ref_context)[
                        (para_context * ref_context) != 0])))
                    distances.append((d, c))

                distances = sorted(distances)
                for dist, c in distances:
                    if dist==np.inf: break
                    if len(cluster_mems[c].intersection(G.nodes[para]['members']))==0:
                        best_cluster = c
                        cluster_mems[c] = cluster_mems[c] | G.nodes[para]['members']
                        break

            tim_cont += time()-tim

            if best_cluster is None:
                # we couldn't merge due to conflict so add additional reference node.
                ref_paralogs.append(para)
                ref_path_lengths.append(single_source_shortest_path_length_mod(G, para, para_nodes[ip:]))
                cluster_mems[len(ref_paralogs)-1] = G.nodes[para]['members']
                best_cluster = len(ref_paralogs)-1

            cluster_dict[best_cluster].add(para)
        
        print("shortest path ", tim_sp)
        print("number shortest path searches ", n_sp_searches)
        print("context ", tim_cont)

        # merge
        print("merging..")
        tim=time()
        for cluster in cluster_dict:
            if len(cluster_dict[cluster]) < 2: continue
            node_count += 1
            G = merge_node_cluster(G,
                    cluster_dict[cluster],
                    node_count)
        print("merge completed in ", time()-tim)


    return (G)


def merge_paralogs(G):

    node_count = max(list(G.nodes())) + 10

    # group paralog nodes by centroid
    paralog_centroids = defaultdict(list)
    for node in G.nodes():
        if G.nodes[node]['paralog']:
            for centroid in G.nodes[node]['centroid']:
                paralog_centroids[centroid].append(node)

    # find nodes that share common centroids
    paralog_centroids = paralog_centroids.values()
    merge_clusters = []
    while len(paralog_centroids) > 0:
        first, *rest = paralog_centroids
        first = set(first)
        lf = -1
        while len(first) > lf:
            lf = len(first)
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2
        merge_clusters.append(first)
        paralog_centroids = rest

    # merge paralog nodes that share the same centroid
    for temp_c in merge_clusters:
        if len(temp_c) > 1:
            node_count += 1
            G = merge_node_cluster(G, temp_c, node_count,
                            check_merge_mems=False)

    return (G)


def clean_misassembly_edges(G, edge_support_threshold):

    bad_edges = set()
    max_weight = 0

    # remove edges with low support near contig ends
    for node in G.nodes():
        max_weight = max(max_weight, G.nodes[node]['size'])
        for neigh in G.neighbors(node):
            if G.nodes[neigh]['hasEnd']:
                if G[node][neigh]['size'] < edge_support_threshold:
                    bad_edges.add((node, neigh))

    # remove edges that have much lower support than the nodes they connect
    for edge in G.edges():
        if float(G.edges[edge]['size']) < (0.05 * min(
                int(G.nodes[edge[0]]['size']), int(G.nodes[edge[1]]['size']))):
            if float(G.edges[edge]['size']) < edge_support_threshold:
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
