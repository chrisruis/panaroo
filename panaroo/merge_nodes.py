import itertools
from collections import Counter
from .isvalid import del_dups

def gen_node_iterables(G, nodes, feature, split=None):
    for n in nodes:
        if split is None:
            yield G.nodes[n][feature]
        else:
            yield G.nodes[n][feature].split(split)

def temp_iter(list_list):
    for n in list_list:
        yield n

def iter_del_dups(iterable):
    seen = {}
    for f in itertools.chain.from_iterable(iterable):
        seen[f] = None
    return (list(seen.keys()))

def merge_node_cluster(G,
                nodes,
                newNode,
                multi_centroid=True,
                check_merge_mems=True):

    if check_merge_mems:
        mem_count = Counter(itertools.chain.from_iterable(gen_node_iterables(G, nodes, 'members')))
        if max(mem_count.values()) > 1:
            raise ValueError("merging nodes with the same genome IDs!")

    # take node with most support as the 'consensus'
    nodes = sorted(nodes, key=lambda x: G.nodes[x]['size'])

    # First create a new node and combine the attributes
    dna = iter_del_dups(gen_node_iterables(G,nodes,'dna'))
    maxLenId = 0
    max_l = 0
    for i, s in enumerate(dna):
        if len(s) >= max_l:
            max_l = len(s)
            maxLenId = i
    members = set(iter_del_dups(gen_node_iterables(G,nodes,'members')))

    if multi_centroid:
        mergedDNA = any(gen_node_iterables(G,nodes,'mergedDNA'))
    else:
        mergedDNA = True

    G.add_node(
        newNode,
        size=len(members),
        centroid=iter_del_dups(gen_node_iterables(G,nodes,'centroid')),
        maxLenId=maxLenId,
        members=members,
        seqIDs=set(iter_del_dups(gen_node_iterables(G,nodes,'seqIDs'))),
        hasEnd=any(gen_node_iterables(G,nodes,'hasEnd')),
        protein=iter_del_dups(gen_node_iterables(G,nodes,'protein')),
        dna=dna,
        annotation=";".join(iter_del_dups(gen_node_iterables(G,nodes,'annotation',split=";"))),
        description=";".join(iter_del_dups(gen_node_iterables(G,nodes,'description',split=";"))),
        lengths=list(itertools.chain.from_iterable(gen_node_iterables(G,nodes,'lengths'))),
        longCentroidID=max(gen_node_iterables(G,nodes,'longCentroidID')),
        paralog=any(gen_node_iterables(G,nodes,'paralog')),
        mergedDNA=mergedDNA
    )
    if "prevCentroids" in G.nodes[nodes[0]]:
        G.nodes[newNode]['prevCentroids'] = ";".join(
            set(iter_del_dups(gen_node_iterables(G,nodes,'prevCentroids',split=";"))))

    # Now iterate through neighbours of each node and add them to the new node
    merge_nodes = set(nodes)
    for node in nodes:
        for neighbour in G.neighbors(node):
            if neighbour in merge_nodes: continue
            if G.has_edge(newNode, neighbour):
                G[newNode][neighbour]['size'] += G[node][neighbour]['size']
                G[newNode][neighbour]['members'] |= G[node][neighbour]['members']
            else:
                G.add_edge(newNode, neighbour,
                       size=G[node][neighbour]['size'],
                       members=G[node][neighbour]['members'])

    # remove old nodes from Graph
    G.remove_nodes_from(nodes)

    return G

def merge_nodes(G,
                nodeA,
                nodeB,
                newNode,
                multi_centroid=True,
                check_merge_mems=True):

    if check_merge_mems:
        if len(G.nodes[nodeA]['members'] & G.nodes[nodeB]['members']) > 0:
            raise ValueError("merging nodes with the same genome IDs!")

    # take node with most support as the 'consensus'
    if G.nodes[nodeA]['size'] < G.nodes[nodeB]['size']:
        nodeB, nodeA = nodeA, nodeB

    # First create a new node and combine the attributes
    dna = del_dups(G.nodes[nodeA]['dna'] + G.nodes[nodeB]['dna'])
    maxLenId = 0
    max_l = 0
    for i, s in enumerate(dna):
        if len(s) >= max_l:
            max_l = len(s)
            maxLenId = i

    if multi_centroid:
        G.add_node(
            newNode,
            size=len(G.nodes[nodeA]['members'] | G.nodes[nodeB]['members']),
            centroid=del_dups(G.nodes[nodeA]['centroid'] +
                              G.nodes[nodeB]['centroid']),
            maxLenId=maxLenId,
            members=G.nodes[nodeA]['members'] | G.nodes[nodeB]['members'],
            seqIDs=G.nodes[nodeA]['seqIDs'] | G.nodes[nodeB]['seqIDs'],
            hasEnd=(G.nodes[nodeA]['hasEnd'] or G.nodes[nodeB]['hasEnd']),
            protein=del_dups(G.nodes[nodeA]['protein'] +
                             G.nodes[nodeB]['protein']),
            dna=dna,
            annotation=";".join(
                del_dups(G.nodes[nodeA]['annotation'].split(";") +
                         G.nodes[nodeB]['annotation'].split(";"))),
            description=";".join(
                del_dups(G.nodes[nodeA]['description'].split(";") +
                         G.nodes[nodeB]['description'].split(";"))),
            lengths=G.nodes[nodeA]['lengths'] + G.nodes[nodeB]['lengths'],
            longCentroidID=max(G.nodes[nodeA]['longCentroidID'],
                               G.nodes[nodeB]['longCentroidID']),
            paralog=(G.nodes[nodeA]['paralog'] or G.nodes[nodeB]['paralog']),
            mergedDNA=(G.nodes[nodeA]['mergedDNA']
                       or G.nodes[nodeB]['mergedDNA']))
        if "prevCentroids" in G.nodes[nodeA]:
            G.nodes[newNode]['prevCentroids'] = ";".join(
                set(G.nodes[nodeA]['prevCentroids'].split(";") +
                    G.nodes[nodeB]['prevCentroids'].split(";")))
    else:
        G.add_node(
            newNode,
            size=len(G.nodes[nodeA]['members'] | G.nodes[nodeB]['members']),
            centroid=del_dups(G.nodes[nodeA]['centroid'] +
                              G.nodes[nodeB]['centroid']),
            maxLenId=maxLenId,
            members=G.nodes[nodeA]['members'] | G.nodes[nodeB]['members'],
            seqIDs=G.nodes[nodeA]['seqIDs'] | G.nodes[nodeB]['seqIDs'],
            hasEnd=(G.nodes[nodeA]['hasEnd'] or G.nodes[nodeB]['hasEnd']),
            protein=del_dups(G.nodes[nodeA]['protein'] +
                             G.nodes[nodeB]['protein']),
            dna=dna,
            annotation=G.nodes[nodeA]['annotation'],
            description=G.nodes[nodeA]['description'],
            paralog=(G.nodes[nodeA]['paralog'] or G.nodes[nodeB]['paralog']),
            lengths=G.nodes[nodeA]['lengths'] + G.nodes[nodeB]['lengths'],
            longCentroidID=max(G.nodes[nodeA]['longCentroidID'],
                               G.nodes[nodeB]['longCentroidID']),
            mergedDNA=True)
        if "prevCentroids" in G.nodes[nodeA]:
            G.nodes[newNode]['prevCentroids'] = ";".join(
                set(G.nodes[nodeA]['prevCentroids'].split(";") +
                    G.nodes[nodeB]['prevCentroids'].split(";")))

    # Now iterate through neighbours of each node and add them to the new node
    neigboursB = list(G.neighbors(nodeB))
    neigboursA = list(G.neighbors(nodeA))
    for neighbor in neigboursA:
        if neighbor in neigboursB:
            G.add_edge(newNode,
                       neighbor,
                       size=G[nodeA][neighbor]['size'] +
                       G[nodeB][neighbor]['size'],
                       members=G[nodeA][neighbor]['members']
                       | G[nodeB][neighbor]['members'])
            neigboursB.remove(neighbor)
        else:
            G.add_edge(newNode,
                       neighbor,
                       size=G[nodeA][neighbor]['size'],
                       members=G[nodeA][neighbor]['members'])

    for neighbor in neigboursB:
        G.add_edge(newNode,
                   neighbor,
                   size=G[nodeB][neighbor]['size'],
                   members=G[nodeB][neighbor]['members'])

    # remove old nodes from Graph
    G.remove_nodes_from([nodeA, nodeB])

    if len(max(G.nodes[newNode]["dna"], key=len)) <= 0:
        print(G.nodes[newNode]["dna"])
        raise NameError("Problem!")

    return G


def delete_node(G, node):
    # add in new edges
    for mem in G.nodes[node]['members']:
        mem_edges = list(
            set([e[1] for e in G.edges(node) if mem in G.edges[e]['members']]))
        if len(mem_edges) < 2: continue
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] |= set([mem])
                G[n1][n2]['size'] += 1
            else:
                G.add_edge(n1, n2, size=1, members=set([mem]))

    # now remove node
    G.remove_node(node)

    return G


def remove_member_from_node(G, node, member):

    # add in replacement edges if required
    mem_edges = list(
        set([e[1] for e in G.edges(node) if member in G.edges[e]['members']]))
    if len(mem_edges) > 1:
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] |= [member]
                G[n1][n2]['size'] += 1
            else:
                G.add_edge(n1, n2, size=1, members=set([member]))

    # remove member from node
    while str(member) in G.nodes[node]['members']:
        G.nodes[node]['members'].remove(str(member))
    G.nodes[node]['seqIDs'] = set([
        sid for sid in G.nodes[node]['seqIDs']
        if sid.split("_")[0] != str(member)
    ])
    G.nodes[node]['size'] -= 1

    return G