import subprocess
import tempfile
import os
import sys
import re
from collections import defaultdict
import networkx as nx
from Bio.Seq import reverse_complement, Seq
import edlib
import itertools
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from joblib import Parallel, delayed
from tqdm import tqdm
import math

def check_cdhit_version(cdhit_exec='cd-hit'):
    """Checks that cd-hit can be run, and returns version.

    Args:
        cdhit_exec (str)
            Location of cd-hit executable

            [default = 'cd-hit']

    Returns:
        version (int)
            Major version of cd-hit
    """
    p = str(
        subprocess.run(cdhit_exec + ' -h', stdout=subprocess.PIPE, shell=True))
    version = False
    find_ver = re.search(r'CD-HIT version \d+\.\d+', p)
    if find_ver:
        version = float(find_ver[0].split()[-1])
    if not version:
        sys.stderr.write("Need cd-hit to be runnable through: " + cdhit_exec +
                         "\n")
        sys.exit(1)

    return (version)


def run_cdhit(
        input_file,
        output_file,
        id=0.95,
        n_cpu=1,
        s=0.0,  # length difference cutoff (%), default 0.0
        aL=0.0,  # alignment coverage for the longer sequence
        AL=99999999,  # alignment coverage control for the longer sequence
        aS=0.0,  # alignment coverage for the shorter sequence
        AS=99999999,  # alignment coverage control for the shorter sequence
        accurate=True,  # use the slower but more accurate options
        use_local=False,  #whether to use local or global sequence alignment
        quiet=False):

    cmd = "cd-hit"
    cmd += " -T " + str(n_cpu)
    cmd += " -i " + input_file
    cmd += " -o " + output_file
    cmd += " -c " + str(id)
    cmd += " -s " + str(s)
    cmd += " -aL " + str(aL)
    cmd += " -AL " + str(AL)
    cmd += " -aS " + str(aS)
    cmd += " -AS " + str(AS)
    cmd += " -M 0 -d 999"

    if use_local:
        cmd += " -G 0"

    if accurate:
        cmd += " -g 1 -n 2"

    if not quiet:
        print("running cmd: " + cmd)
    else:
        cmd += " > /dev/null"

    subprocess.run(cmd, shell=True, check=True)

    return


def run_cdhit_est(
        input_file,
        output_file,
        id=0.99,
        n_cpu=1,
        s=0.0,  # length difference cutoff (%), default 0.0
        aL=0.0,  # alignment coverage for the longer sequence
        AL=99999999,  # alignment coverage control for the longer sequence
        aS=0.0,  # alignment coverage for the shorter sequence
        AS=99999999,  # alignment coverage control for the shorter sequence
        accurate=True,  # use the slower but more accurate options
        use_local=False,  #whether to use local or global sequence alignment
        strand=1,  # default do both +/+ & +/- alignments if set to 0, only +/+
        print_aln=False,  # print alignment overlap in cluster file
        mask=True,
        quiet=False):

    cmd = "cd-hit-est"
    cmd += " -T " + str(n_cpu)
    cmd += " -i " + input_file
    cmd += " -o " + output_file
    cmd += " -c " + str(id)
    cmd += " -s " + str(s)
    cmd += " -aL " + str(aL)
    cmd += " -AL " + str(AL)
    cmd += " -aS " + str(aS)
    cmd += " -AS " + str(AS)
    cmd += " -r " + str(strand)
    cmd += " -M 0 -d 999"

    if mask:
        cmd += " -mask NX"

    if use_local:
        cmd += " -G 0"

    if accurate:
        cmd += " -g 1 -n 6"

    if print_aln:
        cmd += " -p 1"

    if not quiet:
        print("running cmd: " + cmd)
    else:
        cmd += " > /dev/null"

    subprocess.run(cmd, shell=True, check=True)

    return


def cluster_nodes_cdhit(
        G,
        nodes,
        outdir,
        id=0.95,
        dna=False,
        s=0.0,  # length difference cutoff (%), default 0.0
        aL=0.0,  # alignment coverage for the longer sequence
        AL=99999999,  # alignment coverage control for the longer sequence
        aS=0.0,  # alignment coverage for the shorter sequence
        AS=99999999,  # alignment coverage control for the shorter sequence
        accurate=True,  # use the slower but more accurate options
        use_local=False,  #whether to use local or global sequence alignment
        strand=1,  # default do both +/+ & +/- alignments if set to 0, only +/+
        quiet=False,
        prevent_para=True,
        n_cpu=1):

    # create the files we will need
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file.close()

    with open(temp_input_file.name, 'w') as outfile:
        for node in nodes:
            outfile.write(">" + str(node) + "\n")
            if dna:
                outfile.write(
                    max(G.node[node]["dna"].split(";"), key=len) + "\n")
            else:
                outfile.write(
                    max(G.node[node]["protein"].split(";"), key=len) + "\n")

    # run cd-hit
    if dna:
        run_cdhit_est(input_file=temp_input_file.name,
                      output_file=temp_output_file.name,
                      id=id,
                      s=s,
                      aL=aL,
                      AL=AL,
                      aS=AS,
                      accurate=accurate,
                      use_local=use_local,
                      strand=strand,
                      quiet=quiet,
                      n_cpu=n_cpu)
    else:
        run_cdhit(input_file=temp_input_file.name,
                  output_file=temp_output_file.name,
                  id=id,
                  s=s,
                  aL=aL,
                  AL=AL,
                  aS=AS,
                  accurate=accurate,
                  use_local=use_local,
                  quiet=quiet,
                  n_cpu=n_cpu)

    # process the output
    clusters = []
    with open(temp_output_file.name + ".clstr", 'rU') as infile:
        c = []
        for line in infile:
            if line[0] == ">":
                clusters.append(c)
                c = []
            else:
                c.append(int(line.split(">")[1].split("...")[0]))
        clusters.append(c)
    clusters = clusters[1:]

    # optionally split clusters to ensure we don't collapse paralogs
    if prevent_para:
        nodes = list(nodes)
        # set up node to cluster dict
        cluster_dict = {}
        for i, c in enumerate(clusters):
            for n in c:
                cluster_dict[n] = i

        # set up subgraph and new_cluster dict
        sub_G = G.subgraph(nodes)
        if not nx.is_connected(sub_G):
            raise ValueError("Sub graph is not connected!")

        new_clusters = defaultdict(list)

        # ref node with max size and degree > 2
        ref_node = nodes[0]
        for n in nodes[1:]:
            if sub_G.degree[n] > 2:
                if sub_G.node[n]['size'] >= sub_G.node[ref_node]['size']:
                    ref_node = n

        # nodes in Breadth First Search order
        nodes_BFS = [ref_node] + [v for u, v in nx.bfs_edges(sub_G, ref_node)]

        # iterate through making new clusters that satisfy conditions
        for node in nodes_BFS:
            c1 = cluster_dict[node]
            if len(new_clusters[c1]) < 1:
                new_clusters[c1].append([node])
            else:
                # try and add to first valid cluster
                found = False
                for i, c2 in enumerate(new_clusters[c1]):
                    if is_valid(G, node, c2):
                        new_clusters[c1][i].append(node)
                        found = True
                        break
                if not found:
                    # create a new cluster
                    new_clusters[c1].append([node])

        # collapse dictionary into original list format
        clusters = []
        for c1 in new_clusters:
            for c2 in new_clusters[c1]:
                clusters.append(c2)

    # check all nodes are accounted for
    clust_node_set = set([item for sublist in clusters for item in sublist])
    for node in nodes:
        if node not in clust_node_set:
            print("nodes:", nodes)
            print("clust_node_set:", clust_node_set)
            raise ValueError('Clusters are missing a node!')

    # DEBUG: check
    # for c in clusters:
    #     members = [G.node[n]['members'] for n in c]
    #     members = [item for sublist in members for item in sublist]
    #     seen = set()
    #     for m in members:
    #         if m in seen:
    #             raise ValueError("duplicate members!")
    #         seen.add(m)

    # remove temporary files
    os.remove(temp_input_file.name)
    os.remove(temp_output_file.name)
    os.remove(temp_output_file.name + ".clstr")

    return clusters


def is_valid(G, node, cluster):
    found = True
    for n in cluster:
        if len(set(G.node[node]['members']) & set(G.node[n]['members'])) > 0:
            found = False
    return found


def align_dna_cdhit(
        query,
        target,
        temp_dir,
        id=0.99,
        n_cpu=1,
        s=0.0,  # length difference cutoff (%), default 0.0
        aL=0.0,  # alignment coverage for the longer sequence
        AL=99999999,  # alignment coverage control for the longer sequence
        aS=0.0,  # alignment coverage for the shorter sequence
        AS=99999999,  # alignment coverage control for the shorter sequence
        accurate=True,  # use the slower but more accurate options
        use_local=False,  #whether to use local or global sequence alignment
        strand=1,  # default do both +/+ & +/- alignments if set to 0, only +/+
        mask=True,
        quiet=False):

    # create the files we will need
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir)
    temp_output_file.close()

    # prepare files for cdhit
    with open(temp_input_file.name, 'w') as outfile:
        outfile.write(">query\n" + query + "\n")
        outfile.write(">target\n" + target + "\n")

    # run cdhit
    run_cdhit_est(input_file=temp_input_file.name,
                  output_file=temp_output_file.name,
                  id=id,
                  s=s,
                  aL=aL,
                  AL=AL,
                  aS=AS,
                  accurate=accurate,
                  use_local=use_local,
                  print_aln=True,
                  strand=strand,
                  mask=mask,
                  quiet=True)

    # process resulting alignment
    # process the output
    found_seq = ""
    with open(temp_output_file.name + ".clstr", 'rU') as infile:
        rev = False
        for line in infile:
            if "at" in line:
                align = line.split(" at ")[1].split("/")[0].split(":")
                if "query" in line:
                    if int(align[0]) > int(align[1]): rev = True
                    bounds = sorted([int(align[0]), int(align[1])])
                else:
                    if int(align[2]) > int(align[3]): rev = True
                    bounds = sorted([int(align[2]), int(align[3])])
                found_seq = query[bounds[0] - 1:bounds[1]]
                if rev:
                    found_seq = reverse_complement(found_seq)

    # remove temporary files
    os.remove(temp_input_file.name)
    os.remove(temp_output_file.name)
    os.remove(temp_output_file.name + ".clstr")

    return found_seq

def iterative_cdhit(G,
        nodes,
        outdir,
        dna=False,
        s=0.0,  # length difference cutoff (%), default 0.0
        aL=0.0,  # alignment coverage for the longer sequence
        AL=99999999,  # alignment coverage control for the longer sequence
        aS=0.0,  # alignment coverage for the shorter sequence
        AS=99999999,  # alignment coverage control for the shorter sequence
        accurate=True,  # use the slower but more accurate options
        use_local=False,  #whether to use local or global sequence alignment
        strand=1,  # default do both +/+ & +/- alignments if set to 0, only +/+
        quiet=False,
        thresholds=[0.99,0.95,0.90,0.85,0.8,0.75,0.7],
        n_cpu=1):

    # create the files we will need
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file.close()

    clusters = []
    with open(temp_input_file.name, 'w') as outfile:
        for node in nodes:
            clusters.append([node])
            outfile.write(">" + str(node) + "\n")
            if dna:
                outfile.write(
                    max(G.node[node]["dna"].split(";"), key=len) + "\n")
            else:
                outfile.write(
                    max(G.node[node]["protein"].split(";"), key=len) + "\n")
    
    for cid in thresholds:
        # run cd-hit
        if dna:
            run_cdhit_est(input_file=temp_input_file.name,
                        output_file=temp_output_file.name,
                        id=cid,
                        s=s,
                        aL=aL,
                        AL=AL,
                        aS=AS,
                        accurate=accurate,
                        use_local=use_local,
                        strand=strand,
                        quiet=quiet,
                        n_cpu=n_cpu)
        else:
            run_cdhit(input_file=temp_input_file.name,
                    output_file=temp_output_file.name,
                    id=cid,
                    s=s,
                    aL=aL,
                    AL=AL,
                    aS=AS,
                    accurate=accurate,
                    use_local=use_local,
                    quiet=quiet,
                    n_cpu=n_cpu)

        # process the output
        temp_clusters = []
        with open(temp_output_file.name + ".clstr", 'rU') as infile:
            c = []
            for line in infile:
                if line[0] == ">":
                    temp_clusters.append(c)
                    c = []
                else:
                    c.append(int(line.split(">")[1].split("...")[0]))
            temp_clusters.append(c)
        temp_clusters = temp_clusters[1:]

        # collapse previously clustered
        temp_clust_dict = {}
        for c, clust in enumerate(temp_clusters):
            for n in clust:
                temp_clust_dict[n] = c
        clust_dict = defaultdict(list)  
        for clust in clusters:
            c = -1
            for n in clust:
                if n in temp_clust_dict:
                    c = temp_clust_dict[n]
            for n in clust:
                clust_dict[c].append(n)
        clusters = clust_dict.values()

        # cleanup and rename for next round
        os.remove(temp_input_file.name)
        os.remove(temp_output_file.name + ".clstr")
        temp_input_file.name = temp_output_file.name
        temp_output_file.name = temp_output_file.name + "t" + str(cid)

    return(clusters)


def pwdist_edlib(G, cdhit_clusters, threshold, dna=False, n_cpu=1):
    
    # map nodes to centroids
    node_to_centroid = {}
    centroid_to_index = {}
    node_index = 0
    for node in G.nodes():
        centroid = G.node[node]["centroid"].split(";")[0]
        node_to_centroid[node] = centroid
        if centroid not in centroid_to_index:
            centroid_to_index[centroid] = node_index
            node_index+=1

    # Prepare sequences
    seqs = {}
    for node in G.nodes():
        if node_to_centroid[node] in seqs: continue
        if dna:
            seqs[node_to_centroid[node]] = max(G.node[node]["dna"].split(";"), key=len).upper()
        else:
            seqs[node_to_centroid[node]] = max(G.node[node]["protein"].split(";"), key=len).upper()

    # get pairwise id between sequences in the same cdhit clusters
    distances_bwtn_centroids = defaultdict(lambda: 100)

    all_distances = []
    for cluster in tqdm(cdhit_clusters):
        all_distances += Parallel(n_jobs=n_cpu)(
                    delayed(run_pw)(seqs[node_to_centroid[n1]], seqs[node_to_centroid[n1]], n1, n2, dna)
                    for n1, n2 in itertools.combinations(cluster, 2))
    
    data=[]
    row_ind=[]
    col_ind=[]
    for d in all_distances:
        if d[2] >= threshold:
            data.append(1)
            row_ind.append(centroid_to_index[node_to_centroid[d[0]]])
            col_ind.append(centroid_to_index[node_to_centroid[d[1]]])

    distances_bwtn_centroids = csr_matrix((data,
        (row_ind, 
         col_ind)),
         shape=(node_index, node_index))

    return distances_bwtn_centroids, centroid_to_index

def run_pw(seqA, seqB, n1, n2, dna):

    if len(seqA) > len(seqB):
        seqA, seqB = seqB, seqA

    if dna:
        pwid = 0.0
        for sA in [seqA, str(Seq(seqA).reverse_complement())]:
            aln = edlib.align(sA, seqB, mode="HW", task='locations', k=10*len(seqA),
                    additionalEqualities=[('A','N'), ('C','N'), ('G','N'), ('T','N')])
        pwid = max(pwid, 1.0-aln['editDistance']/float(len(seqA)))
    else:
        aln = edlib.align(seqA, seqB, mode="HW", task='locations', k=10*len(seqA),
                    additionalEqualities=[('*', 'X'),('A', 'X'),('C', 'X'),('B', 'X'),('E', 'X'),
                                            ('D', 'X'),('G', 'X'),('F', 'X'),('I', 'X'),('H', 'X'),
                                            ('K', 'X'),('M', 'X'),('L', 'X'),('N', 'X'),('Q', 'X'),
                                            ('P', 'X'),('S', 'X'),('R', 'X'),('T', 'X'),('W', 'X'),
                                            ('V', 'X'),('Y', 'X'),('X', 'X'),('Z', 'X'),('D', 'B'),
                                            ('N', 'B'),('E', 'Z'),('Q', 'Z')])
        pwid = 1.0-aln['editDistance']/float(len(seqA))

    return((n1,n2,pwid))


