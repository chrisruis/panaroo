#Runs PastML to reconstruct gene gain and loss onto a given tree
#To run: python3 pastMLGeneReconstruction.py -g $GenePresenceAbsence.csv -t $Tree

import argparse
from pastml.acr import pastml_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help = "Gene presence absence file produced by convertPanarooPastML.py")
    parser.add_argument("-t", help = "Rooted tree containing all and only samples used in the pangenome reconstruction, with exactly the same sample names")
    args = parser.parse_args()

    genePresenceAbsence = open(args.g).readlines()
    html_compressed = "~/Documents/abscessus/manchester.samples/squeaky/abscessus/cluster3/compressed.html"
    html = "~/Documents/abscessus/manchester.samples/squeaky/abscessus/cluster3/full.html"

    columns = genePresenceAbsence[0].strip().split(",")[1:] #The columns to be reconstructed

    pastml_pipeline(data = args.g,data_sep = ",", columns = columns, tree = args.t, html_compressed = html_compressed, html = html, verbose = True)