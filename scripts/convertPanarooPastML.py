#Converts the gene_presence_absence.Rtab file from Panaroo to a gene presence abscence file that can be used in PastML

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", help = "gene_presence_abscence.Rtab file from Panaroo")
    parser.add_argument("-o", help = "Name of output file")
    args = parser.parse_args()

    genePresenceAbsence = open(args.r).readlines() #Import the Rtab file from Panaroo
    outFile = open(args.o,"w")

    numberOfSamples = len(genePresenceAbsence[0].strip().split("\t"))

    reconstructGenes = [] #Will be filled with the genes that are not present in all samples that will be reconstructed

    for group in genePresenceAbsence: #Iterate through the genes to identify those that will be reconstructed
        geneNumber = "".join(group.strip().split("\t")[1:])
        if len(set(geneNumber)) != 1: #Check if the gene is absent in at least one sample
            reconstructGenes.append(group)

    for sample in range(numberOfSamples):
        sampleName = 0
        for gene in reconstructGenes: #Iterate through the genes
            if sampleName == 0:
                outFile.write(gene.strip().split("\t")[sample])
                sampleName = 1
            else:
                outFile.write("," + gene.strip().split("\t")[sample])
        outFile.write("\n")

    outFile.close()