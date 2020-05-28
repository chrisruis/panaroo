#Extracts a subset of samples from a Panaroo gene_presence_absence.Rtab file
#Uses the column names to identify samples to be extracted, the names must match exactly with those in the given sample names file
#Only keeps genes if they are present in one or more of the samples in the subset

import argparse
from operator import itemgetter

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help = "gene_presence_absence.Rtab file from Panaroo")
    parser.add_argument("-s", help = "Sample names in a single column, no header")
    parser.add_argument("-o", help = "Name of output file")
    args = parser.parse_args()

    genePresenceAbsence = open(args.g).readlines() #Import the gene presence absence information
    sampleNames = open(args.s).readlines() #Import the sample names
    outFile = open(args.o, "w")

    headerLine = genePresenceAbsence[0].strip().split("\t")

    samplePositions = [] #Will be filled with the positions of the samples in the gene presence absence file

    for sampleName in sampleNames: #Iterate through the names, check if they're in the gene presence absence information and get their position
        if sampleName.strip() in headerLine:
            samplePositions.append(headerLine.index(sampleName.strip()))
        else:
            print sampleName.strip() + " is not in the gene presence absence file"
    
    sortedSamplePositions = sorted(samplePositions) #Sort the positions

    for i,gene in enumerate(genePresenceAbsence): #Iterate through the genes, extract their presence absence in the target samples and save if they are present
        if i == 0:
            outFile.write(gene.strip().split("\t")[0] + "\t" + "\t".join(list(itemgetter(*sortedSamplePositions)(gene.strip().split("\t")))) + "\n")
        else:
            geneSamples = itemgetter(*sortedSamplePositions)(gene.strip().split("\t")) #Extract the presence absence of the gene in the samples
            if "1" in set(geneSamples): #Check if the gene is present in the samples
                outFile.write(gene.strip().split("t")[0] + "\t" + "\t".join(geneSamples) + "\n")
    
    outFile.close()