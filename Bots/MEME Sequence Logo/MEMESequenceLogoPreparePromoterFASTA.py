
import numpy as np

import xlrd
import xlsxwriter
import csv

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import GenBank
from Bio import motifs

from timeit import default_timer as timer


startTime = timer()






targetGenesOldLocuses = ["SP_1859", "SP_1810"]

# Get Locus and Genome Names of targetGenes
oldLocusToLocuses = {}
oldLocusToGenomeNames = {}
for targetGenesOldLocus in targetGenesOldLocuses:
    oldLocusToLocuses[targetGenesOldLocus] = []
    oldLocusToGenomeNames[targetGenesOldLocus] = []

lineCount = 0
with open("../../Data/All_RNAseq.csv") as csvFile:
    csvReader = csv.reader(csvFile, delimiter = ",")
    # for main Clusters and Sub_clusters
    for row in csvReader:
        if lineCount == 0:
            print("first Line of All_RNAseq.csv")
        else: 
            for targetGenesOldLocus in targetGenesOldLocuses:
                if row[15] == targetGenesOldLocus:
                    oldLocusToLocuses[targetGenesOldLocus].append(row[4]) 
                    oldLocusToGenomeNames[targetGenesOldLocus].append(row[3]) 
        lineCount += 1



GBFiles = ["CP035234.gbk","CP035235.gbk","CP035236.gbk","CP035237.gbk","CP035238.gbk","CP035239.gbk","CP035240.gbk","CP035241.gbk","CP035242.gbk","CP035243.gbk",
"CP035244.gbk","CP035245.gbk","CP035246.gbk","CP035247.gbk","CP035248.gbk","CP035249.gbk","CP035250.gbk","CP035251.gbk","CP035252.gbk","CP035253.gbk",
"CP035254.gbk","CP035255.gbk","CP035256.gbk","CP035257.gbk","CP035258.gbk","CP035259.gbk","CP035260.gbk","CP035261.gbk","CP035262.gbk","CP035263.gbk",
"CP035264.gbk","CP035265.gbk","CP038251.gbk","CP038252.gbk","CP038253.gbk","SCPQ00000000.gbk","SCPR00000000.gbk","SCPS00000000.gbk","SCPT00000000.gbk"]
GBNames = [GBFile[:len(GBFile) - 4] for GBFile in GBFiles]


# Start all output fasta Files
fastaFiles = {}
for targetGenesOldLocus in targetGenesOldLocuses:

    fastaFiles[targetGenesOldLocus] = open(r"OrganizedFastaFiles/" + targetGenesOldLocus + r".fasta", "w")


# Real Work
progress = 0
for GBFile in GBFiles:
    progress += 1
    print("Progress: " + str(progress) + "/" + str(len(GBFiles)))
    needToOpen = False
    for targetGenesOldLocus in targetGenesOldLocuses:
        if GBFile[:len(GBFile) - 4] in oldLocusToGenomeNames[targetGenesOldLocus]:
            needToOpen = True
            print("1")
        
    if needToOpen:
        print("2")
        for seq_record in SeqIO.parse("../../Data/Genome Data/GBKs/" + GBFile, "genbank"):
            seq = seq_record.seq
            for targetGenesOldLocus in targetGenesOldLocuses:
                if GBFile[:len(GBFile) - 4] in oldLocusToGenomeNames[targetGenesOldLocus]:
                    print("3")
                    for feature in seq_record.features:
                        if feature.type == "gene":
                            
                            locusTag = feature.qualifiers["locus_tag"][0]
                            # print("4")
                            # print(locusTag)
                            if locusTag in oldLocusToLocuses[targetGenesOldLocus]:
                                print("5")
                                if feature.location.start > 100:
                                    print("6")
                                    fastaFiles[targetGenesOldLocus].write(">" + locusTag + "\n")
                                    fastaFiles[targetGenesOldLocus].write(str(seq[feature.location.start - 100 : feature.location.start]) + "\n")
                                else:
                                    print("promoter region too short")



for fastaFile in fastaFiles.values():
    fastaFile.close()

print()
print("total time: ", timer() - startTime)