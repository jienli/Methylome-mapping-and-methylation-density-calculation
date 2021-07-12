### Since v1, added circular and linear graphing/plotting functionality

# import sys
# sys.path.append('/usr/local/lib/python3.4/site-packages')


import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import GenBank
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
import FindMotif
import GraphCircularGenome as Graph

"""
Use GenBank file (.gbk)
Format: ...
"""

###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################



# my_seq = Seq("ATGCATGC")
# my_seq


# motifs = Motif.Motif(alphabet=IUPAC.ambiguous_dna)
# motif

### Seq with IUPAC.ambiguous_dna does not work, it's just a tag with no real functionality:

# # x = Seq("AATTCTCT", IUPAC.Alphabet)
# # y = Seq("TTN", IUPAC.Alphabet)
# # x = Seq("AAT", generic_dna)
# # y = Seq("AAN", generic_dna)
# x = Seq("AAT", IUPAC.unambiguous_dna)
# y = Seq("AAN", IUPAC.ambiguous_dna)
# # x.find(y)
# x == y


### motifs package also does not work because IUPAC.ambiguous_dna is not considered funionally

# z = motifs.create([Seq("AANA", alphabet=IUPAC.ambiguous_dna)])
# z.degenerate_consensus
# z.counts
# # z.instances.search(Seq("GCGCAAN"))
# pos = 0
# for pos, seq in z.instances.search(Seq("GCGCAAGA", IUPAC.ambiguous_dna)):
#     print(pos, seq)
#     pos
# # motifs.write(z, format=clusterbuster)

### have to code a small find function of my own to find motifs with IUPAC.ambiguous_dna  compatability


# x = [0] * 10                    # works
# x[1] = 1
# x
#
# count = [[0, 0],[0,0]] * 10     # does not work
# count[0][0] = 1
# count
# count[1][0]
#
# y = [[0,0] for i in range(10)]  # works!
# y[0][1] = 1
# y


###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################





for seq_record in SeqIO.parse("../Data/Genome Data/CP035265.gbk", "genbank"):
    genome = seq_record.seq
    location_start = []
    location_end = []
    locus_tag = []
    gene = []

    if seq_record.features:
        print("Check: features available")
    else:
        print("Error: features not available")

    # Get all annotated sesquences and thier features
    for feature in seq_record.features:
        if feature.type == "gene":
            location_start.append(feature.location.start)
            location_end.append(feature.location.end)
            locus_tag.append(feature.qualifiers["locus_tag"])
            if "gene" in feature.qualifiers:
                gene.append(feature.qualifiers["gene"])
            else:
                gene.append(None)



    motif = "TCTAGA"
    # motif = "GCATC"
    # motif = "CRAANNNNNNNNNTTC"
    # motif = "CRAANNNNNNNNCTK"
    # motif = "CRAANNNNNNNNDTTC"



    ### FIND motif

    positions = FindMotif.posits_In_Seq(genome, motif)
    feature_i = 0
    ## Each list of 2 corresponds to a annotated sequence/gene
    ## Within each list of 2:
        # The first number represents the number of motifs in the promoter region (within (promotor_length)nt before the start of gene).
        # The second number represents the number of motifs in the gene.
    # count = [[0, 0]] * len(location_start)        # does not work, each elements are pointed to the same list of 2
    count =[[0, 0] for i in range(len(location_start))]
    gene_total = 0
    promoter_total = 0
    off_count = 0           # not in promoters or genes
    promotor_length = 150   # promotor length
    sum = len(positions)    # total number of motifs found in the genome
    for i in positions:
        while location_end[feature_i] < i:
            feature_i += 1
        if location_start[feature_i] < i:
            count[feature_i][1] += 1
            gene_total += 1
        elif location_start[feature_i] < i + promotor_length:
            count[feature_i][0] += 1
            promoter_total += 1
        else:
            off_count += 1
        # print(feature_i, count[feature_i], 1, count[1])
        # print(count)     # crashed the program, too much


    print('\n{:<24s}{:<20s}{:<12s}{:<15s}{:<15s}'.format("location", "locus_tag", "gene", "motifs (pro)", "motifs (gene)"))
    for i in range(0, len(count)):
        if count[i][0] + count[i][1] != 0:
            print('{:<12s}{:<12s}{:<20s}{:<12s}{:<15s}{:<15s}'.format(str(location_start[i]), str(location_end[i]), str(locus_tag[i]), str(gene[i]), str(count[i][0]), str(count[i][1])))

    print("total # of motifs in genes is", gene_total)
    print("total # of motifs in promoters is", promoter_total)
    print("total # of motifs not in gene/promotor is", off_count)
    print("total # of motifs in entire genome is", sum)

    max_i = 0
    for i in range(len(count)):
        if count[i][1] > count[max_i][1]:
            max_i = i
    print("The gene with the most motifs is ", locus_tag[max_i], count[max_i])
    # plt.hist([1,2,3,4],[1,2,3,4],[1,2,3,4])
    # plt.show()

    # plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
    # plt.axis([0, 6, 0, 20])
    # plt.show()
    Graph.graph_Circular_with_Features_and_Methylation(seq_record, positions, motif)

# TCTAGmA
# mAGATCT

# GCATC
# CGTAG

# CRAANNNNNNNNNTTC
# CRAANNNNNNNNCTK
# CRAANNNNNNNNDTTC
