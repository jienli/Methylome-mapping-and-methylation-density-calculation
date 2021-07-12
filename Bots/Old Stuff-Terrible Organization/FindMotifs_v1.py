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

### Seq with IUPAC.ambiguous_dna does not work, it's just a tag with no real functionality

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


###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################





for seq_record in SeqIO.parse("../Data/Genome Data/CP035265.gbk", "genbank"):
    # print(seq_record.id)
    # print(repr(seq_record.seq))
    # str = seq_record.seq
    #
    # # print(str)   # this step takes some time
    # print(len(seq_record))
    # print(seq_record.description)
    # print(seq_record.features)

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
        # print(1)
        if feature.type == "gene":
            # print(feature.qualifiers["locus_tag"])
            # print(feature.location)
            # # print(feature.qualifiers)
            # # print(feature.location.extract(seq_record))
            #
            # if "gene" in feature.qualifiers:
            #     print(feature.qualifiers["gene"])
            # if "locus_tag" in feature.qualifiers:
            #     print(feature.qualifiers["locus_tag"])
            #
            # print()


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

    # FIND motif
    genome_i = 0
    feature_i = 0
    count = [0] * len(location_start)
    off_count = 0
    len(location_start)
    while genome.find(motif, genome_i) != -1:
        genome_i = genome.find(motif, genome_i)
        while location_end[feature_i] < genome_i:
            feature_i += 1
        if location_start[feature_i] < genome_i:
            count[feature_i] += 1
        else:
            off_count += 1
        genome_i += 1

    sum = 0
    for i in count:
        sum += i
    # print(sum)

    print('\n{:<24s}{:<20s}{:<12s}{:<12s}'.format("location", "locus_tag", "gene", "# of motifs"))
    for i in range(0, len(count)):
        if count[i] != 0:
            print('{:<12s}{:<12s}{:<20s}{:<12s}{:<12s}'.format(str(location_start[i]), str(location_end[i]), str(locus_tag[i]), str(gene[i]), str(count[i])))

    print("total # of motifs in genes is", sum)
    print("total # of motifs in entire genome is", sum + off_count)

    max_i = 0
    for i in range(len(count)):
        if count[i] > count[max_i]:
            max_i = i
    print("The gene with the most motifs is ", locus_tag[max_i], count[max_i])
    # plt.hist([1,2,3,4],[1,2,3,4],[1,2,3,4])
    # plt.show()

    # plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
    # plt.axis([0, 6, 0, 20])
    # plt.show()

# TCTAGA
# GCATC
# CRAANNNNNNNNNTTC
# CRAANNNNNNNNCTK
# CRAANNNNNNNNDTTC
