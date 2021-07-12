
# import sys
# sys.path.append('/usr/local/lib/python3.4/site-packages')


import numpy as np
import xlrd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import GenBank
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
import FindMotif
# import GraphCircularGenome as Graph


from timeit import default_timer as timer
import time

### Record program time
start_time_ = time.time()
start_time = timer()

"""
Use GenBank file (.gbk)
Format: ...
"""

###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################



my_seq = Seq("ATGCATGC")
my_seq


# motifs = Motif.Motif(alphabet=IUPAC.ambiguous_dna)
# motif

## Seq with IUPAC.ambiguous_dna does not work, it's just a tag with no real functionality:

# x = Seq("AATTCTCT", IUPAC.Alphabet)
# y = Seq("TTN", IUPAC.Alphabet)
# x = Seq("AAT", generic_dna)
# y = Seq("AAN", generic_dna)
x = Seq("AAT", IUPAC.unambiguous_dna)
y = Seq("AAN", IUPAC.ambiguous_dna)
# x.find(y)
x == y


## motifs package also does not work because IUPAC.ambiguous_dna is not considered funionally

z = motifs.create([Seq("AANA", alphabet=IUPAC.ambiguous_dna)])
z.degenerate_consensus
z.counts
# z.instances.search(Seq("GCGCAAN"))
pos = 0
for pos, seq in z.instances.search(Seq("GCGCAAGA", IUPAC.ambiguous_dna)):
    print(pos, seq)
    pos
# motifs.write(z, format=clusterbuster)

## have to code a small find function of my own to find motifs with IUPAC.ambiguous_dna  compatability


x = [0] * 10                    # works
x[1] = [1,2]
x

count = [[0, 0],[0,0]] * 10     # does not work
count[0][0] = 1
count
count[1][0]

y = [[0,0] for i in range(10)]  # works!
y[0][1] = 1
y


a = [(1, 2), (3, 4), (5, 6)]
print(a[1])
print(a[1][1])
a[1] = (33,44)
print(a)
# a[1][1] = 45


# a = [0,1,1,1,0]
# b = [1,0,0,1,0]
# print(a&b)      # TypeError: unsupported operand type(s) for &: 'list' and 'list'


a = [0]
a.append(1)
print(a)

a = [(1,3),(5,2),(-1,4),(-1,2)]
a = [[1,3],[5,2],[-1,4],[-1,2]]
a.sort()
print(a)



motifs = [
        ("TCTAGA", [5], 1), 
        ("AGATCT", [5], -1), 
        ("GCATC", [2, 3], 1), 
        ("CGTAG", [2, 3], -1), 
        ("CRAANNNNNNNNNTTC", [3, 13], 1), 
        ("GYTTNNNNNNNNNAAG", [3, 13], -1), 
        ("CRAANNNNNNNNCTK", [3, 13], 1), 
        ("GYTTNNNNNNNNGAM", [3, 13], -1), 
        ("CRAANNNNNNNNDTTC", [13], 1),
        ("GYTTNNNNNNNNHAAG", [13], -1)
        ]
print(motifs[1][0])





###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################



a = [1,2,3,4]
a = a * 2
print(a)


location_end = [10,20,30]
location_start = [5,15,25]
count = [[1,2], [1,3], [0,4]]
gene_length = np.array(location_end) - np.array(location_start)
gene_length2 = [[i, i] for i in gene_length]
density = np.array(count) / gene_length2
print(density)



print("time elapsed_: ", time.time() - start_time_, "time elapsed: ", timer() - start_time)
# print("time elapsed: ", timer() - start_time)


a = "asd"


# rnaseq_wb = xlrd.open_workbook("../../Data/PG1_RNAseq.xlsx")
# rnaseq_sheet = rnaseq_wb.sheet_by_index(0)
# # rnaseq_sheet.cell_value(0,0)
# # for i in range(1, rnaseq_sheet.nrows):

# rnaseq_info = [[rnaseq_sheet.cell_value(i, 4), rnaseq_sheet.cell_value(i, 9), rnaseq_sheet.cell_value(i, 10), rnaseq_sheet.cell_value(i, 11), rnaseq_sheet.cell_value(i, 12), ] for i in range(1, rnaseq_sheet.nrows)]
# # print(rnaseq_info)
# a = rnaseq_sheet.cell_value(1,1)
# print(rnaseq_sheet.cell_value(1,1))


a = ["asa"]
b = str(a)
c = a[0]
print(a, b, c)


a = "AA"
b = "A"
print(a[0] == b)

a = np.array([0, 0, 0])
b = [0, 1, 1]
a = a/b
print(a)
a[np.isnan(a)] = 0
print(a)


# a = [123,32,51,23,71,23,21,231]
# plt.plot(a, 'o', color='black')
# plt.show()


a = "aa"
b = 23
c = a + str(b)
print(c)


a = "GAYNNNNNCTRC"
b = Seq(a, generic_dna)
print(Seq.reverse_complement(b))



print()
print()


def findComplementMotif(motif):
    motifSeq = Seq(motif[0], generic_dna)
    complementMotifSeq = str(Seq.reverse_complement(motifSeq))
    complementPosition = [len(motif[0]) - 1 - motif[1][0]]
    complementStrand = - motif[2]
    return (complementMotifSeq, complementPosition, complementStrand)


a = ("GAAHNNNNNNNNTTYG", [2], 1)
b = ("CRAANNNNNNNNDTTC", [13], -1)

print(findComplementMotif(b))


GBFiles = ["CP035234.gbk","CP035235.gbk","CP035236.gbk","CP035237.gbk","CP035238.gbk","CP035239.gbk","CP035240.gbk","CP035241.gbk","CP035242.gbk","CP035243.gbk",
"CP035244.gbk","CP035245.gbk","CP035246.gbk","CP035247.gbk","CP035248.gbk","CP035249.gbk","CP035250.gbk","CP035251.gbk","CP035252.gbk","CP035253.gbk",
"CP035254.gbk","CP035255.gbk","CP035256.gbk","CP035257.gbk","CP035258.gbk","CP035259.gbk","CP035260.gbk","CP035261.gbk","CP035262.gbk","CP035263.gbk",
"CP035264.gbk","CP035265.gbk","CP038251.gbk","CP038252.gbk","CP038253.gbk","SCPQ00000000.gbk","SCPR00000000.gbk","SCPS00000000.gbk","SCPT00000000.gbk"]
GBNames = [GBFile[:len(GBFile) - 4] for GBFile in GBFiles]

for GBFile in GBFiles:
    if GBFile[:len(GBFile) - 4] in ["CP035234"]:
        print("CP035234")