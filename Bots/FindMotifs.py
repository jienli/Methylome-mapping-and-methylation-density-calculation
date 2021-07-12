### Since v2, will add finding, counting, and plotting multiple motifs at the same time. Also include the reserse strand. Also include the specific site of methlyation

import numpy as np

import xlrd
import xlsxwriter

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import GenBank
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('tkagg')


import FindMotif
# import GraphCircularGenome as Graph
import Parameters as params

from timeit import default_timer as timer




### Record program time
start_time = timer()


def findComplementMotif(motif):
    motifSeq = Seq(motif[0], generic_dna)
    complementMotifSeq = str(Seq.reverse_complement(motifSeq))
    complementPosition = [len(motif[0]) - 1 - motif[1][0]]
    complementStrand = - motif[2]
    return (complementMotifSeq, complementPosition, complementStrand)


def completeMotifsWithComplements(motifs):
    completedMotifs = []
    for motif in motifs:
        completedMotifs.append(motif)
        completedMotifs.append(findComplementMotif(motif))
    return completedMotifs




def convert_Motif_Posit_to_Methylation_Positions(motif_positionss, motifs):
    specific_methylation_positions = [(0,0)]
    for i in range(len(motifs)):
        for motif_position in motif_positionss[i]:
            for methylation_position in motifs[i][1]:
                specific_methylation_positions.append((motif_position + methylation_position, motifs[i][2]))
    specific_methylation_positions.pop(0)
    specific_methylation_positions = remove_Duplicate(specific_methylation_positions)
    specific_methylation_positions.sort()
    return specific_methylation_positions


def remove_Duplicate(my_list):
    return list(dict.fromkeys(my_list))




# Moving Averages (of #6mA, #A and 6mA density)

def find_Moving_Average_of_6mA(genome, specific_methylation_positions, box_size):
    number_of_6mA_moving_avg = 0 

    while specific_methylation_positions[number_of_6mA_moving_avg][0] <= box_size:
        number_of_6mA_moving_avg += 1

    first_6mA_in_box = 0
    next_6mA_out_box = number_of_6mA_moving_avg

    number_of_6mA_moving_avg_list = [0] * len(genome)

    for i in range(len(genome)):
        first_bp_in_box = i
        next_bp_out_box = (i + box_size) % len(genome)

        while first_bp_in_box == specific_methylation_positions[first_6mA_in_box][0]:
            number_of_6mA_moving_avg -= 1
            first_6mA_in_box += 1
            first_6mA_in_box = first_6mA_in_box % len(specific_methylation_positions)

        while next_bp_out_box == specific_methylation_positions[next_6mA_out_box][0]:
            number_of_6mA_moving_avg += 1
            next_6mA_out_box += 1
            next_6mA_out_box = next_6mA_out_box % len(specific_methylation_positions)

        number_of_6mA_moving_avg_list[(i + int(box_size / 2)) % len(genome)] = number_of_6mA_moving_avg
    
    return number_of_6mA_moving_avg_list



def find_Moving_Average_of_A(genome, box_size):
    number_of_A_moving_avg = 0 
    i = 0
    for i in range(box_size):
        if genome[i] == "A" or genome[i] == "T":
            number_of_A_moving_avg += 1

    number_of_A_moving_avg_list = [0] * len(genome)

    for i in range(len(genome)):
        first_bp_in_box = i
        next_bp_out_box = (i + box_size) % len(genome)

        if genome[first_bp_in_box] == "A" or genome[first_bp_in_box] == "T":
            number_of_A_moving_avg -= 1

        if genome[next_bp_out_box] == "A" or genome[next_bp_out_box] == "T":
            number_of_A_moving_avg += 1

        number_of_A_moving_avg_list[(i + int(box_size / 2)) % len(genome)] = number_of_A_moving_avg

    return number_of_A_moving_avg_list


## return 0 if division by 0
def safe_Division(a, b):
    if b == 0:
        if a > 0: 
            return -1
        return 0
    return a / b












# Parameters
strainID = 10
strainName = params.strainNameForStrain[strainID]
genomeFileName = params.genomeFileNameForStrain[strainID]
rnaSeqFileName = params.rnaSeqFileNameForStrain[strainID]
motifs = completeMotifsWithComplements(params.motifsForStrain[strainID])

#Constants
promoter_length = 200   # promotor length


for seq_record in SeqIO.parse("../Data/Genome Data/GBKs/" + genomeFileName, "genbank"):

    print("now doing strain " + strainName)
    genome = seq_record.seq
    genomeSize = len(genome)
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



    ### FIND motif
    motif_positionss = [0] * len(motifs)
    specific_methylation_positions = [0]
    for i in range(len(motifs)):
        motif_positionss[i] = FindMotif.posits_In_Seq(genome, motifs[i][0])
    


    ### FIND specific methlyations sites
    specific_methylation_positions = convert_Motif_Posit_to_Methylation_Positions(motif_positionss, motifs)
 

################################################################################################

    # feature_i = 0
    # ## In this list of lists of 2, Each list of 2 corresponds to a annotated sequence/gene
    # ## Within each list of 2:
    #     # The first number represents the number of 6mA in the promoter region (within (promoter_length)nt before the start of gene).
    #     # The second number represents the number of 6mA in the gene.
    # # count = [[0, 0]] * len(location_start)        # does not work, each elements are pointed to the same list of 2
    # count =[[0, 0] for i in range(len(location_start))]
    # gene_total = 0
    # promoter_total = 0
    # off_count = 0           # not in promoters or genes
    # sum = len(specific_methylation_positions)    # total number of 6mA found in the genome
    # for i in specific_methylation_positions:
    #     while (feature_i < len(location_end)) and (location_end[feature_i] < i[0]):
    #         feature_i += 1
    #     if feature_i == len(location_end):
    #         off_count += 1
    #     elif location_start[feature_i] <= i[0]:
    #         count[feature_i][1] += 1
    #         gene_total += 1
    #     elif location_start[feature_i] < i[0] + promoter_length:
    #         count[feature_i][0] += 1
    #         promoter_total += 1
    #     else:
    #         off_count += 1



    # ### Find density of 6mA in each gene
    # # This is a density according to number of Adnines in the gene, better
    # number_of_A = [[0, 0] for i in range(len(location_start))]
    # promoter_actural_lengths = [0] * len(location_start)

    # total_A = 0
    # intragenetic_A = 0
    # promoter_A = 0
    # gene_i = 0
    # for bp_i in range(len(genome)):
    #     if gene_i < len(location_start) - 1 and location_end[gene_i] < bp_i:
    #         gene_i += 1
    #     if genome[bp_i] == "A" or genome[bp_i] == "T":
    #         total_A += 1
    #         if location_start[gene_i] <= bp_i and bp_i <= location_end[gene_i]:
    #             intragenetic_A += 1
    #             number_of_A[gene_i][1] += 1
    #         elif location_start[gene_i] - promoter_length <= bp_i and bp_i < location_start[gene_i]:
    #             promoter_A += 1
    #             number_of_A[gene_i][0] += 1
    #     if location_start[gene_i] - promoter_length <= bp_i and bp_i < location_start[gene_i]:
    #         promoter_actural_lengths[gene_i] += 1
    # density = np.array(count) / number_of_A
    # density[np.isnan(density)] = 0





    # ### Import and organize essential and expresion information from RNAseq Data
    # rnaseq_wb = xlrd.open_workbook("../Data/RNAseq Data/" + rnaSeqFileName)                                                                                                                               # Need change, file names
    # rnaseq_sheet = rnaseq_wb.sheet_by_index(0)
    # rnaseq_info = [
    #     [
    #         rnaseq_sheet.cell_value(i, 4), 
    #         rnaseq_sheet.cell_value(i, 9), 
    #         rnaseq_sheet.cell_value(i, 10), 
    #         rnaseq_sheet.cell_value(i, 11), 
    #         rnaseq_sheet.cell_value(i, 12), 
    #     ] 
    #     for i in range(0, rnaseq_sheet.nrows)
    # ]
    # rnaseq_info_organized = []
    # j = 0
    # for i in locus_tag:
    #     while j < len(rnaseq_info) - 1 and rnaseq_info[j][0] < i[0]:
    #         j += 1
    #     # rnaseq_info_organized.append(["N/A", "N/A", "N/A", "N/A", "N/A"])
    #     if rnaseq_info[j][0] == i[0]:
    #         rnaseq_info_organized.append(rnaseq_info[j])
    #     else:
    #         rnaseq_info_organized.append(["N/A", "N/A", "N/A", "N/A", "N/A"])






    # ### Additional calcualtions to get per uncertain/essential/unessential gene information

    # essential_pro_bp_number =0
    # essential_genetic_bp_number =0
    # essential_pro_A_number =0
    # essential_genetic_A_number =0
    # essential_pro_6mA_number =0
    # essential_genetic_6mA_number =0
    # essential_pro_with_6mA =0
    # essential_pro_without_6mA =0
    # essential_gene_with_6mA =0
    # essential_gene_without_6mA =0
    # number_of_essential_genes = 0

    # nonessential_pro_bp_number =0
    # nonessential_genetic_bp_number =0
    # nonessential_pro_A_number =0
    # nonessential_genetic_A_number =0
    # nonessential_pro_6mA_number =0
    # nonessential_genetic_6mA_number =0
    # nonessential_pro_with_6mA =0
    # nonessential_pro_without_6mA =0
    # nonessential_gene_with_6mA =0
    # nonessential_gene_without_6mA =0
    # number_of_nonessential_genes = 0

    # uncertain_pro_bp_number =0
    # uncertain_genetic_bp_number =0
    # uncertain_pro_A_number =0
    # uncertain_genetic_A_number =0
    # uncertain_pro_6mA_number =0
    # uncertain_genetic_6mA_number =0
    # uncertain_pro_with_6mA =0
    # uncertain_pro_without_6mA =0
    # uncertain_gene_with_6mA =0
    # uncertain_gene_without_6mA =0
    # number_of_uncertain_genes = 0

    # na_pro_bp_number =0
    # na_genetic_bp_number =0
    # na_pro_A_number =0
    # na_genetic_A_number =0
    # na_pro_6mA_number =0
    # na_genetic_6mA_number =0
    # na_pro_with_6mA =0
    # na_pro_without_6mA =0
    # na_gene_with_6mA =0
    # na_gene_without_6mA =0
    # number_of_na_genes = 0

    # for i in range(len(location_start)):
    #     if rnaseq_info_organized[i][4] == "Essential":
    #         number_of_essential_genes += 1
    #         essential_pro_bp_number += promoter_actural_lengths[i]
    #         essential_genetic_bp_number += location_end[i] - location_start[i]
    #         essential_pro_A_number += number_of_A[i][0]
    #         essential_genetic_A_number += number_of_A[i][1]
    #         essential_pro_6mA_number += count[i][0]
    #         essential_genetic_6mA_number += count[i][1]
    #         if count[i][0] > 0:
    #             essential_pro_with_6mA += 1
    #         else: 
    #             essential_pro_without_6mA += 1
    #         if count[i][1] > 0:
    #             essential_gene_with_6mA += 1
    #         else:
    #             essential_gene_without_6mA += 1

    #     elif rnaseq_info_organized[i][4] == "Non-Essential":
    #         number_of_nonessential_genes += 1
    #         nonessential_pro_bp_number += promoter_actural_lengths[i]
    #         nonessential_genetic_bp_number += location_end[i] - location_start[i]
    #         nonessential_pro_A_number += number_of_A[i][0]
    #         nonessential_genetic_A_number += number_of_A[i][1]
    #         nonessential_pro_6mA_number += count[i][0]
    #         nonessential_genetic_6mA_number += count[i][1]
    #         if count[i][0] > 0:
    #             nonessential_pro_with_6mA += 1
    #         else: 
    #             nonessential_pro_without_6mA += 1
    #         if count[i][1] > 0:
    #             nonessential_gene_with_6mA += 1
    #         else:
    #             nonessential_gene_without_6mA += 1

    #     elif rnaseq_info_organized[i][4] == "Uncertain":
    #         number_of_uncertain_genes += 1
    #         uncertain_pro_bp_number += promoter_actural_lengths[i]
    #         uncertain_genetic_bp_number += location_end[i] - location_start[i]
    #         uncertain_pro_A_number += number_of_A[i][0]
    #         uncertain_genetic_A_number += number_of_A[i][1]
    #         uncertain_pro_6mA_number += count[i][0]
    #         uncertain_genetic_6mA_number += count[i][1]
    #         if count[i][0] > 0:
    #             uncertain_pro_with_6mA += 1
    #         else: 
    #             uncertain_pro_without_6mA += 1
    #         if count[i][1] > 0:
    #             uncertain_gene_with_6mA += 1
    #         else:
    #             uncertain_gene_without_6mA += 1
                
    #     elif rnaseq_info_organized[i][4] == "N/A":
    #         number_of_na_genes += 1
    #         na_pro_bp_number += promoter_actural_lengths[i]
    #         na_genetic_bp_number += location_end[i] - location_start[i]
    #         na_pro_A_number += number_of_A[i][0]
    #         na_genetic_A_number += number_of_A[i][1]
    #         na_pro_6mA_number += count[i][0]
    #         na_genetic_6mA_number += count[i][1]
    #         if count[i][0] > 0:
    #             na_pro_with_6mA += 1
    #         else: 
    #             na_pro_without_6mA += 1
    #         if count[i][1] > 0:
    #             na_gene_with_6mA += 1
    #         else:
    #             na_gene_without_6mA += 1




    # ## Some more Macro stats
    
    # num_of_nonmethylated_promsorgenes = 0
    # num_of_essential_in_nonme_promsorgenes = 0
    # num_of_nonessential_in_nonme_promsorgenes = 0

    # num_of_methylated_genes = 0
    # num_of_methylated_proms = 0
    # num_of_methylated_promsorgenes = 0
    # num_of_essential_in_me_promsorgenes = 0
    # num_of_nonessential_in_me_promsorgenes = 0

    # for i in range(0, len(count)):
    #     if count[i][0] + count[i][1] == 0:

    #         num_of_nonmethylated_promsorgenes += 1

    #         if str(rnaseq_info_organized[i][4]) == "Essential":
    #             num_of_essential_in_nonme_promsorgenes += 1
    #         if str(rnaseq_info_organized[i][4]) == "Non-Essential":
    #             num_of_nonessential_in_nonme_promsorgenes += 1


    #     if count[i][0] + count[i][1] != 0:

    #         num_of_methylated_promsorgenes += 1

    #         if str(rnaseq_info_organized[i][4]) == "Essential":
    #             num_of_essential_in_me_promsorgenes += 1
    #         if str(rnaseq_info_organized[i][4]) == "Non-Essential":
    #             num_of_nonessential_in_me_promsorgenes += 1

    #     if count[i][0] != 0:
    #         num_of_methylated_proms += 1
    #     if count[i][1] != 0:
    #         num_of_methylated_genes += 1





    # ### Export Methylation INFO to excel
    # wb = xlsxwriter.Workbook("oragnized Methylation Data by Gene.xlsx")
    # sheet1 = wb.add_worksheet()

    # sheet1.write(0, 0, "location")
    # sheet1.write(0, 2, "locus_tag")
    # sheet1.write(0, 3, "gene")
    # sheet1.write(0, 4, "length (pro)")
    # sheet1.write(0, 5, "length (gene)")
    # sheet1.write(0, 6, "#A (pro)")
    # sheet1.write(0, 7, "#A (gene)")
    # sheet1.write(0, 8, "#6mA (pro)")
    # sheet1.write(0, 9, "#6mA (gene)")
    # sheet1.write(0, 10, "density (pro)")
    # sheet1.write(0, 11, "density (gene)")
    # sheet1.write(0, 12, "Mean")
    # sheet1.write(0, 13, "StDev")
    # sheet1.write(0, 14, "zbar")
    # sheet1.write(0, 15, "Binomial.Call")


    # for i in range(0, len(count)):
    #     sheet1.write(i + 1, 0, str(location_start[i]))
    #     sheet1.write(i + 1, 1, str(location_end[i]))
    #     sheet1.write(i + 1, 2, str(locus_tag[i]))
    #     sheet1.write(i + 1, 3, str(gene[i]))
    #     sheet1.write(i + 1, 4, str(promoter_actural_lengths[i]))
    #     sheet1.write(i + 1, 5, str(location_end[i] - location_start[i]))
    #     sheet1.write(i + 1, 6, str(number_of_A[i][0]))
    #     sheet1.write(i + 1, 7, str(number_of_A[i][1]))
    #     sheet1.write(i + 1, 8, str(count[i][0]))
    #     sheet1.write(i + 1, 9, str(count[i][1]))
    #     sheet1.write(i + 1, 10, str(round(density[i][0], 10)))
    #     sheet1.write(i + 1, 11, str(round(density[i][1], 10)))
    #     sheet1.write(i + 1, 12, str(rnaseq_info_organized[i][1]))
    #     sheet1.write(i + 1, 13, str(rnaseq_info_organized[i][2]))
    #     sheet1.write(i + 1, 14, str(rnaseq_info_organized[i][3]))
    #     sheet1.write(i + 1, 15, str(rnaseq_info_organized[i][4]))

    # wb.close()



    # ### PRINT MACRO INFO into a .txt file

    # macroFile = open(r"macroInformation.txt", "w")

    # macroFile.write(strainName)
    # macroFile.write("\n\n")
    # macroFile.write("Total number of methylated genes(include promoter)      " + str(num_of_methylated_promsorgenes) + "\n")
    # macroFile.write("Total number of non-methylated genes(include promoter)  " + str(num_of_nonmethylated_promsorgenes) + "\n")

    # macroFile.write("\n")
    # macroFile.write("total # of m6A in genes is               " + str(gene_total) + "\n")
    # macroFile.write("total # of m6A in promoters is           " + str(promoter_total) + "\n")
    # macroFile.write("total # of m6A not in gene/promotor is   " + str(off_count) + "\n")
    # macroFile.write("total # of m6A in entire genome is       " + str(sum) + "\n")
   
    # macroFile.write("\n")
    # macroFile.write("# of genes+proms without 6mA:    " + str(num_of_nonmethylated_promsorgenes) + "\n")
    # macroFile.write("# of genes+proms with 6mA:       " + str(num_of_methylated_promsorgenes) + "\n")
    # macroFile.write("# of genes with 6mA:             " + str(num_of_methylated_genes) + "\n")
    # macroFile.write("# of promoters with 6mA:         " + str(num_of_methylated_proms) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Average 6mA density (per bp) is:         " + str(safe_Division(sum, len(genome) * 2)) + "\n")
    # macroFile.write("Average 6mA density (per Adenine) is:    " + str(safe_Division(sum, total_A)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among genes(and pros) that have no methylation: " + "\n")
    # macroFile.write("Essentials:     " + str(num_of_essential_in_nonme_promsorgenes) + " " + str(safe_Division(float(num_of_essential_in_nonme_promsorgenes), num_of_nonmethylated_promsorgenes)) + "\n")
    # macroFile.write("Non-Essentials: " + str(num_of_nonessential_in_nonme_promsorgenes) + " " + str(safe_Division(float(num_of_nonessential_in_nonme_promsorgenes), num_of_nonmethylated_promsorgenes)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among genes(and pros) that have methylation: " + "\n")
    # macroFile.write("Essentials:     " + str(num_of_essential_in_me_promsorgenes) + " " + str(safe_Division(float(num_of_essential_in_me_promsorgenes), num_of_methylated_promsorgenes)) + "\n")
    # macroFile.write("Non-Essentials: " + str(num_of_nonessential_in_me_promsorgenes) + " " + str(safe_Division(float(num_of_nonessential_in_me_promsorgenes), num_of_methylated_promsorgenes)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among essential genes: " + "\n")
    # macroFile.write("The average A   density (per bp)      (pro)  is: " + str(safe_Division(essential_pro_A_number, essential_pro_bp_number * 2)) + "\n")
    # macroFile.write("The average A   density (per bp)      (gene) is: " + str(safe_Division(essential_genetic_A_number, essential_genetic_bp_number * 2)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (pro)  is: " + str(safe_Division(essential_pro_6mA_number, essential_pro_A_number)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (gene) is: " + str(safe_Division(essential_genetic_6mA_number, essential_genetic_A_number)) + "\n")
    # macroFile.write("# promoters that have    6mA: " + str(essential_pro_with_6mA) + " " + str(safe_Division(essential_pro_with_6mA, number_of_essential_genes)) + "\n")
    # macroFile.write("# promoters that have no 6mA: " + str(essential_pro_without_6mA) + " " + str(safe_Division(essential_pro_without_6mA, number_of_essential_genes)) + "\n")
    # macroFile.write("# genes     that have    6mA: " + str(essential_gene_with_6mA) + " " + str(safe_Division(essential_gene_with_6mA, number_of_essential_genes)) + "\n")
    # macroFile.write("# genes     that have no 6mA: " + str(essential_gene_without_6mA) + " " + str(safe_Division(essential_gene_without_6mA, number_of_essential_genes)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among non-essential genes: " + "\n")
    # macroFile.write("The average A   density (per bp)      (pro)  is: " + str(safe_Division(nonessential_pro_A_number, nonessential_pro_bp_number * 2)) + "\n")
    # macroFile.write("The average A   density (per bp)      (gene) is: " + str(safe_Division(nonessential_genetic_A_number, nonessential_genetic_bp_number * 2)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (pro)  is: " + str(safe_Division(nonessential_pro_6mA_number, nonessential_pro_A_number)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (gene) is: " + str(safe_Division(nonessential_genetic_6mA_number, nonessential_genetic_A_number)) + "\n")
    # macroFile.write("# promoters that have    6mA: " + str(nonessential_pro_with_6mA) + " " + str(safe_Division(nonessential_pro_with_6mA, number_of_nonessential_genes)) + "\n")
    # macroFile.write("# promoters that have no 6mA: " + str(nonessential_pro_without_6mA) + " " + str(safe_Division(nonessential_pro_without_6mA, number_of_nonessential_genes)) + "\n")
    # macroFile.write("# genes     that have    6mA: " + str(nonessential_gene_with_6mA) + " " + str(safe_Division(nonessential_gene_with_6mA, number_of_nonessential_genes)) + "\n")
    # macroFile.write("# genes     that have no 6mA: " + str(nonessential_gene_without_6mA) + " " + str(safe_Division(nonessential_gene_without_6mA, number_of_nonessential_genes)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among uncertain genes: " + "\n")
    # macroFile.write("The average A   density (per bp)      (pro)  is: " + str(safe_Division(uncertain_pro_A_number, uncertain_pro_bp_number * 2)) + "\n")
    # macroFile.write("The average A   density (per bp)      (gene) is: " + str(safe_Division(uncertain_genetic_A_number, uncertain_genetic_bp_number * 2)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (pro)  is: " + str(safe_Division(uncertain_pro_6mA_number, uncertain_pro_A_number)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (gene) is: " + str(safe_Division(uncertain_genetic_6mA_number, uncertain_genetic_A_number)) + "\n")
    # macroFile.write("# promoters that have    6mA: " + str(uncertain_pro_with_6mA) + " " + str(safe_Division(uncertain_pro_with_6mA, number_of_uncertain_genes)) + "\n")
    # macroFile.write("# promoters that have no 6mA: " + str(uncertain_pro_without_6mA) + " " + str(safe_Division(uncertain_pro_without_6mA, number_of_uncertain_genes)) + "\n")
    # macroFile.write("# genes     that have    6mA: " + str(uncertain_gene_with_6mA) + " " + str(safe_Division(uncertain_gene_with_6mA, number_of_uncertain_genes)) + "\n")
    # macroFile.write("# genes     that have no 6mA: " + str(uncertain_gene_without_6mA) + " " + str(safe_Division(uncertain_gene_without_6mA, number_of_uncertain_genes)) + "\n")

    # macroFile.write("\n")
    # macroFile.write("Among N/A genes: " + "\n")
    # macroFile.write("The average A   density (per bp)      (pro)  is: " + str(safe_Division(na_pro_A_number, na_pro_bp_number * 2)) + "\n")
    # macroFile.write("The average A   density (per bp)      (gene) is: " + str(safe_Division(na_genetic_A_number, na_genetic_bp_number * 2)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (pro)  is: " + str(safe_Division(na_pro_6mA_number, na_pro_A_number)) + "\n")
    # macroFile.write("The average 6mA density (per Adenine) (gene) is: " + str(safe_Division(na_genetic_6mA_number, na_genetic_A_number)) + "\n")
    # macroFile.write("# promoters that have    6mA: " + str(na_pro_with_6mA) + " " + str(safe_Division(na_pro_with_6mA, number_of_na_genes)) + "\n")
    # macroFile.write("# promoters that have no 6mA: " + str(na_pro_without_6mA) + " " + str(safe_Division(na_pro_without_6mA, number_of_na_genes)) + "\n")
    # macroFile.write("# genes     that have    6mA: " + str(na_gene_with_6mA) + " " + str(safe_Division(na_gene_with_6mA, number_of_na_genes)) + "\n")
    # macroFile.write("# genes     that have no 6mA: " + str(na_gene_without_6mA) + " " + str(safe_Division(na_gene_without_6mA, number_of_na_genes)) + "\n")



    # macroFile.write("\n")
    # max_i = 0
    # for i in range(len(count)):
    #     if count[i][1] > count[max_i][1]:
    #         max_i = i
    # macroFile.write("The gene with the most m6A is " + str(locus_tag[max_i]) + str(count[max_i]) + "\n")
    
    
    # macroFile.close()
##################################################################################################





    ### python graphing with Bio

    # Graph.graph_Circular_with_Features_and_Methylation(seq_record, specific_methylation_positions, motifs)






    ### into Circos

    # ## moving #6mA, #A, 6mA density averages

    # box_size = 1000000
    # number_of_6mA_moving_avg_list = find_Moving_Average_of_6mA(genome, specific_methylation_positions, box_size)
    # # plt.plot(number_of_6mA_moving_avg_list, 'o', color='black')
    # # plt.show()
    # number_of_A_moving_avg_list = find_Moving_Average_of_A(genome, box_size)
    # # plt.plot(number_of_A_moving_avg_list, 'o', color='red')
    # # plt.show()
    # density_of_6mA_per_A_moving_avg_list = np.array(number_of_6mA_moving_avg_list) / np.array(number_of_A_moving_avg_list)
    # # plt.plot(density_of_6mA_per_A_moving_avg_list, 'o', color='blue')
    # # plt.show()

    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/movingAve_1000000.txt", "w")
    # for i in range(len(density_of_6mA_per_A_moving_avg_list)):
    #     if i % 1000 == 0 and i < genomeSize:
    #         circosFile.write("hs1" + " " + str(i) + " " + str(i + 1000) + " " + str(density_of_6mA_per_A_moving_avg_list[i]) + "\n")
    # circosFile.close()



    box_size = 100000
    number_of_6mA_moving_avg_list = find_Moving_Average_of_6mA(genome, specific_methylation_positions, box_size)

    number_of_A_moving_avg_list = find_Moving_Average_of_A(genome, box_size)

    density_of_6mA_per_A_moving_avg_list = np.array(number_of_6mA_moving_avg_list) / np.array(number_of_A_moving_avg_list)

    plt.plot(number_of_6mA_moving_avg_list, 'o', color='black')
    # plt.show()
    plt.savefig('/Users/jienli/Documents/BC/Tim_van_Opijnen/7_2021_Spring/new methylation Figures/' + strainName + '_6mA.png')
    plt.clf()
    plt.plot(number_of_A_moving_avg_list, 'o', color='red')
    # plt.show()
    plt.savefig('/Users/jienli/Documents/BC/Tim_van_Opijnen/7_2021_Spring/new methylation Figures/' + strainName + '_A.png')
    plt.clf()
    plt.plot(density_of_6mA_per_A_moving_avg_list, 'o', color='blue')
    # plt.show()
    plt.savefig('/Users/jienli/Documents/BC/Tim_van_Opijnen/7_2021_Spring/new methylation Figures/' + strainName + '_6mA_per_A.png')
    plt.clf()

    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/movingAve_100000.txt", "w")
    # for i in range(len(density_of_6mA_per_A_moving_avg_list)):
    #     if i % 1000 == 0 and i < genomeSize:
    #         circosFile.write("hs1" + " " + str(i) + " " + str(i + 1000) + " " + str(density_of_6mA_per_A_moving_avg_list[i]) + "\n")
    # circosFile.close()



    # box_size = 10000
    # number_of_6mA_moving_avg_list = find_Moving_Average_of_6mA(genome, specific_methylation_positions, box_size)
    # number_of_A_moving_avg_list = find_Moving_Average_of_A(genome, box_size)
    # density_of_6mA_per_A_moving_avg_list = np.array(number_of_6mA_moving_avg_list) / np.array(number_of_A_moving_avg_list)

    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/movingAve_10000.txt", "w")
    # for i in range(len(density_of_6mA_per_A_moving_avg_list)):
    #     if i % 1000 == 0 and i < genomeSize:
    #         circosFile.write("hs1" + " " + str(i) + " " + str(i + 1000) + " " + str(density_of_6mA_per_A_moving_avg_list[i]) + "\n")
    # circosFile.close()




    # ## genes
    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/karyotype/genes.txt", "w")
    # circosFile.write("chr - hs1 " + strainName + " 0 " + str(genomeSize) + " chr1\n")
    # for i in range(0, len(count)):
    #     circosFile.write("band hs1 " + str(locus_tag[i]) + " " + str(locus_tag[i]) + " " + str(location_start[i]) + " " + str(location_end[i]))
    #     if i % 2 == 0:
    #         circosFile.write(" blue")
    #     else:
    #         circosFile.write(" lblue")
    #     circosFile.write("\n") 
    # circosFile.close()



    # ## highlights (6mAs)
    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/highlights/methylationPlus.txt", "w")
    # for i in specific_methylation_positions:
    #     if i[1] == 1:
    #         circosFile.write("hs1 " + str(i[0]) + " " + str(i[0] + 50) + "\n")
    # circosFile.close()


    # circosFile = open(r"../Circos/Circos 2020 Summer Methylation graphs/data/highlights/methylationMinus.txt", "w")
    # for i in specific_methylation_positions:
    #     if i[1] == -1:
    #         circosFile.write("hs1 " + str(i[0]) + " " + str(i[0] + 50) + "\n")
    # circosFile.close()





print()
print("time elapsed: ", timer() - start_time)

