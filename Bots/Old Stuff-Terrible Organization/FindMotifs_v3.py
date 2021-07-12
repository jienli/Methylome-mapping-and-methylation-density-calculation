### Since v2, will add finding, counting, and plotting multiple motifs at the same time. Also include the reserse strand. Also include the specific site of methlyation

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
import GraphCircularGenome as Graph

from timeit import default_timer as timer
# import time

### Record program time
# start_time_ = time.time()
start_time = timer()

"""
Use GenBank file (.gbk)
Format: ...
"""


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


for seq_record in SeqIO.parse("../Data/Genome Data/CP035265.gbk", "genbank"):
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



    ### Define the motifs of interest as well as methlayion locations

    ## list of touples of 3. each touple include the sequence of the motif, the locations of methylation, and the strand (+1 or -1)
    # motif = "TCTAGA"
    # motif = "GCATC"
    # motif = "CRAANNNNNNNNNTTC"
    # motif = "CRAANNNNNNNNCTK"
    # motif = "CRAANNNNNNNNDTTC"


    # motifs = [
    #     ("TCTAGA", [5]), 
    #     ("GCATC", [2, 3]), 
    #     ("CRAANNNNNNNNNTTC", [3, 13]), 
    #     ("CRAANNNNNNNNCTK", [3, 13]), 
    #     ("CRAANNNNNNNNDTTC", [13]),
    #     ("GYTTNNNNNNNNHAAG", [13])
    #     ]
    

    ## Based on the information in the methylation organization data I created, there are multiple methylations on one motif (some not on A? maybe counted the other strand too?)
    ## OK this is wrong, use the one below this one
    # motifs = [
    #     ("TCTAGA", [5], 1), 
    #     ("AGATCT", [5], -1), 
    #     ("GCATC", [2, 3], 1), 
    #     ("CGTAG", [2, 3], -1), 
    #     ("CRAANNNNNNNNNTTC", [3, 13], 1), 
    #     ("GYTTNNNNNNNNNAAG", [3, 13], -1), 
    #     ("CRAANNNNNNNNCTK", [3, 13], 1), 
    #     ("GYTTNNNNNNNNGAM", [3, 13], -1), 
    #     ("CRAANNNNNNNNDTTC", [13], 1),
    #     ("GYTTNNNNNNNNHAAG", [13], -1)
    #     ]

    ## Based on 2018.11.29_SC FRosconi_Spneumoniae.xlsx. Only one methylation per motif
    motifs = [
        ("GCATC", [2], 1), 
        ("GATGC", [1], 1), 
        ("GCATC", [3], -1), 
        ("GATGC", [2], -1), 

        ("TCTAGA", [5], 1), 
        ("TCTAGA", [0], -1), 

        ("CRAANNNNNNNNNTTC", [3], 1), 
        ("GAANNNNNNNNNTTYG", [12], -1), 

        ("CRAANNNNNNNNCTK", [3], 1), 
        ("MAGNNNNNNNNTTYG", [1], 1), 
        ("CRAANNNNNNNNCTK", [13], -1), 
        ("MAGNNNNNNNNTTYG", [11], 1), 

        ("GAAHNNNNNNNNTTYG", [2], 1),
        ("CRAANNNNNNNNDTTC", [13], -1)
        ]


    ### FIND motif
    motif_positionss = [0] * len(motifs)
    specific_methylation_positions = [0]
    for i in range(len(motifs)):
        motif_positionss[i] = FindMotif.posits_In_Seq(genome, motifs[i][0])
    


    ### FIND specific methlyations sites
    specific_methylation_positions = convert_Motif_Posit_to_Methylation_Positions(motif_positionss, motifs)
    # positions = motif_positionss[0]

    # Constants
    promoter_length = 200   # promotor length

    feature_i = 0
    ## In this list of lists of 2, Each list of 2 corresponds to a annotated sequence/gene
    ## Within each list of 2:
        # The first number represents the number of motifs in the promoter region (within (promoter_length)nt before the start of gene).
        # The second number represents the number of motifs in the gene.
    # count = [[0, 0]] * len(location_start)        # does not work, each elements are pointed to the same list of 2
    count =[[0, 0] for i in range(len(location_start))]
    gene_total = 0
    promoter_total = 0
    off_count = 0           # not in promoters or genes
    sum = len(specific_methylation_positions)    # total number of motifs found in the genome
    for i in specific_methylation_positions:
        while (feature_i < len(location_end)) and (location_end[feature_i] < i[0]):
            feature_i += 1
        if feature_i == len(location_end):
            off_count += 1
        elif location_start[feature_i] < i[0]:
            count[feature_i][1] += 1
            gene_total += 1
        elif location_start[feature_i] < i[0] + promoter_length:
            count[feature_i][0] += 1
            promoter_total += 1
        else:
            off_count += 1
        # print(feature_i, count[feature_i], 1, count[1])
        # print(count)     # crashed the program, too much



    ### Find density and more

    # This is density according to gene length, not a good measure
    # gene_length_ = np.array(location_end) - np.array(location_start)
    # gene_length = [[i, i] for i in gene_length_]
    # density = np.array(count) / gene_length

    # This is a density according to number of Adnines in the gene, better
    number_of_A = [[0, 0] for i in range(len(location_start))]
    promoter_actural_lengths = [0] * len(location_start)

    total_A = 0
    intragenetic_A = 0
    promoter_A = 0
    gene_i = 0
    for bp_i in range(len(genome)):
        if gene_i < len(location_start) - 1 and location_end[gene_i] < bp_i:
            gene_i += 1
        if genome[bp_i] == "A":
            total_A += 1
            if location_start[gene_i] <= bp_i and bp_i <= location_end[gene_i]:
                intragenetic_A += 1
                number_of_A[gene_i][1] += 1
            elif location_start[gene_i] - promoter_length <= bp_i and bp_i < location_start[gene_i]:
                promoter_A += 1
                number_of_A[gene_i][0] += 1
        if location_start[gene_i] - promoter_length <= bp_i and bp_i < location_start[gene_i]:
            promoter_actural_lengths[gene_i] += 1
    density = np.array(count) / number_of_A
    density[np.isnan(density)] = 0





    ### Import and organize essential and expresion information from RNAseq Data
    rnaseq_wb = xlrd.open_workbook("../Data/PG1_RNAseq.xlsx")
    rnaseq_sheet = rnaseq_wb.sheet_by_index(0)
    rnaseq_info = [
        [
            rnaseq_sheet.cell_value(i, 4), 
            rnaseq_sheet.cell_value(i, 9), 
            rnaseq_sheet.cell_value(i, 10), 
            rnaseq_sheet.cell_value(i, 11), 
            rnaseq_sheet.cell_value(i, 12), 
        ] 
        for i in range(1, rnaseq_sheet.nrows)
    ]
    rnaseq_info_organized = []
    j = 0
    for i in locus_tag:
        while j < len(rnaseq_info) - 1 and rnaseq_info[j][0] < i[0]:
            j += 1
        # print(i[0], "    ", rnaseq_info[j][0] == i, "    ",rnaseq_info[j][0])
        if rnaseq_info[j][0] == i[0]:
            rnaseq_info_organized.append(rnaseq_info[j])
        else:
            rnaseq_info_organized.append(["N/A", "N/A", "N/A", "N/A", "N/A"])

    print(rnaseq_info_organized)



    ### Additional calcualtions to get per uncertain/essential/unessential gene information


    # print("The average A   density (per bp)      (pro)  is: ", essential_pro_A_number / essential_pro_bp_number)
    # print("The average A   density (per bp)      (gene) is: ", essential_genetic_A_number / essential_genetic_bp_number)
    # print("The average 6mA density (per Adenine) (pro)  is: ", essential_pro_6mA_number / essential_pro_A_number)
    # print("The average 6mA density (per Adenine) (gene) is: ", essential_genetic_6mA_number / essential_genetic_A_number)
    # print("# promoters that have    6mA: ", essential_pro_with_6mA, essential_pro_with_6mA / number_of_essential_genes)
    # print("# promoters that have no 6mA: ", essential_pro_without_6mA, essential_pro_without_6mA / number_of_essential_genes)
    # print("# genes     that have    6mA: ", essential_gene_with_6mA, essential_gene_with_6mA / number_of_essential_genes)
    # print("# genes     that have no 6mA: ", essential_gene_without_6mA, essential_gene_without_6mA / number_of_essential_genes)

    essential_pro_bp_number =0
    essential_genetic_bp_number =0
    essential_pro_A_number =0
    essential_genetic_A_number =0
    essential_pro_6mA_number =0
    essential_genetic_6mA_number =0
    essential_pro_with_6mA =0
    essential_pro_without_6mA =0
    essential_gene_with_6mA =0
    essential_gene_without_6mA =0
    number_of_essential_genes = 0

    nonessential_pro_bp_number =0
    nonessential_genetic_bp_number =0
    nonessential_pro_A_number =0
    nonessential_genetic_A_number =0
    nonessential_pro_6mA_number =0
    nonessential_genetic_6mA_number =0
    nonessential_pro_with_6mA =0
    nonessential_pro_without_6mA =0
    nonessential_gene_with_6mA =0
    nonessential_gene_without_6mA =0
    number_of_nonessential_genes = 0

    uncertain_pro_bp_number =0
    uncertain_genetic_bp_number =0
    uncertain_pro_A_number =0
    uncertain_genetic_A_number =0
    uncertain_pro_6mA_number =0
    uncertain_genetic_6mA_number =0
    uncertain_pro_with_6mA =0
    uncertain_pro_without_6mA =0
    uncertain_gene_with_6mA =0
    uncertain_gene_without_6mA =0
    number_of_uncertain_genes = 0

    na_pro_bp_number =0
    na_genetic_bp_number =0
    na_pro_A_number =0
    na_genetic_A_number =0
    na_pro_6mA_number =0
    na_genetic_6mA_number =0
    na_pro_with_6mA =0
    na_pro_without_6mA =0
    na_gene_with_6mA =0
    na_gene_without_6mA =0
    number_of_na_genes = 0

    for i in range(len(location_start)):
        if rnaseq_info_organized[i][4] == "Essential":
            number_of_essential_genes += 1
            essential_pro_bp_number += promoter_actural_lengths[i]
            essential_genetic_bp_number += location_end[i] - location_start[i]
            essential_pro_A_number += number_of_A[i][0]
            essential_genetic_A_number += number_of_A[i][1]
            essential_pro_6mA_number += count[i][0]
            essential_genetic_6mA_number += count[i][1]
            if count[i][0] > 0:
                essential_pro_with_6mA += 1
            else: 
                essential_pro_without_6mA += 1
            if count[i][1] > 0:
                essential_gene_with_6mA += 1
            else:
                essential_gene_without_6mA += 1

        elif rnaseq_info_organized[i][4] == "Non-Essential":
            number_of_nonessential_genes += 1
            nonessential_pro_bp_number += promoter_actural_lengths[i]
            nonessential_genetic_bp_number += location_end[i] - location_start[i]
            nonessential_pro_A_number += number_of_A[i][0]
            nonessential_genetic_A_number += number_of_A[i][1]
            nonessential_pro_6mA_number += count[i][0]
            nonessential_genetic_6mA_number += count[i][1]
            if count[i][0] > 0:
                nonessential_pro_with_6mA += 1
            else: 
                nonessential_pro_without_6mA += 1
            if count[i][1] > 0:
                nonessential_gene_with_6mA += 1
            else:
                nonessential_gene_without_6mA += 1

        elif rnaseq_info_organized[i][4] == "Uncertain":
            number_of_uncertain_genes += 1
            uncertain_pro_bp_number += promoter_actural_lengths[i]
            uncertain_genetic_bp_number += location_end[i] - location_start[i]
            uncertain_pro_A_number += number_of_A[i][0]
            uncertain_genetic_A_number += number_of_A[i][1]
            uncertain_pro_6mA_number += count[i][0]
            uncertain_genetic_6mA_number += count[i][1]
            if count[i][0] > 0:
                uncertain_pro_with_6mA += 1
            else: 
                uncertain_pro_without_6mA += 1
            if count[i][1] > 0:
                uncertain_gene_with_6mA += 1
            else:
                uncertain_gene_without_6mA += 1
                
        elif rnaseq_info_organized[i][4] == "N/A":
            number_of_na_genes += 1
            na_pro_bp_number += promoter_actural_lengths[i]
            na_genetic_bp_number += location_end[i] - location_start[i]
            na_pro_A_number += number_of_A[i][0]
            na_genetic_A_number += number_of_A[i][1]
            na_pro_6mA_number += count[i][0]
            na_genetic_6mA_number += count[i][1]
            if count[i][0] > 0:
                na_pro_with_6mA += 1
            else: 
                na_pro_without_6mA += 1
            if count[i][1] > 0:
                na_gene_with_6mA += 1
            else:
                na_gene_without_6mA += 1






    ### PRINT genes with no Methylation

    # print("These are the genes with no methylation in promotor or coding region")
    # print('\n{:<24s}{:<20s}{:<12s}'.format("location", "locus_tag", "gene"))
    # num_of_nonmethylated_promsorgenes = 0
    # for i in range(0, len(count)):
    #     if count[i][0] + count[i][1] == 0:
    #         print('{:<12s}{:<12s}{:<20s}{:<12s}'.format(str(location_start[i]), str(location_end[i]), str(locus_tag[i]), str(gene[i])))
    #         num_of_nonmethylated_promsorgenes += 1
    # print("For a total of ", num_of_nonmethylated_promsorgenes)

    print("These are the genes with no methylation in promotor or coding region")
    print(
        '\n{:<24s}{:<20s}{:<12s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}'.
        format(
            "location", 
            "locus_tag", 
            "gene", 
            "length (pro)", 
            "length (gene)", 
            "#A (pro)", 
            "#A (gene)", 
            "#6mA (pro)", 
            "#6mA (gene)", 
            "density (pro)", 
            "density (gene)", 
            "Mean", 
            "StDev", 
            "zbar", 
            "Binomial.Call"
        )
    )
    num_of_nonmethylated_promsorgenes = 0
    num_of_essential_in_nonme_promsorgenes = 0
    num_of_nonessential_in_nonme_promsorgenes = 0
    for i in range(0, len(count)):
        if count[i][0] + count[i][1] == 0:
            print(
                '{:<12s}{:<12s}{:<20s}{:<12s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}'.
                format(
                    str(location_start[i]), 
                    str(location_end[i]), 
                    str(locus_tag[i]), 
                    str(gene[i]), 
                    str(promoter_actural_lengths[i]),
                    str(location_end[i] - location_start[i]),
                    str(number_of_A[i][0]),
                    str(number_of_A[i][1]),
                    str(count[i][0]), 
                    str(count[i][1]), 
                    str(round(density[i][0], 10)), 
                    str(round(density[i][1], 10)), 
                    str(rnaseq_info_organized[i][1]), 
                    str(rnaseq_info_organized[i][2]), 
                    str(rnaseq_info_organized[i][3]), 
                    str(rnaseq_info_organized[i][4])
                )
            )

            num_of_nonmethylated_promsorgenes += 1

            if str(rnaseq_info_organized[i][4]) == "Essential":
                num_of_essential_in_nonme_promsorgenes += 1
            if str(rnaseq_info_organized[i][4]) == "Non-Essential":
                num_of_nonessential_in_nonme_promsorgenes += 1

    print("For a total of ", num_of_nonmethylated_promsorgenes)




    ### PRINT Methylation INFO
    print()
    print("Now these are genes with methylaytion in either promoter or coding region")
    print(
        '\n{:<24s}{:<20s}{:<12s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}'.
        format(
            "location", 
            "locus_tag", 
            "gene", 
            "length (pro)", 
            "length (gene)", 
            "#A (pro)", 
            "#A (gene)", 
            "#6mA (pro)", 
            "#6mA (gene)", 
            "density (pro)", 
            "density (gene)", 
            "Mean", 
            "StDev", 
            "zbar", 
            "Binomial.Call"
        )
    )
    num_of_methylated_genes = 0
    num_of_methylated_proms = 0
    num_of_methylated_promsorgenes = 0
    num_of_essential_in_me_promsorgenes = 0
    num_of_nonessential_in_me_promsorgenes = 0
    for i in range(0, len(count)):
        if count[i][0] + count[i][1] != 0:
            print(
                '{:<12s}{:<12s}{:<20s}{:<12s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}{:<15s}'.
                format(
                    str(location_start[i]), 
                    str(location_end[i]), 
                    str(locus_tag[i]), 
                    str(gene[i]), 
                    str(promoter_actural_lengths[i]),
                    str(location_end[i] - location_start[i]),
                    str(number_of_A[i][0]),
                    str(number_of_A[i][1]),
                    str(count[i][0]), 
                    str(count[i][1]), 
                    str(round(density[i][0], 10)), 
                    str(round(density[i][1], 10)), 
                    str(rnaseq_info_organized[i][1]), 
                    str(rnaseq_info_organized[i][2]), 
                    str(rnaseq_info_organized[i][3]), 
                    str(rnaseq_info_organized[i][4])
                )
            )

            num_of_methylated_promsorgenes += 1

            if str(rnaseq_info_organized[i][4]) == "Essential":
                num_of_essential_in_me_promsorgenes += 1
            if str(rnaseq_info_organized[i][4]) == "Non-Essential":
                num_of_nonessential_in_me_promsorgenes += 1

        if count[i][0] != 0:
            num_of_methylated_proms += 1
        if count[i][1] != 0:
            num_of_methylated_genes += 1
        

    ### PRINT MACRO INFO
    print()
    print("total # of m6A in genes is", gene_total)
    print("total # of m6A in promoters is", promoter_total)
    print("total # of m6A not in gene/promotor is", off_count)
    print("total # of m6A in entire genome is", sum)
    
    print()
    print("# of genes+proms without 6mA: ", num_of_nonmethylated_promsorgenes)
    print("# of genes+proms with 6mA: ", num_of_methylated_promsorgenes)
    print("# of genes with 6mA: ", num_of_methylated_genes)
    print("# of promoters with 6mA: ", num_of_methylated_proms)

    print()
    print("Average 6mA density (per bp) is: ", sum / len(genome))
    print("Average 6mA density (per Adenine) is: ", sum / total_A)


    print()
    print("Among genes(and pros) that have no methylation: ")
    print("Essentials:     ", num_of_essential_in_nonme_promsorgenes, float(num_of_essential_in_nonme_promsorgenes) / num_of_nonmethylated_promsorgenes)
    print("Non-Essentials: ", num_of_nonessential_in_nonme_promsorgenes, float(num_of_nonessential_in_nonme_promsorgenes) / num_of_nonmethylated_promsorgenes)

    print()
    print("Among genes(and pros) that have methylation: ")
    print("Essentials:     ", num_of_essential_in_me_promsorgenes, float(num_of_essential_in_me_promsorgenes) / num_of_methylated_promsorgenes)
    print("Non-Essentials: ", num_of_nonessential_in_me_promsorgenes, float(num_of_nonessential_in_me_promsorgenes) / num_of_methylated_promsorgenes)


    print()
    print("Among essential genes: ")
    print("The average A   density (per bp)      (pro)  is: ", essential_pro_A_number / essential_pro_bp_number)
    print("The average A   density (per bp)      (gene) is: ", essential_genetic_A_number / essential_genetic_bp_number)
    print("The average 6mA density (per Adenine) (pro)  is: ", essential_pro_6mA_number / essential_pro_A_number)
    print("The average 6mA density (per Adenine) (gene) is: ", essential_genetic_6mA_number / essential_genetic_A_number)
    print("# promoters that have    6mA: ", essential_pro_with_6mA, essential_pro_with_6mA / number_of_essential_genes)
    print("# promoters that have no 6mA: ", essential_pro_without_6mA, essential_pro_without_6mA / number_of_essential_genes)
    print("# genes     that have    6mA: ", essential_gene_with_6mA, essential_gene_with_6mA / number_of_essential_genes)
    print("# genes     that have no 6mA: ", essential_gene_without_6mA, essential_gene_without_6mA / number_of_essential_genes)

    print()
    print("Among non-essential genes: ")
    print("The average A   density (per bp)      (pro)  is: ", nonessential_pro_A_number / nonessential_pro_bp_number)
    print("The average A   density (per bp)      (gene) is: ", nonessential_genetic_A_number / nonessential_genetic_bp_number)
    print("The average 6mA density (per Adenine) (pro)  is: ", nonessential_pro_6mA_number / nonessential_pro_A_number)
    print("The average 6mA density (per Adenine) (gene) is: ", nonessential_genetic_6mA_number / nonessential_genetic_A_number)
    print("# promoters that have    6mA: ", nonessential_pro_with_6mA, nonessential_pro_with_6mA / number_of_nonessential_genes)
    print("# promoters that have no 6mA: ", nonessential_pro_without_6mA, nonessential_pro_without_6mA / number_of_nonessential_genes)
    print("# genes     that have    6mA: ", nonessential_gene_with_6mA, nonessential_gene_with_6mA / number_of_nonessential_genes)
    print("# genes     that have no 6mA: ", nonessential_gene_without_6mA, nonessential_gene_without_6mA / number_of_nonessential_genes)

    print()
    print("Among uncertain genes: ")
    print("The average A   density (per bp)      (pro)  is: ", uncertain_pro_A_number / uncertain_pro_bp_number)
    print("The average A   density (per bp)      (gene) is: ", uncertain_genetic_A_number / uncertain_genetic_bp_number)
    print("The average 6mA density (per Adenine) (pro)  is: ", uncertain_pro_6mA_number / uncertain_pro_A_number)
    print("The average 6mA density (per Adenine) (gene) is: ", uncertain_genetic_6mA_number / uncertain_genetic_A_number)
    print("# promoters that have    6mA: ", uncertain_pro_with_6mA, uncertain_pro_with_6mA / number_of_uncertain_genes)
    print("# promoters that have no 6mA: ", uncertain_pro_without_6mA, uncertain_pro_without_6mA / number_of_uncertain_genes)
    print("# genes     that have    6mA: ", uncertain_gene_with_6mA, uncertain_gene_with_6mA / number_of_uncertain_genes)
    print("# genes     that have no 6mA: ", uncertain_gene_without_6mA, uncertain_gene_without_6mA / number_of_uncertain_genes)

    print()
    print("Among N/A genes: ")
    print("The average A   density (per bp)      (pro)  is: ", na_pro_A_number / na_pro_bp_number)
    print("The average A   density (per bp)      (gene) is: ", na_genetic_A_number / na_genetic_bp_number)
    print("The average 6mA density (per Adenine) (pro)  is: ", na_pro_6mA_number / na_pro_A_number)
    print("The average 6mA density (per Adenine) (gene) is: ", na_genetic_6mA_number / na_genetic_A_number)
    print("# promoters that have    6mA: ", na_pro_with_6mA, na_pro_with_6mA / number_of_na_genes)
    print("# promoters that have no 6mA: ", na_pro_without_6mA, na_pro_without_6mA / number_of_na_genes)
    print("# genes     that have    6mA: ", na_gene_with_6mA, na_gene_with_6mA / number_of_na_genes)
    print("# genes     that have no 6mA: ", na_gene_without_6mA, na_gene_without_6mA / number_of_na_genes)



    print()
    max_i = 0
    for i in range(len(count)):
        if count[i][1] > count[max_i][1]:
            max_i = i
    print("The gene with the most m6A is ", locus_tag[max_i], count[max_i])
    # plt.hist([1,2,3,4],[1,2,3,4],[1,2,3,4])
    # plt.show()

    # plt.plot([1, 2, 3, 4], [1, 4, 9, 16], 'ro')
    # plt.axis([0, 6, 0, 20])
    # plt.show()
    Graph.graph_Circular_with_Features_and_Methylation(seq_record, specific_methylation_positions, motifs)




    ### moving #6mA average, per box_size = 1000 bp

    box_size = 10000
    number_of_6mA_moving_avg = 0 

    while specific_methylation_positions[number_of_6mA_moving_avg][0] <= box_size:
        number_of_6mA_moving_avg += 1

    first_6mA_in_box = 0
    next_6mA_out_box = number_of_6mA_moving_avg



    # genome += genome
    # specific_methylation_positions += specific_methylation_positions


    number_of_6mA_moving_avg_list = [0] * (len(genome) - box_size)
    for i in range(len(genome) - box_size):
        if first_6mA_in_box < len(specific_methylation_positions) and i > specific_methylation_positions[first_6mA_in_box][0]:
            number_of_6mA_moving_avg -= 1
            first_6mA_in_box += 1
        if next_6mA_out_box < len(specific_methylation_positions) and i + box_size > specific_methylation_positions[next_6mA_out_box][0]:
            number_of_6mA_moving_avg += 1
            next_6mA_out_box += 1
        number_of_6mA_moving_avg_list[i] = number_of_6mA_moving_avg


    # Unsuccessful Attempt

    # number_of_6mA_moving_avg_list = [0] * (len(genome))
    # for i in range(len(genome)):
    #     if i % genomeSize > specific_methylation_positions[first_6mA_in_box % len(specific_methylation_positions)][0]:
    #         number_of_6mA_moving_avg -= 1
    #         first_6mA_in_box += 1
    #     if (i + box_size) % genomeSize > specific_methylation_positions[next_6mA_out_box % len(specific_methylation_positions)][0]:
    #         number_of_6mA_moving_avg += 1
    #         next_6mA_out_box += 1
    #     number_of_6mA_moving_avg_list[i] = number_of_6mA_moving_avg
    

    # ## Space out a little (actually no, this is useless and made it worse)
    # space_out = 1
    # for i in range(len(genome) - box_size):
    #     number_of_6mA_moving_avg_list[i] = number_of_6mA_moving_avg_list[i - i % space_out]

    plt.plot(number_of_6mA_moving_avg_list, 'o', color='black')
    plt.show()







    ################################### for circos ########################################
    print()
    print()
    print()
    print()
    

    # for i in specific_methylation_positions: 
    #     if i[1] == 1:
    #         print("hs1", i[0], i[0] + 1)
    
    # print()
    # print()
    # print("################################################################")
    # print()
    # print()

    # for i in specific_methylation_positions: 
    #     if i[1] == -1:
    #         print("hs1", i[0], i[0] + 1)


    for i in range(len(number_of_6mA_moving_avg_list)):
        if i % 1000 == 0 and i < genomeSize:
            print("hs1", i, i + 1000, number_of_6mA_moving_avg_list[i])


    print()
    print()
    print()
    print()
    ################################### for circos ########################################





print()
# print("time elapsed_: ", time.time() - start_time_)
print("time elapsed: ", timer() - start_time)


# TCTAGmA
# mAGATCT

# GCATC
# CGTAG

# CRAANNNNNNNNNTTC
# CRAANNNNNNNNCTK
# CRAANNNNNNNNDTTC
