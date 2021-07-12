from Bio.Data.IUPACData import ambiguous_dna_values


###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################


# # from Bio.Data.IUPACData import ambiguous_dna_values
# ambiguous_dna_values["A"]
# "A" not in ambiguous_dna_values["N"]



# i = -1
# for i in range(10):
#     print(i)
#     # break
# else:                   ## Only will execute "else" if ended loop normally
#     print("normal", i)
# print("all the time", i)

###############################################################################################
######################################## TESTING AREA #########################################
###############################################################################################



def posits_In_Seq(seq, motif):
    positions = []
    for i in range(len(seq) - len(motif)):
        for j in range(len(motif)):
            if seq[i+j] not in ambiguous_dna_values[motif[j]]:
                break
        else:       ## only if the loop ended normally without breaks, meaning match found
            positions.append(i)
    return positions

posits_In_Seq("AAAAACCCGGTTCGATCGCATGGA","CNCN")
