
import numpy as np

import xlrd
import xlsxwriter
import csv

from timeit import default_timer as timer



startTime = timer()


# clusterIDs = set()
# lineCount = 0
# with open("../Data/All_RNAseq.csv") as csvFile:
#     csvReader = csv.reader(csvFile, delimiter = ",")
#     for row in csvReader:
#         if lineCount == 0:
#             print("first Line")
#         else:
#             clusterIDs.add(int(row[0]))
#         lineCount += 1

# print(sorted(clusterIDs))
# print(len(clusterIDs))

### result: there are 2973 Cluster ids ranging from 1 to 4103



locusTagToCluster = {}
lineCount = 0
with open("../Data/All_RNAseq.csv") as csvFile:
    csvReader = csv.reader(csvFile, delimiter = ",")
    for row in csvReader:
        if lineCount == 0:
            print("first Line")
        else: 
            locusTagToCluster["['" + row[4] + "']"] = row[0]
        lineCount += 1




clusters = []
# for i in range(0, 4104):
#     clusters.append({"strainsWith": 0, "methylationDensities": []})       # index number is cluster id;  "strainsWith" are strains that has this cluster;   "methylationDensities": only includes strains that has this cluster and has 6mA methylation in this cluster

clusters = [{"strainsWith": 0, "methylationDensitiesG": [], "methylationDensitiesP": []} for _ in range(0, 4104)]

filesForStrain = ["BHN97 (No Motifs!)", "CT22F", "D39", "GA22F", "Taiwan19F", "TIGR4", "PG1", "PG2", "PG3", "PG4", "PG5", "PG6", "PG7", "PG8", "PG9", "PG10", "PG11", "PG12", "PG13", "PG14", "PG15", "PG16", "PG17", "PG18", "PG19", "PG20", "PG21", "PG22", "PG23??!!?!!", "PG24", "PG25", "PG26", "PG27 (No Motifs!)", "PG28", "PG29", "PG30"]

foundCluster = 0
cannotFindCluster = 0

for filepaths in filesForStrain:
    print("working on " + filepaths)
    wb = xlrd.open_workbook("../Combined Result/" + filepaths + "/oragnized Methylation Data by Gene.xlsx")   
    sheet = wb.sheet_by_index(0)
    for i in range(1, sheet.nrows):
        locusTag = sheet.cell_value(i, 2)
        methylationDensityP = float(sheet.cell_value(i, 10))
        methylationDensityG = float(sheet.cell_value(i, 11))
        if locusTag in locusTagToCluster:
            foundCluster += 1
            clusterID = int(locusTagToCluster[locusTag])
            clusters[clusterID]["strainsWith"] += 1
            if methylationDensityP != 0:
                if methylationDensityP == np.inf:
                    print(locusTag)
                clusters[clusterID]["methylationDensitiesP"].append(methylationDensityP)
            if methylationDensityG != 0:
                clusters[clusterID]["methylationDensitiesG"].append(methylationDensityG)
        else:
            cannotFindCluster += 1


wb = xlsxwriter.Workbook("MethylationByCluster.xlsx")
ws = wb.add_worksheet()

ws.write(0, 0, "Cluster")
ws.write(0, 1, "Present in")
ws.write(0, 2, "Methylated in")
ws.write(0, 3, "Avg Methylation Density")
ws.write(0, 4, "Methylation Std")
ws.write(0, 5, "Methylated in (pro)")
ws.write(0, 6, "Avg Methylation Density (pro)")
ws.write(0, 7, "Methylation Std (pro)")

excelColumn = 0
for i in range(1, len(clusters)):
    # print()
    # print(clusters[i]["strainsWith"])
    # print(clusters[i]["methylationDensitiesG"])
    # print(clusters[i]["methylationDensitiesP"])
    if clusters[i]["strainsWith"] != 0:
        strainsWithCluster = clusters[i]["strainsWith"]

        methylationDensityCountG = len(clusters[i]["methylationDensitiesG"])
        if methylationDensityCountG != 0:
            methylationDensityMeanG = np.mean(clusters[i]["methylationDensitiesG"])
            methylationDensityStdG = np.std(clusters[i]["methylationDensitiesG"])
        else:
            methylationDensityMeanG = "N/A"
            methylationDensityStdG = "N/A"

        methylationDensityCountP = len(clusters[i]["methylationDensitiesP"])
        if methylationDensityCountP != 0:
            methylationDensityMeanP = np.mean(clusters[i]["methylationDensitiesP"])
            methylationDensityStdP = np.std(clusters[i]["methylationDensitiesP"])
        else:
            methylationDensityMeanP = "N/A"
            methylationDensityStdP = "N/A"

        excelColumn += 1
        ws.write(excelColumn, 0, i)
        ws.write(excelColumn, 1, strainsWithCluster)
        ws.write(excelColumn, 2, methylationDensityCountG)
        ws.write(excelColumn, 3, methylationDensityMeanG)
        ws.write(excelColumn, 4, methylationDensityStdG)
        ws.write(excelColumn, 5, methylationDensityCountP)
        ws.write(excelColumn, 6, methylationDensityMeanP)
        ws.write(excelColumn, 7, methylationDensityStdP)

wb.close()

print()
print("Found " + str(foundCluster) + " genes in Clusters")
print("There are " + str(cannotFindCluster) + " genes without associated known cluster")


print()
print("total time: ", timer() - startTime)