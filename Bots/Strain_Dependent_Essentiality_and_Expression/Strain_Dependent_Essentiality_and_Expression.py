

from scipy.stats import pearsonr
from scipy.stats import spearmanr
import numpy as np

import xlrd
import xlsxwriter
import csv

from timeit import default_timer as timer




startTime = timer()

print()
print()
print()

# Generate a dictionary (key: locus tag, value: cluster) and (key: locus tag, value: SubCluster)
locusTagToCluster = {}
locusTagToSubCluster = {}
lineCount = 0
with open("../../Data/All_RNAseq.csv") as csvFile:
    csvReader = csv.reader(csvFile, delimiter = ",")
    # for main Clusters and Sub_clusters
    for row in csvReader:
        if lineCount == 0:
            print("first Line of All_RNAseq.csv")
        else: 
            locusTagToCluster["['" + row[4] + "']"] = row[0]
            locusTagToSubCluster["['" + row[4] + "']"] = row[1]
        lineCount += 1




# get the target strain dependent essential clusters in 101.csv
targetClusters = []
targetClusterDescriptions = []
essentialityExpressionByCluster = {}
lineCount = 0
with open("101.csv") as csvFile:
    csvReader = csv.reader(csvFile, delimiter = ",")
    for row in csvReader:
        if lineCount == 0:
            print("first Line of 101.csv")
        else: 
            targetClusters.append(row[0])
            targetClusterDescriptions.append(row[1])
            essentialityExpressionByCluster[row[0]] = [[], []]          # for each cluster, has a array of essentiality and corresponding expression
        lineCount += 1




filesForStrain = ["BHN97 (No Motifs!)", "CT22F", "D39", "GA22F", "Taiwan19F", "TIGR4", "PG1", "PG2", "PG3", "PG4", "PG5", "PG6", "PG7", "PG8", "PG9", "PG10", "PG11", "PG12", "PG13", "PG14", "PG15", "PG16", "PG17", "PG18", "PG19", "PG20", "PG21", "PG22", "PG23??!!?!!", "PG24", "PG25", "PG26", "PG27 (No Motifs!)", "PG28", "PG29", "PG30"]


foundCluster = 0
cannotFindCluster = 0

NANumber = 0
for filepaths in filesForStrain:
    print("working on " + filepaths)
    wb = xlrd.open_workbook("../../Combined Result/" + filepaths + "/oragnized Methylation Data by Gene.xlsx")   
    sheet = wb.sheet_by_index(0)


    for i in range(1, sheet.nrows):
        locusTag = sheet.cell_value(i, 2)

        validNumber = True
        try:
            expressionMean = float(sheet.cell_value(i, 12))
            essentialityZbar = float(sheet.cell_value(i, 14))
        except ValueError:
            validNumber = False

        if locusTag in locusTagToCluster:
            foundCluster += 1
            if locusTagToCluster[locusTag] in targetClusters:
                if validNumber:
                    essentialityExpressionByCluster[locusTagToCluster[locusTag]][0].append(essentialityZbar)
                    essentialityExpressionByCluster[locusTagToCluster[locusTag]][1].append(expressionMean)
                else:
                    NANumber += 1
            elif locusTagToSubCluster[locusTag] in targetClusters:
                if validNumber:
                    essentialityExpressionByCluster[locusTagToSubCluster[locusTag]][0].append(essentialityZbar)
                    essentialityExpressionByCluster[locusTagToSubCluster[locusTag]][1].append(expressionMean)
                else:
                    NANumber += 1
        else:
            cannotFindCluster += 1


# Writing 

wb = xlsxwriter.Workbook("101_Strain_Dependent_Essentiality_and_Expression.xlsx")
ws = wb.add_worksheet()

ws.write(0, 0, "Cluster")
ws.write(0, 1, "Description")
ws.write(0, 2, "Pearson's Correlation")
ws.write(0, 3, "Spearman's Correlation")
ws.write(0, 5, "Essentialities and Expressions")

for i in range(1, len(targetClusters)):
    currentCluster = targetClusters[i - 1]
    currentClusterDescription = targetClusterDescriptions[i - 1]

    essentialities = essentialityExpressionByCluster[currentCluster][0]
    expressions = essentialityExpressionByCluster[currentCluster][1]
    corrP, pValueP = pearsonr(essentialities, expressions)
    corrS, pValueS = spearmanr(essentialities, expressions)
    if np.isnan(pValueP):
        pValueP = "N/A"
    if np.isnan(pValueS):
        pValueS = "N/A"

    ws.write(i * 2 - 1, 0, currentCluster)
    ws.write(i * 2 - 1, 1, currentClusterDescription)
    ws.write(i * 2 - 1, 2, corrP)
    ws.write(i * 2, 2, pValueP)
    ws.write(i * 2 - 1, 3, corrS)
    ws.write(i * 2, 3, pValueS)

    for j in range(0, len(essentialities)):
        ws.write(i * 2 - 1, j + 5, essentialities[j])
        ws.write(i * 2, j + 5, expressions[j])



ws2 = wb.add_worksheet()

ws2.write(0, 0, "Cluster")
ws2.write(0, 1, "Description")
ws2.write(0, 2, "Pearson's Correlation")
ws2.write(0, 3, "p value")
ws2.write(0, 4, "Spearman's Correlation")
ws2.write(0, 5, "p value")

for i in range(1, len(targetClusters)):
    currentCluster = targetClusters[i - 1]
    currentClusterDescription = targetClusterDescriptions[i - 1]

    essentialities = essentialityExpressionByCluster[currentCluster][0]
    expressions = essentialityExpressionByCluster[currentCluster][1]
    corrP, pValueP = pearsonr(essentialities, expressions)
    corrS, pValueS = spearmanr(essentialities, expressions)
    if np.isnan(pValueP):
        pValueP = "N/A"
    if np.isnan(pValueS):
        pValueS = "N/A"

    ws2.write(i , 0, currentCluster)
    ws2.write(i , 1, currentClusterDescription)
    ws2.write(i , 2, corrP)
    ws2.write(i , 3, pValueP)
    ws2.write(i , 4, corrS)
    ws2.write(i , 5, pValueS)



wb.close()



print()
print("total time: ", timer() - startTime)