# No this doesn't work for some reasons
import xlrd
import xlsxwriter

import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile



# df = pd.read_excel("All_RNAseq copy2.xlsx")

rnaseq_wb = xlrd.open_workbook(r"All_RNAseq copy4.xlsx") 
# rnaseq_wb = xlrd.open_workbook("PG3_RNAseq.xlsx")
print()  
print()  
print()  
print(rnaseq_wb.nsheets)
for name in rnaseq_wb.sheets():
    print(name)
print("=========")
for name in rnaseq_wb.sheet_names():
    print(name)


# rnaseq_sheet = rnaseq_wb.sheet_by_name("All_RNAseq copy")
rnaseq_sheet = rnaseq_wb.sheet_by_index(0)
print(rnaseq_wb._all_sheets_count)
print("===================================================================")

for i in range(7, 31):
    print("now on PG" + str(i))
    wb = xlsxwriter.Workbook("PG" + str(i) + "_RNAseq.xlsx")
    sheet1 = wb.add_worksheet()
    newRow = 0
    for row in range(0, rnaseq_sheet.nrows):
        if rnaseq_sheet.cell_value(row, 2) == "PG0" + str(i) or rnaseq_sheet.cell_value(row, 2) == "PG" + str(i):
            for column in range(0, 17):
                sheet1.write(newRow, column, rnaseq_sheet.cell_value(row, column))
            newRow += 1
    wb.close()


