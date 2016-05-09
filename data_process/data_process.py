#This is an initial data process step for rice data
#Written by Chuqiao Ren

import math
import numpy as np

#Open raw data file (SOFT)
file = open("/Users/rachelren/Google Drive/old mac/CBMF4761/project/data/BP/GDS2365_full.soft")

#Read file once and get all the useful information
toRead = False
data = []
overall_data = []
for line in file:
    if not toRead:
        if "!dataset_table_begin":
            toRead = True
    else:
        line = line.strip()
        if line[0:6] == "ID_REF":
            title = line
        else:
            if "general secretory pathway protein" in line:
                data.append(line)
            elif "Chromosome 1" in line:
                overall_data.append(line)

file.close()

#Process raw data
title_list = title.split("\t")
ID = []
for title in title_list:
    if "GSM" in title:
        ID.append(title)

data_nr_list = []
data_name_list = []

for d in data:
    data_nr = d.split("\t")[2:49]
    data_name = d.split("\t")[50]
    if "null" not in data_nr:
        data_nr_list.append(data_nr)
        data_name_list.append(data_name)

overall_nr_list = []
overall_name_list = []
for m in overall_data:
    overall_nr = m.split("\t")[2:49]
    overall_name = m.split("\t")[1]
    if "null" not in overall_nr:
        overall_nr_list.append(overall_nr)
        overall_name_list.append(overall_name)


#Initialize the data matrix
data_matrix = []
for row in range(len(ID) + 1):
    temp = []
    for col in range(len(data_name_list) + 1):
        temp.append(0)
    data_matrix.append(temp)

overall_data_matrix = []
for row in range(len(ID) + 1):
    temp = []
    for col in range(len(overall_name_list) + 1):
        temp.append(0)
    overall_data_matrix.append(temp)

overall_data_matrix_T = []
for row in range(len(overall_name_list) + 1):
    temp = []
    for col in range(len(ID) + 1):
        temp.append(0)
    overall_data_matrix_T.append(temp)

#Construct data matrix
for row in range(len(data_matrix)):
    for col in range(len(data_matrix[0])):
        if row == 0 :
            if col == 0:
                data_matrix[row][col] = "data"
            else:
                data_matrix[row][col] = data_name_list[col - 1]
        else:
            if col == 0:
                data_matrix[row][col] = ID[row - 1]
            else:
                data_matrix[row][col] = 2**(float(data_nr_list[col - 1][row -1]))

for row in range(len(overall_data_matrix)):
    for col in range(len(overall_data_matrix[0])):
        if row == 0 :
            if col == 0:
                overall_data_matrix[row][col] = "data"
            else:
                overall_data_matrix[row][col] = overall_name_list[col - 1]
        else:
            if col == 0:
                overall_data_matrix[row][col] = ID[row - 1]
            else:
                overall_data_matrix[row][col] = float(overall_nr_list[col - 1][row -1])

for row in range(len(overall_data_matrix_T)):
    for col in range(len(overall_data_matrix_T[0])):
        if row == 0 :
            if col == 0:
                overall_data_matrix_T[row][col] = "data"
            else:
                overall_data_matrix_T[row][col] = ID[row - 1]
        else:
            if col == 0:
                overall_data_matrix_T[row][col] = overall_name_list[col - 1]
            else:
                overall_data_matrix_T[row][col] = float(overall_nr_list[row - 1][col -1])

#Overwrite to file
file_w = open("./processed_data.csv", 'w')
file_w_overall = open("./processed_data_chromosome_1_T.txt", 'w')

for row in range(len(data_matrix)):
    line = ''
    for col in range(len(data_matrix[0])):
        line += str(data_matrix[row][col])
        line += ","
    line = line[:-1] + '\n'
    file_w.write(line)

for row in range(len(overall_data_matrix_T)):
    line = ''
    for col in range(len(overall_data_matrix_T[0])):
        line += str(overall_data_matrix_T[row][col])
        line += "\t"
    line = line[:-1] + '\n'
    file_w_overall.write(line)

file_w.close()
file_w_overall.close()

#End of file
print "Program finished without interruptions..."


