#This file is used to simulate data for validating our dynamic linear model
#Based on work from http://bioinformatics.rutgers.edu/Supplements/ExpAna/
#Written by Chuqiao Ren

import numpy as np
import math

num = 500

def c1(x):
    return 0.15 * x - 0.7 + sample_normal(1, 0.3)

def c2(x):
    return sample_normal(1, 0.6)

def c3(x):
    return (0-0.3) * x - 0.3 + sample_normal(1,0.3)

def c4(x):
    return sample_normal(1,0.1) * math.sin(1.2 * sample_normal(1, 0.05) * x + 0.8 * 2 * math.pi) + sample_normal(0,0.4)

def sample_normal(mean, sigma):
    return np.random.normal(loc=mean, scale=sigma)

def simulation(num):
    print "start simulation"
    step = (6 * math.pi - 2)/num
    all = []
    for i in range(5):
        c1_data = []
        c2_data = []
        c3_data = []
        c4_data = []
        for t in range(num):
            x = t * step + 2
            c1_data.append(c1(x))
            c2_data.append(c2(x))
            c3_data.append(c3(x))
            c4_data.append(c4(x))
        all.append(c1_data)
        all.append(c2_data)
        all.append(c3_data)
        all.append(c4_data)

    #Extra sample for c2
    c2_data = []
    for t in range(num):
        x = t * step + 2
        c2_data.append(c2(x))
    all.append(c2_data)

    return all

def main():
    global num

    raw_data = simulation(num)
    #Initialize the data matrix
    data_matrix = []
    for row in range(num + 1):
        temp = []
        for col in range(21 + 1):
            temp.append(0)
        data_matrix.append(temp)

    #Construct data matrix
    for row in range(len(data_matrix)):
        for col in range(len(data_matrix[0])):
            if row == 0 :
                if col == 0:
                    data_matrix[row][col] = "data"
                else:
                    data_matrix[row][col] = col
            else:
                if col == 0:
                    data_matrix[row][col] = row
                else:
                    data_matrix[row][col] = raw_data[col - 1][row - 1]

    #Write to file
    file_w = open("./simulated_data.csv", 'w')

    for row in range(len(data_matrix)):
        line = ''
        for col in range(len(data_matrix[0])):
            line += str(data_matrix[row][col])
            line += ","
        line = line[:-1] + '\n'
        file_w.write(line)

    file_w.close()

main()
