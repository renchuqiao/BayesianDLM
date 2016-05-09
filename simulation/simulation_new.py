##===============
##==== libraries
##===============
import numpy as np
from numpy.linalg import inv
from scipy.stats import wishart
from scipy.stats import bernoulli
import math
from numpy import linalg as LA
# import matplotlib.pyplot as plt
# import seaborn as sns

##=====================
##==== global variables
##=====================
n_time = 0 #Note: Don't change value here; change below
n_gene = 0

# n_tissue = 0
#The following parameters need to be determined by test-and-trials
#According to Barbara, they used alpha=beta=1 for the uniform on sparsity
#alpha = 1 beta = 2 is a line of y = -2x + 2
beta = [1,1]                #parameters of beta[alpha, beta]
gamma = [1,2,1,2]     #parameters of NG[mu, kappa, alpha, beta]
normalWishart = [[2,2],2,[[10,5],[5,10]],3]   #parameters of NW[mu, kappa, Lambda, v]

##=================================================
##===The Following code is adapted from Shuo's code
##=================================================

##=====================
##==== Sampling
##=====================

##==== sampling from Gaussian (with mean and std)
def sampler_Normal(mu, sigma):
    sample = np.random.normal(mu, sigma)
    return sample



##==== sampling from Wishart
def sampler_W(df, scale):
    #
    sample = wishart.rvs(df, scale, size=1, random_state=None)
#    matrix = sample[0]
    return sample


##==== sampling from Gamma
def sampler_Gamma(para1, para2):
    para2 = 1.0/para2
    x = np.random.gamma(para1, para2, 1)
    return x[0]


## ==== End of adaptation

## ==== sampling Beta
def sampler_beta(a, b):
    return np.random.beta(a, b)



## ==== Start to simulate

def simulation():
    global n_gene
    global n_time
    global beta
    global gamma
    global normalWishart

    print "start to simulate theta_0"
    theta_0 = []

    for sample in range(n_gene):
        precisionMatrix = sampler_W(normalWishart[3], normalWishart[2])
        precisionMatrix_scaled = []
        for i in range(len(precisionMatrix)):
            temp = []
            for j in range(len(precisionMatrix[0])):
                temp.append(precisionMatrix[i][j] / normalWishart[1])
            precisionMatrix_scaled.append(temp)
        mu = np.random.multivariate_normal(normalWishart[0], precisionMatrix_scaled)
        theta_0.append(np.random.multivariate_normal(mu, precisionMatrix))
    print "finish simulating theta_0"

    print "start to simulate "

    print "start to simulate all time points"
    print "finish simulating all time points"




if __name__ == '__main__':



    # DEBUG
    print "enter program..."


    # DEBUG
    print "now start preparing the data..."



    ##==================================
    ##==== loading and preparing dataset
    ##==================================
    # data_prepare()			# prepare the "dataset" and "markerset"



    # DEBUG
    print "finish data preparation..."


    # DEBUG
    print "now initializing all the variables..."

    ##================================
    ##==== initialize global variables
    ##================================
    n_time = 500			# TODO: this is tunable, and the number 400 comes from results of other more complex methods
    n_gene = 10

    #initialize normal wishart parameter
    mu = []
    precision = []
    for i in range(n_gene):
        mu.append(2)

    for i in range(n_gene):
        temp = []
        for j in range(n_gene):
            if (i == j):
                temp.append(1)
            else:
                temp.append(0)
        precision.append(temp)

    normalWishart[0] = mu
    normalWishart[2] = precision
    normalWishart[3] = n_gene + 1


    #initialize sparsity parameter
    beta = [2, 2]  # parameters of beta[alpha, beta]
    gamma = [1, 2]  # parameters of Gamma[alpha, beta]

    #DEBUG
    print "finish initializing all the variables..."

    #DEBUG
    print "now start simulation..."

    simulation()



