# Bayesian Dynamic Linear Model
__This is an implementation of Bayesian Dynamic Linear Model__     
Author: Chuqiao Ren and Ruilin Zhong          
@Columbia University      
Final Project for CBMF W4761 Computational Genomics Spring 2016     
Special thanks to  Dr. Itsik Pe'er and Shuo Yang  
   
This repo has the following folders:  
- data_process: this folder contains all the python scripts that process B.Pseudomallei data and the processed B.Pseudomallei data. 
  - data_process.py: This is the python script that process B.Pseudomallei data  
  - processed_data.csv: This csv file contains all genes  
  - processed_data_chromosome_1.csv: This csv file contains all genes from chromosome 1  
  - processed_data_chromosome_2.csv: This csv file contains all genes from chromosome 2  

- DLM: This folder contains the Matlab implementation of DLM. main.py is the main script that would call the function ltpdf.m  

- simulation: This folder contains the python script that simulate DLM data  
  - simulation_new.py: This is the python script that simulate data based on graphical model  
  - simulated_data.xlsx: this is the simulated data: 10 genes with 47 data points. This is the input file for our Matlab code.  
  - simulation.py: This is the python script that simulate data based on other paper. This is not used any more.  
  - simulated_data.csv: This is the simulated data raw data from python script. We used Number on Mac Book Pro to transform this file to xlsx file.  

- LTVAR: This folder contains LT-VAR model. It is implemented in ox programming language. In order to run the files, you have to first install the ox programming language. Ltvar_ex.ox is the main function and it will call LTVAR.ox. Note that this code is adapt from the following paper: Nakajima, Jouchi, and Mike West. "Bayesian analysis of latent threshold dynamic models." Journal of Business & Economic Statistics 31.2 (2013): 151-164.  

- results: This folder contains all the figures that we referred to in our report. Please refer to the README in this folder for more details.   

- SSClust: This folder contains all the R script for SSClust method. Note that this code is adapt from the following paper:  Ma, Ping, et al. "A data-driven clustering method for time course gene expression data." Nucleic Acids Research 34.4 (2006): 1261-1269.  

- Example: This folder contains the sample input file along with the sample output for our implementation of DLM.   

You can find the processed data in the data_process folder. However, if you want to download original data, please go to http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2365 and the dataSet record is GDS2365.  
Furthermore, our supplementary figures are in the results folder.  
### Instruction on how to run SSClust
In order to run the SSClust, first, download and install R for your operating system from: http://cran.us.r-project.org/ Then, open R in the terminal and type in   
```R
>chooseCRANmirror()
>install.packages(c("mvtnorm", "gss")) 
```
This will install all the dependencies. Then run SSClust.R as the main script.   

Note: in order to perform well, you need to specify the correct path to the file in line 28 in the SSClust.R script
```R
my.data = read.table("processed_data_chromosome_1_T.txt", header=T, na.strings =" ", sep="\t")
``` 

### Instruction on how to run LT-TVR
In order to run LT-TVR, first, download and install ox programming language from http://www.doornik.com/ox/  
The main script is Ltvar_ex.ox  
To run, please open Ltvar_ex.ox and then specify the data path in line 22 
```
my = loadmat("usdata.xls");	//data
``` 
Then specify the gene name in line 26
```
asvar = {"p", "x", "i"};
```
And lastly, specify the maximum number of iterations for MCMC in line 45
```
Ltvar.MCMC(50000);			//MCMC estimation
```

### Instruaction on how to run DLM
In order to run DLM, you have to first download and install Matlab from https://www.mathworks.com/campaigns/products/ppc/google/matlab-trial-request.html?s_eid=ppc_5852767762&q=download%20matlab&refresh=true  

The input file should be a Excel file containing a matrix of data with rows being time points and columns being genes.  

The example folder contains the test input file. It contains 25 genes from real data with 47 time points.   
In order to run DLM
* Specify the correct path in main.m file in the DLM folder. (line 5 `[Ylog,Ynames,time]=xlsread('25toy.xlsx');`)  
* Follow the instructions in the console  
* Output prediction figures  
* Output prediction error figures  
* Output correlation figures
Note that it may take a while for the program to run, and the outputs are predictions. The sampler output is in the output folder in Example folder.  


