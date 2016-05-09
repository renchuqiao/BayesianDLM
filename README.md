# Bayesian Dynamic Linear Model
This is an implementation of Bayesian Dynamic Linear Model by Chuqiao Ren and Ruilin Zhong at Columbia University  
This repo has the following folders:  
- data_process: this folder contains all the python scripts that process B.Pseudomallei data and the processed B.Pseudomallei data. 
  - data_process.py: This is the python script that process B.Pseudomallei data  
  - processed_data.csv: This csv file contains all genes  
  - processed_data_chromosome_1.csv: This csv file contains all genes from chromosome 1  
  - processed_data_chromosome_2.csv: This csv file contains all genes from chromosome 2  
- DLM: This folder contains the Matlab implementation of DLM. main.py is the main script that would call the function ltpdf.m
- simulation: This folder contains the python script that simulate DLM data  
  - simulation_new.py: This is the python script that simulate data based on graphical model  
  - simulation_10genes.xlsx: this is the simulated data: 10 genes with 47 data points. This is the small input file that can be tested on our DLM code.  
  - simulation.py: This is the python script that simulate data based on other paper. This is not used any more.  
- LTVAR: This folder contains LT-VAR model. It is implemented in ox programming language. In order to run the files, you have to first install the ox programming language. Ltvar_ex.ox is the main function and it will call LTVAR.ox. Note that this code is adapt from the following paper: Nakajima, Jouchi, and Mike West. "Bayesian analysis of latent threshold dynamic models." Journal of Business & Economic Statistics 31.2 (2013): 151-164.   
- results: This folder contains all the figures that we referred to in our report. Please refer to the README in this folder for more details.  
- SSClust: This folder contains all the R script for SSClust method. Note that this code is adapt from the following paper:  Ma, Ping, et al. "A data-driven clustering method for time course gene expression data." Nucleic Acids Research 34.4 (2006): 1261-1269.  

You can find the processed data in the data_process folder. However, if you want to download original data, please go to http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2365 and the dataSet record is GDS2365.  
Furthermore, our supplementary figures are in the results folder. 
