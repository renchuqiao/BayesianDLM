# Bayesian Dynamic Linear Model
This is an implementation of Bayesian Dynamic Linear Model by Chuqiao Ren and Ruilin Zhong at Columbia University  
There are three folders:  
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
