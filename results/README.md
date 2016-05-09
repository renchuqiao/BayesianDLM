#Results and Figures
This is the results and supplementary figures for our project. It has two parts:  
- simulated_data: This folder contains DLM results for simulated data.  
- real_data: This folder contains DLM results for B. Pseudomallei data chromosome 1 first cluster.    
- GM.png: This is the overall graphical model  

###Results for B.Pseudomallei data
Although we performed SSClust to cluster the whole genome from chromosome 1 to over 100 clusters, we will show the results for the first cluster, in order to demonstrate that our model indeed perform well. 

In the real_data folder, it has the following folder:

- covariance_between_genes: This folder contains figures to show correlations between genes (300 figures).  
- prediction: This folder contains figures to show prediction for next half hour for each gene (25 figures).  
- prediction_error: This folder contains figures to show the errors for prediction of next half hour for each gene (25 figures).  

__Prediction__: The figures in this folder shows the forecast data along with the real data point over time. Also, we report the 95% confidence credible band. The results are great; the forecast data points follows the trend of the real data point. The 95% confidence interval successfully capture the variation of the prediction.  
__Prediction Error__: The figures shows the error between predicted value and the true value at each time for each gene. Suppose `$\hat{y_t}$` is the predicted value at time `$t$` using data at or before time `$t-1$`. And `$y_t$` is the real number. Now, we can calculate the error as the difference between `$y_t$` and `$\hat{y_t}$`: `$\hat{y_t} - y_t$`. For almost all figures, the 1-step errors are small around 0. This is what we expected.  
__Covariance Between Genes__: The figures in this folder shows the correlation between each pairs of genes over time. That is, for every time t, we can capture the correlations between each pairs of genes. For most of the gene pairs, there is no correlations over time; the covariance plot is fluctuating around 0. However, this is not true for some pairs of gene. For example, in figure 181.png, it captures the covariance between gene BPSL0025 and gene BPSL0050. The correlation fluctuating around 0.3 to 0.6 over 47 time points. Another example is in figure 179.png, where the correlation between gene BPSL0025 and gene BPSL0048 fluctuates between -0.7 and -0.5. So, we may say that gene BPSL0025 is positively correlated with gene BPSL0050 but negatively correlated to gene BPSL0048.  

###Results for Simulated data:
We did simulation followed the graphical model. The hyper-parameters are specified in our report. The simulated data contains 10 genes with 47 time points.  

In the simulated folder, it has the following folder:
- covariance_between_genes (45 figures)  
- prediction (10 figures)  
- prediction_error (10 figures)  

The folder contents are the same for simulated data and real data. Please refer above for the description of each folder.

__Prediction__: The figures shows that our model works great for each gene. This validate our model.  
__Prediction Error__: The figures further indicate that our model performs well.  
__Covariance Between Genes__: The figures capture the correlation between each simulated gene. For most of the gene pairs, the covariance is around 0, indicating no correlation. However, for some gene pairs, we can see positive correlation and for some gene pairs we can see negative correlation. We will examine these correlations more in the real data set.



