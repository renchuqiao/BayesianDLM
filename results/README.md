#Results and Figures
This is the results and supplementary figures for our project. It has two parts:  
- simulation: This folder contains DLM results for simulated data.  
- real_data: This folder contains DLM results for B. Pseudomallei data chromosome 1 first cluster.    
- GM.png: This is the overall graphical model  

###Results for B.Pseudomallei data
Although we performed SSClust to cluster the whole genome from chromosome 1 to over 100 clusters, we will show the results for the first cluster, in order to demonstrate that our model indeed perform well. 

In the real_data folder, it has the following folder:

- covariance_between_genes: This folder contains figures to show correlations between genes.  
- prediction: This folder contains figures to show prediction for next half hour for each gene.  
- prediction_error: This folder contains figures to show the errors for prediction of next half hour for each gene.  Suppose $\hat{y_t}$ is the predicted value at time $t$ using data at or before time $t-1$. And $y_t$ is the real number. Now, we can calculate the error as the difference between $y_t$ and $\hat{y_t}$: $\hat{y_t} - y_t$. 

__Prediction__: The figures in this folder shows the forecast data along with the real data point over time. Also, we report the 95% confidence credible band. The results are great; the forecast data points follows the trend of the real data point. The 95% confidence interval successfully capture the variation of the prediction.  
__Prediction Error__: The figures  
