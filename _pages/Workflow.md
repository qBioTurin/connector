---
layout: single
title: "CONNECTOR - Workflow"
permalink: /workflow/
--- 

# Workflow

The CONNECTOR package is developed to fit and cluster longitudinal data (e.g., time course data), according to the model proposed by James and Sugar (2003).  This method is particularly effective when the observations are sparse, irregularly spaced, and curves are sampled at different times. Their functional clustering is obtained through two phases. Firstly, the dimension of the input data is reduced. The original infinite-dimensional problem is converted into a finite-dimensional problem using basis functions (natural cubic splines) with a random-effects model for the coefficients.
The unknown cluster memberships are treated as missing data and included into the Expectation and Maximization algorithm for the estimation of the functional clustering model parameters. To suggest the proper number of clusters k, the fDB index and the Elbow plot are returned.

The workflow implemented is reported in Figure.
![](/assets/images/Workflow.png)

In the first step the input data, i.e. the growth data and an annotation file in which each sample is associated with relevant information, are loaded. Once a CONNECTOR list is created, this data structure will be updated with the new results at each step of the analysis in order to easily browse all data and results obtained from the further analysis.  Graphical visualization of the main results is given at each step of the analysis to easily mine the results. A pre-processing step to verify the sparsity of the time points and to possibly define a time at which truncate the growth data is implemented.



