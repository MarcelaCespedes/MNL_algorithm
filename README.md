# MNL_algorithm
Repository for the Maximisation of Network Likelihood (MNL) algorithm, as discussed in manuscript titled _An efficient algorithm for estimating brain covariance networks_ published in PLOS One, available [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0198583).

## Reproduce the simulation study 
The **Simulation_study** folder contains all the necessary R code to generate simulated data from binary matrices <a href="http://www.codecogs.com/eqnedit.php?latex=S_1" target="_blank"><img src="http://latex.codecogs.com/gif.latex?S_1" title="S_1" /></a> and <a href="http://www.codecogs.com/eqnedit.php?latex=S_2" target="_blank"><img src="http://latex.codecogs.com/gif.latex?S_2" title="S_2" /></a> and then perform the simulated study for the MNL algorithm, as well as for the graphical LASSO (gLASSO) and Pearson's pairwise correlation (PPC) analyses and reproduce the plots in the manuscript.

Please note, that in order to re-generate the plots the code must be run in the following sequential order
1. Generate simulated data for sample sizes <a href="http://www.codecogs.com/eqnedit.php?latex=N&space;=&space;100,&space;250,&space;500,&space;1000" target="_blank"><img src="http://latex.codecogs.com/gif.latex?N&space;=&space;100,&space;250,&space;500,&space;1000" title="N = 100, 250, 500, 1000" /></a> from binary solution matrices. Each sample size consists of 10 replicates. The simulated data will be saved as an .Rdata file.
2. Run the simulation study for the MNL algorithm. As this is performed on all sample sizes, for each replicate and on both networks (4 x 10 x 2= 80 times), this may take a considerable amount of time. We recommend running the MNL algorithm in parallel for each sample size via four R sessions as each simulation run is independent of each other. The output of the simulation study will be saved as .Rdata files. 
3. Run the simulation study for the gLASSO and PPC algorithms. The results of this final stage will reproduce the plots in the paper. 

This folder includes additional functions required for the simulation study to run, such as the C++ code to evaluate the MVN density, generation of 1st, 2nd and random binary symmetric matrices, evaluate sensitivity and specificity as well as several matrix visualisatin functions.

## Supporting Information
The **Supporting_Information** folder contains all the supplementary material for the manuscript in pdf form.

## Comments or feedback
Please contact Marcela.Cespedes@csiro.au
