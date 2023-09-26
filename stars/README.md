Companion material for the paper "Stars In Their Constellations: Great person or great team?"
by Mindruta, Bercovitz, Mares, and Feldman.

We provide the following information to replicate the results in the paper: 
1.	Data files

stars_replication_step1	: Anonymized data on “PI-constellation” collaborations (matches). An observation is a tuple of {market identifier, PI identifier, Constellation identifier, Match indicator}. PIs are denoted as upstream agents; constellations as donwstream agents. Data is described in Sections 3.2 and 4.1 of the paper. 

stars_replication_step2: This file is an extension of stars_replication_step1. Relative to that file, it contains two additional columns (removeup and removedown). These variables take the value of 1 for the upstream, respectively, downstream agents whose contribution intervals need to be calculated. 

savedGroupsReplication.m Contains the list of 2500 subsamples used for generating the 95% confidence intervals in the paper. 

These input files are necessary to replicate the results in the paper, but the code allows users to work with their own data. 

2. Code
   The main analyses in the paper were done in Mathematica (12.3). We used Stata for descriptive statistics (Tables 1, 3a, 3b, & 4). Here, we provide the Mathematica code for the maximum score estimator (Equation 1) and the code for Step 2 in the paper (calculation of contribution intervals).
   In addition, we provide a code in R for the calculation of the contribution intervals. The code in R takes as input the estimates of the matching function obtained in Mathematica 12.3 and produces the contribution intervals.
   
The credit for the code goes to Theodore Chronis and Panaghis Mavrokefalos. A very early version of the R code was created by Christina Tatli. The code was created in close consultation with Denisa Mindruta.

mse_replication_point_estimates_and_CI.nb (Mathematica): This file provides the code for estimating the matching coefficients (Equation 1). The code also generates the 95% confidence intervals. See Step 1 in the paper. 
modifyR-remove_replication_code.nb : This file provides the code for calculating the contribution intervals (See Step 2 in the paper). The matching estimates are calculated internally. 
modify_workbook.R : This code calculates the contribution intervals (Step 2 in the paper). The matching estimates need to be provided by the user. 
