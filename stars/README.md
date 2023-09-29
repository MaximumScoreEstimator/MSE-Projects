# Companion material for the paper "Stars In Their Constellations: Great person or great team?"
## by Mindruta, Bercovitz, Mares, and Feldman.

We provide the following information to replicate the results in the paper: 
1.	Data files

`stars_replication_step1`: Anonymized data on “PI-constellation” collaborations (matches). An observation is a tuple of {market identifier, PI identifier, Constellation identifier, Match indicator}. PIs are denoted as upstream agents; constellations as donwstream agents. Data (variables, structure) is described in Sections 3.2 and 4.1 of the paper. 

`stars_replication_step2`: This file is an extension of `stars_replication_step1`. Relative to that file, it contains two additional columns (removeup and removedown). These variables take the value of 1 for the upstream, respectively, downstream agents whose contribution intervals need to be calculated. 

`savedGroupsReplication.m`: Contains the list of 2500 subsamples used for generating the 95% confidence intervals in the paper. 

These input files are necessary to replicate the results in the paper, but the code allows users to work with their own data, as long as input files follow the same format.  

2. Code
   
   The main analyses in the paper were done in Mathematica (version 12.3). We used Stata for descriptive statistics (Tables 1, 3a, 3b, & 4). Since the Stata code for generating summary statistics is trivial, we do not include it here but it is available upon request. We provide the Mathematica code for the maximum score estimator (Equation 1, Step 1) and for the calculation of contribution intervals: min and max (Step 2).

In addition, we provide a code in R for the calculation of the contribution intervals (Step 2 only). The code in R takes as input the estimates of the matching function obtained in Mathematica 12.3 and produces the contribution intervals (min and max). Note also that the code for Step 2 (both the one in Mathematica and the one in R) allows users to export the configuration of new matchings after an agent is "removed" from a market. The code produces this output for each market separately. We used this functionality for our discussion in Section 5.3 of the paper. 
      
The credit for the code goes to Theodore Chronis ([tchronis](https://github.com/tchronis)) and Panaghis Mavrokefalos ([ConsideredHarmless](https://github.com/ConsideredHarmless)). The code was created in close consultation with Denisa Mindruta.

`mse_replication_point_estimates_and_CI.nb` (written in Mathematica): This file provides the code for estimating the matching coefficients (Equation 1) and the 95% confidence intervals. See Step 1 in the paper. 

`modifyR-remove_replication_code.nb` (written in Mathematica): This file provides the code for calculating min and max of contribution intervals (See Step 2 in the paper). The code starts by calculating the matching estimates internally and then proceeds to calculate the min and max of value capture intervals. The matching estimates are identical to the ones yielded by `mse_replication_point_estimates_and_CI.nb`.

`modify_workbook.R` (written in R): This code in R is the counterpart of the `modifyR-remove_replication_code.nb` in Mathematica. It calculates the contribution intervals (Step 2 in the paper), but the matching estimates need to be provided by the user as an input in the code. 

Disclaimer: The code in Mathematica was generated with version 12.3 of the software. We obtained identical results in Mathematica 13. We used Differential Evolution (DE) Method to solve the discrete optimization problem that produces the matching estimates. However, differences may exist between versions of Mathematica (and between Mathematica and R) regarding the variant of the DE method provided by default by the software producers, and these differences may impact the empirical results. Please use Mathematica versions 12.3 or 13 to ensure the results are the same.
