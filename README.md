
# Spatial Models for the Cell Specific Analysis of DNAm in Whole Blood

Supplementary simulation study and analysis R scripts accompanying the manuscript "A Bayesian spatial correlation model improves detection of cell specific effects of DNA methylation in whole blood"

## Usage
Key tasks are accomplished by the main scripts located in the project directory. To fit the Bayesian models OpenBUGS has to first be installed by the user. Within each script, key parameters that have to be modified are: "inpath" which should point to the project directory and "outpath" which specifies were created files will be saved. If there are additional things that need to be specified a script will note this with comments. Commands that install required R packages are provided at the top of each script.

## Main Scripts

### processGSE147430.R
Performs simple analyses based on linear models to identify differences in methylation between smokers and non-smokers in isolated CD8 T-cells.

### fit_Bayes_models_to_GSE53045.R
Fits Bayesian spatial and non-spatial models to the GSE53045 dataset in order to obtain estimates of cell specific differences between smokers and non-smokers. This analysis exclusively targets chromosome 1 and was presented in the manuscript.

### perform_tests_GSE53045.R
Performs statistical tests for a difference in methylation between smokers and non-smokers in CD8 T-cells, and compares results derived from bulk samples with results obtained from isolated cells. For the Bayesian models, estimates obtained from "fit_Bayes_models_to_GSE53045.R" are loaded to construct test statistics. For TCA and CellDMC, model fits and statistical tests based on GSE53045 are performed in one go, since their run-times are much shorter. After all hypothesis tests are conducted, the script loads analysis results for CD8 T-cell isolates obtained from the "processGSE147430.R" script, and performs comparisons. 

### marginal_simulation_A.R
Runs simulation scenarios of type A, in which there is only a single study arm.

### two_arm_simulation_B2.R
Runs simulation scenarios of type B, in which there are two study arms, cases and controls, for which methylation levels differ between certain blood cell types.

## The "utility_scripts" Folder

### marginal_model_functions.R
Contains functions that define marginal models and their MCMC estimation procedure utilized in simulation study A.

### two_arm_model_functions_sim.R
Contains functions that define two-arm models and their MCMC estimation procedure utilized in simulation study B.

### two_arm_model_functions_analysis.R
Contains functions that define two-arm models and their MCMC estimation procedure utilized in the real data analysis. These models are equivalent to those presented in "two_arm_model_functions_sim.R", but they have been stripped down to not sample RMSPE values for the sake of processing speed.

### utility_functions.R
Contains utility functions required to simulate data and apply the approach. Includes both: the cluster generating algorithm ("construct.cluster.structures") and the scaling factor prediction function ("prfun").

## The "deconvolution" Folder
Contains R objects and functions required to perform reference based deconvolution following the approach of Salas and Koestler (PMID: 29843789).

## The "data" Folder

### DNA_methylation_data.RData
Contains DNA methylation data (PMID: 29843789) of isolated blood leukocytes, also available under the GEO accession GSE110554. This data is used to simulate artificial bulk sample DNAm data in whole blood.
