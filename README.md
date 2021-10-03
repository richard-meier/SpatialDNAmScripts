
# Spatial Models for the Cell Specific Analysis of DNAm in Whole Blood

Supplementary simulation study and analysis R scripts accompanying the manuscript "A Bayesian spatial correlation model improves detection of cell specific effects of DNA methylation in whole blood"

## Usage
Key tasks are accomplished by the main scripts located in the project directory. To fit the Bayesian models OpenBUGS has to first be installed by the user. Within each script, key parameters that have to be modified are: "inpath" which should point to the project directory and "outpath" which specifies were created files will be saved. It is important that paths are specified in Unix style. If there are additional things that need to be specified a script will note this with comments. Commands that install required R packages are provided at the top of each script.

## Main Scripts

### analysis_example.R
Example script that demonstrates how to perform cell specific analysis of a small example dataset based on the SCM2 and TM models in conjunction with the cluster generating algorithm.

### process_GSE42861.R
Fits Bayesian spatial and non-spatial models to the GSE42861 dataset in order to obtain estimates of cell specific differences between smokers and non-smokers. Only super clusters containing CpGs included in "combined_set.csv" will be utilized.

### perform_tests_GSE42861.R
Performs statistical tests for a difference in methylation between ever smokers and never smokers. For the Bayesian models, estimates obtained from "process_GSE42861.R" are loaded to construct test statistics. For TCA and CellDMC, model fits and statistical tests based on GSE42861 are performed in one go, since their run-times are much shorter.

### marginal_simulation_A.R
Runs simulation scenarios of type A, in which there is only a single study arm.

### two_arm_simulation_B2.R
Runs simulation scenarios of type B, in which there are two study arms, cases and controls, for which methylation levels differ between certain blood cell types.

## The "utility_scripts" Folder

### marginal_model_functions.R
Contains functions that define marginal models and their MCMC estimation procedure utilized in simulation study A.

### two_arm_model_functions_sim.R
Contains functions that define two-arm models and their MCMC estimation procedure utilized in simulation study B.

### two_arm_model_functions_three_variables.R
Contains functions that define two-arm models and their MCMC estimation procedure utilized in the real data analysis. These models adjust for two additional covariates and have been stripped down to not sample RMSPE values for the sake of processing speed.

### utility_functions.R
Contains utility functions required to simulate data and apply the approach. Includes both: the cluster generating algorithm ("construct.cluster.structures") and the scaling factor prediction function ("prfun").

## The "deconvolution" Folder
Contains R objects and functions required to perform reference based deconvolution following the approach of Salas and Koestler (PMID: 29843789).

## The "data" Folder

### DNA_methylation_data.RData
Contains DNA methylation data (PMID: 29843789) of isolated blood leukocytes, also available under the GEO accession GSE110554. This data is used to simulate artificial bulk sample DNAm data in whole blood.

### example_data.RData
Contains a small example dataset that is used in the "analysis_example.R" script.

### combined_set.csv
Contains gold-standard, smoking-associated CpGs as defined identified by Joehanes et al., 2016 (PMID: 27651444).
