# biomass_loss

Code and data to reproduce analyses from "**Accelerating biomass loss from forest disturbances across Europe**"

All R code used for the analyses is provided in the folder *lib*. Intermediate datasets are provided as zip archives in the subdirectories of *data* and should be unpacked directly in their location for re-creating the analyses.

To carry out the analysis steps, the code should be run in the provided order in *lib*. Please note that for some parts of the analyses require parallel processing to achieve reasonable processing times. For model fitting (*lib/06_sensitivity_analysis.R*) a BLAS/OpenBLAS installation is needed.

R version: 4.5.0
