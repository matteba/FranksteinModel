# Frankenstein SDM: A Composite Species Distribution Modelling Framework

This repository contains the R code used to implement the **Frankenstein Species Distribution Model (SDM)**, a composite modelling framework that integrates multiple spatio-temporal structures to account for environmental heterogeneity across ecological subregions.

The approach was developed and tested using both simulated and real-world data (see manuscript for details). This repository includes the necessary scripts to **simulate spatio-temporal abundance data** and **fit different SDM configurations**, including the Frankenstein model.

## Repository Structure

├── 01_simulate_dataset_2reg.R # Script to generate simulated spatio-temporal dataset composed of 2 regions

├── 02_FIT&CV-FRANKENSTEIN.R # Script to fit SDM configurations and Frankenstein model for 2 regions

├── 03_COVARIATE SIM&FIT&CV.R # Script to simulate covariate + fit SDM configurations and Frankenstein model + CV for 2 regions

├── 04_simulate_dataset_4reg.R # Script to generate simulated spatio-temporal dataset composed of 4 regions

├── 05_FIT&CV-FRANKENSTEIN.R # Script to fit SDM configurations and Frankenstein model for 4 regions

├── 06_Projection.R # Script to project the three SDM configurations and Frankenstein model and create the contour plot

├── 07_REAL-DATASET-FIT&CV-FRANKENSTEIN.R # Script to fit SDM configurations and Frankenstein model for the real dataset

├── /data/ # Folder to store output datasets 

├── /figure/ # Folder to store plots or results


## Description of Scripts

- **01_simulate_dataset.R**  
  This script generates a synthetic dataset with two ecological regions (A and B), each exhibiting different spatio-temporal dynamics. Region A is modelled with a progressive AR(1) process, while Region B follows a persistent spatial structure with a linear temporal trend. Abundance is simulated from a Negative Binomial distribution with observational error.

- **02_fit_SDM_Frankenstein.R**  
  This script fits four individual spatio-temporal SDM configurations (*only-intercept, persistent, opportunistic, progressive*) to each subregion and then constructs the **Frankenstein SDM**, which integrates the best configuration per region. The output includes model comparisons (e.g., WAIC scores) and visual diagnostics.

  -**03_Projection.R**
This script is to project the three single configurations and the Frankestein SDM using the first simulated dataset of the list object that can be found in the repository /data/ . Finally the contour plot is created for the observed data and projected output of the three configuration and the Frankestein model. This is to create figure 3 in the manuscript 
