# Frankenstein SDM: A Composite Species Distribution Modelling Framework

This repository contains the R code used to implement the **Frankenstein Species Distribution Model (SDM)**, a composite modelling framework that integrates multiple spatio-temporal structures to account for environmental heterogeneity across ecological subregions.

The approach was developed and tested using both simulated and real-world data (see manuscript for details). This repository includes the necessary scripts to **simulate spatio-temporal abundance data** and **fit different SDM configurations**, including the Frankenstein model.

## Repository Structure

├── 01_simulate_dataset.R # Script to generate simulated spatio-temporal dataset
├── 02_fit_SDM_Frankenstein.R # Script to fit SDM configurations and Frankenstein model
├── /data/ # Folder to store output datasets (optional)
├── /figures/ # Folder to store plots or results (optional)


## Description of Scripts

- **01_simulate_dataset.R**  
  This script generates a synthetic dataset with two ecological regions (A and B), each exhibiting different spatio-temporal dynamics. Region A is modelled with a progressive AR(1) process, while Region B follows a persistent spatial structure with a linear temporal trend. Abundance is simulated from a Negative Binomial distribution with observational error.

- **02_fit_SDM_Frankenstein.R**  
  This script fits four individual spatio-temporal SDM configurations (*only-intercept, persistent, opportunistic, progressive*) to each subregion and then constructs the **Frankenstein SDM**, which integrates the best configuration per region. The output includes model comparisons (e.g., WAIC scores) and visual diagnostics.
