# Historical Functional Linear Models for Human Movement Data

### Authors: Arlo Hook, Mark Watsford, Paul Wu, Edward Gunning, Giles Hooker, and John Warmenhoven.

## About 

This repository contains code for the paper *Historical Functional Linear Models for Human Movement Data*. 
The repository is strutured as follows:

#### Functions

This folder contains functions for implementing the finite element and point wise HFLM as well as a data generation function used to generate the data for the simulation part of the study.

#### Simulation Study

This folder contains all required scripts and associated files to run the simulations presented in the study.

- **Data:** Generated data will be stored here once created using the `Data Generation.R` script which leverages the `Data Gen.rds` file that contains the relevant fPCA and generative model coefficients required.
- **Simulation:** The Scripts required to run the simulations.
- **Results:** The processed result files and where raw results from each simulation will be saved when run.
- **Analysis:** Scripts for processing raw results from the simulations and scripts for generating the figures in the paper.

#### CRP and Vector Coding Comparison

This folder contains all required scripts to reproduce the continuous relative phase and vector coding analysis and estimate the HFLM for comparison.
