# A Monotone Single-Index Accelerated Failure Time Model for Correlated Heavily Right-Censored Data

## Description

This repository provides the R, Stan, and supporting code needed to implement the proposed method from the titled paper, along with several comparison methods used in the simulation studies.

The repository includes:

1. R and Stan code for fitting the proposed monotone single-index accelerated failure time model.
2. Code for reproducing the simulation results, including the figures and tables reported in the paper.
3. Simulation datasets, provided under the `Data/` directory.
4. Code for implementing four comparison methods considered in the paper.

Before running the analysis files, please run:

```r
source("Directory.R")
```
to generate the necessary folders for saving the results. 

## Repository Structure

The repository contains four main folders, each corresponding to a candidate model. These models are defined by two modeling choices:

The error distribution:

(a) Normal <br>
(b) Generalized Extreme Value (GEV)

The basis used to model the single-index function:

(a) Bernstein Polynomial (BP) <br>
(b) Gaussian Process (GP) 

Each folder is named according to the corresponding method and contains the files needed to reproduce the simulation results.

## Usage 
Each method-specific folder contains a sequence of R and Stan files. To reproduce the simulation results from the paper, run the R files in the specified order.

The user should specify <code>sim_id</code>, which indicates the simulated dataset used under the two simulation setups described in the paper. By default, 

```r
sim_id = 1
```
The files are organized as follows:

### 1. Stan Model File

```r
*.stan
```
This file contains the Stan model used for the behind-the-scenes Hamiltonian Monte Carlo algorithm.
No edits are necessary for this file.

### 2. Initial Analysis Using All Variables

```r
*_main_simulation.R
```

This script performs the initial model fitting using all available variables.

### 3. Variable Selection Procedure

```r
analyze_*_simulation.R
```

This script implements the variable selection procedure using the sequential 2-means algorithm of Li and Pati (2017)

#### Reference: <br>
Li H. &amp; Pati D. (2017). Variable selection using shrinkage priors,  <em>Computational Statistics &amp; Data Analysis</em>, 107, 107--119.

### 4. Secondary Analysis Using Selected Variables

```r
*_main_simulation_SelectedVariables.R
```

This script refits the model using only the variables selected in the previous step.

### 5. Final Analysis and Result Generation

```r
analyze_*_simulation_SelectedVariables.R
```

This script produces the final simulation summaries, figures, and tables reported in the paper.


## Recommended Workflow

A typical workflow is:

```r
# Step 1: Create required directories
source("Directory.R")

# Step 2: Move to the folder corresponding to the selected method

# Step 3: Run the initial analysis
source("*_main_simulation.R")

# Step 4: Run the variable selection analysis
source("analyze_*_simulation.R")

# Step 5: Run the secondary analysis using selected variables
source("*_main_simulation_SelectedVariables.R")

# Step 6: Generate final simulation results
source("analyze_*_simulation_SelectedVariables.R")
```

Please replace * with the appropriate method-specific file name.

## Data
The simulation datasets are provided under the <code>Data/</code> directory. These datasets are used to reproduce the simulation studies reported in the paper.

## License

All code in this repository is free to use with proper citation of the paper: 

Su Y., Pati D  &amp; Bandyopadhyay D. (2026+). A monotone single-index accelerated failure time model for correlated heavily right-censored data,  <em>Under Review</em>
