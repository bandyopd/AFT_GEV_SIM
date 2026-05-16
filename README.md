## A monotone single-index accelerated failure time model for correlated highly right-censored data
### Description
We provide R files and the supporting codes (stan and R) to deliver the proposed method and some comparison methods. In addition, we provide the 
simulation data sets (under directory Data), and the R file to perform analysis. Run Directory.R first to generate the necessary folders to save 
the results.
This repository contains the R codes to 
> 1. facilitate the application of the developed method in the titled paper. 
> 2. reproduce the simulation results including Figures and Tables in the titled paper.
> 3. this repo consists of four comparison methods in the paper
### Usage
Four individual folders are provided under four candidate models corresponding to choices for the error distribution (Normal or GEV) and the 
basis used to model the single index function (BP or GP). The folders are named according to the methods. 
Each folder contains the codes for reproducing the simulation results in the paper. R files need to be run in sequential order for reproduction 
purposes with user specification of sim_id (indicating the simulated data under the two setups in the paper, 1 by default). They are organized 
as follows:
•	Stan file *.stan for the behind the scene HMC algorithm (no edits necessary)
•	R file *_main_simulation.R for the initial analysis using all variables provided
•	R file analyze_*_simulation.R for the variable selection procedure using the algorithm sequential 2-means (Li and Debdeep (2017))
•	R file *_main_simulation_SelectedVariables.R for the secondary analysis using only the variables selected in the previous step
•	R file analyze_*_simulation_SelectedVariables.R for the final analysis for producing the simulation results in the paper.
### License
All the codes are free to use with a proper citation of the paper.

