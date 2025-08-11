### e-Photosynthesis-rice
This repository contains the code and data necessary to reproduce the results presented in the paper by Vijayakumar, S., Wang, Y.,Lin, H.C., Carmo-Silva, E., Long, S.P. and Taylor, S. H. entitled "Tailoring a dynamic model of photosynthetic metabolism towards greater carbon assimilation in rice".

## Overview
This repository provides a generalisable workflow for calibration, parameter optimisation, and input scaling that can be applied to any plant species with available gas exchange and Rubisco kinetic and activity data.
Our paper describes the steps for optimising resource investment among photosynthetic enzymes to maximise carbon assimilation in a model of C3 photosynthetic metabolism calibrated to rice, summarised in `pipeline_simple_flowchart.png` (see below). 

This version of e-Photosynthesis integrates species-specific temperature dependencies for the balance of oxygenation and carboxylation (Vomax/Vcmax) and ribulose-1,6-bisphosphate carboxylase/oxygenase (Rubisco) activation, as well as accounting for tight binding Rubisco inhibitors.
The workflow begins by combining species-specific equations for calculating temperature dependencies of Rubisco catalytic properties with leaf-level gas exchange measurements for Oryza sativa cv. IR64 to derive photosynthetic parameters describing Calvin-Benson-Bassham (CBB) cycle activity, i.e. Vcmax and J. 
These photosynthetic parameters are used to re-scale enzyme activities in e-Photosynthesis before running the model to identify redistributions that are optimal for CO2 assimilation at different [CO2] levels. Following this, additional analyses are applied to evaluate strategies for increasing photosynthesis through overexpressing subsets of 2-6 enzymes to engineer improved photosynthesis under drought stress, current ambient [CO2] and future elevated [CO2].

## Structure
The files within the main repository contains scripts, variables and data adapted from the version of e-Photosynthesis stored in https://github.com/cropsinsilico/C3-metabolic-and-leaf-model.
All code licensing information is detailed in `LICENSE`.

Files are stored in various directories defining their function:

1) CO_2 response - for fitting leaf-level gas exchange measurements to the Farquhar-vonCaemmerer-Berry (FvCB) model

The Data sub-folder contains:
   - `Gas_exchange_measurement_WT_plants.csv` - A `.csv` file containing 8 sets of leaf level gas exchange measurements for wild type Oryza sativa cv. IR64
   - Filenames beginning with `IR64-A009-07-33-05` contain subsets of data for each curve (Photo, Ci, Tleaf, Press)
 
The Scripts sub-folder contains:
   - `msuRACiFit_params.R` - R script containing fitting procedure for obtaining photosynthetic parameters

2) Parameterisation - describes modification of model inputs and parameters prior to model optimisation.
   
   - Rubisco Parameters:
      - Two `.R` scripts (`Temp-resp-Rubisco.R` and `Vmax-temp.R`) determine the temperature dependences of Vomax/Vcmax and Rubisco activation state. 
     The values/formulas obtained are then used to modify `CM_Rate.m`.

   - Input Scaling
     - Prior to running the scripts, the MATLAB working directory can be set interactively using `selpath = uigetdir()` 
and subsequently added to the MATLAB path using `addpath(genpath(selpath))`. This only needs to be completed *once*.

   - Model Comparison
     - Running `Farq_ePhoto_comparison_new.m` saves assimilation rates under Rubisco, RuBP regeneration and TPU limitation states for the FvCB and e-Photosynthesis models, enabling the comparison of untuned e-and tuned e-Photosynthesis model (with scaling of Vmax).

3) Optimisation - these are the files used to run the optimisation of e-Photosynthesis

- Shell_Scripts - batch job scripts (.sh) for running `gpmain` simulations on the High-End Computing (HEC) cluster: https://lancaster-hec.readthedocs.io/
  - simulations can be run by submitting single scripts or arrays (to run the same program multiple times)
  - running `job_array_gpmain_rice_129.sh` creates the photorespiratory constraints for all other simulations

- gpmain scripts - called by batch job scripts to run model optimisation simulations at Cc = 129 umol mol-1, then at Cc = 130-380 umol mol-1

- `CalculateGrossAssimilation` creates stacked enzyme strategies using average optimized Vmax (1-10) modelled under Cc = 130, 250 and 360 and calculates gross assimilation 

4) Sensitivity Analysis
   - a single script used to evaluate the effect of variability in protein expression on predicted improvements to CO2 assimilation

## Running the Analysis
The code is run in the following order, where the final script of each analysis plots the results:

1. CO2 response fitting 
   - `msuRACiFit_params.R`
      
2. Parameterisation
   - (a) Rubisco_Parameters
      - `Temp-resp-Rubisco.R` -> `Vmax-temp.R` - > modify `CM_Rate.m` as required
   - (b) Input_Scaling
      - `Jmax_adj_simple.m` -> `Vcmax_adj_simple.m`
   - (c) Model_Comparison
      - `Farq_ePhoto_comparison_new.m` 
   
3) Optimisation
   - `job_gpmain_rice_129_new.m` -> `job_gpmain_rice_130...380_new.m` -> `CalculateGrossAssimilation.m` 
   
5) Sensitivity Analysis
   - `enzyme_adjustment_test_2000.m` 
