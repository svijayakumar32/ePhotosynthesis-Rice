## Structure

This repository contains code for fitting CO₂ response data using the `msuRACiFit` package.

The analysis is fully scripted within a single .R file for portability and reproducibility.

All core logic is contained in `msuRACiFit_params.R`, located in the `Scripts/` folder.

The `Data/' subdirectory contains the gas exchange data files used in the analysis.

- Raw data:
`Gas_exchange_measurement_WT_plants.csv`

- Cleaned input files:
These are named using the pattern:
`IR64-A009-07-33-05-0x_Wildtypex.csv`

## Requirements
Running the code requires installation of the following R packages:
- `devtools` - for installing other packages.
- `msuRACiFit` - for fitting photosynthetic CO2 response curves to assimilation. 
- `here` - for constructing paths to project files.
- `readxl` - to import Excel files into R.

## Running the Analysis

To run the full analysis:

1. Open a terminal or command prompt.

2. Navigate to the repository root directory.

3. Run the `msuRACiFit_params.R` script. 

     (a) The script installs the `rChoiceDialogs` package to interactively set the working directory (`CO2_response_fitting`) using `setwd(rchoose.dir(default = getwd(), caption = "Select Directory"))` to match your local path to the repository.
   
     (b) After setting the location of the new `CO2_response_fitting directory`, you can also source the script by running:
     `source("Scripts/msuRACiFit_params.R")`

This allows the script to run from any environment within the directory.
