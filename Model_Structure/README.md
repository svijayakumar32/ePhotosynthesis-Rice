This repository contains scripts associated with the structure of the e-Photosynthesis model:

- `cdn.m` defines the environmental conditions, such as CO2 concentration and photon flux density.
- `SYSInitial.m` defines the simulation time.

For the model files, the first part of the filename describes the portion of the e-Photosynthesis model being simulated:
- `CM_` describes carbon metabolism (CBB cycle).
- `PS_` describes the CBB cycle in addition to the starch synthesis and triose phosphate export. 
- `PR_` describes the photorespiration model. 
- `PS_PR_` describes the CBB cycle, starch synthesis, triose phosphate export, and photorespiration process. 
- `FI_` describes the light energy absorption, transfer, primary charge separation, and electron transfer around PSII. 
- `BF_` describes the electron transfer from reduced plastiquinone until the generation of ATP and NADPH, including the ion transfer through thylakoid membrane and ATP synthesis process. 
- `FIBF_` describes reactions covered by FI_Drive, and BF_Drive. 
- `RuACT_` describes reactions of Rubisco activation process. 
- `XanCycle_` describes reaction of the xanthophyll cycle. 
- `EPS_` describes reactions covered by PS_PRDrive and FIBF_Drive.
- `RA_` describes reactions covered by EPS_Drive and RA_Drive. 
- `RedoxReg_` describes reactions covered by EPS_Drive and the redox regulation of enzyme activities. 
- `DynaPS_` describes reactions covered by EPS_Drive, RuACT_Drive and XanCycle_Drive.
- `tr DynaPS_` describes reactions covered by EPS_Drive, RuACT_Drive, XanCycle_Drive and RROEA_Drive.

The second part of these filenames use suffixes to describe the purpose of the script:
- `_AddTitle` defines titles for concentration/time plots plots.
- `_Ini` defines initial values and parameters.
- `_Drive` defines model simulation settings.
- `_Graph` defines plotting commands.
- `_Rate` defines rate equations.
- `_mb` defines differential equations.

Additional Scripts
- `IModelCom.m`  is called by CM_Drive to initialise the structure of the model using different components of the full photosynthesis model
- `IniModelCom.m` is called by EPS_Drive to to initialise the structure of the model using different components of the full photosynthesis model
