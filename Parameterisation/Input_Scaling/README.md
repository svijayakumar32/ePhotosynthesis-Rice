The original enzyme activities in `Einput7.txt` are not species-specific and must therefore be scaled using photosynthetic parameters obtained from `CO2_response_fitting`.
These parameters are then integrated into the following scripts to get optimal α scaling factors which are multiplied by the enzyme Vmax list in `Einput7.txt`.
To obtain the adjusted `Einput` values, the Vmax of Rubisco (`Einput(1)`) must be multiplied by αRubisco, whereas the Vmax of all remaining enzymes must be multiplied by αEnzymes.
This calculation is included in `gpmain` and other scripts so does not need to be added explicitly.

- `Jmax_adj_simple.m` runs model calibration to find optimal αEnzymes by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the RuBP-regeneration limited range of [CO2]

- `Vcmax_adj_simple.m` - runs model calibration to find optimal αRubisco by minimising SSR between assimilation rates of FvCB and e-Photosynthesis models in the Rubisco limited range of [CO2]

- Through running `Vcmax_adj_simple.m`, `WeatherTemp.mat` is generated for use in downstream analysis.
