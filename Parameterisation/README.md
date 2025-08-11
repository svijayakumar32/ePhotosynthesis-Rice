This folder contains scripts used to scale enzyme inputs and incorporate additional parameters into the e-Photosynthesis model.
The scripts in the sub-directories must be run in the following order:

- `Rubisco Parameters` - calculates species-specific Rubisco parameters prior to model optimisation.
- `Input Scaling` - rescales input enzyme activities to the model must be rescaled to match the photosynthetic response of the species being modelled.
- `Model Comparison` - minimises differences between FvCB and e-Photosynthesis model assimilation rates across Rubisco- and RuBP regeneration- limiting ranges of CO2 concentrations.
