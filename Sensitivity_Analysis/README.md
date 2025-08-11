The sensitivity analysis can be run using a single script - as in the case of the optimisation routine, the use of high performance computing is recommended.

## Analysis
- `enzyme_adjustment_test_new_2000.m` runs the sensitivity analysis to determine assimilation rates associated for random fold changes (n = 2000) of various enzyme combinations drawn between FC = 1 (no change) and FC = 1.25 for Rubisco and FC = 2.8 for SBPase, aldolase, PRK, FBPase and TK at Cc = 130 umol mol-1, 250 umol mol-1 or Cc = 360 umol mol-1
- `job_enzyme_adjustment_new.sh` is a job script for submitting `enzyme_adjustment_test_new_2000.m` from the login node to the SLURM scheduling system.  
