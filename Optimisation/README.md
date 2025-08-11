This folder contains scripts associated with the optimisation of the e-Photosynthesis model, described briefly in `Optimisation_Procedure.txt`.

## Inputs
- `Einput7.txt` - a list of the original enzyme activities for e-Photosynthesis model
- `MW&Kcat.txt` - a list of molecular weights and kcat values for photosynthetic enzymes

## Evolutionary Algorithm
The evolutionary algorithm is run using the following scripts:
- `average.m` calculates the average Vmax in each generation
- `mutate.m` generates variation in Vmax populations by modifying mutatePercentage
- `rankPop.m` ranks populations based on their CO2 uptake values
- `resizePop.m` resizes the population after each generation
- `stdev.m` calculates the standard deviation of Vmax in each generation
- `twoPoint.m` is a variation of `resizePop.m`, but it is not used in the current simulations

## High Performance Computing
These files are called by scripts in the `gpmain` folder, which optimise the model at a specific CO2 concentration (Cc = 130-380 umol mol^-1 at intervals of 10 umol mol^-1).
- Simulations can be run by submitting single scripts or arrays (to run the same program multiple times). For each Cc input, an array of 15 batch jobs is submitted in parallel. While only 10 successful replicates are required for downstream analysis, the additional 5 jobs are included as a buffer. This overprovisioning accounts for potential failures in optimisation due to convergence issues or other runtime errors that may occur in a subset of jobs.

- Running `job_array_gpmain_rice_129.sh` creates the photorespiratory constraints for all other simulations, so it must be run *first*.
- `extract_PR_constraints.m` imports the results of optimisation at Cc = 129 umol mol^-1, then averages and exports PR enzyme Vmax values to serve as constraints for full optimisations

As the optimisation process is computationally intensive, it is recommended that users run the simulations using a high performance computing (HPC) cluster, allotting sufficient memory of 8GB per job.

For the purposes of this work, all optimisations were run on the Lancaster University High End Computing (HEC) cluster using basic Linux commands.
- For further documentation on the Lancaster University HEC, please refer to https://lancaster-hec.readthedocs.io/en/latest/index.html. 
- For general documentation about SLURM, please see: https://slurm.schedmd.com/.

The MATLAB scripts in gpmain (`gpmain_rice_xxx_new.m`) require corresponding job scripts (`job_gpmain_rice_xxx.sh`.sh) that are submitted from the login node to the SLURM scheduling system.
In our analysis, the batch job scripts are run as arrays - enabling the same program to be run multiple times.

Each job script is a text file containing the commands to be run along with the compute resources required to run them. 
* Prior to running a script, it may be necessary to convert line endings of the text files from DOS (Windows) to Unix, e.g.

`dos2unix job_array_gpmain_rice_xxx.sh`

Batch jobs are run on the HEC by creating a batch job script and submitting it to the system using the command `sbatch`. 
For example, a job script named `job_array_gpmain_rice_xxx.sh` can be submitted with the command:

`sbatch job_array_gpmain_rice_xxx.sh`

After completing the optimisations, `CalculateGrossAssimilation.m` can be used to calculate gross CO2 assimilation as a result of using enzyme stacking strategies (Stress/Current/Future, termed here as Low/Ambient/Elevated) by combining non-optimised and optimised Vmax data of targeted enzymes at Cc = 130, 250 and 360 umol mol^-1.
