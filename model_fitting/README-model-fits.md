# Model-fitting for GhV

This is how the models were fit simultaneously to three different datasets for GhV: age-seroprevalence data, seasonal serological data across a single year, organized by biweek based on the birth-pulse for the site in question, and seasonal viral shedding data organized by biweek based on the birth pulse for the site in question.


1. First, we used the file "matrix-pars.csv" to specify initial conditions for the 20 models tested. These were 5 different model structures, each tested under 4 different conditions: the base model, the base model + elevated mortality on the infectious class ("mu-sick"), the base model + elevated mortality on the recovered class ("mu-rec"), and the base model + elevated mortality on both the infectious and recovered classes ("mu-sick-mu-rec"). 

2. We then run the script "parse-parameter-sets-matrix.R" to resample the fitted parameters for all 20 models and produce 100 distinct parameter sets for each model that get saved as an individual csv, half to the folder "matrix-fit-1" and half to the folder "matrix-fit-2".

3. In addition to the csv files of parameter sets, the "matrix-fit-1" and "matrix-fit-2" folders also contain all three datasets, the R script for model fitting (fit-models.R) and the SLURM script used to run these array jobs on Savio.

4. That "matrix-fit-1" or "matrix-fit-2" folder then gets uploaded to Savio, and we run an array job to fit each model to the three datasets, using each initial parameter set.


## pause

5. The output from that fitting (.Rdata files with the fitted parameters, likelihood scores, and AIC for each model) are saved in the sub-folder "out-1". 

6. From there, we run the 'process-fits.R' script to consolidate and parse the initial parameter fits and plot them. A series of output plots are saved in the subfolder "out-1-plots". 

7. We then run the script 'parse-parameter-sets-refit.R' to parse those initial parameter fits, decide on the top 20 models (5 best base models + all four iterations of each) and generate new starting values (62 unique parameter sets per all 20 models - or 248 unique parameter sets for all 5 base models) to fit again. 

8. These then get added as individual .csv files to the 'refit-1' subdirectory, along with the fitting script ('fit-models.R'), the three data files, and a new SLURM script. 

9. When these finish, we download the resulting model fits into "refit-out", process and plot them with the script 'process-refits.R' with plots stored in the 'refit-out-plots' folder. The top 50 overall model fits are saved in the files 'top-refits-AIC.csv' and 'top-refits-llik.csv'. The best fit parameter combinations for each of the 20 model strucutres are listed in 'top_refits_per_modeltype.csv'. All of these top fit files are stored in the broader 'GhV-fits' subfolder.

10. In general, the best fitting models did not produce reliable confidence intervals around each parameter fit, so we will need to add these later. NEXT STEP!

