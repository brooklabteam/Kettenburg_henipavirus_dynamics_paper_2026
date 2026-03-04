# Model-fitting for GhV and AngV

This is how the models were fit to the two different datasets for GhV and AngV.
For both pathogens, models weere fit simultaneously to three different datasets: age-seroprevalence data, seasonal serological data across a single year, organized by biweek based on the birth-pulse for the site in question, and seasonal viral shedding data organized by biweek based on the birth pulse for the site in question.

## GhV fits. 
1. First, we used the file "fit-pars-init.csv" to specify initial conditions for 52 models tested. These were 13 different model structures, each tested under 4 different conditions: the base model, the base model + elevated mortality on the infectious class ("mu-sick"), the base model + elevated mortality on the recovered class ("mu-rec"), and the base model + elevated mortality on both the infectious and recovered classes ("mu-sick-mu-rec"). 
2. We then run the script "parse-parameter-sets-1000-total.R" to resample the fitted parameters for all 52 models and produce 19 distinct parameter sets for each model that get saved as an individual csv in the folder "fit-1".
3. In addition to the csv files of parameter sets, the "fit-1" folder also contains all three datasets, the R script for model fitting (fit_models.R) and the SLURM script used to run this array job on Savio.
3. That "fit-1" folder then gets uploaded to Savio, and we run an array job to fit each model to the three datasets, using each initial parameter set.
4. The output from that fitting (.Rdata files with the fitted parameters, likelihood scores, and AIC for each model) are saved in the sub-folder "out-1".
5. From there, we run the 'process-fits.R' script to consolidate and parse the initial parameter fits and plot them. A series of output plots are saved in the subfolder "out-1-plots".
6. We then run the script 'parse-parameter-sets-refit.R' to parse those initial parameter fits, decide on the top 20 models (5 best base models + all four iterations of each) and generate new starting values (XXXX per model) to fit again.










