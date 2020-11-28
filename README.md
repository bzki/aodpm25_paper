# Project for gap-filling AOD to improve PM25 prediction models

I provide a description here of the main files and folders, and the process to run the analyses in question. 

This project was run on an Ubuntu 18.04 machine with 16GB RAM running R version 3.6. Certain files were run on a cluster, also using R 3.6. The key packages I used are `LatticeKrig` (version 8.4) and `ranger` (version 0.12.1).

No data is presently provided in this repository, as Github is not meant to distribute data. I provide a link to the relevant data for reproducing analyses [here](http://www.mediafire.com/file/hx9c8tjobk7vtzj/aodpm25_data.zip/file). This was prepared by Xuefei Hu; details are provided in the paper as well as Hu et al. 2017. The sole purpose of the provided data is to replicate the analyses in the paper. 

Unzip the downloaded `.zip` file, which consists of another `.zip` file and other `.rds` files. The `n_aod_pred_2011_new8_.zip` compressed directory consists of `.rds` files that should be placed in the `rawrds/` subdirectory of this project before running the scripts. So, for example, we would have `rawrds/n_aod_pred_2011_new8_182.rds`. The other file in the downloaded `.zip` file is `region_assignment.rds`, which should be placed in the `cmaq/` subdirectory. In addition, I use a U.S. 2010 TIGER shapefile for U.S. states taken from IPUMS NHGIS -- this is used for creating state overlays on the maps, and should be downloaded directly from [IPUMS NHGIS](https://www.nhgis.org/) and unzipped/renamed so that we have the following folder: `cmaq/shapefile_tl2010_us_state_2010/` which contains the shapefile and related files. The scripts produce `holdout_4_circleB.rds`, but if desired, this should be placed in the `aod_holdouts/` folder. 

In the directory `functions/`, I include some helper functions that will be sourced for some of the analysis scripts. In the directory `tests/`, I include a few unit tests to ensure that my functions behave as expected with some simplified example datasets.

In the main directory, there are three main scripts. These R files process and prepare the main data:

- `00_csv_to_rds.R`: This reads in the raw .csv files and outputs .rds files provided in the downloaded link (so you would skip this script, since I provide the `.rds` files).
- `01_fill_nldas.R`: *(Start here)*. This fills in 3 missing values for NLDAS variables using nearest neighbor and nearest 4 neighbors. Missing values, infinite, and extreme values are also replaced for NARR variables, but these are not used in any of the analyses.  
- `02_aod_observed`: This creates subsets of the data, keeping the observations where AOD is observed and omitting the others. 
- `bash_commands.sh`: This contains shell script code for running the above files. The scripts throughout this project can be run in this fashion: `R CMD BATCH --vanilla filetorun.R Rout/filetorun.Rout`, where I designate `Rout/` as a folder that contains the logs. 

The `aod_holdouts/` directory contains the files pertaining to AOD gap-filling. A spread sheet `aod_holdouts/code/aod_descriptions.ods` contains the file names and descriptions. Briefly, you may find in `aod_holdouts/code/`:

- `create_holdouts.R` creates the testing data for AOD daily experiments. This is done several different ways, but only the last method (using large circle areas as holdouts) is used in the analyses. `make_split.R` creates easy to access splits of the training/testing/validation (for the most part, I refer to validation + training as training.). 
- `make_splitcv.R` and `make_splitcv_full.R` create cross-validation splits in the training AOD data and the observed AOD data, respectively, using the `blockCV` package. `cv_blockmaps.R` creates a nice daily visualization of this. `overall_holdout_map.R` creates additionals maps. 
- `cv_lk.R`, `cv_rf.R`, `cv_lk_full.R`, `cv_rf_full.R`, are fitting models to the 10-fold cross-validation made by the previous scripts.
- `test_lk.R` and `test_rf.R` train on the training data and predict on the test data. 
- `fit_lk_full.R` and `fit_ranger_full.R` train to the full observed data and make predictions where not observed.
- `combine_predictions_cv.R` and `combine_predictions_cv_full.R` stack the cross-validation predictions on the training and observed AOD data, respectively, and determine how they ought to be combined using the methods mentioned in the paper.
- `assess_results.R` looks at performance in the training/testing split. 
- `full_maps.R` creates full maps of predictions and other maps used in the paper. 


The `pm25/` directory contains the files pertaining to PM2.5 prediction. A spread sheet `pm25/code/pm25_descriptions.ods` contains the relevant file names and descriptions. Briefly you may find in the `pm25/code/` folder:

- `make_cvfold.R` creates the cross-validation folds for PM2.5 analyses (again, using `blockCV` for the spatially clustered CVs). `plot_cvfold.R` produces accompanying visualizations.
- `cv_rf_daily.R` and `cv_rf_daily_m5.R` fits the daily random forest models M1 through M5 to the cross-validation folds. 
- `cv_rf_st.R` and `cv_rf_st_m5.R` similarly fits the spatio-temporal random forest models M1 through M5 to the cross-validation folds.
- `daily_results.R`, `st_results.R`, and `combined_tables.R` summarizes the predictions. 
- `full_rf.R`, `full_rf_m5.R` fit the full spatio-temporal models to all of the observed data and generate predictions where not observed. `maps.R` creates full maps and tables based on these results.


