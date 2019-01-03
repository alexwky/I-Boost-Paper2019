# I-Boost-Paper2018

This directory contains the codes to perform all the analyses and reproduce all the figures presented in the paper: Wong KY, Fan C, Tanioka M, Parker JS, Nobel AB, Zeng D, Lin DY, Perou CM. I-Boost: an integrative boosting approach for predicting survival time with multiple genomics platforms. 2019.

The R-package **IBoost** can be found at [https://github.com/alexwky/I-Boost](https://github.com/alexwky/I-Boost).

Before running the codes, download the processed TCGA pan-cancer data set `TCGA_8cancer_rmmis.csv` from [https://doi.org/10.5281/zenodo.2530387] (https://doi.org/10.5281/zenodo.2530387) and store it in the directory `Data/`.

## Simulation Data Sets

The simulation data sets used in the paper are stored in the directory `SimulationData`. The zip files in the directory contain the file `allPredictorValues.csv` and 3,000 files with names in the form of `Data-setting[setting number]-[replication number].csv` that consist of rows in the form:
```
0.240290380066643,0,276
0.578869428061267,1,531
2.14527656943497,1,837
0.272058507083588,1,1241
.
.
.
```
Each of the 3,000 files contains the simulated data for 500 subjects for a specific setting and simulation replicate. Each row in a file contains data for a subject. The first element on each row is the survival or censoring time, the second element is the event indicator that equals 1 or 0 if the event is observed or right-censored, respectively, and the third element is the index of the row in the file `allPredictorValues.csv` that contains the values of the clinical/genomic predictors of this subject. For instance, a value of `1` for the third element represents that the values of the clinical/genomic predictors of this subject are stored in the first row of `allPredictorValues.csv`.

## Simulation Studies

Run the bash code `simulations.sh` in the home directory. It creates the directories `SimulationData`, `SimulationResults`, `SimulationSettings`, and R codes in the `Programs` directory. Run the program `sim_setting.R` in the `Programs` directory; it will generate files that contain parameter values for simulating data. Run the program `gendata.R` in the `Programs` directory; it will generate 1,000 simulated data sets for each simulation setting in the directory `SimulationData` (and the file `allPredictorValues.csv` described above). (This step is not necessary for the simulation studies, because the (same) data sets are generated internally in the simulation programs.)

Run the R programs `sim-[method name]-s[setting].R` in the `Programs` directory. Each program performs analysis on 1,000 simulated data sets and output the results in the directory `SimulationResults`. Each row of the output files represent:
```
replication number, number of variables selected, number of true signal variables selected, risk correlation, number of clinical variables selected, number of gene modules selected,  number of protein expressions selected, number of miRNA expressions selected, number of mutations selected, number of copy number variations selected, MSE for clinical variables, MSE for gene modules,  MSE for protein expressions, MSE for miRNA expressions, MSE for mutations, MSE for copy number variations
```
Some of the programs may take very long time to run. To save time, one may run the analysis for separate replicates in separate programs and combine the results in output files in the format given above.

## Analysis of TCGA data

Run the bash codes `analysis_splits.sh`, `analysis_splits_cox.sh`, and `analysis_wholedata.sh` in the home directory. They will create all programs necessary to perform all analyses on the TCGA data sets presented in the paper. Run the created programs in the `Programs` directory. The programs with prefix `DataAnalysis-` in the name will perform analyses over 30 training/testing data splits using maximum likelihood estimation, LASSO, elastic net, I-Boost-CV, or I-Boost-Permutation. The results will be written in files in corresponding directories in `Results`. Files with prefix `Model_` in the name contain the selected predictors and the estimated regression parameters, and files with prefix `summary_` in the name contain summaries of the analysis results, including the C-index in the testing set. Upon completion of all analyses, run the program `gather.R` in the `Programs` directory, which will combine the `summary_` files.

The programs with prefix `WholeDataAnalysis-` in the name perform the analyses on the whole LUAD, KIRC, or pan-cancer data sets (without setting aside subjects for testing). The selected predictors and estimated regression parameters are written in corresponding directories in `Results/WholeData`.

## Generation of Figures

Upon completion of all analyses, run the programs `plotFigure1.R`, `plotFigure2and3andS1.R`, `plotFigure4.R`, `plotFigure5.R`, and `plotFigure6.R` in `Programs/`. They will generate Figures 1-6 and Figure S1 in the directory `Plots`.
