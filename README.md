FT-RSGT: Data Analysis
=========
The analyses associated with the TMS-EEG experiment are available in `tms_analyses`. The folder should not be renamed.

## Preprocessing

`data_preprocessing.R` preprocesses raw .csv files: it computes approximate entropy (AE) and behavioral variability (BV) for each probe trial based on the preceding 25 trials. It also changes some variable types and saves the preprocessed dataset as an `.Rdata` file. The script also returns a long data frame and outputs some graphs.

## Modeling functions
Before proceeding to model fitting, make sure the necessary functions are loaded by starting `model_fitting_functions.R`.

## Model Fitting
The script does model fitting on preprocessed data and plots parameter coefficients.

## Results

Hypothesis testing is done in markdown for better visualization. It includes brms non-linear hypothesis testing, Kruskal-Wallis tests and post-hoc Wilcoxon multiple comparisons.