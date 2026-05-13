# R scripts

Replication scripts for "Threshold Without Trigger". Run in numbered order. 
See the main README in the repository root for installation and data access.

## Scripts

- `01_dag.R` - Generates the directed acyclic graph for the Methods chapter. 
  Uses the `dagitty` package to specify the causal structure and `ggplot2` 
  with `ggforce` to render the figure. Independent of the cleaned data and 
  can be run at any time.

- `02_clean_data.R` - Loads the DCID dataset from `../data/`, filters to 
  state-initiated incidents between 2000 and 2020, constructs the NATO-target 
  and non-NATO comparison series, and aggregates to monthly, quarterly, 
  half-yearly, and yearly counts. Writes the derived series to 
  `../data/derived_series.rds`.

- `03_modeltest.R` - Specification diagnostics on the monthly series. 
  Compares Poisson against negative binomial, zero-inflated alternatives, 
  spline degrees of freedom, linear against spline trends, post-cut-off 
  trend breaks, residual autocorrelation, temporal granularity, and prior 
  sensitivity for the quarterly Bayesian model. Tables print to the console 
  in LaTeX form; figures save to `../Figures/`.

- `04_regression.R` - Bayesian negative binomial regressions for the main 
  analysis. Estimates the baseline model, short-exposure-window models 
  (0 to 3, 6, and 12 months), a dynamic-exposure-window model, comparative 
  regression models (NATO target versus non-NATO target), and appendix 
  models for extended short windows and pre-recognition windows.

- `05_descriptive.R` - Compositional figures and tables for the Results 
  chapter. Produces pre/post bar charts for method and severity, half-yearly 
  composition time series for target type, coercive objective, damage type, 
  and critical infrastructure sector, and comparative pre/post tables.

## Output

Figures save to a `Figures/` directory in the repository root, which is created automatically when the scripts run. Tables print to the console as LaTeX code.
