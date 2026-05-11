# Threshold Without Trigger

Replication code for the 2026 MA thesis "Threshold Without Trigger: On NATO's 
2014 Recognition of Cyber Attacks Under Article 5", Department of Political 
Science, University of Oslo.

## Author

Marina Johansen  
GitHub: [@marinajoh](https://github.com/marinajoh)

## Overview

This thesis examines whether NATO's 5 September 2014 recognition that cyber 
attacks may lead to the invocation of Article 5 is associated with changes in 
the frequency or composition of observed cyber incidents directed at NATO 
member states. The analysis combines a Bayesian negative binomial regression of 
incident counts with a descriptive analysis of compositional change across six 
dimensions, using the Dyadic Cyber Incident and Campaign Dataset (DCID) version 
2.0 for the period 2000 to 2020.

## Repository structure

```text
threshold-without-trigger/
├── R/                  R scripts for cleaning, modelling, and figures
├── data/               Data access instructions (data not redistributed)
├── README.md           This file
├── LICENSE             MIT License for the code in this repository
└── .gitignore          Standard R gitignore
```

## Requirements

The analysis was run on R version 4.5.1 (2025-06-13). The following R packages 
are required:

- Data handling: `readxl`, `dplyr`, `tidyr`, `tibble`, `janitor`, `lubridate`
- Modelling: `rstanarm`, `MASS`, `pscl`, `splines`, `rethinking`
- Tables and output: `knitr`, `kableExtra`, `modelsummary`
- Figures and DAG: `ggplot2`, `scales`, `ggforce`, `dagitty`

Install the required CRAN packages with:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "tibble", "janitor", "lubridate",
  "rstanarm", "MASS", "pscl",
  "knitr", "kableExtra", "modelsummary",
  "ggplot2", "scales", "ggforce", "dagitty"
))
```

The `rethinking` package is installed separately. See installation instructions 
at https://github.com/rmcelreath/rethinking.

## Data

This analysis uses the February 2023 update of the Dyadic Cyber Incident and 
Campaign Dataset (DCID) version 2.0. See `data/README.md` for download 
instructions. The dataset is not included in this repository.

## Usage

Run the R scripts in `R/` in numbered order:

1. `01_clean_data.R` - Loads and cleans the DCID dataset
2. `02_dag.R` - Generates the directed acyclic graph for the analysis
3. `03_modeltest.R` - Model comparison and prior diagnostics
4. `04_regression.R` - Bayesian negative binomial regression analyses
5. `05_descriptive.R` - Descriptive compositional analyses

## Citation

If you use this code, please cite the thesis:

Johansen, M. (2026). *Threshold without trigger: On NATO's 2014 recognition of cyber attacks under Article 5* [Master's thesis, Department of Political Science, University of Oslo].

## License

The code in this repository is licensed under the MIT License (see `LICENSE`). 
The DCID dataset is published separately by its authors and is not covered by 
this license.
