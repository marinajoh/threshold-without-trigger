# Data

This thesis uses the February 2023 update of the Dyadic Cyber Incident and
Campaign Dataset (DCID) version 2.0.

## Citation

Maness, Ryan C., Brandon Valeriano, Kathryn Hedgecock, Benjamin M. Jensen, and
Jose M. Macias. 2022. *The Dyadic Cyber Incident and Campaign Dataset, version 2.0*.
Available at: https://drryanmaness.wixsite.com/cyberconflict/cyber-conflict-dataset

## Access

The dataset file is not included in this repository and must be downloaded
separately from the DCID project website.

Direct link: https://drryanmaness.wixsite.com/cyberconflict/cyber-conflict-dataset

The file used in this analysis is `DCID_2.0_Release_update_February_2023.xlsx`.
The R scripts read the first sheet of the workbook, which contains the standard
DCID 2.0 release used in this analysis.

An earlier version of DCID 2.0 (September 2022) is also archived on Harvard
Dataverse under a Creative Commons CC0 1.0 Universal Public Domain Dedication:
https://doi.org/10.7910/DVN/CQOMYV. This thesis uses the February 2023 update
rather than the September 2022 Dataverse release.

## Replication instructions

1. Download `DCID_2.0_Release_update_February_2023.xlsx` from the DCID project
   website above.
2. Place the file in this `data/` directory.
3. Run the R scripts in `R/` in numbered order, starting with `01_clean_data.R`.
