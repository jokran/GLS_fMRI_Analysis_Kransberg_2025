# Analysis Script for Kransberg et al. (2025): Grid-Like Signal fMRI Study

## Repository Contents

This repository contains the R script (`GLS_analyses.R`) used for the primary statistical analyses reported in the paper:

* **[Failure to Detect Entorhinal Grid-Like Signals in a Pre-registered Human fMRI Study]**
* Kransberg, J., Bråthen, A. C. S., Falch, E. S., Øverbye, K. E. Ø., Garrido, P. F., Fjell, A. M., Stangl, M., Wolbers, T., Sneve, M. H., & Walhovd, K. B. (Expected Year/Journal Info).
* [Optional: Add Link to Pre-print or Published Paper Here when available]

**Contact:** Jonas Kransberg - Jonas.kransberg@psykologi.uio.no 

## Script Description (`GLS_analyses.R`)

This script performs the analyses investigating grid-like signals (GLS) in human entorhinal cortex fMRI data, including:
* Tests for 6-fold, 5-fold, and 7-fold GLS magnitude against zero.
* Analyses stratified by age group.
* Linear mixed-effects models examining effects of age, sex, and hemisphere.
* Temporal and spatial stability analyses.
* Control analyses for spatial smoothing, segmentation method, high task performance, and potential 1-fold signal interference.

## IMPORTANT: Data Note

* **This script uses simulated dummy data.** The structure mirrors the actual data used in the study, but the values are randomly generated for demonstration purposes.
* The original dataset used in the study has restricted access due to privacy regulations concerning sensitive participant information.
* Requests for access to the original data can be directed to the corresponding author as detailed in the manuscript, subject to ethical and data protection approvals.

## Requirements

* **R:** Developed using R version [Specify R version, e.g., 4.x.x]. Compatibility with other versions is likely but not guaranteed.
* **R Packages:**
    * `tidyverse` (for data manipulation and plotting utilities)
    * `lmerTest` (for linear mixed-effects models)

You can install these packages in R using:
```R
install.packages(c("tidyverse", "lmerTest"))