[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14279450.svg)](https://doi.org/10.5281/zenodo.14279450)

This repository contains data, code, and statistical output for **Goodman, Reum, Barnes, Punt, Ianelli, McHuron, De Leo, & Holsman (2024), "Climate Covariate Choice and Uncertainty in Projecting Species Range Shifts: A Case Study in the Eastern Bering Sea"**, in the journal *Fish and Fisheries*.

**Directory structure:**

:file_folder: **R**: Functions used in this project

:file_folder: **analysis**: R scripts for running project analyses

:file_folder: **data**: All ROMS-NPZ, trawl survey, and other datasets used in this project

:file_folder: **inst**: Files needed for functions to work, copied to package directory

**Use:**

Not all data needed to run this repository are publicly available at the time of publication. Specifically:

1. Estimates derived from the Bering10K ROMS-NPZ model on the standard survey grid ([K20P19_CMIP6.zip](https://drive.google.com/file/d/1qth5mKYP_voKaskRZyAndGN0jZMFsdMo/view?usp=share_link) and [K20P19_CMIP5.zip](https://drive.google.com/drive/folders/1t_JqDBQU-Fyy5nvIYRAmVcqzWi4mq7mk?usp=share_link)). The `BC_ACLIMsurveyrep` folder from the CMIP6 zip file must be copied into in the `data` directory, and the contents from the CMIP5 `BC_ACLIMsurveyrep` folder copied into it. Message the lead author for access if needed.
2. The EBS size-binned weight CPUE [data files](https://drive.google.com/drive/folders/189RnXV_PZUpKei1BfKFkZBcyvM6X4P8r) must be placed in `data/trawl_surveys_size_binned`. Message the lead author for access if needed.

All analyses can then be run by sourcing `make.R`. Downloading and bias-correcting the ROMSNPZ output, running the SDMs, and producing range maps and other derived quantitates may take up to two days depending on your computing resources. Because this project make's use of the `rstudioapi` package to run models for each species as RStudio jobs, the makefile must be sourced from RStudio.
