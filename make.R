
if (!require("aclim2sdms")) {
  devtools::install()
  require("aclim2sdms")
}

pkgs <- c("here", "dplyr", "tidyr", "purrr", "ggplot2", "sf", "stars", "ncmeta", "mgcv", "aclim2sdms", "foreach", "doParallel")

sapply(pkgs, require, character.only = TRUE)

# Download & bias-correct ROMS level 2 data ---------------------------------------------

## This takes up to a couple days - do not re-run if it can be avoided

## Download and bias-correct ROMS outputs
source(here("analysis", "get_roms_level2.R"))

## Extract ROMS outputs for each summer, stack variables for prediction
source(here("analysis", "join_roms_level2.R"))

## Extract covariates corresponding to survey locations and dates 1982-2022
source(here("analysis", "hindcast_extract.R"))

# Parameters ----------------------------------------------------------------------------

## Model formulas
mod_forms <- list(
  ~ s(temp_bottom5m, k = 3) + s(oxygen_bottom5m, k = 3) + s(depth_m, k = 3) + s(phi, k = 3), 
  ~ s(temp_bottom5m, k = 3) + s(pH_bottom5m, k = 3) + s(depth_m, k = 3) + s(phi, k = 3), 
  ~ s(temp_bottom5m, k = 3) + s(oxygen_bottom5m, k = 3) + s(X, Y, k = 20) + s(X, Y, by = cold_pool_2C, k = 20), 
  ~ s(temp_bottom5m, k = 3) + s(pH_bottom5m, k = 3) + s(X, Y, k = 20) + s(X, Y, by = cold_pool_2C, k = 20)
)

## Specifications for each species
## Threshold is size threshold for classifying individuals as adults, in cm
## select is whether to choose best model based on time-series ("TSCV") or 5-fold ("CV") cross validation
specs <- list(
  
  species = c("arrowtooth flounder", "Pacific halibut", "Pacific cod", 
              "walleye pollock", "yellowfin sole", "northern rock sole", 
              "snow crab", "red king crab"), 
  
  survey_file = c(
    here("data", "trawl_surveys_size_binned", c(
      "ebs.srvy98.atf.cpue_data.Rdata", "ebs.srvy98.halibut.cpue_data.Rdata", 
      "ebs.srvy98.pcod.cpue_data.Rdata", "ebs.srvy98.plk.cpue_data.Rdata",
      "ebs.srvy98.yfs.cpue_data.Rdata", "ebs.srvy98.nrs.cpue_data.Rdata"
    )), 
    here("data", "trawl_surveys", c(
      "snow_crab.csv", "red_king_crab.csv"
    ))
  ), 
  
  threshold = c(48, 50, 58, 38, 30, 31, NA, NA)
  
)

# Run Models ----------------------------------------------------------------------------

dir.create(here("output"))

## Read in ROMS-NPZ data
source(here("analysis", "load_ROMS.R"))
all_vars <- c("sim", "year", "station_id", "latitude", "longitude", "area_swept_km2", unique(unlist(lapply(mod_forms, all.vars))))
ROMS_full <- drop_na(ROMS_full[,all_vars])

## Save ROMS data
write.csv(ROMS_full, here("data", "ROMS_surveyrep_joined.csv"), row.names = FALSE)
write.csv(ebs_stations, here("data", "ebs_stations.csv"), row.names = FALSE)

## Expand specs data frame by adding rows for juveniles and adults
specs <- as.data.frame(specs, row.names = NULL)
specs_rep <- rep(seq_len(dim(specs)[1]), 1 + !is.na(specs$threshold))
specs <- specs[specs_rep,]
specs$length_bin <- c("juvenile", "adult")[c(1, 2 - diff(specs_rep))]
specs$length_bin[is.na(specs$threshold)] <- NA
rm(specs_rep)

## Number of simultaneous jobs to run
cores <- 6

## Launch each species / length bin on new R processes as they become available
for (i in 1:nrow(specs)) {
  
  species <- specs$species[i]
  survey_file <- specs$survey_file[i]
  threshold <- specs$threshold[i]
  length_bin <- specs$length_bin[i]
  
  if (!is.na(length_bin)) {
   
    source(here("analysis", "load_trawl_size_binned.R"))
    
    model_data <- joined_data |> filter(bin == length_bin)
    val_data <- val_data |> filter(bin == length_bin)
     
  } else {
    
    source(here("analysis", "load_trawl.R"))
    
    model_data <- joined_data
    
  }
  
  rstudioapi::jobRunScript(
    here("analysis", "gam_averaging.R"),
    name = ifelse(is.na(length_bin), species, paste0(species, " (", length_bin, ")")),
    workingDir = here(),
    importEnv = TRUE
  )
  
  Sys.sleep(10)
  n_running <- length(list.files(here("output"), pattern = "running", recursive = TRUE))
  
  while(n_running >= cores) {
    Sys.sleep(10)
    n_running <- length(list.files(here("output"), pattern = "running", recursive = TRUE))
  }
  
}

## Wait until all jobs are done to continue
n_complete <- length(list.files(here("output"), pattern = "complete", recursive = TRUE))
while(n_complete < nrow(specs)) {
  Sys.sleep(10)
  n_complete <- length(list.files(here("output"), pattern = "complete", recursive = TRUE))
}

# Derived outputs -----------------------------------------------------------------------

## Compute overlap for all species / size bin pairs
source(here("analysis", "compute_overlap.R"))

## Compute empirical range metrics (centroids and area occupied)
source(here("analysis", "range_metrics_empirical.R"))

## Projected and fitted vs. observed maps 
source(here("analysis", "distribution_maps.R"))

## Plots of the cold pool SVCs
source(here("analysis", "cold_pool_svc_plots.R"))

## Compare fits and range shifts projected by different candidate models
source(here("analysis", "model_comparison.R"))

## Write outputs as netcdf
source(here("analysis", "save_ncdf.R"))