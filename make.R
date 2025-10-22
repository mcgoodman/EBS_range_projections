
if (!require("aclim2sdms")) {
  devtools::install()
  require("aclim2sdms")
}

if (!require("Bering10KThredds")) {
  devtools::install_github("mcgoodman/Bering10KThredds")
}

pkgs <- c("here", "dplyr", "tidyr", "purrr", "ggplot2", "sf", "stars", "ncmeta", "mgcv", "aclim2sdms", "Bering10KThredds", "foreach", "doParallel")

sapply(pkgs, require, character.only = TRUE)

# Download & bias-correct ROMS level 2 data ---------------------------------------------

# Years to download data / fit model for
years <- c(1993:2019, 2021:2022)

## This takes up to a couple days - do not re-run if it can be avoided
process_roms <- FALSE

if (process_roms) {
  
  ## Download and bias-correct ROMS outputs
  source(here("analysis", "get_mom6_level2.R"))
  
  ## Download ROMS level 2 forecasts
  source(here("analysis", "mom6_forecast.R"))
  
}

## Read in ROMS-NPZ data
ROMS_data <- read.csv(here("data", "surveyrep_observed_1982-2022.csv")) |> 
  group_by(year) |> 
  mutate(cold_pool_2C = sum(temp_bottom5m < 2)/n()) |> 
  ungroup()

## Read in unique EBS stations
ebs_stations <- read.csv(here("data", "ebs_stations.csv"))

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

## Pre-estimated model weights
weights <- read.csv(here("data", "model_weights.csv"))

# Run Models ----------------------------------------------------------------------------

dir.create(here("output"))

## Expand specs data frame by adding rows for juveniles and adults
specs <- as.data.frame(specs, row.names = NULL)
specs_rep <- rep(seq_len(dim(specs)[1]), 1 + !is.na(specs$threshold))
specs <- specs[specs_rep,]
specs$length_bin <- c("juvenile", "adult")[c(1, 2 - diff(specs_rep))]
specs$length_bin[is.na(specs$threshold)] <- NA
rm(specs_rep)

## Number of simultaneous jobs to run
cores <- 8

## Launch each species / length bin on new R processes as they become available
for (i in 1:nrow(specs)) {
  
  species <- specs$species[i]
  survey_file <- specs$survey_file[i]
  threshold <- specs$threshold[i]
  length_bin <- specs$length_bin[i]
  
  if (!is.na(length_bin)) {
   
    source(here("analysis", "load_trawl_size_binned.R"))
    
    model_data <- model_data |> filter(bin == length_bin)
    
    w_binom <- weights$weight[weights$species == species & weights$bin == length_bin & weights$component == "binomial"]
    w_tw <- weights$weight[weights$species == species & weights$bin == length_bin & weights$component == "tweedie"]
     
  } else {
    
    source(here("analysis", "load_trawl.R"))
    
    w_binom <- weights$weight[weights$species == species & weights$component == "binomial"]
    w_tw <- weights$weight[weights$species == species & weights$component == "tweedie"]
    
  }
  
  rstudioapi::jobRunScript(
    here("analysis", "gam_predictions.R"),
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

# Obtain forecast maps
source(here("analysis", "forecast_maps.R"))