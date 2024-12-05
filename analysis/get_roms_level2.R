
pkgs <- c("here", "stars", "dplyr", "curl", "aclim2sdms")
sapply(pkgs, require, character.only = TRUE)

# Path to write files to
dir.create(roms_dir <- here("data", "roms_level2"), recursive = TRUE)
dir.create(roms_bc_dir <- here("data", "roms_level2_bc"), recursive = TRUE)

# List all variables, climate scenarios, and earth models to download data for
vars <- c("temp_bottom5m", "oxygen_bottom5m", "pH_bottom5m")
scenarios <- c("SSP126", "SSP585")
models <- c("GFDL", "CESM", "MIROC")

# All possible combinations
specs <- expand.grid(var = vars, scenario = scenarios, earth_model = models, stringsAsFactors = FALSE)
specs$version <- ifelse(grepl("SSP", specs$scenario), "CMIP6", "CMIP5")
specs <- specs |> arrange(version, var)

# Run
for (i in 1:nrow(specs)) {
  
  print(paste0(i, " (", round((i/nrow(specs))*100), ")%"))
  roms_name <- ifelse(grepl("SSP", specs$scenario[i]), "B10K-K20P19_", "B10K-H16_")
  
  # Read in hindacst and obtain weekly weighted means
  if (i == 1 || !exists("hind") || !(specs$var[i] == specs$var[i - 1] & specs$version[i] == specs$version[i -1])) {
    hind <- get_level2(specs$var[i], "hindcast", version = specs$version[i], write_dir = NA, start = 1970, end = 2022)
    saveRDS(hind, paste0(roms_dir, "/", roms_name, "CORECFS_", specs$var[i], ".rds"))
    hind_wkly <- weight_weeks(hind, start = 1985, end = 2014)
  }
  
  # Download historical and forecast runs
  hist <- get_level2(specs$var[i], "historical", scenario = specs$scenario[i], 
                     earth_model = specs$earth_model[i], version = specs$version[i], 
                     write_dir = NA, start = 1985, end = 2014)
  fcst <- get_level2(specs$var[i], "forecast", scenario = specs$scenario[i], 
                     earth_model = specs$earth_model[i], version = specs$version[i], 
                     write_dir = NA, start = 2020)
  
  # Write out non-bias-corrected forecast
  fcst_path <- paste0(roms_dir, "/", roms_name, specs$version[i], "_", tolower(specs$earth_model[i]), "_", tolower(specs$scenario[i]), "_", specs$var[i], ".rds")
  saveRDS(fcst, fcst_path)
  
  # Obtain weekly weighted means for historical
  hist_wkly <- weight_weeks(hist, start = 1985, end = 2014); rm(hist)
  
  # Bias correct
  roms_bc <- bias_correct(fcst, hind_wkly, hist_wkly); rm(fcst, hist_wkly)
  
  # Write out
  bc_path <- paste0(roms_bc_dir, "/", roms_name, specs$version[i], "_", tolower(specs$earth_model[i]), "_", tolower(specs$scenario[i]), "_", specs$var[i], "_bc.rds" )
  saveRDS(roms_bc, bc_path)
  
  rm(roms_bc); gc()
  
}

rm(list = ls())
