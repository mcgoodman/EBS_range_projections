
pkgs <- c("here", "stars", "sf", "dplyr", "aclim2sdms")
sapply(pkgs, require, character.only = TRUE)

# This script joins bias-corrected ROMS outputs pertaining to the same scenario / earth model
# together, adding in sediment grain size, depth, and UTM coordinates, so that the file
# can later be used to obtain predictions from fitted GAM models

# List ROMS files and associated variable names, scenarios, and earth models
roms_files <- list.files(here("data", "roms_level2_bc"), full.names = TRUE)
roms_files <- c(roms_files, list.files(here("data", "roms_level2"), pattern = "CORECFS", full.names = TRUE))
roms_base <- basename(roms_files)
roms_info <- strsplit(gsub("bc.rds|.rds", "", roms_base), "_")
roms_var <- vapply(roms_info, \(x) paste(tail(x, 2), collapse = "_"), "character")
roms_version <- vapply(roms_info, \(x) x[2], "character") 
roms_model <- ifelse(roms_version == "CORECFS", "CORECFS", vapply(roms_info, \(x) x[3], "character"))
roms_scenario <- ifelse(roms_version == "CORECFS", "hindcast", toupper(vapply(roms_info, \(x) x[4], "character")))

# Unique ROMS sims
roms_sims <- as.data.frame(unique(cbind(model = roms_model, scenario = roms_scenario)))

dir.create(save_dir <- here("data", "roms_level2_bc_annual"))

for (i in 1:nrow(roms_sims)) {
 
  vars_i <- roms_var[roms_model == roms_sims$model[i] & roms_scenario == roms_sims$scenario[i]]
  
  cat(paste0(paste(roms_sims[i,], collapse = " "), " (", i, "/", nrow(roms_sims), ")\n"))
  
  for (j in 1:length(vars_i)) {
    
    cat(paste("   ", vars_i[j], "\n"))
    
    roms_bc <- readRDS(roms_files[roms_var == vars_i[j] & roms_model == roms_sims$model[i] & roms_scenario == roms_sims$scenario[i]])
    
    # Find dates in each year which are closest to July 1st
    dates <- as.Date(st_get_dimension_values(roms_bc, "ocean_time"))
    years <- lubridate::year(dates)
    date_sub <- aggregate(dates, by = list(years), \(x) x[which.min(abs(x - as.Date(paste0(lubridate::year(x), "-07-01"))))])$x
    if(roms_sims$scenario[i] == "hindcast") {
      date_sub <- which(years <= 2022 & dates %in% date_sub)
    } else {
      date_sub <- which(years >= 2023 & dates %in% date_sub) 
    }
    
    # Subset to those dates
    roms_bc <- roms_bc |> slice(date_sub, along = "ocean_time")
    
    # Rename, strip units info (causes problems during model prediction later)
    names(roms_bc) <- vars_i[j]
    roms_bc[[vars_i[j]]] <- as.numeric(roms_bc[[vars_i[j]]])
    
    # For first ROMS file, add in sediment grain size, depth, and UTM coordinates as attributes
    # For ROMS files after that, append data to first ROMS file as a new attribute
    if (j == 1) {
      
      roms_grd <- st_redimension(roms_bc[,,,1], new_dims = st_dimensions(roms_bc)[1:2])
      
      phi <- st_warp(get_sediment(), roms_grd)
      phi_rep <- roms_bc |> mutate(phi = rep(phi$phi.grd, dim(roms_bc)[3])) |> dplyr::select(phi)
      
      depth <- st_warp(get_bathymetry(), roms_grd)
      depth_rep <- roms_bc |> mutate(depth_m = -rep(depth$etopo_bedrock_15arcsecond.tif, dim(roms_bc)[3])) |> select(depth_m)
      
      coords <- roms_bc |> 
        st_coordinates() |> 
        mutate(xi_rho = rotate_lon(xi_rho)) |> 
        add_utm(c("xi_rho", "eta_rho"), utm_crs = "+proj=utm +zone=2 +datum=WGS84")
      
      roms_full <- roms_bc |> mutate(X = coords$X, Y = coords$Y) 
      
      roms_full <- c(roms_full, phi_rep, depth_rep)
      
    } else {
      
      roms_full <- c(roms_full, roms_bc)
      
    }
    
  }
  
  # Subset the ROMS grid to get rid of whitespace
  roms_full <- roms_full |> slice(xi_rho, 70:170) |> slice(eta_rho, 35:155)
  
  saveRDS(roms_full, here(save_dir, paste0(paste(roms_sims[i,], collapse = "_"), ".rds")))
  
}

rm(list = ls())
