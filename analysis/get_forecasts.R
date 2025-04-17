
pkgs <- c("here", "stars", "dplyr", "curl", "aclim2sdms", "Bering10KThredds")
sapply(pkgs, require, character.only = TRUE)

# Setup -----------------------------------------------------------------------

dir.create(tmp <- tempdir())

url_base <- "https://data.pmel.noaa.gov/aclim/thredds/fileServer/files/B10K-K20_CORECFS/Level3/"

# Variables to download and corresponding links
fcst_urls <- list(
  temp_bottom5m = paste0(url_base, "B10K-K20_CORECFS_2025forecast_temp_bottom5m.nc"), 
  pH_bottom5m = paste0(url_base, "B10K-K20_CORECFS_2025forecast_pH_bottom5m.nc"),
  oxygen_bottom5m = paste0(url_base, "B10K-K20_CORECFS_2025forecast_oxygen_bottom5m.nc")
)

file_names <- paste0(tmp, "/", names(fcst_urls), ".nc")

roms_fcst <- setNames(vector("list", length(fcst_urls)), names(fcst_urls))

# Download --------------------------------------------------------------------

for (i in 1:length(fcst_urls)) {
  
  curl::curl_download(fcst_urls[[i]], file_names[i])
  
  roms_i <- read_ncdf(file_names[i], var = paste0("forecast_", names(fcst_urls)[i]))
  st_crs(roms_i) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  # Choose date closest to July 1st
  date_sel <- which.min(abs(as.Date(st_get_dimension_values(roms_i, "time")) - as.Date("2025-07-01")))
  roms_i <- roms_i |> slice(date_sel, along = "time")
  
  # Crop to EBS
  ebs <- get_ebs_shapefile() |> st_transform("+proj=longlat +datum=WGS84") |> st_shift_longitude()
  roms_i <- roms_i |> st_crop(ebs)
  
  # Rename, strip units
  names(roms_i) <- names(fcst_urls)[i]
  roms_i[[names(fcst_urls)[i]]] <- as.numeric(roms_i[[names(fcst_urls)[i]]])
  
  roms_fcst[[i]] <- roms_i
  
}

# Join, add other variables ---------------------------------------------------

roms_fcst <- Reduce("c", roms_fcst)

phi <- setNames(st_warp(get_sediment(), roms_fcst), "phi")
depth <- setNames(st_warp(get_bathymetry(), roms_fcst), "depth_m")

coords <- roms_fcst |> 
  st_coordinates() |> 
  mutate(xi_rho = rotate_lon(xi_rho)) |> 
  add_utm(c("xi_rho", "eta_rho"), utm_crs = "+proj=utm +zone=2 +datum=WGS84")

roms_fcst <- roms_fcst |> c(phi, depth) |> mutate(X = coords$X, Y = coords$Y)

roms_fcst$cold_pool_2C <- sum(roms_fcst$temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(roms_fcst$temp_bottom5m))

saveRDS(roms_fcst, here("data", "roms_forecast.rds"))