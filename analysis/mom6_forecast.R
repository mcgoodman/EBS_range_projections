
library("tidyverse")
library("stars")
library("units")
library("here")
library("aclim2sdms")
library("BeringSeaData")

slice_proj <- function(x, start = 1993, end = 2022) {
  
  time <- st_get_dimension_values(x, "ocean_time")
  year <- lubridate::year(time)
  year_slice <- which(year >= start & year <= end)
  x |> slice(year_slice, along = "ocean_time")
  
}

summarize_proj <- function(x, start = 1993, end = 2022, var = p_occurrence, f = mean) {
  
  var <- enquo(var)
  x <- x |> slice_proj(start, end) |> dplyr::select(!!var)
  x <- st_apply(x, 1:2, f)
  names(x) <- quo_name(var)
  x
  
}

# Read in MOM6 forecast -----------------------------------

# Reference grid with latitude and longitude
mom6_grid <- read_stars(
  here("data", "mom6", "mom6nep_hc202507_ocean_static_ak.nc"),
  sub = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Forecast bottom temperature, o2, and H+ concentration
mom6_fcst <- read_stars(
  here("data", "mom6", "mom6nep_hc202507_daily_anomaly_20250701.nc"), 
  sub = c("tob", "btm_o2", "btm_htotal")
)

# Subset to July 1st (middle of survey)
july1 <- which(as.Date(st_get_dimension_values(mom6_fcst, "time")) == as.Date("2025-07-01"))
mom6_fcst <- mom6_fcst |> slice(july1, along = "time")

# Append forecasted variables to curvilinear grid,
# and transform O2 from mol/kg --> (mmol/m3)
# Leave htotal untransformed to apply as anomaly before transforming later
mom6_fcst <- mom6_grid |> mutate(
  temp_bottom5m = drop_units(c(mom6_fcst$tob)), 
  oxygen_bottom5m = as.numeric(c(mom6_fcst$btm_o2)) * 1e3 * 1.025e3,
  htotal_bottom5m = drop_units(c(mom6_fcst$btm_htotal))
)

# Warp to ROMS grid ---------------------------------------

# Hacky: Interpolate latitude and longitudes for grid cells over land
# which are missing coordinates. Needed for warping to ROMS grid
# Should not affect outcome as these cells are well outside EBS survey

# Polynomial interpolation of longitude by row
x <- drop_units(st_get_dimension_values(mom6_fcst, "x"))
mmx <- cbind(1, poly(I(1:ncol(x)), 3))
for (i in 1:nrow(x)) x[i, is.na(x[i,])] <- (mmx %*% coef(lm(I(x[i,]) ~ 0 + I(mmx))))[is.na(x[i,])]

# Polynomial interpolation of latitude by row
y <- drop_units(st_get_dimension_values(mom6_fcst, "y"))
mmy <- cbind(1, poly(I(1:ncol(y)), 4))
for (i in 1:nrow(y)) y[i, is.na(y[i,])] <- (mmy %*% coef(lm(I(y[i,]) ~ 0 + I(mmy))))[is.na(y[i,])]

# Overwrite coordinates
ll_units <- units(st_dimensions(mom6_fcst)$x$values)
st_dimensions(mom6_fcst)$x$values <- as_units(x, ll_units)
st_dimensions(mom6_fcst)$y$values <- as_units(y, ll_units)

# Transform to UTM
ak_coast <- BeringSeaData::get_ak_coast()
mom6_fcst <- st_transform(mom6_fcst, st_crs(ak_coast))

# Warp to ROMS grid
roms_grid <- readRDS(here("data", "roms_level2_bc_annual", "CORECFS_hindcast.rds"))
roms_grid <- roms_grid |> select() |> slice(1, along = "ocean_time")
mom6_fcst <- mom6_fcst |> st_warp(roms_grid)

# Crop to EBS survey region
ebs <- get_ebs_shapefile() |> st_transform("+proj=longlat +datum=WGS84") |> st_shift_longitude()
mom6_fcst <- mom6_fcst[ebs]

# Append additional variables -----------------------------

phi <- setNames(st_warp(get_sediment(), mom6_fcst), "phi")
depth <- setNames(st_warp(get_bathymetry(), mom6_fcst), "depth_m")

coords <- mom6_fcst |> 
  st_coordinates() |> 
  mutate(xi_rho = rotate_lon(xi_rho)) |> 
  add_utm(c("xi_rho", "eta_rho"), utm_crs = "+proj=utm +zone=2 +datum=WGS84")

mom6_fcst <- mom6_fcst |> c(phi, depth) |> mutate(X = coords$X, Y = coords$Y)

mom6_fcst$cold_pool_2C <- sum(mom6_fcst$temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(mom6_fcst$temp_bottom5m))

# Apply MOM6 anomaly to ROMS climatology ------------------

# H+ (mol/kg) --> pH (with 1.025 adjustment for seawater density)
# pH = -log10(H+ * 1.025)
# H+ = (10^(-pH))/1.025

roms_hindcast <- readRDS(here("data", "roms_level2_bc_annual", "CORECFS_hindcast.rds"))

roms_hindcast <- roms_hindcast |> 
  mutate(htotal_bottom5m = (10^(-pH_bottom5m))/1.025)

roms_clim <- c(
  summarize_proj(roms_hindcast, var = temp_bottom5m), 
  summarize_proj(roms_hindcast, var = oxygen_bottom5m), 
  summarize_proj(roms_hindcast, var = htotal_bottom5m)
)

# Add MOM6 anomaly to ROMS climatology
mom6_fcst_adj <- mom6_fcst |> mutate(
  temp_bottom5m = c(roms_clim$temp_bottom5m) + temp_bottom5m, 
  oxygen_bottom5m = pmax(c(roms_clim$oxygen_bottom5m) + oxygen_bottom5m, min(roms_hindcast$oxygen_bottom5m, na.rm = TRUE)),
  htotal_bottom5m = c(roms_clim$htotal_bottom5m) + htotal_bottom5m, 
  pH_bottom5m = -log10(htotal_bottom5m * 1.025)
)

mom6_fcst_adj$cold_pool_2C <- sum(mom6_fcst_adj$temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(mom6_fcst_adj$temp_bottom5m))

saveRDS(mom6_fcst_adj, here("data", "mom6", "mom6_forecast_adj.rds"))
