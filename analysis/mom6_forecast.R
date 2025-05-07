
library("tidyverse")
library("stars")
library("units")
library("here")
library("aclim2sdms")
library("Bering10KThredds")

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
  here("data", "mom6", "ocean_static_ak.nc"),
  sub = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Forecast bottom temperature, o2, and H+ concentration
mom6_fcst <- read_stars(
  here("data", "mom6", "mom6nep_hc202411_forecast_2025.nc"), 
  sub = c("tob", "btm_o2", "btm_htotal")
)

# Subset to July 1st (middle of survey)
july1 <- which(as.Date(st_get_dimension_values(mom6_fcst, "time")) == as.Date("2022-07-01"))
mom6_fcst <- mom6_fcst |> slice(july1, along = "time")

# Append forecasted variables to curvilinear grid,
# and transform O2 from mol/kg --> (mmol/m3),
# H+ (mol/kg) --> pH (with 1.025 adjustment for seawater density)
mom6_fcst <- mom6_grid |> mutate(
  temp_bottom5m = drop_units(c(mom6_fcst$tob)), 
  oxygen_bottom5m = as.numeric(c(mom6_fcst$btm_o2)) * 1e3 * 1.025e3,
  pH_bottom5m = -drop_units(log10(c(mom6_fcst$btm_htotal) * 1.025))
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
ak_coast <- Bering10KThredds::get_ak_coast()
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

saveRDS(mom6_fcst, here("data", "mom6", "mom6_forecast.rds"))

# MOM6 climatology ----------------------------------------

mom6_clim <- read_stars(
  here("data", "mom6", "mom6nep_hc202411_daily_clim_1993-2022.nc"), 
  sub = c("tob", "btm_o2", "btm_htotal")
)

mom6_clim <- mom6_clim |> slice(
  which(as.Date(st_get_dimension_values(mom6_clim, "time")) == as.Date("2022-07-01")), 
  along = "time"
)

mom6_clim <- mom6_grid |> mutate(
  temp_bottom5m = drop_units(c(mom6_clim$tob)), 
  oxygen_bottom5m = as.numeric(c(mom6_clim$btm_o2)) * 1e3 * 1.025e3,
  pH_bottom5m = -drop_units(log10(c(mom6_clim$btm_htotal) * 1.025))
)

# Overwrite coordinates
st_dimensions(mom6_clim)$x$values <- as_units(x, ll_units)
st_dimensions(mom6_clim)$y$values <- as_units(y, ll_units)

# Transform to UTM
mom6_clim <- st_transform(mom6_clim, st_crs(ak_coast))

# Warp to ROMS grid
mom6_clim <- mom6_clim |> st_warp(roms_grid)

# Crop to EBS survey region
mom6_clim <- mom6_clim[ebs]

saveRDS(mom6_clim, here("data", "mom6", "mom6_climatology.rds"))

# "Bias-corrected" MOM6 forecast --------------------------

roms_hindcast <- readRDS(here("data", "roms_level2_bc_annual", "CORECFS_hindcast.rds"))

roms_clim <- c(
  summarize_proj(roms_hindcast, var = temp_bottom5m), 
  summarize_proj(roms_hindcast, var = oxygen_bottom5m), 
  summarize_proj(roms_hindcast, var = pH_bottom5m)
)

# Difference between ROMS and MOM6 climatology
roms_mom6_diff <- mom6_clim |> 
  select(-geolon, -geolat) |> 
  mutate(
    temp_diff = temp_bottom5m - c(roms_clim$temp_bottom5m), 
    oxygen_diff = oxygen_bottom5m - c(roms_clim$oxygen_bottom5m),
    pH_diff = pH_bottom5m - c(roms_clim$pH_bottom5m)
  )

mom6_fcst_adj <- mom6_fcst |> mutate(
  temp_bottom5m = temp_bottom5m - c(roms_mom6_diff$temp_diff), 
  oxygen_bottom5m = pmax(oxygen_bottom5m - c(roms_mom6_diff$oxygen_diff), 0),
  pH_bottom5m = pH_bottom5m - c(roms_mom6_diff$pH_diff)
)


mom6_fcst_adj$cold_pool_2C <- sum(mom6_fcst_adj$temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(mom6_fcst_adj$temp_bottom5m))

saveRDS(mom6_fcst_adj, here("data", "mom6", "mom6_forecast_adj.rds"))