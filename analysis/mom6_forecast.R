
library("tidyverse")
library("stars")
library("units")
library("here")
library("aclim2sdms")
library("BeringSeaData")

# Read in MOM6 forecast -----------------------------------

# Reference grid with latitude and longitude
mom6_grid <- read_stars(
  here("data", "mom6", "mom6nep_hc202507_ocean_static_ak.nc"),
  sub = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Forecast bottom temperature, o2, and H+ concentration
mom6_fcst <- read_stars(
  here("data", "mom6", "mom6nep_hc202507_selected_daily_20250701.nc"), 
  sub = c("tob", "btm_o2", "btm_htotal")
)

# Subset to July 1st (middle of survey)
july1 <- which(as.Date(st_get_dimension_values(mom6_fcst, "time")) == as.Date("2025-07-01"))
mom6_fcst <- mom6_fcst |> slice(july1, along = "time")

# Append forecasted variables to curvilinear grid,
# and transform O2 from mol/kg --> (mmol/m3)
# H+ (mol/kg) --> pH (with 1.025 adjustment for seawater density)
mom6_fcst <- mom6_grid |> mutate(
  temp_bottom5m = drop_units(c(mom6_fcst$tob)), 
  oxygen_bottom5m = as.numeric(c(mom6_fcst$btm_o2)) * 1e3 * 1.025e3,
  pH_bottom5m = drop_units(-log10(c(mom6_fcst$btm_htotal) * 1.025))
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

# Warp to grid returned by BeringSeaData::get_mom6_nep
cefi_grid <- readRDS(here("data", "mom6", "mom6_hindcast.rds"))
cefi_grid <- cefi_grid |> select() |> slice(1, along = "ocean_time")
mom6_fcst <- mom6_fcst |> st_warp(cefi_grid)

# Crop to EBS survey region
ebs <- get_ebs_shapefile()
mom6_fcst <- mom6_fcst[ebs]

# Append additional variables -----------------------------

phi <- setNames(st_warp(get_sediment(), mom6_fcst), "phi")
depth <- setNames(st_warp(get_bathymetry(), mom6_fcst), "depth_m")

coords <- mom6_fcst |> st_coordinates()
mom6_fcst <- mom6_fcst |> c(phi, depth) |> mutate(X = coords$x/1000, Y = coords$y/1000)

mom6_fcst$cold_pool_2C <- sum(mom6_fcst$temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(mom6_fcst$temp_bottom5m))

saveRDS(mom6_fcst, here("data", "mom6", "mom6_forecast_adj.rds"))
