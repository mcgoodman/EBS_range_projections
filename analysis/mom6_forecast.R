
library("tidyverse")
library("stars")
library("units")
library("here")
library("aclim2sdms")
library("BeringSeaData")

# Read in MOM6 forecast -----------------------------------

# Reference grid with latitude and longitude
mom6_grid <- read_ncdf(
  here("data", "mom6_forecast", "ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc"),
  var = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Forecast bottom temperature, o2, and H+ concentration
mom6_fcst <- list(
  tob = read_ncdf(here("data", "mom6_forecast", "tob.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "tob"), 
  btm_o2 = read_ncdf(here("data", "mom6_forecast", "btm_o2.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "btm_o2"), 
  btm_htotal = read_ncdf(here("data", "mom6_forecast", "btm_htotal.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "btm_htotal")
)

# Join forecast variables
mom6_fcst <- mom6_fcst |> purrr::map2(names(mom6_fcst), setNames)
mom6_fcst <- Reduce("c", mom6_fcst)

# Subset to July 1st (middle of survey)
july1 <- which(as.Date(st_get_dimension_values(mom6_fcst, "time")) == as.Date("2026-04-30"))
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
x <- drop_units(mom6_grid$geolon)
mmx <- cbind(1, poly(I(1:ncol(x)), 3))
for (i in 1:nrow(x)) x[i, is.na(x[i,])] <- (mmx %*% coef(lm(I(x[i,]) ~ 0 + I(mmx))))[is.na(x[i,])]

# Polynomial interpolation of latitude by row
y <- drop_units(mom6_grid$geolat)
mmy <- cbind(1, poly(I(1:ncol(y)), 3))
for (i in 1:nrow(y)) y[i, is.na(y[i,])] <- (mmy %*% coef(lm(I(y[i,]) ~ 0 + I(mmy))))[is.na(y[i,])]

# Set dimensions of forecast object using interpolated coordinates
coords <- list(ih = rotate_lon(x, from = "0/360"), jh = y)
mom6_fcst <- st_as_stars(mom6_fcst, curvilinear = coords)

# Warp to grid returned by BeringSeaData::get_mom6_nep
mom6_hindcast <- readRDS(here("data", "mom6_hindcast", "mom6_hindcast.rds"))
cefi_grid <- mom6_hindcast |> select() |> slice(1, along = "time")
mom6_fcst <- mom6_fcst |> st_warp(cefi_grid, threshold = 0.5)

# Append additional variables -----------------------------

# Add UTM coordinates, depth, and sediment grain size
coords <- mom6_fcst |> st_coordinates()
phi <- setNames(st_warp(get_sediment(), mom6_fcst), "phi")
depth <- setNames(st_warp(get_bathymetry(), mom6_fcst), "depth_m")
mom6_fcst <- mom6_fcst |> c(phi, depth) |> mutate(X = coords$x/1000, Y = coords$y/1000)

# Crop to EBS survey region
ebs <- get_ebs_shapefile()
mom6_fcst <- mom6_fcst[ebs]

# Add cold pool, truncate depth, remove unecessary variables
mom6_fcst <- mom6_fcst |> 
  mutate(
    cold_pool_2C = sum(temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(temp_bottom5m)),
    depth_m = pmin(-pmin(0, depth_m), max(c(mom6_hindcast$depth_m), na.rm = TRUE))
  ) |> 
  select(temp_bottom5m, oxygen_bottom5m, pH_bottom5m, phi, depth_m, cold_pool_2C, X, Y)

saveRDS(mom6_fcst, here("data", "mom6_forecast", "mom6_forecast.rds"))
