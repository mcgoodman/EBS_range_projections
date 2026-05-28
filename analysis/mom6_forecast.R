
library("tidyverse")
library("stars")
library("units")
library("here")
library("aclim2sdms")
library("BeringSeaData")

# Script: 
# Reads in and formats MOM6 April 30th 2026 conditions ("nowcast")
# Reads in MOM6 April 30th and July 1st hindcast and computes April 30th anomaly
# Applies April 30th anomaly to July 1st climatology to derive forecast

# MOM6 April 30th -----------------------------------------

## Read in ------------------------------------------------

# Reference grid with latitude and longitude
mom6_grid <- read_ncdf(
  here("data", "mom6_forecast", "ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc"),
  var = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Forecast bottom temperature, o2, and H+ concentration
mom6_nowcast <- list(
  tob = read_ncdf(here("data", "mom6_forecast", "tob.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "tob"), 
  btm_o2 = read_ncdf(here("data", "mom6_forecast", "btm_o2.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "btm_o2"), 
  btm_htotal = read_ncdf(here("data", "mom6_forecast", "btm_htotal.nep.iq0-342jq446-743.hcast.daily.e202604.20260101.nc"), var = "btm_htotal")
)

# Join forecast variables
mom6_nowcast <- mom6_nowcast |> purrr::map2(names(mom6_nowcast), setNames)
mom6_nowcast <- Reduce("c", mom6_nowcast)

# Subset to July 1st (middle of survey)
july1 <- which(as.Date(st_get_dimension_values(mom6_nowcast, "time")) == as.Date("2026-04-30"))
mom6_nowcast <- mom6_nowcast |> slice(july1, along = "time")

# Append forecasted variables to curvilinear grid,
# and transform O2 from mol/kg --> (mmol/m3)
# H+ (mol/kg) --> pH (with 1.025 adjustment for seawater density)
mom6_nowcast <- mom6_grid |> mutate(
  temp_bottom5m = drop_units(c(mom6_nowcast$tob)), 
  oxygen_bottom5m = as.numeric(c(mom6_nowcast$btm_o2)) * 1e3 * 1.025e3,
  pH_bottom5m = drop_units(-log10(c(mom6_nowcast$btm_htotal) * 1.025))
)

## Warp to ROMS grid --------------------------------------

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
mom6_nowcast <- st_as_stars(mom6_nowcast, curvilinear = coords)

# Warp to grid returned by BeringSeaData::get_mom6_nep
mom6_apr30 <- readRDS(here("data", "mom6_hindcast", "mom6_hindcast_apr30.rds"))
cefi_grid <- mom6_apr30 |> select() |> slice(1, along = "time")
mom6_nowcast <- mom6_nowcast |> st_warp(cefi_grid, threshold = 0.5)

## Append additional variables ----------------------------

# Add UTM coordinates, depth, and sediment grain size
coords <- mom6_nowcast |> st_coordinates()
phi <- setNames(st_warp(get_sediment(), mom6_nowcast), "phi")
depth <- setNames(st_warp(get_bathymetry(), mom6_nowcast), "depth_m")
mom6_nowcast <- mom6_nowcast |> c(phi, depth) |> mutate(X = coords$x/1000, Y = coords$y/1000)

# Crop to EBS survey region
ebs <- get_ebs_shapefile()
mom6_nowcast <- mom6_nowcast[ebs]

# Add cold pool, truncate depth, remove unecessary variables
mom6_nowcast <- mom6_nowcast |> 
  mutate(
    cold_pool_2C = sum(temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(temp_bottom5m)),
    depth_m = pmin(-pmin(0, depth_m), max(c(mom6_apr30$depth_m), na.rm = TRUE))
  ) |> 
  select(temp_bottom5m, oxygen_bottom5m, pH_bottom5m, phi, depth_m, cold_pool_2C, X, Y)

# Compute April 30th anomaly ------------------------------

# Climatology for temperature and oxygen (standard mean)
mom6_clim_apr30 <- st_apply(mom6_apr30, 1:2, mean, na.rm = TRUE)

# Years for trend projection
years <- lubridate::year(st_get_dimension_values(mom6_apr30, "time"))
forecast_year <- 2026

# Linear trend projection function
project_trend <- function(y, years, forecast_year = 2026) {
  if (sum(!is.na(y)) < 2) return(NA)
  fit <- lm(y ~ years)
  return(as.numeric(predict(fit, newdata = data.frame(years = forecast_year))))
}

# 1. Project April 30th pH trend to 2026
pH_trend_apr30 <- st_apply(mom6_apr30["pH_bottom5m"], 1:2, project_trend, years = years, forecast_year = forecast_year)

# 2. Compute anomalies
# Temperature and oxygen: standard climatological anomalies
# pH: detrended anomaly (relative to the projected 2026 trend baseline)
mom6_nowcast <- mom6_nowcast |> mutate(
  temp_anom = c(temp_bottom5m) - c(mom6_clim_apr30[['temp_bottom5m']]), 
  oxygen_anom = c(oxygen_bottom5m) - c(mom6_clim_apr30[['oxygen_bottom5m']]), 
  pH_anom = c(pH_bottom5m) - c(pH_trend_apr30[[1]])
)

# Apply to July climatology -------------------------------

mom6_july1 <- readRDS(here("data", "mom6_hindcast", "mom6_hindcast_july1.rds"))

# 3. Climatology/Trend projection for July 1st
mom6_clim_july1 <- st_apply(mom6_july1, 1:2, mean, na.rm = TRUE)
pH_trend_july1 <- st_apply(mom6_july1["pH_bottom5m"], 1:2, project_trend, years = years, forecast_year = forecast_year)

# 4. Construct persistence forecast for July 1st
# Start from the spatial grid structure of the nowcast
mom6_fcst <- mom6_nowcast |> mutate(
  temp_bottom5m = c(mom6_clim_july1[['temp_bottom5m']]) + temp_anom,
  oxygen_bottom5m = c(mom6_clim_july1[['oxygen_bottom5m']]) + oxygen_anom,
  pH_bottom5m = c(pH_trend_july1[[1]]) + pH_anom
)

# 5. Update derived variables (Cold Pool Index)
mom6_fcst <- mom6_fcst |> mutate(
  cold_pool_2C = sum(temp_bottom5m < 2, na.rm = TRUE)/sum(!is.na(temp_bottom5m))
)

# Keep only the columns expected by the downstream SDM models
mom6_fcst <- mom6_fcst |> select(temp_bottom5m, oxygen_bottom5m, pH_bottom5m, phi, depth_m, cold_pool_2C, X, Y)

# Save the final forecast
saveRDS(mom6_fcst, here("data", "mom6_forecast", "mom6_forecast.rds"))

ggplot() + 
  geom_stars(aes(fill = temp_anom), data = mom6_fcst) + 
  scale_fill_gradientn(
    colors = c("#151a44", "#1156bf", "#69a5bd", "white", "#d08c75", "#ab2923", "#400813"),
    rescaler = ~ scales::rescale_mid(.x, mid = 0), 
    na.value = NA
  ) + 
  coord_sf()
