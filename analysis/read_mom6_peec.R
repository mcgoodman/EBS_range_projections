
# Script to read Kelly's outputs from 2026 PEEC and extract at survey locations / times

library("here")
library("stars")
library("dplyr")
library("curl")
library("aclim2sdms")
library("BeringSeaData")
library("units")

# hauls
if(!exists("years")) years <- c(1993:2019, 2021:2022)
hauldata <- get_hauldata("EBS", years = years)
hauldata$date <- as.Date(as.POSIXct(hauldata$date_time))
hauldata$year <- lubridate::year(hauldata$date)

# Directory containing downloaded MOM6 AK hindcast for 2026 PEEC
dir.create(mom6_dir <- here("data", "mom6_hindcast"), recursive = TRUE)

# List all variables to read in data for
vars <- c("tob", "btm_o2", "btm_htotal")

mom6 <- setNames(vector("list", length(vars)), vars)

# Read in MOM6 hindcast -------------------------------------------------------

for (i in seq_along(vars)) {
  
  print(paste0(i, " (", round((i/length(vars))*100), ")%"))
  
  mom6_i <- vector("list", length(years))

  files_i <- list.files(file.path(mom6_dir, vars[i]), pattern = "*.nc")
  file_years <- as.integer(substr(gsub(".*e202604\\.(.*)\\.nc.*", "\\1", files_i), 1, 4))
  
  for (y in seq_along(years)) {
    
    mom6_i[[y]] <- suppressWarnings(suppressMessages(
      read_ncdf(
        file.path(mom6_dir, vars[i], files_i[file_years == years[y]]), 
        var = vars[i]
      )
    ))
    
    start_date <- min(c(hauldata$date[hauldata$year == years[y]], as.Date(paste0(years[y], "-04-30")))) - 1
    end_date <- max(hauldata$date[hauldata$year == years[y]]) + 1
    dates_iy <- st_get_dimension_values(mom6_i[[y]], "time")

    mom6_i[[y]] <- mom6_i[[y]] |> 
      slice(which(dates_iy >= start_date & dates_iy <= end_date), along = "time")
    
  }
  
  mom6[[vars[i]]] <- do.call("c", mom6_i)
  mom6[[vars[i]]] <- mom6[[vars[i]]] |> st_set_dimensions(
    values = do.call("c", lapply(mom6_i, st_get_dimension_values, "time")),
    which = 3, names = "time"
  )
  
}

# Join
mom6 <- do.call("c", lapply(mom6, drop_units))
mom6_dates <- as.Date(st_get_dimension_values(mom6, "time"), origin = as.Date("1993-01-01"))
mom6 <- st_set_dimensions(mom6, "time", values = mom6_dates)

# Format output for extraction ------------------------------------------------

# Reference grid with latitude and longitude
mom6_grid <- read_ncdf(
  here("data", "mom6_forecast", "ocean_static.nep.iq0-342jq446-743.hcast.static.e202604.20240101.nc"),
  var = c("geolon", "geolat"),
  curvilinear = c("geolon", "geolat")
)

# Polynomial interpolation of longitude by row
x <- drop_units(mom6_grid$geolon)
mmx <- cbind(1, poly(I(1:ncol(x)), 3))
for (i in 1:nrow(x)) x[i, is.na(x[i,])] <- (mmx %*% coef(lm(I(x[i,]) ~ 0 + I(mmx))))[is.na(x[i,])]

# Polynomial interpolation of latitude by row
y <- drop_units(mom6_grid$geolat)
mmy <- cbind(1, poly(I(1:ncol(y)), 3))
for (i in 1:nrow(y)) y[i, is.na(y[i,])] <- (mmy %*% coef(lm(I(y[i,]) ~ 0 + I(mmy))))[is.na(y[i,])]

# Overwrite coordinates
mom6_grid$geolon <- as_units(x, units(mom6_grid$geolon))
mom6_grid$geolat <- as_units(y, units(mom6_grid$geolat))

# Replicate across time
mom6_grid <- st_replicate(mom6_grid, "time", mom6_dates)

# Append forecasted variables to curvilinear grid,
# and transform O2 from mol/kg --> (mmol/m3)
# H+ (mol/kg) --> pH (with 1.025 adjustment for seawater density)
mom6 <- mom6_grid |> mutate(
  temp_bottom5m = c(mom6$tob), 
  oxygen_bottom5m = as.numeric(c(mom6$btm_o2)) * 1e3 * 1.025e3,
  pH_bottom5m = -log10(c(mom6$btm_htotal) * 1.025)
)

# Transform to UTM
ak_coast <- BeringSeaData::get_ak_coast()
mom6 <- st_transform(mom6, st_crs(ak_coast))

# Extract for hauls -----------------------------------------------------------

# Convert haul data to sf for extracting from MOM6 outputs
hauldata <- hauldata |> 
  select(year, station_id = station, date, lon = longitude_dd_start, 
         lat = latitude_dd_start, area_swept_km2) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  st_transform(st_crs(mom6)) |> 
  mutate(X = st_coordinates(geometry)[,1], Y = st_coordinates(geometry)[,2])

# Nearest-neighbor extraction
mom6_survey <- mom6 |> st_extract(hauldata, time_column = "date")

hauldata <- mom6_survey |>
  select(temp_bottom5m, oxygen_bottom5m, pH_bottom5m) |> 
  cbind(st_drop_geometry(hauldata))

# Add static variables
bathy <- st_warp(get_bathymetry(), slice(mom6, 1, along = "time"))
phi <- st_warp(get_sediment(), slice(mom6, 1, along = "time"))

hauldata <- hauldata |> 
  mutate(
    phi = do.call("c", lapply(split(hauldata, hauldata$year), \(x) st_extract(phi, x)$phi)), 
    depth_m = -do.call("c", lapply(split(hauldata, hauldata$year), \(x) st_extract(bathy, x)$depth_m))
  ) |> 
  rename(latitude = lat, longitude = lon)

hauldata <- st_drop_geometry(hauldata)

write.csv(hauldata, here("data", "surveyrep_observed_1982-2022.csv"), row.names = FALSE)

# Create annual snapshot (April 30) comparable to persistence forecast --------

# Subset to April 30th
apr_dates <- as.Date(paste0(years, "-07-01"))
mom6_apr30 <- mom6 |> slice(which(mom6_dates %in% apr_dates), along = "time") |> units::drop_units()
mom6_apr30 <- mom6_apr30 |> st_set_dimensions(3, values = apr_dates, names = "time")

# Add bathymetry and sediment grain size
phi_rep <- st_replicate(phi, "time", values = apr_dates)
bathy_rep <- st_replicate(bathy, "time", values = apr_dates)
mom6_apr30 <- c(mom6_apr30, phi_rep)
mom6_apr30 <- c(mom6_apr30, bathy_rep)

# Warp to regrid for cropping
ebs <- get_ebs_shapefile("EBS", type = "boundary")
regrid <- get_mom6_nep(start_date = as.Date("2000-01-01"), end_date =  as.Date("2000-02-01"), extent = st_bbox(ebs))
regrid <- slice(regrid, 1, along = "time")
coords <- list(ih = rotate_lon(x, from = "0/360"), jh = y)
mom6_apr30 <- mom6_apr30 |> st_as_stars(curvilinear = coords) |> st_warp(regrid, threshold = 0.5)

# Crop to EBS
mom6_apr30 <- mom6_apr30[ebs]

# Assign UTM coordinates as attribute
coords <- st_coordinates(mom6_apr30)
mom6_apr30 <- mom6_apr30 |> mutate(X = coords$x/1000, Y = coords$y/1000) 

# Remove unecessary variables, clip bathymetry
mom6_apr30 <- mom6_apr30 |> 
  mutate(depth_m = pmin(-pmin(0, depth_m), max(hauldata$depth_m))) |> 
  select(temp_bottom5m, oxygen_bottom5m, pH_bottom5m, phi, depth_m, X, Y)

saveRDS(mom6_apr30, here("data", "mom6_hindcast", "mom6_hindcast.rds"))
