
pkgs <- c("here", "stars", "dplyr", "curl", "aclim2sdms", "BeringSeaData")
sapply(pkgs, require, character.only = TRUE)

# hauls
if(!exists("years")) years <- c(1993:2019, 2021:2022)
hauldata <- get_hauldata("EBS", years = years)
hauldata$date <- as.Date(as.POSIXct(hauldata$date_time))
hauldata$year <- lubridate::year(hauldata$date)

# Path to write files to
dir.create(mom6_dir <- here("data", "mom6_level2"), recursive = TRUE)

# List all variables, climate scenarios, and earth models to download data for
specs <- data.frame(
  cefi_name = c("tob", "btm_o2", "btm_htotal"),
  cefi_category = c("ocean_daily", "ocean_cobalt_daily_2d", "ocean_cobalt_daily_2d")
)

mom6 <- setNames(vector("list", nrow(specs)), specs$cefi_name)

# Download MOM6 hindcast ------------------------------------------------------

for (i in 1:nrow(specs)) {
  
  print(paste0(i, " (", round((i/nrow(specs))*100), ")%"))
  
  mom6_i <- vector("list", length(years))
  
  for (y in seq_along(years)) {
    
    start_date <- min(hauldata$date[hauldata$year == years[y]]) - 1
    end_date <- max(hauldata$date[hauldata$year == years[y]]) + 1
    
    mom6_i[[y]] <- get_mom6_nep(
      specs$cefi_name[i], freq = "daily", category = specs$cefi_category[i], 
      release = "r20250818", start_date = start_date, end_date = end_date
    )
    
  }
  
  mom6[[specs$cefi_name[i]]] <- do.call("c", mom6_i)
  
}

mom6 <- do.call("c", mom6)

# Extract for hauls -----------------------------------------------------------

# Convert haul data to sf for extracting from MOM6 outputs
hauldata <- hauldata |> 
  select(year, station_id = station, date, lon = longitude_dd_start,
         lat = latitude_dd_end, area_swept_km2) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |> 
  st_transform(st_crs(mom6)) |> 
  mutate(X = st_coordinates(geometry)[,1], Y = st_coordinates(geometry)[,2])

# Extract
mom6_survey <- mom6 |> st_extract(hauldata, time_column = "date")

mom6_dates <- st_get_dimension_values(mom6, "time")

# Extract nearest for coordinates outside raster
for (i in seq_len(nrow(mom6_survey))) {
  
  if (any(is.na(mom6_survey[i, names(mom6)]))) {
    
    mom6_i <- slice(mom6, which(mom6_dates == hauldata$date[i]), along = "time")
    centroids <- suppressWarnings(st_centroid(st_as_sf(mom6_i)))
    nearest <- centroids[which.min(st_distance(mom6_survey[i,], centroids)),]
    mom6_survey[i,names(mom6)] <- st_drop_geometry(st_extract(mom6_i, nearest)[,names(mom6)])
    
  }
  
}

# Convert units
mom6_survey <- mom6_survey |> 
  select(names(mom6)) |> 
  st_drop_geometry() |> 
  mutate(across(everything(), units::drop_units)) |> 
  mutate(
    oxygen_bottom5m = btm_o2 * 1e3 * 1.025e3, # mol/kg --> (mmol/m3)
    pH_bottom5m = -log10(btm_htotal * 1.025)
  ) |> 
  select(temp_bottom5m = tob, oxygen_bottom5m, pH_bottom5m)

hauldata <- cbind(hauldata, mom6_survey)

# Add static variables
bathy <- st_warp(get_bathymetry(), slice(mom6, 1, along = "time"))
phi <- st_warp(get_sediment(), slice(mom6, 1, along = "time"))

hauldata <- hauldata |> 
  cbind(st_drop_geometry(st_extract(phi, hauldata))) |> 
  cbind(st_drop_geometry(st_extract(bathy, hauldata))) |>
  mutate(depth_m = -depth_m) |> 
  rename(latitude = lat, longitude = lon)

hauldata <- st_drop_geometry(hauldata)

write.csv(hauldata, here("data", "surveyrep_observed_1982-2022.csv"), row.names = FALSE)

# Create annual summer snapshot (July 1) dataset ------------------------------

july_dates <- as.Date(paste0(years, "-07-01"))
mom6_july1 <- mom6 |> slice(which(mom6_dates %in% july_dates), along = "time") |> units::drop_units()

phi_rep <- st_replicate(phi, "ocean_time", values = july_dates)
bathy_rep <- st_replicate(bathy, "ocean_time", values = july_dates)

mom6_july1 <- mom6_july1 |> st_set_dimensions(3, values = july_dates, names = "ocean_time")
mom6_july1 <- c(mom6_july1, phi_rep)
mom6_july1 <- c(mom6_july1, bathy_rep)

mom6_july1 <- mom6_july1 |> 
  mutate(
    oxygen_bottom5m = btm_o2 * 1e3 * 1.025e3,
    pH_bottom5m = -log10(btm_htotal * 1.025), 
    depth_m = -depth_m
  ) |> 
  select(temp_bottom5m = tob, oxygen_bottom5m, pH_bottom5m, phi, depth_m)

coords <- st_coordinates(mom6_july1)
mom6_july1 <- mom6_july1 |> mutate(X = coords$x/1000, Y = coords$y/1000) 

saveRDS(mom6_july1, here("data", "mom6", "mom6_hindcast.rds"))
