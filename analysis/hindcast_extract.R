
pkgs <- c("aclim2sdms", "Bering10KThredds", "stars", "sf", "dplyr", "here")
sapply(pkgs, require, character.only = TRUE)

# Read in ROMS level 2 hindcast
hind_lvl2 <- lapply(list.files(here("data", "roms_level2"), pattern = "B10K-K20P19_CORECFS", full.names = TRUE), readRDS)
names(hind_lvl2) <- vapply(hind_lvl2, names, character(1))
roms_dates <- as.Date(st_get_dimension_values(hind_lvl2[[1]], "ocean_time"))

# Read in bathymetry and sediment
bathy <- get_bathymetry()
phi <- get_sediment()

# Read in sampled survey stations by year
srvy <- read.csv(here("data", "trawl_surveys", "ebs_stations_by_year.csv"))
srvy <- srvy |> mutate(date = as.Date(as.POSIXct(DATETIME, format = "%d-%b-%Y %H:%M%OS")))
names(srvy) <- tolower(names(srvy))
srvy <- srvy |> 
  group_by(year, station_id = station) |> 
  summarize(
    latitude = mean(latitude), 
    longitude = mean(longitude), 
    date = unique(date),
    area_swept_km2 = mean(area_swept_ha) / 100, 
    .groups = "drop"
  ) |> 
  mutate(sampled = TRUE)

# Read in average location and date of sampling, expand to all years to create "survey-replicated" dataset
ebs_stns <- read.csv(here("data", "ebs_stations.csv"))
srvy_rep <- do.call("rbind", lapply(1982:2024, \(x) mutate(ebs_stns, year = x)))
srvy_rep$date <- as.Date(paste0(srvy_rep$year, "-01-01")) + srvy_rep$day_of_year
srvy_rep <- srvy_rep |> 
  dplyr::select(year, station_id, latitude, longitude, date) |> 
  mutate(
    area_swept_km2 = round(mean(srvy$area_swept_km2), 4), 
    sampled = FALSE
  )
                       
# Join, overwrite survey-replicated locations / dates with actual, where relevant
srvy <- srvy_rep |> 
  rows_update(srvy, by = c("year", "station_id"), unmatched = "ignore") |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |> 
  add_utm(utm_crs = "+proj=utm +zone=2 +datum=WGS84")

srvy[,c(names(hind_lvl2), "depth_m", "phi")] <- NA

# Loop over survey observations, extract covariates
for (i in 1:nrow(srvy)) {
  
  cat(paste0("\r", i, "/", nrow(srvy), " (", round(100*(i/nrow(srvy))), "%)"))
  
  roms_date <- which.min(abs(roms_dates - srvy$date[i]))
  
  srvy[i, names(hind_lvl2)] <- vapply(hind_lvl2, \(x) st_extract(x[,,,roms_date], srvy[i,])[1, 1, drop = TRUE], numeric(1))
  srvy$depth_m[i] <- -as.numeric(st_extract(bathy, srvy[i,])[1, 1, drop = TRUE])
  srvy$phi[i] <- as.numeric(st_extract(phi, srvy[i,])[1, 1, drop = TRUE])
  
  if (any(is.na(srvy[i, names(hind_lvl2)]))) {
    
    centroids <- st_centroid(st_as_sf(hind_lvl2$temp[,,,roms_date]))
    nearest <- centroids[which.min(st_distance(srvy[i,], centroids)),]
    srvy[i, names(hind_lvl2)] <- vapply(hind_lvl2, \(x) st_extract(x[,,,roms_date], nearest)[1, 1, drop = TRUE], numeric(1))
    
  }
  
}

write.csv(st_drop_geometry(srvy), "data/surveyrep_observed_1982-2024.csv", row.names = FALSE)

rm(list = ls())
