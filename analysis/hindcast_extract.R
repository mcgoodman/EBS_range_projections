
pkgs <- c("stars", "sf", "dplyr")
sapply(pkgs, require, character.only = TRUE)

# Read in ROMS level 2 hindcast
hind_lvl2 <- lapply(list.files(here("data", "roms_level2"), pattern = "B10K-K20P19_CORECFS", full.names = TRUE), readRDS)
names(hind_lvl2) <- vapply(hind_lvl2, names, character(1))
roms_dates <- st_get_dimension_values(hind_lvl2[[1]], "ocean_time")

# Read in bathymetry and sediment
bathy <- get_bathymetry()
phi <- get_sediment()

# Read in survey stations by year
srvy <- read.csv(here("data", "trawl_surveys", "ebs_stations_by_year.csv"))
srvy <- srvy |> filter(YEAR <= 2022) |> mutate(DATETIME = as.POSIXct(DATETIME, format = "%m/%d/%Y %H:%M:%OS")) |> arrange(DATETIME)
srvy <- st_as_sf(srvy, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
srvy[,c(names(hind_lvl2), "depth_m", "phi")] <- NA

# Loop over survey observations, extract covariates
for (i in 1:nrow(srvy)) {
  
  cat(paste0("\r", i, "/", nrow(srvy), " (", round(100*(i/nrow(srvy))), "%)"))
  
  roms_date <- which.min(abs(roms_dates - srvy$DATETIME[i]))
  
  srvy[i, names(hind_lvl2)] <- vapply(hind_lvl2, \(x) st_extract(x[,,,roms_date], srvy[i,])[1, 1, drop = TRUE], numeric(1))
  srvy$depth_m[i] <- -as.numeric(st_extract(bathy, srvy[i,])[1, 1, drop = TRUE])
  srvy$phi[i] <- as.numeric(st_extract(phi, srvy[i,])[1, 1, drop = TRUE])
  
}

write.csv(st_drop_geometry(srvy), "data/surveyrep_observed_1982-2022.csv", row.names = FALSE)

rm(list = ls())
