
pkgs <- c("here", "stars", "sf", "dplyr", "tidyr", "aclim2sdms")
sapply(pkgs, require, character.only = TRUE)

ak_coast <- get_ak_coast() ## Alaska coastline map

# Hindcast ------------------------------------------------------------------------------

## Read in survey-replicated CMIP6 ROMS-NPZ hindcast
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CORECFS_BC_hind.Rdata"))

# Read in temperature, pH, and oxygen corresponding to surveys
hind_srvy <- read.csv(here("data", "surveyrep_observed_1982-2022.csv"))

# Survey-replicated datasets contain some surveys which are never actually sampled
stns_keep <- trimws(unique(hind_srvy$STATION[hind_srvy$STATION %in% trimws(hind$station_id)]))

# If the same station is surveyed twice in one year, take mean
hind_srvy <- hind_srvy |>
  group_by(year = YEAR, station_id = STATION) |> 
  summarize(
    temp_bottom5m = mean(temp, na.rm = TRUE), 
    oxygen_bottom5m = mean(oxygen, na.rm = TRUE), 
    pH_bottom5m = mean(pH, na.rm = TRUE),
    area_swept_km2 = mean(AREA_SWEPT_HA) / 100,
    .groups = "drop"
  )

# Average area swept for use in predictions
area_avg <- round(mean(hind_srvy$area_swept_km2), 5)

## Function to clean up station names, put parameter values in columns
clean_roms <- function(data, col) {
  
  col <- enquo(col)
  
  data |> 
    dplyr::select(station_id, year, latitude, longitude, var, !!col) |> 
    mutate(station_id = trimws(station_id)) |> 
    filter(station_id %in% stns_keep) |> 
    pivot_wider(names_from = "var", values_from = quo_name(col)) |> 
    mutate(longitude = rotate_lon(longitude, from = "0/360"), area_swept_km2 = area_avg) |> 
    add_utm(utm_crs = "+proj=utm +zone=2 +datum=WGS84") |> 
    group_by(year) |> 
    mutate(cold_pool_2C = (sum(temp_bottom5m <= 2)/n())) |> 
    ungroup()
  
}

## Clean up station names, put parameter values in columns
ROMS_data <- hind |> clean_roms(val_raw) |> mutate(sim = "hindcast") |> filter(year <= 2020)

## Replace survey-replicated values with temp, oxygen + pH extracted separately
ROMS_data <- ROMS_data |> rows_update(drop_na(hind_srvy), by = c("year", "station_id"), unmatched = "ignore")

## Data for 2021-2022 for model validation
ROMS_new <- expand_grid(unique(ROMS_data[,c("station_id", "latitude", "longitude", "X", "Y")]), year = 2021:2022)
ROMS_new <- ROMS_new |> left_join(hind_srvy, by = c("year", "station_id")) |> mutate(sim = "hindcast")
ROMS_new <- ROMS_new |> drop_na() |> group_by(year) |> mutate(cold_pool_2C = sum(temp_bottom5m <= 2)/n())

ROMS_data <- bind_rows(ROMS_data, ROMS_new)

## Use actual survey lat/longs for surveyed stations
srvy <- read.csv(here("data", "trawl_surveys", "ebs_stations_by_year.csv"))
srvy <- srvy |> mutate(station_id = trimws(STATION)) |> group_by(year = YEAR, station_id) |> 
  summarize(latitude = mean(LATITUDE), longitude = mean(LONGITUDE), .groups = "drop")

ROMS_data <- ROMS_data |> rows_update(srvy, by = c("year", "station_id"), unmatched = "ignore") |> 
  add_utm(utm_crs = "+proj=utm +zone=2 +datum=WGS84")

# MIROC ---------------------------------------------------------------------------------

sims <- list()

## Read in survey-replicated MIROC SSP126 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_miroc_ssp126_BC_fut.Rdata"))

sims$SSP126_MIROC <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP126 MIROC")

## Read in survey-replicated MIROC SSP585 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_miroc_ssp585_BC_fut.Rdata"))

sims$SSP585_MIROC <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP585 MIROC")

## Read in survey-replicated MIROC RCP4.5 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP5_MIROC_rcp45_BC_fut.Rdata"))

sims$RCP45_MIROC <- fut |> clean_roms(val_delta) |> mutate(sim = "RCP45 MIROC")

# CESM ----------------------------------------------------------------------------------

## Read in survey-replicated CESM SSP126 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_cesm_ssp126_BC_fut.Rdata"))

sims$SSP126_CESM <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP126 CESM")

## Read in survey-replicated CESM SSP585 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_cesm_ssp585_BC_fut.Rdata"))

sims$SSP585_CESM <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP585 CESM")

## Read in survey-replicated MIROC RCP4.5 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP5_CESM_rcp45_BC_fut.Rdata"))

sims$RCP45_CESM <- fut |> clean_roms(val_delta) |> mutate(sim = "RCP45 CESM") |> filter(year < 2080)

# GFDL ----------------------------------------------------------------------------------

## Read in survey-replicated GFDL SSP126 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_gfdl_ssp126_BC_fut.Rdata"))

sims$SSP126_GFDL <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP126 GFDL")

## Read in survey-replicated GFDL SSP585 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP6_gfdl_ssp585_BC_fut.Rdata"))

sims$SSP585_GFDL <- fut |> clean_roms(val_delta) |> mutate(sim = "SSP585 GFDL")

## Read in survey-replicated MIROC RCP4.5 Projection
load(here("data", "BC_ACLIMsurveyrep", "ACLIMsurveyrep_B10K-K20P19_CMIP5_GFDL_rcp45_BC_fut.Rdata"))

sims$RCP45_GFDL <- fut |> clean_roms(val_delta) |> mutate(sim = "RCP45 GFDL")

# Join data -----------------------------------------------------------------------------

ebs_stations <- ROMS_data |> group_by(station_id) |> select(latitude, longitude, X, Y) |> summarise_all(mean) |> ungroup()

sims <- filter(do.call("rbind", sims), year >= 2023 & year <= 2094)

ROMS_full <- ROMS_data |> bind_rows(sims) |> arrange(year, sim, station_id)

# Depth data ----------------------------------------------------------------------------

## Read NOAA digital elevation model
bedrock <- get_bathymetry()

## Extract depth of digital elevation model at survey locations
ROMS_full$depth_m <- -round(st_extract(bedrock, as.matrix(ROMS_full[,c("longitude", "latitude")]))[[1]])

# Sediment grain size -------------------------------------------------------------------

## Read sediment grain size raster
phi <- get_sediment()

## Extract sediment grain size at survey locations
ROMS_full$phi <- st_extract(phi, as.matrix(ROMS_full[,c("longitude", "latitude")]))[[1]]

rm(hind, fut, bedrock, ROMS_data, ROMS_new, phi, sims, hind_srvy, srvy, stns_keep, area_avg)
