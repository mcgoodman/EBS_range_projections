
pkgs <- c("here", "dplyr", "tidyr")
sapply(pkgs, require, character.only = TRUE)

## Read in species data for EBS
cpue_data <- read.csv(survey_file)
cpue_data <- cpue_data |> dplyr::select(YEAR, STATION, HAUL, WTCPUE)

## Read in data containing stations surveyed in each year, to add missing zeroes and area swept
srvy <- read.csv("data/surveyrep_observed_1982-2022.csv")
srvy <- srvy |> mutate(area_swept_km2 = AREA_SWEPT_HA / 100) |> dplyr::select(YEAR, STATION, HAUL, area_swept_km2)

## Bin data into juvenile and adult, aggregate catches by bin
cpue_data <- cpue_data |>
  right_join(srvy, by = c("YEAR", "STATION", "HAUL")) |> 
  mutate(WTCPUE = ifelse(is.na(WTCPUE), 0, WTCPUE * 100)) |> 
  mutate(station_id = factor(STATION, levels = ebs_stations$station_id)) |> 
  filter(!is.na(station_id)) |> 
  rename(year = YEAR, haul = HAUL, cpue_kgkm2 = WTCPUE) |>
  mutate(year_chr = factor(year, levels = sort(unique(year))))

## Remove years with no presences
drop_yrs <- cpue_data |> group_by(year) |> summarize(p = sum(cpue_kgkm2 == 0)/n())
drop_yrs <- drop_yrs$year[drop_yrs$p == 1]

## Merge survey data with ROMS-NPZ data
joined_data <- cpue_data |>
  filter(!(year %in% drop_yrs)) |> 
  left_join(filter(dplyr::select(ROMS_full, -area_swept_km2), sim == "hindcast"), by = c("station_id", "year")) |> 
  mutate(present = as.numeric(cpue_kgkm2 > 0)) |> 
  filter(year <= 2019)

## Merge 2021-2022 survey data with ROMS-NPZ data
val_data <- cpue_data |> 
  filter(year %in% c(2021, 2022)) |> 
  left_join(dplyr::select(filter(ROMS_full, sim == "hindcast" & year %in% 2021:2022), -area_swept_km2), by = c("station_id", "year")) |> 
  mutate(present = as.numeric(cpue_kgkm2 > 0)) |> 
  drop_na()

rm(cpue_data, srvy, drop_yrs)
