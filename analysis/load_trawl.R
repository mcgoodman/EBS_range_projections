
pkgs <- c("here", "dplyr", "tidyr")
sapply(pkgs, require, character.only = TRUE)

## Read in species data for EBS
cpue_data <- read.csv(survey_file)
cpue_data <- cpue_data |> dplyr::select(YEAR, STATION, WTCPUE)

## Read in data containing stations surveyed in each year, to add missing zeroes and area swept
srvy <- read.csv(here("data", "trawl_surveys", "ebs_stations_by_year.csv"))
srvy <- srvy |> group_by(YEAR, STATION) |> summarize(area_swept_km2 = mean(AREA_SWEPT_HA) / 100, .groups = "drop")

## Bin data into juvenile and adult, aggregate catches by bin
cpue_data <- cpue_data |>
  right_join(srvy, by = c("YEAR", "STATION")) |> 
  mutate(WTCPUE = ifelse(is.na(WTCPUE), 0, WTCPUE * 100)) |> 
  mutate(station_id = factor(STATION, levels = ebs_stations$station_id)) |> 
  filter(!is.na(station_id)) |> 
  rename(year = YEAR, cpue_kgkm2 = WTCPUE) |>
  mutate(year_chr = factor(year, levels = sort(unique(year))))

## Remove years with no presences
drop_yrs <- cpue_data |> group_by(year) |> summarize(p = sum(cpue_kgkm2 == 0)/n())
drop_yrs <- drop_yrs$year[drop_yrs$p == 1]

## Merge survey data with ROMS-NPZ data
model_data <- cpue_data |>
  filter(!(year %in% drop_yrs)) |> 
  left_join(dplyr::select(ROMS_data, -area_swept_km2, -sampled), by = c("station_id", "year")) |> 
  mutate(present = as.numeric(cpue_kgkm2 > 0))

## Merge 2021-2022 survey data with ROMS-NPZ data

rm(cpue_data, srvy, drop_yrs)
