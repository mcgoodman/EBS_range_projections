
pkgs <- c("here", "dplyr", "tidyr")
sapply(pkgs, require, character.only = TRUE)

load(survey_file)

# Full list of standard hauls
srvy <- read.csv(here("data", "trawl_surveys", "ebs_stations_by_year.csv"))
srvy <- srvy |> dplyr::select(STATIONID = STATION, YEAR, HAUL)

# Drop hauls where species was observed but no length data was recorded (common earlier in the survey)
cpue_data <- cpue_data$CPUE_station_bin_yr |> 
  dplyr::select(YEAR, STATIONID, BIN, HAUL = HAULJOIN, haul_CPUE_KGKM2, bin_CPUE_KGKM2) |> 
  group_by(YEAR, STATIONID, HAUL) |> 
  filter(!all(is.na(bin_CPUE_KGKM2) & haul_CPUE_KGKM2 > 0)) |> 
  mutate(bin_CPUE_KGKM2 = ifelse(haul_CPUE_KGKM2 == 0, 0, bin_CPUE_KGKM2)) |> 
  right_join(srvy, by = c("STATIONID", "YEAR", "HAUL"))
  
## Bin data into juvenile and adult, aggregate catches by bin
## Repeat surveys of the same station in the same year are averaged
cpue_data <- cpue_data |> 
  mutate(
    station_id = factor(STATIONID, levels = ebs_stations$station_id), 
    bin = ifelse(BIN < (threshold * 10), "juvenile", "adult")
  ) |> 
  filter(!is.na(station_id)) |> 
  group_by(year = YEAR, station_id, bin, haul = HAUL) |> 
  summarize(
    cpue_kgkm2 = sum(bin_CPUE_KGKM2), 
    cpue_tot = unique(haul_CPUE_KGKM2), 
    .groups = "drop"
  ) |> 
  mutate(year_chr = factor(year, levels = sort(unique(year)))) |> 
  filter(!is.na(cpue_kgkm2))

## Remove years with no presences
drop_yrs <- cpue_data |> group_by(year) |> summarize(p = sum(cpue_kgkm2 == 0)/n())
drop_yrs <- drop_yrs$year[drop_yrs$p == 1]

## Merge survey data with hindcast ROMS-NPZ data
model_data <- cpue_data |>
  filter(!(year %in% drop_yrs) & year %in% years) |> 
  left_join(ROMS_data, by = c("station_id", "year")) |> 
  mutate(present = as.numeric(cpue_kgkm2 > 0)) |> 
  drop_na()

rm(cpue_data, srvy, drop_yrs)
