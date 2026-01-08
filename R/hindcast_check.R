
library("tidyverse")
library("stars")
library("BeringSeaData")
library("units")

sp_bin <- "pacific_cod-adult"
species <- gsub("_", " ", strsplit(sp_bin, "-")[[1]][1])
bin <- strsplit(sp_bin, "-")[[1]][2]
if (is.na(bin)) bin <- "" else bin <- paste0(bin, " ")

loc <- "St. Paul"
# center <- c(-165.44, 64.5) # Nome
center <- c(-170.28, 57.16)

hindcast <- readRDS(file.path("output", sp_bin, "level2_projections/CORECFS_hindcast.rds"))

annual_means <- hindcast |> 
  st_apply(3, mean, na.rm = TRUE) |> 
  as.data.frame() |> 
  mutate(year = lubridate::year(ocean_time))

annual_means |> 
  ggplot(aes(year, biomass_fit)) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(1982, 2022, 2))

ebs <- get_ebs_shapefile("EBS", type = "boundary")

center <- st_sfc(st_point(center, dim = "XY"), crs = 4326) # Nome
center <- st_transform(center, st_crs(ebs))

buffers_nm <- as_units(c(25, 50, 75, seq(100, 500, 100)), "nautical_mile")
buffers_chr <- paste0(c(0, lag(drop_units(buffers_nm))[-1]), "-", drop_units(buffers_nm), "nm")

buffers <- buffers_nm |> 
  lapply(\(x) st_as_sf(st_buffer(center, x))) |> 
  bind_rows() |> 
  st_difference() |> 
  st_intersection(ebs) |> 
  mutate(dist = buffers_chr)

hindcast <- st_transform(hindcast, st_crs(ebs))

years <- lubridate::year(st_get_dimension_values(hindcast, "ocean_time"))
buffer_means <- setNames(vector("list", length(years)), years)

for (i in seq_along(years)) {
  
  buffer_means[[i]] <- hindcast |> 
    slice(i, along = "ocean_time") |> 
    aggregate(by = buffers, mean, na.rm = TRUE) |> 
    as.data.frame() |> 
    mutate(year = years[i])
  
}

buffer_means <- buffer_means |> 
  bind_rows() |> 
  left_join(mutate(buffers, geometry = x)) |> 
  mutate(dist = factor(dist, levels = buffers_chr))

buffer_means |> 
  ggplot(aes(year, biomass_fit)) + 
  geom_line(aes(color = dist)) + 
  scale_color_viridis_d(option = "turbo") + 
  scale_x_continuous(breaks = seq(1982, 2022, 4)) + 
  labs(title = paste0(bin, species, ", ", loc), y = "mean biomass")

ggsave("~/Downloads/snow_crab_hindcast.png", height = 4, width = 8, units = "in", dpi = 500)
