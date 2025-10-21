
# Setup -------------------------------------------------------------

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "sf", "stars", "cowplot", "BeringSeaData")
sapply(pkgs, require, character.only = TRUE)

save_dir <- here("output", "figures")

slice_proj <- function(x, start = 1993, end = 2022) {
  
  time <- st_get_dimension_values(x, "ocean_time")
  year <- lubridate::year(time)
  year_slice <- which(year >= start & year <= end)
  x |> slice(year_slice, along = "ocean_time")
  
}

summarize_proj <- function(x, start = 1993, end = 2022, var = p_occurrence, f = mean) {
  
  var <- enquo(var)
  x <- x |> slice_proj(start, end) |> dplyr::select(!!var)
  x <- st_apply(x, 1:2, f)
  names(x) <- quo_name(var)
  x
  
}

# Probability of occurrence forecasts -------------------------------

## Read in hindcasts and forecasts, compute difference --------------

sp_dirs <- list.dirs(here("output"), full.names = TRUE)
sp_dirs <- sp_dirs[!(sp_dirs %in% c(save_dir, here("output"), here("output", "ESR_plots")))]
sp_bin <- lapply(gsub("_", " ", basename(sp_dirs)), \(x) strsplit(x, split = "-")[[1]])

hindcasts <- forecasts <- anomalies <- setNames(vector("list", length(sp_dirs)), basename(sp_dirs))

for (i in seq_along(sp_dirs)) {
  
  # Read in hindcast, compute average probability of occurrence for reference period
  hindcasts[[i]] <- summarize_proj(readRDS(file.path(sp_dirs[i], "hindcast_level2.rds")))
  
  # Read in forecast
  forecasts[[i]] <- readRDS(file.path(sp_dirs[i], "forecast_level2_adj.rds"))
  
  # Compute anomaly
  anomalies[[i]] <- forecasts[[i]] |> 
    mutate(p_occurrence_diff = pmin(pmax(p_occurrence - c(hindcasts[[i]]$p_occurrence), -0.33), 0.33)) |> 
    select(p_occurrence_diff)
  
}

# Read in Alaska coast, transform maps to UTM for plotting
ak_coast <- get_ak_coast()
hindcasts <- lapply(hindcasts, st_transform, crs = st_crs(ak_coast))
forecasts <- lapply(forecasts, st_transform, crs = st_crs(ak_coast))
anomalies <- lapply(anomalies, st_transform, crs = st_crs(ak_coast))
ak_coast <- filter(st_crop(ak_coast, hindcasts$`walleye_pollock-adult`), label != "Russia")
ebs <- get_ebs_shapefile()

saveRDS(list(hindcasts = hindcasts, forecasts = forecasts, anomalies = anomalies), here("output", "occurrence_mapdata.rds"))

## Plot -------------------------------------------------------------

for (i in seq_along(sp_dirs)) {
  
  hindcast_plot <- ggplot() + 
    geom_stars(aes(fill = p_occurrence, color = p_occurrence), data = hindcasts[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_viridis_c(option = "magma", na.value = NA, limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
    scale_fill_viridis_c(option = "magma", na.value = NA, limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar("1993-2022 P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar("1993-2022 P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  forecast_plot <- ggplot() + 
    geom_stars(aes(fill = p_occurrence, color = p_occurrence), data = forecasts[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_viridis_c(option = "magma", na.value = NA, limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
    scale_fill_viridis_c(option = "magma", na.value = NA, limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar("Forecast P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar("Forecast P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  anomaly_plot <- ggplot() + 
    geom_stars(aes(fill = p_occurrence_diff, color = p_occurrence_diff), data = anomalies[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_gradient2(low = "#08306B", high = "#67000D", na.value = NA, limits = c(-0.33, 0.33), breaks = seq(-0.33, 0.33, 0.11)) + 
    scale_fill_gradient2(low = "#08306B", high = "#67000D", na.value = NA, limits = c(-0.33, 0.33), breaks = seq(-0.33, 0.33, 0.11)) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar(expression(Delta~"P(occurrence)"), title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar(expression(Delta~"P(occurrence)"), title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  joined <- cowplot::plot_grid(
    hindcast_plot, forecast_plot, anomaly_plot, nrow = 1
  ) + theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(sp_dirs[i], "occurrence_forecast.png"), joined, height = 4, width = 10, units = "in", dpi = 300)
  
}

# Relative biomass forecasts ----------------------------------------

## Read in hindcasts and forecasts, compute difference --------------

hindcasts <- forecasts <- anomalies <- setNames(vector("list", length(sp_dirs)), basename(sp_dirs))

for (i in seq_along(sp_dirs)) {
  
  # Read in hindcast, compute average probability of occurrence for reference period
  hindcasts[[i]] <- readRDS(file.path(sp_dirs[i], "hindcast_level2.rds"))
  hindcasts[[i]] <- c(
    summarize_proj(mutate(hindcasts[[i]], log_biomass = log(biomass_fit)), var = log_biomass), 
    summarize_proj(hindcasts[[i]], var = biomass_fit)
  )
  
  # Read in forecast
  forecasts[[i]] <- mutate(readRDS(file.path(sp_dirs[i], "forecast_level2_adj.rds")), log_biomass = log(biomass_fit))
  
  # Center both, to compute relative biomass distribution anomaly
  hindcasts[[i]]$log_biomass_ctr <- hindcasts[[i]]$log_biomass - mean(hindcasts[[i]]$log_biomass, na.rm = TRUE)
  forecasts[[i]]$log_biomass_ctr <- forecasts[[i]]$log_biomass - mean(forecasts[[i]]$log_biomass, na.rm = TRUE)
  
  # Biomass proportion
  hindcasts[[i]]$biomass_prop <- hindcasts[[i]]$biomass_fit/sum(hindcasts[[i]]$biomass_fit, na.rm = TRUE)
  forecasts[[i]]$biomass_prop <- forecasts[[i]]$biomass_fit/sum(forecasts[[i]]$biomass_fit, na.rm = TRUE)
  
  # Compute anomaly and anomaly Z-score
  anomalies[[i]] <- forecasts[[i]] |> mutate(
    #log_biomass_diff = pmax(pmin(log_biomass_ctr - c(hindcasts[[i]]$log_biomass_ctr), 3), -3),
    log_biomass_diff = log_biomass_ctr - c(hindcasts[[i]]$log_biomass_ctr),
    log_biomass_zscore = pmax(pmin(log_biomass_diff / biomass_se, 3), -3), 
    biomass_prop_diff = biomass_prop - c(hindcasts[[i]]$biomass_prop)
  )
  
}

# Transform maps to UTM for plotting
hindcasts <- lapply(hindcasts, st_transform, crs = st_crs(ak_coast))
forecasts <- lapply(forecasts, st_transform, crs = st_crs(ak_coast))
anomalies <- lapply(anomalies, st_transform, crs = st_crs(ak_coast))

saveRDS(list(hindcasts = hindcasts, forecasts = forecasts, anomalies = anomalies), here("output", "biomass_mapdata.rds"))

## Plot -------------------------------------------------------------

for (i in seq_along(sp_dirs)) {
  
  plot_range <- c(0, max(c(hindcasts[[i]]$biomass_prop), c(forecasts[[i]]$biomass_prop)))
  
  hindcast_plot <- ggplot() + 
    geom_stars(aes(fill = biomass_prop, color = biomass_prop), data = hindcasts[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_viridis_c(option = "magma", na.value = NA, limits = plot_range) + 
    scale_fill_viridis_c(option = "magma", na.value = NA, limits = plot_range) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar("1993-2022 biomass proportion", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar("1993-2022 biomass proportion", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  forecast_plot <- ggplot() + 
    geom_stars(aes(fill = biomass_prop, color = biomass_prop), data = forecasts[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_viridis_c(option = "magma", na.value = NA, limits = plot_range) + 
    scale_fill_viridis_c(option = "magma", na.value = NA, limits = plot_range) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar("Forecast biomass proportion", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar("Forecast biomass proportion", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  anomaly_plot <- ggplot() + 
    geom_stars(aes(fill = biomass_prop_diff, color = biomass_prop_diff), data = anomalies[[i]]) +
    geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
    geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) +
    scale_color_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
    scale_fill_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    coord_sf(expand = FALSE) + 
    theme_void() + 
    theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
    guides(
      color = guide_colorbar(expression(Delta~"proportion biomass"), title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar(expression(Delta~"proportion biomass"), title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
    )
  
  joined <- cowplot::plot_grid(
    hindcast_plot, forecast_plot, anomaly_plot, nrow = 1
  ) + theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(sp_dirs[i], "biomass_forecast.png"), joined, height = 4, width = 10, units = "in", dpi = 300)
  
}

# MOM6 / ROMS anomaly -----------------------------------------------

roms_hindcast <- readRDS(here("data", "roms_level2_bc_annual", "CORECFS_hindcast.rds"))
mom6_forecast <- readRDS(here("data", "mom6", "mom6_forecast_adj.rds"))

roms_hindcast <- c(
  summarize_proj(roms_hindcast, var = temp_bottom5m), 
  summarize_proj(roms_hindcast, var = oxygen_bottom5m), 
  summarize_proj(roms_hindcast, var = pH_bottom5m)
)

mom6_anomaly <- mom6_forecast |> mutate(
  temp_diff = temp_bottom5m - c(roms_hindcast$temp_bottom5m),
  oxygen_diff = oxygen_bottom5m - c(roms_hindcast$oxygen_bottom5m), 
  pH_diff = pH_bottom5m - c(roms_hindcast$pH_bottom5m)
)

mom6_anomaly <- mom6_anomaly |> st_transform(crs = st_crs(ak_coast))

temp_plot <- ggplot() + 
  geom_stars(aes(fill = temp_diff, color = temp_diff), data = mom6_anomaly) + 
  geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
  geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) + 
  scale_color_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  scale_fill_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  coord_sf(expand = FALSE) + 
  theme_void() + 
  theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
  guides(
    color = guide_colorbar("temperature anomaly", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
    fill = guide_colorbar("temperature anomaly", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
  )

oxygen_lab <- expression("oxygen anomaly mmol"~m^{-3})

oxygen_plot <- ggplot() + 
  geom_stars(aes(fill = oxygen_diff, color = oxygen_diff), data = mom6_anomaly) + 
  geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
  geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) + 
  scale_color_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  scale_fill_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  coord_sf(expand = FALSE) + 
  theme_void() + 
  theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
  guides(
    color = guide_colorbar(oxygen_lab, title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
    fill = guide_colorbar(oxygen_lab, title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
  )

pH_plot <- ggplot() + 
  geom_stars(aes(fill = pH_diff, color = pH_diff), data = mom6_anomaly) + 
  geom_sf(data = ak_coast, fill = "grey80", color = NA) + 
  geom_sf(data = ebs, fill = NA, color = "black", linewidth = 0.5) + 
  scale_color_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  scale_fill_gradient2(low = "#08306B", high = "#67000D", na.value = NA) + 
  coord_sf(expand = FALSE) + 
  theme_void() + 
  theme(legend.position = "top", plot.background = element_rect(fill = "white", color = NA)) + 
  guides(
    color = guide_colorbar("pH anomaly", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
    fill = guide_colorbar("pH anomaly", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
  )

joined <- cowplot::plot_grid(
  temp_plot, oxygen_plot, pH_plot, nrow = 1
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here("output", "roms_anomalies.png"), joined, height = 4, width = 10, units = "in", dpi = 300)
