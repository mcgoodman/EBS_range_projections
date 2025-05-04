
# Setup -------------------------------------------------------------

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "sf", "stars", "cowplot", "Bering10KThredds")
sapply(pkgs, require, character.only = TRUE)

save_dir <- here("output", "figures")

slice_proj <- function(x, start = 1995, end = 2014) {
  
  time <- st_get_dimension_values(x, "ocean_time")
  year <- lubridate::year(time)
  year_slice <- which(year >= start & year <= end)
  x |> slice(year_slice, along = "ocean_time")
  
}

summarize_proj <- function(x, start = 1995, end = 2014, var = p_occurrence, f = mean) {
  
  var <- enquo(var)
  x <- x |> slice_proj(start, end) |> dplyr::select(!!var)
  x <- st_apply(x, 1:2, f)
  names(x) <- quo_name(var)
  x
  
}

# Read in hindcasts and forecasts, compute difference ---------------

sp_dirs <- list.files(here("output"), full.names = TRUE)
sp_dirs <- sp_dirs[sp_dirs != save_dir]
sp_bin <- lapply(gsub("_", " ", basename(sp_dirs)), \(x) strsplit(x, split = "-")[[1]])

hindcasts <- forecasts <- anomalies <- setNames(vector("list", length(sp_dirs)), basename(sp_dirs))

for (i in seq_along(sp_dirs)) {
  
  # Read in hindcast, compute average probability of occurrence for reference period
  hindcasts[[i]] <- summarize_proj(readRDS(file.path(sp_dirs[i], "hindcast_level2.rds")))
  
  # Read in forecast
  forecasts[[i]] <- readRDS(file.path(sp_dirs[i], "forecast_level2.rds"))
  
  # Compute anomaly and anomaly Z-score
  anomalies[[i]] <- forecasts[[i]] |> 
    mutate(p_occurrence_diff = pmin(pmax(p_occurrence - c(hindcasts[[i]]$p_occurrence), -0.33), 0.33)) |> 
    select(p_occurrence_diff)
  
}

# Plot --------------------------------------------------------------

# Read in Alaska coast, transform maps to UTM for plotting
ak_coast <- get_ak_coast()
hindcasts <- lapply(hindcasts, st_transform, crs = st_crs(ak_coast))
forecasts <- lapply(forecasts, st_transform, crs = st_crs(ak_coast))
anomalies <- lapply(anomalies, st_transform, crs = st_crs(ak_coast))
ak_coast <- filter(st_crop(ak_coast, hindcasts$`walleye_pollock-adult`), is.na(DESC_))
ebs <- get_ebs_shapefile()

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
      color = guide_colorbar("1995-2014 P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"), 
      fill = guide_colorbar("1995-2014 P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
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
