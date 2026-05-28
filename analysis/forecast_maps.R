
# Setup -------------------------------------------------------------

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "sf", "stars", "cowplot", "BeringSeaData")
sapply(pkgs, require, character.only = TRUE)

save_dir <- here("output", "figures")

slice_proj <- function(x, start = 1993, end = 2022) {
  
  time <- st_get_dimension_values(x, "time")
  year <- lubridate::year(time)
  year_slice <- which(year >= start & year <= end)
  x |> slice(year_slice, along = "time")
  
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
  forecasts[[i]] <- readRDS(file.path(sp_dirs[i], "forecast.rds"))
  
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
  forecasts[[i]] <- mutate(readRDS(file.path(sp_dirs[i], "forecast.rds")), log_biomass = log(biomass_fit))
  
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

# Derived quantities ------------------------------------------------

sp_levels <- c(
  "walleye pollock", "Pacific cod", "Pacific halibut", "arrowtooth flounder", 
  "northern rock sole", "red king crab", "snow crab"
)

grid <- readRDS(here("data", "mom6_forecast", "mom6_forecast.rds"))

derived <- list.files(file.path("output"), recursive = TRUE, full.names = TRUE, pattern = "derived_quantities.rds")
derived <- derived |> lapply(readRDS) |> setNames(basename(dirname(derived)))

# Center of gravity samples - forecast
cog_fcst <- derived |> lapply(\(x) x$cog) |> bind_rows(.id = "sp_bin")

# Center of gravity mean - hindcast
hindcasts_occ <- readRDS(here("output", "occurrence_mapdata.rds"))$hindcasts
hindcasts_biomass <- readRDS(here("output", "biomass_mapdata.rds"))$hindcasts
hindcasts <- mapply(stars:::c.stars, hindcasts_biomass, hindcasts_occ, SIMPLIFY = FALSE)

hindcast_df <- hindcasts |> 
  lapply(\(x) mutate(x, X = c(grid$X), Y = c(grid$Y))) |> 
  lapply(as.data.frame, add_coordinates = FALSE) |> 
  bind_rows(.id = "sp_bin") |> 
  drop_na()

cog_hindcast <- hindcast_df |> 
  select(sp_bin, X, Y, occurrence = p_occurrence, biomass = biomass_fit) |> 
  pivot_longer(-c(sp_bin, X, Y), names_to = "model", values_to = "value") |> 
  group_by(sp_bin, model) |>
  summarize(cog_E_km = weighted.mean(X, value), cog_N_km = weighted.mean(Y, value), .groups = "drop") |> 
  pivot_longer(c(cog_E_km, cog_N_km), names_to = "coord", values_to = "cog_mean", names_prefix = "cog_")

cog_fcst <- cog_fcst |> 
  left_join(cog_hindcast, by = c("sp_bin", "model", "coord")) |> 
  mutate(cog_anom = cog - cog_mean) |> 
  separate(sp_bin, into = c("species", "bin"), sep = "-") |> 
  mutate(
    species = factor(gsub("_", " ", species), levels = sp_levels), 
    bin = ifelse(is.na(bin), "all", bin), 
    coord = coord |> replace_values(
      "E_km" ~ "COG Eastings anomaly (km)", 
      "N_km" ~ "COG Northings anomaly (km)"
    )
  )

cog_fcst_summary <- cog_fcst |> 
  group_by(species, bin, model, coord) |> 
  summarize(
    est = median(cog_anom), 
    sd = sd(cog_anom), 
    q2.5 = quantile(cog_anom, 0.025), 
    q10 = quantile(cog_anom, 0.25),
    q90 = quantile(cog_anom, 0.75),
    q97.5 = quantile(cog_anom, 0.975), 
    .groups = "drop"
  )

# area occupied - forecast
area_fcst <- derived |> 
  lapply(\(x) x$area_occupied) |> 
  lapply(\(x) data.frame(sample = seq_along(x), area = x)) |> 
  bind_rows(.id = "sp_bin") |> 
  separate(sp_bin, into = c("species", "bin"), sep = "-") |> 
  mutate(
    species = factor(gsub("_", " ", species), levels = sp_levels), 
    bin = ifelse(is.na(bin), "all", bin), 
    model = "occurrence", 
    coord = "area occupied"
  )

area_fcst_summary <- area_fcst |> 
  group_by(species, bin, model, coord) |> 
  summarize(
    est = median(area), 
    q2.5 = quantile(area, 0.025), 
    q10 = quantile(area, 0.25),
    q90 = quantile(area, 0.75),
    q97.5 = quantile(area, 0.975), 
    .groups = "drop"
  )

# area occupied - hindcast
# since area occupied is derived from the posterior predictive distribution, 
# the mean is just the mean probability of occurrence
area_hindcast <- hindcast_df |> 
  summarize(est = mean(p_occurrence), .by = "sp_bin") |> 
  separate(sp_bin, into = c("species", "bin"), sep = "-") |> 
  mutate(species = factor(gsub("_", " ", species), levels = sp_levels), bin = ifelse(is.na(bin), "all", bin)) |> 
  mutate(model = "occurrence hindcast", coord = "area occupied")

derived_quant_plot <- cog_fcst_summary |> 
  bind_rows(area_fcst_summary) |> 
  ggplot(aes(est, bin, color = model)) + 
  geom_vline(aes(xintercept = 0), data = data.frame(coord = c("COG Eastings anomaly (km)", "COG Northings anomaly (km)"))) + 
  geom_violin(aes(x = cog_anom, fill = model), data = cog_fcst, scale = "area", show.legend = FALSE,
             color = NA, position = position_dodge(width = 0.8), alpha = 0.5) + 
  geom_violin(aes(x = area, fill = model), data = area_fcst, scale = "area", show.legend = FALSE,
             color = NA, position = position_dodge(width = 0.8), alpha = 0.5) + 
  geom_linerange(aes(xmin = q2.5, xmax = q97.5), linewidth = 0.5, lineend = "round",
                position = position_dodge(width = 0.8)) + 
  geom_linerange(aes(xmin = q10, xmax = q90), linewidth = 0.75, lineend = "round",
                position = position_dodge(width = 0.8)) + 
  geom_point(position = position_dodge(width = 0.8), shape = 21, size = 1, fill = "white", stroke = 1) + 
  geom_point(data = area_hindcast, position = position_nudge(y = -0.25)) +
  facet_grid(species~coord, scales = "free", space = "free_y", switch = "both") + 
  scale_color_manual(values = c("biomass" = "tomato2", "occurrence" = "dodgerblue", "occurrence hindcast" = "dodgerblue4")) +
  scale_fill_manual(values = c("biomass" = "tomato2", "occurrence" = "dodgerblue", "occurrence hindcast" = "dodgerblue4")) +
  expand_limits(x = c(0, 1)) +
  theme(
    strip.placement = "outside", 
    strip.background = element_blank(),
    axis.title = element_blank(), 
    strip.text.y.left = element_text(angle = 0, hjust = 1), 
    legend.position = "top", 
    panel.grid.minor = element_blank()
  )

ggsave(
  here("output", "derived_quantities.png"), derived_quant_plot,
  height = 5, width = 8, units = "in", dpi = 500
)

# Bhattacharyya similarity

forecasts <- readRDS(here("output", "biomass_mapdata.rds"))$forecasts

hb_arrays <- hindcasts |> lapply(\(x) c(x$biomass_prop))
fb_arrays <- forecasts |> lapply(\(x) c(x$biomass_prop))
bhatta_biomass <- mapply(\(x, y) sum(sqrt(x * y), na.rm = TRUE), hb_arrays, fb_arrays)

hp_arrays <- hindcasts |> lapply(\(x) c(x$p_occurrence)/sum(c(x$p_occurrence), na.rm = TRUE))
fp_arrays <- forecasts |> lapply(\(x) c(x$p_occurrence)/sum(c(x$p_occurrence), na.rm = TRUE))
bhatta_occurrence <- mapply(\(x, y) sum(sqrt(x * y), na.rm = TRUE), hp_arrays, fp_arrays)

bhatta <- data.frame(
  sp_bin = rep(names(forecasts), times = 2), 
  model = rep(c("occurrence", "biomass"), each = length(forecasts)), 
  coef = c(bhatta_occurrence, bhatta_biomass)
)

bhatta <- bhatta |> 
  separate(sp_bin, into = c("species", "bin"), sep = "-") |> 
  mutate(
    species = factor(gsub("_", " ", species), levels = sp_levels), 
    bin = ifelse(is.na(bin), "all", bin)
  ) 

bhatta_plot <- bhatta |> 
  ggplot(aes(coef, bin, fill = model, color = model)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.1) + 
  geom_point(position = position_dodge(width = 0.8), shape = 21, size = 1.5, stroke = 1.5, fill = "white") +
  facet_grid(species~., scales = "free_y", space = "free_y", switch = "y") + 
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_cartesian(clip = "off") +
  xlab("Bhattacharyya coefficient") +
  theme_minimal() +
  theme(
    strip.placement = "outside", 
    strip.background = element_blank(),
    axis.title.y = element_blank(), 
    axis.title.x = element_text(margin = margin(1, 0, 0, 0, "lines")),
    strip.text.y.left = element_text(angle = 0, hjust = 1), 
    legend.position = "top", 
    panel.grid.minor = element_blank()
  )

ggsave(
  here("output", "bhattacharyya.png"), bhatta_plot,
  height = 5, width = 5, units = "in", dpi = 500
)
