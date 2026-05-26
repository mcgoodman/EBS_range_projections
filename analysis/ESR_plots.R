
# Setup -------------------------------------------------------------

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "sf", "stars", "cowplot", "BeringSeaData")
sapply(pkgs, require, character.only = TRUE)

save_dir <- here("output", "figures")

hex_ramp <- function(x, diverging = FALSE) {
  
  cols <- rep(NA_character_, length(x))
  
  if(!diverging) {
    xt <- (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    cols[!is.na(xt)] <- rgb(colorRamp(viridis::magma(20))(xt[!is.na(xt)]), maxColorValue = 255)
  } else {
    xt <- abs(x)/(max(abs(x), na.rm = TRUE))
    cols[!is.na(xt) & x > 0] <- rgb(colorRamp(RColorBrewer::brewer.pal(9, "Reds"))(xt[!is.na(xt) & x > 0]), maxColorValue = 255)
    cols[!is.na(xt) & x <= 0] <- rgb(colorRamp(RColorBrewer::brewer.pal(9, "Blues"))(xt[!is.na(xt) & x <= 0]), maxColorValue = 255)
  }
  
  cols
}

# Probability of occurrence maps ------------------------------------

## Reformat ---------------------------------------------------------

occurrence <- readRDS(here("output", "occurrence_mapdata.rds"))
occurrence$forecasts <- occurrence$forecasts |> lapply(\(x) dplyr::select(x, p_occurrence))

occurrence <- lapply(occurrence, \(x) {
  out <- lapply(names(x), \(y) {names(x[[y]]) <- y; x[[y]]})
  out <- merge(Reduce("c", out), name = "species_bin")
  names(out) <- names(x[[1]]) 
  out$color <- hex_ramp(out[[1]], diverging = grepl("diff", names(out)))
  out
})

occurrence <- merge(Reduce("c", lapply(occurrence, \(x) select(x, color))), name = "category")
names(occurrence) <- "color"
occurrence <- st_set_dimensions(occurrence, 4, values = c("Hindcast (1993-2022)", "Forecast", "Anomaly"))

sp_bins <- st_get_dimension_values(occurrence, "species_bin")
sp_formatted <- gsub(" juvenile", "\n(juvenile)", gsub(" adult", "\n(adult)", gsub("_|-", " ", sp_bins)))
occurrence <- st_set_dimensions(occurrence, 3, values = sp_formatted, name = "species_bin")

## Legend -----------------------------------------------------------

occurrence_seq <- seq(0, 1, 0.1)
occurrence_legend_plot <- data.frame(x = occurrence_seq) |> 
  ggplot(aes(x, x, color = x)) + 
  geom_point() + 
  scale_color_viridis_c(option = "magma") + 
  guides(
    color = guide_colorbar("P(occurrence)", title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
  ) + 
  theme(legend.position = "top")

occurrence_legend <- cowplot::get_legend(occurrence_legend_plot)

anomaly_seq <- seq(-0.33, 0.33, 0.01)
anomaly_legend_plot <- data.frame(x = anomaly_seq) |> 
  ggplot(aes(x, x, color = x)) + 
  geom_point() + 
  scale_color_gradientn(
    colors = hex_ramp(anomaly_seq, diverging = TRUE), 
    values = scales::rescale(anomaly_seq, from = c(-0.33, 0.33), to = c(0, 1))
  ) + 
  guides(
    color = guide_colorbar(expression(Delta~"P(occurrence)"), title.position = "top", title.hjust = 0.5, barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black")
  ) + 
  theme(legend.position = "top")

anomaly_legend <- cowplot::get_legend(anomaly_legend_plot)

blank <- ggplot() + geom_blank() + theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA))

joined_legend <- cowplot::plot_grid(blank, occurrence_legend, anomaly_legend, blank, nrow = 1, rel_widths = c(0.2, 0.3, 0.3, 0.2))

## Map plot ---------------------------------------------------------

ebs <- get_ebs_shapefile()
coast <- st_crop(get_ak_coast(), st_buffer(ebs, 1e5))

occurrence_maps <- ggplot() + 
  geom_sf(data = coast, fill = "grey80", color = "grey60") + 
  geom_stars(
    aes(fill = color, color = color), 
    data = occurrence |> 
      slice(
        which(sp_formatted %in% c("arrowtooth flounder\n(adult)", "red king crab", "snow crab")), 
        along = "species_bin"
      )
  ) + 
  geom_sf(data = ebs, fill = NA, color = "grey40", linewidth = 0.5) +
  facet_grid(species_bin ~ category, switch = "y") + 
  scale_fill_identity() + 
  scale_color_identity() + 
  theme_minimal() +
  coord_sf(expand = FALSE, clip = "off") + 
  scale_x_continuous(breaks = seq(-180, -150, 10)) + 
  theme(panel.border = element_rect(fill = NA, color = "grey60"), 
        panel.spacing = unit(0, "in"), 
        strip.placement = "outside", 
        strip.text = element_text(size = 12, face = "bold"))

occurrence_maps <- cowplot::plot_grid(
  occurrence_maps, joined_legend, ncol = 1, rel_heights = c(0.9, 0.1), rel_widths = c(1, 0.5)
) + theme(plot.background = element_rect(fill = "white", color = NA))

dir.create(here("output", "ESR_plots"))

ggsave(here("output", "ESR_plots", "occurrence_forecasts.png"), occurrence_maps, height = 9, width = 9, units = "in", dpi = 300)

# Biomass maps ------------------------------------------------------

## Reformat ---------------------------------------------------------

biomass <- readRDS(here("output", "biomass_mapdata.rds"))
biomass$hindcasts <- biomass$hindcasts |> lapply(\(x) dplyr::select(x, biomass_prop))
biomass$forecasts <- biomass$forecasts |> lapply(\(x) dplyr::select(x, biomass_prop))
biomass$anomalies <- biomass$anomalies |> lapply(\(x) dplyr::select(x, biomass_prop_diff))

sp_bins <- c("walleye_pollock-adult", "Pacific_cod-adult")
sp_formatted <- gsub(" juvenile", "\n(juvenile)", gsub(" adult", "\n(adult)", gsub("_|-", " ", sp_bins)))

biomass <- lapply(biomass, \(x) {
  out <- x[sp_bins]
  out <- lapply(names(out), \(y) {names(out[[y]]) <- y; out[[y]]})
  out <- merge(Reduce("c", out), name = "species_bin")
  names(out) <- names(x[[1]]) 
  out$color <- hex_ramp(out[[1]], diverging = grepl("diff", names(out)))
  out <- st_set_dimensions(out, which = 3, values = sp_formatted, names = "species_bin")
  out
})

prop_range <- range(c(biomass$hindcasts$biomass_prop, biomass$forecasts$biomass_prop), na.rm = TRUE)
diff_range <- range(biomass$anomalies$biomass_prop_diff, na.rm = TRUE)

biomass <- merge(Reduce("c", lapply(biomass, \(x) select(x, color))), name = "category")
names(biomass) <- "color"
biomass <- st_set_dimensions(biomass, 4, values = c("Hindcast (1993-2022)", "Forecast", "Anomaly"))

## Legend -----------------------------------------------------------

biomass_seq <- seq(0, prop_range[2], length.out = 100)
biomass_legend_plot <- data.frame(x = biomass_seq) |> 
  ggplot(aes(x, x, color = x)) + 
  geom_point() + 
  scale_color_viridis_c(option = "magma") + 
  guides(
    color = guide_colorbar(
      "biomass proportion", title.position = "top", title.hjust = 0.5, 
      barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"
    )
  ) + 
  theme(legend.position = "top")

biomass_legend <- cowplot::get_legend(biomass_legend_plot)

anomaly_seq <- seq(diff_range[1], diff_range[2], length.out = 100)
anomaly_legend_plot <- data.frame(x = anomaly_seq) |> 
  ggplot(aes(x, x, color = x)) + 
  geom_point() + 
  scale_color_gradientn(
    colors = hex_ramp(anomaly_seq, diverging = TRUE), 
    values = scales::rescale(anomaly_seq, from = diff_range, to = c(0, 1))
  ) + 
  guides(
    color = guide_colorbar(
      expression(Delta~"biomass proportion"), title.position = "top", title.hjust = 0.5, 
      barwidth = unit(12, "lines"), ticks.colour = "black", frame.colour = "black"
    )
  ) + 
  theme(legend.position = "top")

anomaly_legend <- cowplot::get_legend(anomaly_legend_plot)

blank <- ggplot() + geom_blank() + theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA))

joined_legend <- cowplot::plot_grid(blank, biomass_legend, anomaly_legend, blank, nrow = 1, rel_widths = c(0.2, 0.3, 0.3, 0.2))

## Map plot ---------------------------------------------------------

biomass_maps <- ggplot() + 
  geom_sf(data = coast, fill = "grey80", color = "grey60") + 
  geom_stars(aes(fill = color, color = color), data = biomass) + 
  geom_sf(data = ebs, fill = NA, color = "grey40", linewidth = 0.5) +
  facet_grid(species_bin ~ category, switch = "y") + 
  scale_fill_identity() + 
  scale_color_identity() +
  theme_minimal() +
  coord_sf(expand = FALSE, clip = "off") + 
  scale_x_continuous(breaks = seq(-180, -150, 10)) + 
  theme(panel.border = element_rect(fill = NA, color = "grey60"), 
        panel.spacing = unit(0, "in"), 
        strip.placement = "outside", 
        strip.text = element_text(size = 12, face = "bold"))

biomass_maps <- cowplot::plot_grid(
  biomass_maps, joined_legend, ncol = 1, rel_heights = c(0.85, 0.15), rel_widths = c(1, 0.5)
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(here("output", "ESR_plots", "biomass_forecasts.png"), biomass_maps, height = 6.5, width = 9, units = "in", dpi = 300)

# Time series -------------------------------------------------------

sp_bins <- c(
  "walleye_pollock-adult", "walleye_pollock-juvenile", 
  "Pacific_cod-adult", "Pacific_cod-juvenile",
  "arrowtooth_flounder-adult", "arrowtooth_flounder-juvenile", 
  "snow_crab", "red_king_crab"
)
sp_formatted <- gsub(" juvenile", "\n(juvenile)", gsub(" adult", "\n(adult)", gsub("_|-", " ", sp_bins)))

hindcasts <- lapply(sp_bins, \(x) readRDS(
  here("output", x, "hindcast_level2.rds")
))

forecasts <- lapply(sp_bins, \(x) readRDS(
  here("output", x, "forecast_level2_adj.rds")
))

forecast_means <- forecasts |> 
  setNames(sp_bins) |> 
  lapply(as.data.frame) |> 
  bind_rows(.id = "species_bin") |> 
  group_by(species_bin) |> 
  summarize(across(everything(), \(x) sum(x > 0.5, na.rm = TRUE)/sum(!is.na(x)))) |> 
  mutate(year = 2025) |> 
  select(-eta_rho, -xi_rho)

hindcast_ts <- hindcasts |> 
  setNames(sp_bins) |> 
  lapply(st_apply, 3, \(x) sum(x > 0.5, na.rm = TRUE)/sum(!is.na(x))) |> 
  lapply(as.data.frame) |> 
  bind_rows(.id = "species_bin") |> 
  mutate(year = lubridate::year(time)) |> 
  filter(as.character(time) != "1976-07-04 11:56:00") |> 
  select(-time) |> 
  bind_rows(forecast_means) |> 
  separate("species_bin", into = c("species", "bin"), remove = FALSE, sep = "-", fill = "right") |> 
  mutate(
    species = gsub("_", " ", species), bin = ifelse(!is.na(bin), bin, "all"), 
    species = factor(species, levels = c(
      "walleye pollock", "Pacific cod", "arrowtooth flounder", 
      "red king crab", "snow crab"
    )), 
    bin = factor(bin, levels = c("adult", "juvenile", "all"))
  )

ts_plot <- hindcast_ts |> 
  filter(year != 2025) |> 
  ggplot(aes(year, p_occurrence, color = bin)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~species, ncol = 1, scales = "free_y") + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 50, bs = "gp", m = c(3, 1)), 
              method.args = list(family = gaussian(link = "logit")), 
              se = FALSE, n = 300, alpha = 0.9) + 
  geom_point(data = hindcast_ts |> filter(year == 2025), size = 2) + 
  scale_color_manual(values = c("dodgerblue3", "firebrick3", "black")) + 
  scale_x_continuous(breaks = seq(1970, 2025, 5)) +
  labs(y = "mean probability of occurrence", color = "life stage") +
  theme_bw() + 
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(size = 12, face = "bold", hjust = 0)
  )
  
ggsave(here("output", "ESR_plots", "time_series_plot.png"), ts_plot, height = 6.5, width = 9, units = "in", dpi = 300)
