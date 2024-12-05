
pkgs <- c("here", "dplyr", "ggplot2", "sf", "aclim2sdms")

sapply(pkgs, require, character.only = TRUE)

theme_set(
  theme_void() +
    theme(
      strip.background = element_blank(), 
      legend.position = "bottom", 
      strip.text = element_text(size = 14), 
      plot.background = element_rect(color = NA, fill = "white")
    )
)

slice_proj <- function(x, start = 2080, end = 2099) {
  
  time <- st_get_dimension_values(x, "ocean_time")
  year <- lubridate::year(time)
  year_slice <- which(year >= start & year <= end)
  x |> slice(year_slice, along = "ocean_time")
  
}

summarize_proj <- function(x, start = 2080, end = 2099, var = p_occurrence, f = mean) {
  
  var <- enquo(var)
  x <- x |> slice_proj(start, end) |> dplyr::select(!!var)
  x <- st_apply(x, 1:2, f)
  names(x) <- quo_name(var)
  x
  
}

# Read in spatial data --------------------------------------------------------

ebs_grid <- get_ebs_shapefile(region = "EBS", type = "grid") |> dplyr::select(station_id, geometry)

coast <- get_ak_coast() |> st_crop(ebs_grid)

# Functions to plot fitted values ---------------------------------------------

## Maps of fitted vs. observed values
fitted_obs_map <- function(plot_data, type = c("biomass", "occurrence")) {
  
  type <- match.arg(type) 
  
  type |> switch(
    biomass = plot_data |> 
      filter(sim == "hindcast") |> 
      pivot_longer(c(fitted, observed)) |> 
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = coast, inherit.aes = FALSE, fill = "grey80", color = "grey80") +
      geom_sf(aes(fill = log(value + 1), color = log(value + 1))) +
      facet_grid(name~year, switch = "y") + 
      scale_color_viridis_c(option = "magma", na.value = "grey60") +
      scale_fill_viridis_c(option = "magma", na.value = "grey60") +
      guides(
        color = guide_colorbar("log(Biomass CPUE + 1)", title.hjust = 0.5, title.position = "top", barwidth = unit(15, "lines")), 
        fill = guide_colorbar("log(Biomass CPUE + 1)")
      ), 
    occurrence = plot_data |> 
      filter(sim == "hindcast") |> 
      pivot_longer(c(fitted, observed)) |> 
      ggplot(aes(geometry = geometry)) +
      geom_sf(data = coast, inherit.aes = FALSE, fill = "grey80", color = "grey80") +
      geom_sf(aes(fill = value, color = value)) +
      facet_grid(name~year, switch = "y") + 
      scale_color_viridis_c(option = "magma", na.value = "grey60") +
      scale_fill_viridis_c(option = "magma", na.value = "grey60") +
      guides(
        color = guide_colorbar("probability of occurrence", title.hjust = 0.5, title.position = "top", barwidth = unit(15, "lines")), 
        fill = guide_colorbar("probability of occurrence")
      )
  )
  
}

# Summer probability of occurrence maps -------------------------------------------------

yr_start <- c(2040, 2060, 2080)
yr_end <- c(2059, 2079, 2099)

sp_paths <- list.dirs(here("output"), full.names = FALSE, recursive = FALSE)

for (i in 1:length(sp_paths)) {
  
  print(sp_paths[i])
  
  sim_files <- list.files(here("output", sp_paths[1], "/level2_projections"))
  
  # Summarize hindcast 1995-2015
  hindcast <- readRDS(here("output", sp_paths[i], "/level2_projections/CORECFS_hindcast.rds"))
  hindcast <- c(summarize_proj(hindcast, start = 1995, end = 2015), summarize_proj(hindcast, var = biomass_fit, start = 1995, end = 2015))
  
  # Read in projections
  sim_files <- sim_files[grepl("SSP", sim_files)]
  sim_names <- gsub(".rds", "", sim_files)
  sims <- lapply(sim_files, \(x) readRDS(here("output", sp_paths[i], "/level2_projections/", x)))
  
  # Variables for use in redimensioning projections
  index <- 1
  template <- sims[[1]][,,,1] |> dplyr::select(p_occurrence) |> st_set_dimensions(which = 3, names = "id", values = 1)
  
  # List to store projections for each time period
  sims_mean <- vector("list", length(yr_start))
  
  # Loop over time periods, store mean for each SSP & ESM
  for (j in 1:length(yr_start)) {
    
    sims_jp <- sims |> lapply(summarize_proj, start = yr_start[j], end = yr_end[j])
    sims_jb <- sims |> lapply(summarize_proj, var = biomass_fit, start = yr_start[j], end = yr_end[j])
    sims_j <- purrr::map2(sims_jb, sims_jp, c)
    
    sims_j <- sims_j |> purrr::map2(sim_names, \(x, y) {
      template <- template |> mutate(
        p_occurrence = c(x$p_occurrence), 
        biomass_fit = c(x$biomass_fit),
        ESM = toupper(vapply(y, \(z) strsplit(z, "_")[[1]][1], character(1))), 
        SSP = vapply(y, \(z) strsplit(z, "_")[[1]][2], character(1)), 
        Period = paste(yr_start[j], yr_end[j], sep = "-")
      )
      template <- template |> st_set_dimensions(which = 3, names = "id", values = index)
      index <<- index + 1
      template
    })
    
    sims_mean[[j]] <- Reduce(\(x, y) c(x, y, along = 3), sims_j)
    
  }
  
  # Join projections and hindcast
  sims_mean <- Reduce(\(x, y) c(x, y, along = 3), sims_mean)
  
  hindcast <- template |> 
    st_set_dimensions(which = 3, names = "id", values = index) |> 
    mutate(p_occurrence = c(hindcast$p_occurrence), biomass_fit = c(hindcast$biomass_fit), SSP = "hindcast", ESM = "", Period = "")
  
  sims <- c(sims_mean, hindcast, along = 3)
  
  sims <- sims |> st_transform(st_crs(ebs_grid))
  
  sims <- sims |> mutate(SSP = case_when(SSP == "SSP126" ~ "SSP 1-2.6", SSP == "SSP585" ~ "SSP 5-8.5", .default = SSP))
  
  # Probability of occurrence plots
  
  ebs_boundary <- st_buffer(st_union(st_as_sf(dplyr::select(sims[,,,1], p_occurrence))), 2500)
  
  hindcast_plot <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = p_occurrence), data = slice(sims, along = "id", 19)) + 
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), breaks = seq(0.25, 0.75, 0.25), na.value = NA) + 
    theme(legend.position = "right",  legend.direction = "horizontal", plot.margin = margin(1, 2, 1, 1, "lines")) + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"))) + 
    labs(fill = "Probability of Occurrence") + 
    ggtitle("  hindcast 1995-2015")
    
  plot_2040 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = p_occurrence), data = slice(sims, along = "id", 1:6)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2040-2059")
  
  plot_2060 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = p_occurrence), data = slice(sims, along = "id", 7:12)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2060-2079")
  
  plot_2080 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = p_occurrence), data = slice(sims, along = "id", 13:18)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2080-2099")
  
  joined_plot <- cowplot::plot_grid(hindcast_plot, plot_2040, plot_2060, plot_2080, ncol = 2) + theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(here("output", sp_paths[i], "projection_occurrence.png"), joined_plot, height = 8, width = 12, units = "in", dpi = 500)
  
  # Biomass plots
  
  sims$log_biomass <- pmax(-10, log(sims$biomass_fit))
  
  biomass_lim <- range(sims$log_biomass, na.rm = TRUE)
  
  hindcast_plot <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = log_biomass), data = slice(sims, along = "id", 19)) + 
    scale_fill_viridis_c(option = "magma", limits = biomass_lim, na.value = NA) + 
    theme(legend.position = "right",  legend.direction = "horizontal", plot.margin = margin(1, 2, 1, 1, "lines")) + 
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"))) + 
    labs(fill = expression(paste("Log(CPUE kg"~"km"^{2}, ")"))) + 
    ggtitle("  hindcast 1995-2015")
  
  plot_2040 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = log_biomass), data = slice(sims, along = "id", 1:6)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = biomass_lim, na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2040-2059")
  
  plot_2060 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = log_biomass), data = slice(sims, along = "id", 7:12)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = biomass_lim, na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2060-2079")
  
  plot_2080 <- ggplot() + 
    geom_sf(data = ebs_boundary, fill = "grey60", color = "grey60", linewidth = 1) + 
    geom_stars(aes(fill = log_biomass), data = slice(sims, along = "id", 13:18)) + 
    facet_grid(SSP ~ ESM, switch = "y") + 
    scale_fill_viridis_c(option = "magma", limits = biomass_lim, na.value = NA) + 
    theme(
      legend.position = "none", 
      panel.spacing = unit(-3, "lines"),
      strip.text = element_text(size = 12, margin = margin(0.4,0,0,0, "lines"), face = "bold", color = "grey60"),
      strip.text.y = element_text(size = 12, angle = 90)
    ) + 
    ggtitle("2080-2099")
  
  joined_plot <- cowplot::plot_grid(hindcast_plot, plot_2040, plot_2060, plot_2080, ncol = 2) + theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(here("output", sp_paths[i], "projection_biomass.png"), joined_plot, height = 8, width = 12, units = "in", dpi = 500)
  
  # Fitted vs. observed map 
  
  ## Match directory with species info
  species <- gsub("_", " ", strsplit(sp_paths[i], "-")[[1]][1])
  length_bin <- strsplit(sp_paths[i], "-")[[1]][2]
  survey_file <- unique(specs$survey_file[species == specs$species])
  threshold <- unique(specs$threshold[species == specs$species])
  
  ## Read in observed data
  if (!is.na(length_bin)) {
    source(here("analysis", "load_trawl_size_binned.R"))
    model_data <- joined_data |> filter(bin == length_bin) |> dplyr::select(year, station_id, observed = cpue_kgkm2)
    model_data <- val_data |> dplyr::select(year, station_id, observed = cpue_kgkm2) |> bind_rows(model_data)
  } else {
    source(here("analysis", "load_trawl.R"))
    model_data <- joined_data |> dplyr::select(year, station_id, observed = cpue_kgkm2)
    model_data <- val_data |> dplyr::select(year, station_id, observed = cpue_kgkm2) |> bind_rows(model_data)
  }
  
  ## Read in model projections
  fit <- readRDS(here("output", sp_paths[i], "projection_surveyrep.rds"))
  fit <- fit |> filter(sim == "hindcast")
  
  ## Join model fits with observed data
  plot_data <- fit |> left_join(model_data, by = c("year", "station_id"))
  
  ## PLot fitted vs. observed
  fit_obs_plot <- plot_data |> 
    filter(year %in% c(2010, 2015, 2019, 2021, 2022)) |> 
    left_join(ebs_grid, by = "station_id") |> 
    mutate(fitted = biomass_fit) |> 
    fitted_obs_map()
  
  ggsave(here("output", sp_paths[i], "fitted_observed_map.png"), fit_obs_plot,
         height = 5, width = 10, units = "in", dpi = 500)
  
}
