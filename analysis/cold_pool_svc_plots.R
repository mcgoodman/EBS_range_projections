
pkgs <- c("here", "dplyr", "stars", "sf", "tidyr", "ggplot2")
sapply(pkgs, require, character.only = TRUE)

theme_set(
  theme_void() + theme(
    strip.placement = "outside", 
    plot.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(size = 16, margin = margin(0.4,0,0.4,0, "lines")), 
    strip.text.y = element_text(size = 16, angle = 90),
    plot.title = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 16, color = "black"), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
)

ebs <- get_ebs_shapefile()

cold_pool_files <- list.files(here("output"), pattern = "cold_pool_SVC.rds", recursive = TRUE)
sp_dirs <- dirname(cold_pool_files)

for (i in seq_along(sp_dirs)) {
  
  svc_data <- readRDS(here("output", cold_pool_files[i]))
  
  svc_data$cloglog_occ <- pmax(-10, cloglog(svc_data$p_occurrence, inverse = TRUE))
  svc_data$log_biomass_fit <- pmax(-10, log(svc_data$biomass_fit))
  
  binom_plot <- svc_data |> 
    mutate(extent = paste0(extent*100, "%")) |> 
    ggplot() + 
    geom_sf(data = ebs, fill = NA, color = "grey60", linewidth = 2) +
    geom_sf(aes(fill = cloglog_occ, color = cloglog_occ)) + 
    facet_wrap(~extent, nrow = 1) + 
    scale_fill_viridis_c(option = "magma") + 
    scale_color_viridis_c(option = "magma") + 
    theme(legend.position = "bottom") + 
    guides(
      fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "grey60", frame.linewidth = 0.8), 
      color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "grey60", frame.linewidth = 0.8)
    ) + 
    labs(fill = "s(log-odds occurrence)", color = "s(log-odds occurrence)")
  
  ggsave(paste0(here("output"), "/", sp_dirs[i], "/binomial_SVC.png"), binom_plot, height = 4, width = 12, units = "in", dpi = 500)
  
  tw_plot <- svc_data |> 
    mutate(extent = paste0(extent*100, "%")) |> 
    ggplot() + 
    geom_sf(data = ebs, fill = NA, color = "grey60", linewidth = 2) +
    geom_sf(aes(fill = log_biomass_fit, color = log_biomass_fit)) + 
    facet_wrap(~extent, nrow = 1) + 
    scale_fill_viridis_c(option = "magma") + 
    scale_color_viridis_c(option = "magma") + 
    theme(legend.position = "bottom") + 
    guides(
      fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "grey60", frame.linewidth = 0.8), 
      color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "grey60", frame.linewidth = 0.8)
    ) + 
    labs(fill = "s(log CPUE)", color = "s(log CPUE)")
  
  ggsave(paste0(here("output"), "/", sp_dirs[i], "/tweedie_SVC.png"), tw_plot, height = 4, width = 12, units = "in", dpi = 500)
  
  
}