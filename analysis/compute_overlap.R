
pkgs <- c("here", "dplyr", "tidyr", "stars", "sf", "aclim2sdms")
sapply(pkgs, require, character.only = TRUE)

# Combinations of species / bins to compute overlap for -----------------------------------------

output_dirs <- list.dirs(here("output"), full.names = FALSE, recursive = FALSE)

predators <- output_dirs[grepl("adult", output_dirs)]
prey <- output_dirs[grepl("juvenile", output_dirs) | grepl("crab", output_dirs)]

sp_grid <- arrange(expand.grid(pred_dir = predators, prey_dir = prey), pred_dir)

sp_grid <- sp_grid |> 
  separate(pred_dir, into = c("pred", "pred_bin"), sep = "-", fill = "right", remove = FALSE) |> 
  separate(prey_dir, into = c("prey", "prey_bin"), sep = "-", fill = "right", remove = FALSE) |> 
  mutate_at(c("pred", "prey"), \(x) gsub("_", " ", x))

sp_overlap <- vector("list", nrow(sp_grid))

# Area of EBS station cells -----------------------------------------------------------------------

stn_area <- get_ebs_shapefile(type = "grid") |> dplyr::select(station_id, AREA_KM2) |> st_drop_geometry()

# Compute overlap ---------------------------------------------------------------------------------

for(i in seq(nrow(sp_grid))) {
  
  print(paste0(i, " (~" ,round(100 * i / nrow(sp_grid)), "%)"))
  
  pred <- sp_grid$pred[i]
  prey <- sp_grid$prey[i]
  
  pred_bin <- sp_grid$pred_bin[i]
  prey_bin <- sp_grid$prey_bin[i]
  
  pred_dir <- paste0(here("output", sp_grid$pred_dir[i]), "/")
  prey_dir <- paste0(here("output", sp_grid$prey_dir[i]), "/")
  
  ## Read in fitted means for each species
  pred_fit <- readRDS(paste0(pred_dir, "projection_surveyrep.rds"))
  prey_fit <- readRDS(paste0(prey_dir, "projection_surveyrep.rds"))
  stopifnot(identical(pred_fit[,c("year", "sim", "longitude", "latitude")], prey_fit[,c("year", "sim", "longitude", "latitude")]))
  
  ## Read in posterior simulations for each species
  pred_sim <- readRDS(paste0(pred_dir, "posterior_predictive.rds"))
  prey_sim <- readRDS(paste0(prey_dir, "posterior_predictive.rds"))
  
  ## Years and climate projections corresponding to each row
  fit_yr <- pred_fit$year
  fit_sim <- pred_fit$sim
  fit_x <- pred_fit$longitude
  fit_y <- pred_fit$latitude
  
  # Area vector
  area <- stn_area$AREA_KM2[match(pred_fit$station_id, stn_area$station_id)]
  area <- ifelse(is.na(area), 0, area)
  
  ## Loop over posterior simulations, store overlap
  overlap <- vector("list", ncol(pred_sim$tweedie))
  
  for (j in 1:length(overlap)) {
    
    overlap[[j]] <- 
      data.frame(
        year = fit_yr, sim = fit_sim, sample = j, 
        pred_biomass = pred_sim$tweedie[,j], prey_biomass = prey_sim$tweedie[,j],
        pred_present = pred_sim$binomial[,j], prey_present = prey_sim$binomial[,j],
        X = fit_x, Y = fit_y, area = area
      ) |> 
      group_by(year, sim, sample) |> 
      summarize(
        loc_colloc = loc_collocfn(prey_biomass, pred_biomass),
        bhattacharyya = bhatta_coeffn(prey_biomass, pred_biomass), 
        gbl_colloc = glob_collocfn(X, Y, prey_biomass, X, Y, pred_biomass),
        area_overlap = area_overlapfn(prey_present, pred_present, area),
        range_overlap = range_overlapfn(prey_present, pred_present, area),
        bhattacharyya_encounter = bhatta_coeffn(prey_present, pred_present),
        gbl_colloc_encounter = glob_collocfn(X, Y, prey_present, X, Y, pred_present),
        AB_overlap = AB_overlapfn(prey_biomass, pred_biomass),
        .groups = "drop"
      ) |> 
      pivot_longer(cols = c("loc_colloc", "bhattacharyya", "gbl_colloc", "area_overlap", "range_overlap",
                            "bhattacharyya_encounter", "gbl_colloc_encounter", "AB_overlap"), 
                   names_to = "index", values_to = "overlap")
    
  }
  
  overlap <- do.call("rbind", overlap)
  
  sp_overlap[[i]] <- overlap |> 
    group_by(year, sim, index) |> 
    summarise(
      mean = mean(overlap), 
      se = sd(overlap),
      `2.5%` = quantile(overlap, 0.025), 
      `10%` = quantile(overlap, 0.1), 
      `90%` = quantile(overlap, 0.9), 
      `97.5%` = quantile(overlap, 0.975), 
      .groups = "drop"
    ) |> 
    mutate(
      pred = pred, prey = prey, 
      pred_bin = pred_bin, prey_bin = prey_bin, 
      .before = year
    )
  
}

overlap_full <- do.call("rbind", sp_overlap)

write.csv(overlap_full, here("output", "overlap_summary.csv"), row.names = FALSE)
