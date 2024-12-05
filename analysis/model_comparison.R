
pkgs <- c("here", "dplyr", "tidyr", "aclim2sdms", "mgcv", "foreach", "doParallel")
sapply(pkgs, require, character.only = TRUE)

# Fit models for each species, evaluate fit metrics and projections ---------------------

reduced_forms <- list(
  ~ s(temp_bottom5m, k = 3) + s(depth_m, k = 3) + s(phi, k = 3), 
  ~ s(temp_bottom5m, k = 3) + s(X, Y, k = 20) + s(X, Y, by = cold_pool_2C, k = 20)
)

all_forms <- c(mod_forms, reduced_forms)

ROMS_sub <- ROMS_full |> filter((year %in% 1995:2014) | ((year %in% 2075:2094) & (sim == "SSP585 GFDL")))

cl <- makeCluster(cores)
registerDoParallel(cl)

model_comparison <- foreach (
  i = 1:nrow(specs), .combine = "rbind", 
  .packages = c("dplyr", "aclim2sdms", "mgcv", "here"), 
  .export = c("specs", "ROMS_full", "ROMS_sub", "all_forms", "ebs_stations")
) %dopar% {
  
  species <- specs$species[i]
  survey_file <- specs$survey_file[i]
  threshold <- specs$threshold[i]
  length_bin <- specs$length_bin[i]
  
  if (!is.na(length_bin)) {
    
    source(here("analysis", "load_trawl_size_binned.R"), local = TRUE)
    
    model_data <- joined_data |> filter(bin == length_bin)
    val_data <- val_data |> filter(bin == length_bin)
    
  } else {
    
    source(here("analysis", "load_trawl.R"), local = TRUE)
    
    model_data <- joined_data
    
  }
  
  ## List to store model metrics (e.g. AIC, RMSE)
  metrics <- binom_models <- tw_models <- vector("list", length(all_forms))
  
  ## Fit candidate models, store forecast AUC / RMSE / AIC / Deviance explained
  for (j in seq_along(all_forms)) {
    
    print(paste0("model ", j, "/", length(all_forms), ": ", Reduce(paste, trimws(deparse(all_forms[[j]][[2]])))))
    
    ## Run model on full data
    binom_models[[j]] <- tryCatch(gam(update(all_forms[[j]], present ~ . + offset(log(area_swept_km2))), data = model_data, family = binomial(link = "cloglog"), select = TRUE, method = "REML", optimizer = "efs"), error = \(x) NA)
    tw_models[[j]] <- tryCatch(gam(update(all_forms[[j]], cpue_kgkm2 ~ . + s(year_chr, bs = "re")), data = model_data, family = tw(link = "log"), select = TRUE, method = "REML", optimizer = "efs"), error = \(x) NA)
    
    if (("gam" %in% class(binom_models[[j]])) & ("gam" %in% class(tw_models[[j]]))) {
      
      ## Change in mean prevalence and range centroid
      ROMS_change <- ROMS_sub |> mutate(
        p_occurrence = predict(binom_models[[j]], newdata = ROMS_sub, type = "response"), 
        biomass = predict(tw_models[[j]], newdata = ROMS_sub, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)
      ) |> group_by(sim) |> 
        summarize(prevalence = mean(p_occurrence), northings_biomass = weighted.mean(Y, biomass), northings_occurrence = weighted.mean(Y, p_occurrence))
      
      ## Fit, hindcast, and forecast stats
      metrics[[j]] <- data.frame(
        id = j,
        species = rep(species, 2),
        bin = rep(length_bin, 2),
        component = c("binomial", "tweedie"),
        formula = Reduce(paste, deparse(all_forms[[j]])),
        AIC = c(AIC(binom_models[[j]]), AIC(tw_models[[j]])),
        dev.expl = c(summary(binom_models[[j]])$dev.expl, summary(tw_models[[j]])$dev.expl), 
        AUC_train = c(AUC(model_data$present, predict(binom_models[[j]], type = "response")), NA),
        AUC_val = c(AUC(val_data$present, predict(binom_models[[j]], newdata = val_data, type = "response")), NA),
        biomass_spearman_train = c(NA, cor(model_data$cpue_kgkm2, predict(tw_models[[j]], type = "response", use = "complete.obs"), method = "spearman")), 
        biomass_spearman_val = c(NA, cor(val_data$cpue_kgkm2, predict(tw_models[[j]],  newdata = val_data, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE), method = "spearman")),
        prevalence_change = c(ROMS_change$prevalence[ROMS_change$sim == "SSP585 GFDL"]/ROMS_change$prevalence[ROMS_change$sim == "hindcast"] - 1, NA), 
        northings_change_biomass = c(NA, ROMS_change$northings_biomass[ROMS_change$sim == "SSP585 GFDL"] - ROMS_change$northings_biomass[ROMS_change$sim == "hindcast"] - 1), 
        northings_change_occurrence = c(NA, ROMS_change$northings_occurrence[ROMS_change$sim == "SSP585 GFDL"] - ROMS_change$northings_occurrence[ROMS_change$sim == "hindcast"] - 1)
      )
      
    }
    
  }
  
  metrics <- do.call("rbind", metrics)
  
  return(metrics)
  
}

stopCluster(cl)

write.csv(model_comparison, here("output", "model_comparison.csv"), row.names = FALSE)
