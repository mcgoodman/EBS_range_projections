
pkgs <- c("here", "dplyr", "tidyr", "purrr", "ggplot2", "sf", "stars", "mgcv", "aclim2sdms", "dismo")
sapply(pkgs, require, character.only = TRUE)

save_dir <- paste0(here("output", paste0(gsub(" ", "_", species), ifelse(is.na(length_bin), "", paste0("-", length_bin)))), "/")
dir.create(save_dir)

file.create(paste0(save_dir, "running"))

# Run candidate models ------------------------------------------------------------------

## List to store model metrics (e.g. AIC, RMSE)
metrics <- binom_models <- tw_models <- vector("list", length(mod_forms))

## Fit candidate models, store forecast AUC / RMSE / AIC / Deviance explained
for (i in seq_along(mod_forms)) {

  print(paste0("model ", i, "/", length(mod_forms), ": ", Reduce(paste, trimws(deparse(mod_forms[[i]][[2]])))))

  ## Run model on full data
  binom_models[[i]] <- gam(update(mod_forms[[i]], present ~ . + offset(log(area_swept_km2))), data = model_data, 
                           family = binomial(link = "cloglog"), select = TRUE, method = "REML", optimizer = "efs")
  tw_models[[i]] <- gam(update(mod_forms[[i]], cpue_kgkm2 ~ . + s(year_chr, bs = "re")), data = model_data, 
                        family = tw(link = "log"), select = TRUE, method = "REML", optimizer = "efs")

  metrics[[i]] <- data.frame(
    id = i,
    species = rep(species, 2),
    component = c("binomial", "tweedie"),
    formula = Reduce(paste, deparse(mod_forms[[i]])),
    AIC = c(AIC(binom_models[[i]]), AIC(tw_models[[i]])),
    dev.expl = c(summary(binom_models[[i]])$dev.expl, summary(tw_models[[i]])$dev.expl), 
    AUC = c(AUC(model_data$present, predict(binom_models[[i]], type = "response")), NA),
    biomass_pearson = c(NA, cor(model_data$cpue_kgkm2, predict(tw_models[[i]], type = "response", use = "complete.obs"), method = "pearson")),
    biomass_spearman = c(NA, cor(model_data$cpue_kgkm2, predict(tw_models[[i]], type = "response", use = "complete.obs"), method = "spearman"))
  )

}

metrics <- do.call("rbind", metrics)

# Predict on model data -----------------------------------------------------------------

fit_obs <- model_data |> 
  dplyr::select(year, station_id, cpue_kgkm2) |> 
  mutate(
    p_occurrence = apply(simplify2array(lapply(binom_models, \(x) predict(x, type = "response"))), 1, weighted.mean, w = w_binom),
    fitted_cpue_kgkm2 = apply(simplify2array(lapply(tw_models, \(x) predict(x, type = "response"))), 1, weighted.mean, w = w_tw)
  )

write.csv(fit_obs, paste0(save_dir, "fitted_observed.csv"), row.names = FALSE)

# Save model estimates ------------------------------------------------------------------

## Model predictions & SE - survey-replicated scale -------------------------------------

ROMS_fit <- ROMS_data
ROMS_data <- ROMS_data |> mutate(year_chr = as.character(year))

# Binomial model estimates and standard errors
binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = ROMS_data, type = "response")))
binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = ROMS_data, type = "response", se.fit = TRUE)$se.fit))

# Binomial ensemble estimates and standard errors
binom_se <- weighted_se(binom_fit, binom_se, w_binom)
binom_fit <- apply(binom_fit, 1, weighted.mean, w = w_binom)

# Tweedie model estimates and standard errors
tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = ROMS_data, type = "response")))
tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = ROMS_data, type = "response", se.fit = TRUE)$se.fit))

# Tweedie ensemble estimates and standard errors
tw_se <- weighted_se(tw_fit, tw_se, w_tw)
tw_fit <- apply(tw_fit, 1, weighted.mean, w = w_tw)

ROMS_fit <- cbind(
  ROMS_fit[,c("year", "station_id", "longitude", "latitude")], 
  data.frame(p_occurrence = binom_fit, p_occurrence_se = binom_se, biomass_fit = tw_fit, biomass_se = tw_se)
)

write.csv(ROMS_fit, paste0(save_dir, "hindcast_surveyrep_fit.csv"))

## Model predictions & SE - MOM6 level 2 hindcast ---------------------------------------

# Average area swept in km2
area_avg <- round(mean(ROMS_data$area_swept_km2), 5)

roms <- readRDS(here("data", "mom6", "mom6_hindcast.rds"))

## Two-degree cold pool extent
cold_pool_2C <- roms |> dplyr::select(temp_bottom5m) |> 
  st_apply(3, \(x) {x <- x[!is.na(x)]; sum(x < 2)/length(x)})
cold_pool_2C <- cold_pool_2C$temp_bottom5m

## Add cold pool extent to ROMS raster 
roms$cold_pool_2C <- rep(cold_pool_2C, each = dim(roms)[1] * dim(roms)[2])

roms_yrs <- st_get_dimension_values(roms, "ocean_time")
p_occ <- se_p_occ <- cpue <- se_cpue <- vector("list", length(roms_yrs))

for (j in 1:length(roms_yrs)) {
  
  roms_yr <- slice(roms, along = "ocean_time", j)
  
  ## Need to work with data frame for predicting ignoring random effects
  roms_yr_df <- as.data.frame(roms_yr)
  roms_yr_df$area_swept_km2 <- area_avg
  roms_keep <- which(complete.cases(as.data.frame(roms_yr)))
  fit_vec <- rep(NA, nrow(roms_yr_df))
  
  ## Predict binomial model average on ROMS grid
  binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = roms_yr_df[roms_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = roms_yr_df[roms_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  fit_vec[roms_keep] <- c(apply(binom_fit, 1, weighted.mean, w = w_binom)); p_occ[[j]] <- fit_vec
  fit_vec[roms_keep] <- c(weighted_se(binom_fit, binom_se, w_binom)); se_p_occ[[j]] <- fit_vec
  
  ## Predict Tweedie model average on ROMS grid
  tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = roms_yr_df[roms_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = roms_yr_df[roms_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  fit_vec[roms_keep] <- c(apply(tw_fit, 1, weighted.mean, w = w_tw)); cpue[[j]] <- fit_vec
  fit_vec[roms_keep] <- c(weighted_se(tw_fit, tw_se, w_tw)); se_cpue[[j]] <- fit_vec
  
}

roms <- roms |> 
  mutate(
    p_occurrence = c(unlist(p_occ)), 
    p_occurrence_se = c(unlist(se_p_occ)), 
    biomass_fit = c(unlist(cpue)),
    biomass_se = c(unlist(se_cpue))
  ) |> 
  dplyr::select(
    p_occurrence, p_occurrence_se, biomass_fit, biomass_se
  )

saveRDS(roms, paste0(save_dir, "hindcast_level2.rds"))

## Forecast, level-2, ROMS/MOM6 climatology adjusted ------------------------------------

mom6_fcst <- readRDS(here("data", "mom6", "mom6_forecast_adj.rds"))

fcst_df <- as.data.frame(mom6_fcst)
fcst_df$area_swept_km2 <- area_avg
mom6_keep <- which(complete.cases(as.data.frame(fcst_df)))
fit_vec <- rep(NA, nrow(fcst_df))

## Predict binomial model average for MOM6 forecast on ROMS grid
binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
fit_vec[mom6_keep] <- c(apply(binom_fit, 1, weighted.mean, w = w_binom)); p_occ <- fit_vec
fit_vec[mom6_keep] <- c(weighted_se(binom_fit, binom_se, w_binom)); se_p_occ <- fit_vec

## Predict Tweedie model average for MOM6 forecast on ROMS grid
tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
fit_vec[mom6_keep] <- c(apply(tw_fit, 1, weighted.mean, w = w_tw)); cpue <- fit_vec
fit_vec[mom6_keep] <- c(weighted_se(tw_fit, tw_se, w_tw)); se_cpue <- fit_vec

mom6_fcst <- mom6_fcst |> 
  mutate(
    p_occurrence = c(unlist(p_occ)), 
    p_occurrence_se = c(unlist(se_p_occ)), 
    biomass_fit = c(unlist(cpue)),
    biomass_se = c(unlist(se_cpue))
  ) |> 
  dplyr::select(
    p_occurrence, p_occurrence_se, biomass_fit, biomass_se
  )

saveRDS(mom6_fcst, paste0(save_dir, "forecast_level2_adj.rds"))

# Exit ----------------------------------------------------------------------------------

file.remove(paste0(save_dir, "running"))
file.create(paste0(save_dir, "complete"))
