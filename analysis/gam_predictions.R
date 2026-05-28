
library("here")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("sf")
library("stars")
library("mgcv")
library("aclim2sdms")
library("dismo")

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

MOM6_fit <- MOM6_data
MOM6_data <- MOM6_data |> mutate(year_chr = as.character(year))

# Binomial model estimates and standard errors
binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = MOM6_data, type = "response")))
binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = MOM6_data, type = "response", se.fit = TRUE)$se.fit))

# Binomial ensemble estimates and standard errors
binom_se <- weighted_se(binom_fit, binom_se, w_binom)
binom_fit <- apply(binom_fit, 1, weighted.mean, w = w_binom)

# Tweedie model estimates and standard errors
tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = MOM6_data, type = "response")))
tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = MOM6_data, type = "response", se.fit = TRUE)$se.fit))

# Tweedie ensemble estimates and standard errors
tw_se <- weighted_se(tw_fit, tw_se, w_tw)
tw_fit <- apply(tw_fit, 1, weighted.mean, w = w_tw)

MOM6_fit <- cbind(
  MOM6_fit[,c("year", "station_id", "longitude", "latitude")], 
  data.frame(p_occurrence = binom_fit, p_occurrence_se = binom_se, biomass_fit = tw_fit, biomass_se = tw_se)
)

write.csv(MOM6_fit, paste0(save_dir, "hindcast_surveyrep_fit.csv"))

## Model predictions & SE - MOM6 level 2 hindcast ---------------------------------------

# Average area swept in km2
area_avg <- round(mean(MOM6_data$area_swept_km2), 5)

MOM6 <- readRDS(here("data", "mom6_hindcast", "mom6_hindcast_july1.rds"))

## Two-degree cold pool extent
cold_pool_2C <- MOM6 |> dplyr::select(temp_bottom5m) |> 
  st_apply(3, \(x) {x <- x[!is.na(x)]; sum(x < 2)/length(x)})
cold_pool_2C <- cold_pool_2C$temp_bottom5m

## Add cold pool extent to MOM6 raster 
MOM6$cold_pool_2C <- rep(cold_pool_2C, each = dim(MOM6)[1] * dim(MOM6)[2])

MOM6_yrs <- lubridate::year(st_get_dimension_values(MOM6, "time"))
p_occ <- se_p_occ <- cpue <- se_cpue <- vector("list", length(MOM6_yrs))

for (j in seq_along(MOM6_yrs)) {
  
  MOM6_yr <- slice(MOM6, along = "time", j)
  
  ## Need to work with data frame for predicting ignoring random effects
  MOM6_yr_df <- as.data.frame(MOM6_yr, add_coordinates = FALSE)
  MOM6_yr_df$area_swept_km2 <- area_avg
  MOM6_keep <- which(complete.cases(as.data.frame(MOM6_yr)))
  fit_vec <- rep(NA, nrow(MOM6_yr_df))
  
  ## Predict binomial model average on MOM6 grid
  binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = MOM6_yr_df[MOM6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = MOM6_yr_df[MOM6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  fit_vec[MOM6_keep] <- c(apply(binom_fit, 1, weighted.mean, w = w_binom)); p_occ[[j]] <- fit_vec
  fit_vec[MOM6_keep] <- c(weighted_se(binom_fit, binom_se, w_binom)); se_p_occ[[j]] <- fit_vec
  
  ## Predict Tweedie model average on MOM6 grid
  tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = MOM6_yr_df[MOM6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = MOM6_yr_df[MOM6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  fit_vec[MOM6_keep] <- c(apply(tw_fit, 1, weighted.mean, w = w_tw)); cpue[[j]] <- fit_vec
  fit_vec[MOM6_keep] <- c(weighted_se(tw_fit, tw_se, w_tw)); se_cpue[[j]] <- fit_vec
  
}

MOM6 <- MOM6 |> 
  mutate(
    p_occurrence = c(unlist(p_occ)), 
    p_occurrence_se = c(unlist(se_p_occ)), 
    biomass_fit = c(unlist(cpue)),
    biomass_se = c(unlist(se_cpue))
  ) |> 
  dplyr::select(
    p_occurrence, p_occurrence_se, biomass_fit, biomass_se
  )

saveRDS(MOM6, paste0(save_dir, "hindcast_level2.rds"))

## MOM6 forecast --------------------------------------------------------------

mom6_fcst <- readRDS(here("data", "mom6_forecast", "mom6_forecast.rds"))

fcst_df <- as.data.frame(mom6_fcst, add_coordinates = FALSE)
fcst_df$area_swept_km2 <- area_avg
mom6_keep <- which(complete.cases(as.data.frame(fcst_df)))
fit_vec <- rep(NA, nrow(fcst_df))

## Predict binomial model average for MOM6 forecast on MOM6 grid
binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
fit_vec[mom6_keep] <- c(apply(binom_fit, 1, weighted.mean, w = w_binom)); p_occ <- fit_vec
fit_vec[mom6_keep] <- c(weighted_se(binom_fit, binom_se, w_binom)); se_p_occ <- fit_vec

## Predict Tweedie model average for MOM6 forecast on MOM6 grid
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

saveRDS(mom6_fcst, paste0(save_dir, "forecast.rds"))

# Derived quantities ----------------------------------------------------------

# Number of samples from each model
nsim <- 1000
binom_n <- table(sort(sample(rep(1:length(binom_models), round(w_binom * 1010)), nsim, replace = FALSE)))
tw_n <- table(sort(sample(rep(1:length(tw_models), round(w_tw * 1010)), nsim, replace = FALSE)))

## Predict on link scale for binomial models
binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
binom_fit <- binom_sim <- binom_fit[,rep(seq_along(binom_n), times = binom_n)]
binom_se <- binom_se[,rep(seq_along(binom_n), times = binom_n)]
binom_sim[] <- cloglog(rnorm(prod(dim(binom_fit)), c(binom_fit), c(binom_se)))
binom_sim[] <- rbinom(prod(dim(binom_sim)), size = 1, prob = c(binom_sim))

# Predict for tweedie models
tw_n <- rep(seq_along(tw_n), times = tw_n)
tw_power <- vapply(tw_models, \(model) as.numeric(stringr::str_extract_all(unclass(model$family)$family, "\\d+([.,]\\d+)?")[[1]]), numeric(1))
tw_scale <- vapply(tw_models, \(model) model$scale, numeric(1))
tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = fcst_df[mom6_keep,], type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
tw_fit <- tw_sim <- tw_fit[,tw_n]; tw_se <- tw_se[,tw_n]
tw_sim[] <- exp(rnorm(prod(dim(tw_fit)), c(tw_fit), c(tw_se)))
for (i in seq_len(ncol(tw_sim))) {
  tw_sim[,i] <- tweedie::rtweedie(nrow(tw_fit), mu = c(tw_sim[,i]), power = tw_power[tw_n[i]], phi = tw_scale[tw_n[i]])
}

# Center of gravity
cog <- data.frame(
  model = rep(c("occurrence", "biomass"), each = 2 * nsim), 
  coord = rep(rep(c("E_km", "N_km"), each = nsim), times = 2), 
  sim = rep(seq_len(nsim), times = 2 * nsim),
  cog = c(
    apply(binom_sim, 2, \(x) weighted.mean(fcst_df$X[mom6_keep], x)), 
    apply(binom_sim, 2, \(x) weighted.mean(fcst_df$Y[mom6_keep], x)), 
    apply(tw_sim, 2, \(x) weighted.mean(fcst_df$X[mom6_keep], x)), 
    apply(tw_sim, 2, \(x) weighted.mean(fcst_df$Y[mom6_keep], x))
  )
)

# area occupied
area <- apply(binom_sim, 2, \(x) sum(x)/length(x))

saveRDS(list(cog = cog, area_occupied = area), paste0(save_dir, "derived_quantities.rds"))

# Exit ----------------------------------------------------------------------------------

file.remove(paste0(save_dir, "running"))
file.create(paste0(save_dir, "complete"))
