
pkgs <- c("here", "dplyr", "tidyr", "purrr", "ggplot2", "sf", "stars", "mgcv", "aclim2sdms", "dismo", "tweedie")
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

# Obtain model weights ------------------------------------------------------------------

nrep <- 100

binom_weights <- tw_weights <- matrix(NA_real_, nrow = length(mod_forms), ncol = nrep)

for (i in 1:nrep) {
  
  cat(paste0("obtaining weights (rep ", i, "/", nrep, ")\n"))
    
  model_data$fold <- dismo::kfold(model_data, k = 5, by = model_data$year)
  
  train <- model_data[!(model_data$fold == 5), ]
  test <- model_data[model_data$fold == 5, ]
  
  tw_i <- binom_i <- vector("list", length(mod_forms))
  for (j in seq_along(mod_forms)) {
    tw_i[[j]] <- gam(update(mod_forms[[j]], cpue_kgkm2 ~ . + s(year_chr, bs = "re")), data = train, 
                     family = tw(link = "log"), select = TRUE, method = "REML", optimizer = "efs")
    binom_i[[j]] <- gam(update(mod_forms[[j]], present ~ . + offset(log(area_swept_km2))), data = train, 
                        family = binomial(link = "cloglog"), select = TRUE, method = "REML", optimizer = "efs")
  }
  
  yhat <- simplify2array(lapply(tw_i, \(x) predict(x, test, type = "response")))
  tw_weights[,i] <- solve_weights(test$cpue_kgkm2, yhat)
  
  yhat <- simplify2array(lapply(binom_i, \(x) predict(x, test, type = "response")))
  binom_weights[,i] <- solve_weights(test$present, yhat)
  
}

saveRDS(list(binomial = binom_weights, tweedie = tw_weights), paste0(save_dir, "weight_stacking.rds"))

metrics$weight[metrics$component == "tweedie"] <- w_tw <- rowMeans(tw_weights)/sum(rowMeans(tw_weights))
metrics$weight[metrics$component == "binomial"] <- w_binom <- rowMeans(binom_weights)/sum(rowMeans(binom_weights))

write.csv(metrics, paste0(save_dir, "candidate_models.csv"), row.names = FALSE)

# Predict on validation data ------------------------------------------------------------

# Ensemble mean - probability of occurrence
val_binom <- lapply(binom_models, \(x) predict(x, newdata = val_data, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE))
val_binom <- simplify2array(val_binom)
val_data$p_occurrence <- apply(val_binom, 1, weighted.mean, w = w_binom)

# Ensemble mean - biomass (carrying last year's random effect estimate forward)
val_tw <- lapply(tw_models, \(x) predict(x, newdata = val_data, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE))
val_tw <- simplify2array(val_tw)
val_data$cpue_fitted <- apply(val_tw, 1, weighted.mean, w = w_tw)

# Compute validation stats for 2021-2022
val_stats <- val_data |> 
  group_by(year = as.character(year)) |>
  summarize(
    AUC = AUC(present, p_occurrence), 
    biomass_pearson = cor(cpue_kgkm2, cpue_fitted, method = "pearson"), 
    biomass_spearman = cor(cpue_kgkm2, cpue_fitted, method = "spearman") 
  ) |> 
  add_row(
    year = "all", 
    AUC = AUC(val_data$present, val_data$p_occurrence), 
    biomass_pearson = cor(val_data$cpue_kgkm2, val_data$cpue_fitted), 
    biomass_spearman = cor(val_data$cpue_kgkm2, val_data$cpue_fitted, method = "spearman") 
  )

write.csv(val_stats, paste0(save_dir, "validation_stats.csv"), row.names = FALSE)

# Re-run models with all data -----------------------------------------------------------

model_data <- bind_rows(model_data, val_data)

## List to store model metrics (e.g. AIC, RMSE)
binom_models <- tw_models <- vector("list", length(mod_forms))

## Re-fit candidate models
for (i in seq_along(mod_forms)) {
  
  print(paste0("model ", i, "/", length(mod_forms), ": ", Reduce(paste, trimws(deparse(mod_forms[[i]][[2]])))))
  
  ## Run model on full data
  binom_models[[i]] <- gam(update(mod_forms[[i]], present ~ . + offset(log(area_swept_km2))), data = model_data, 
                           family = binomial(link = "cloglog"), select = TRUE, method = "REML", optimizer = "efs")
  tw_models[[i]] <- gam(update(mod_forms[[i]], cpue_kgkm2 ~ . + s(year_chr, bs = "re")), data = model_data, 
                        family = tw(link = "log"), select = TRUE, method = "REML", optimizer = "efs")
  
}

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

years <- sort(unique(ROMS_full$year))

ROMS_fit <- ROMS_full |> group_by(year, sim) |> nest()

ROMS_fit_all <- list(binomial = vector("list", nrow(ROMS_fit)), tweedie = vector("list", nrow(ROMS_fit)))

for (i in 1:nrow(ROMS_fit)) {
  
  message(paste0("obtaining ensemble predictions (", ROMS_fit$sim[i], " ", ROMS_fit$year[i], ")\r"), appendLF = FALSE)

  # Binomial model estimates and standard errors
  binom_fit <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = ROMS_fit$data[[i]], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  binom_se <- simplify2array(lapply(binom_models, \(x) predict(x, newdata = ROMS_fit$data[[i]], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  ROMS_fit_all$binomial[[i]] <- binom_fit
  
  # Binomial ensemble estimates and standard errors
  binom_se <- weighted_se(binom_fit, binom_se, w_binom)
  binom_fit <- apply(binom_fit, 1, weighted.mean, w = w_binom)
  
  # Tweedie model estimates and standard errors
  tw_fit <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = ROMS_fit$data[[i]], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
  tw_se <- simplify2array(lapply(tw_models, \(x) predict(x, newdata = ROMS_fit$data[[i]], type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
  ROMS_fit_all$tweedie[[i]] <- tw_fit
  
  # Tweedie ensemble estimates and standard errors
  tw_se <- weighted_se(tw_fit, tw_se, w_tw)
  tw_fit <- apply(tw_fit, 1, weighted.mean, w = w_tw)
  
  ROMS_fit$data[[i]] <- cbind(
    ROMS_fit$data[[i]][,c("station_id", "longitude", "latitude")], 
    data.frame(p_occurrence = binom_fit, p_occurrence_se = binom_se, biomass_fit = tw_fit, biomass_se = tw_se)
  )
  
}

cat("\n")

names(ROMS_fit_all$binomial) <- names(ROMS_fit_all$tweedie) <- paste(ROMS_fit$sim, ROMS_fit$year, sep = "_")
saveRDS(ROMS_fit_all, paste0(save_dir, "projection_surveyrep_allmodels.rds"))

ROMS_fit <- ROMS_fit |> unnest(cols = c(data)) |> ungroup()

## Environmental extrapolation (Mahalanobis distance) -----------------------------------

extrap <- matrix(NA_real_, nrow = nrow(ROMS_full), ncol = length(mod_forms))

for (i in 1:length(mod_forms)) {
  
  ## Centroid and covariance matrix for hindcast covariates in binomial model
  vars <- all.vars(mod_forms[[i]])
  vars <- vars[vars != "X" & vars != "Y"]
  
  ## Compute Mahalanobis distance for binomial model
  if (length(vars) > 0) {
    hindcast_center <- colMeans(ROMS_full[ROMS_full$sim == "hindcast" & ROMS_full$year >= 1995 & ROMS_full$year <= 2015, vars])
    hindcast_cov <- cov(ROMS_full[ROMS_full$sim == "hindcast" & ROMS_full$year >= 1995 & ROMS_full$year <= 2015, vars])
    extrap[,i] <- c(mahalanobis(ROMS_full[, vars], hindcast_center, hindcast_cov))
  } else {
    extrap[,i] <- NA
  }
  
}

ROMS_fit$mahalanobis_binom <- apply(extrap, 1, weighted.mean, w = w_binom)
ROMS_fit$mahalanobis_tweedie <- apply(extrap, 1, weighted.mean, w = w_tw)

saveRDS(ROMS_fit, paste0(save_dir, "projection_surveyrep.rds"))

## Posterior samples - survey-replicated scale ------------------------------------------

set.seed(100)

# Extract Tweedie power and dispersion parameters
tw_power <- vapply(tw_models, \(model) as.numeric(stringr::str_extract_all(unclass(model$family)$family, "\\d+([.,]\\d+)?")[[1]]), numeric(1))
tw_scale <- vapply(tw_models, \(model) model$scale, numeric(1))

# Method to portion out 1000 samples from candidate models with probability very close to weights
binom_n <- table(sort(sample(rep(1:length(binom_models), round(w_binom * 1010)), 1000, replace = FALSE)))
tw_n <- table(sort(sample(rep(1:length(tw_models), round(w_tw * 1010)), 1000, replace = FALSE)))

# Sample from binomial candidate models
binom_post <- binom_samples <- vector("list", length(binom_n))
for (i in 1:length(binom_n)) {
  mod_i <- as.integer(names(binom_n))[i]
  message(paste0("sampling from binomial model ", mod_i, "\r"), appendLF = FALSE)
  fit <- predict(binom_models[[mod_i]], newdata = ROMS_full, type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)
  bp_i <- t(simplify2array(purrr::map2(fit$fit, fit$se.fit, ~cloglog(rnorm(binom_n[i], .x, .y)))))
  bs_i <- simplify2array(apply(bp_i, 2, \(x) rbinom(length(x), size = 1, prob = x)))
  if (binom_n[i] > 1) {
    binom_post[[i]] <- bp_i
    binom_samples[[i]] <- bs_i 
  } else { 
    binom_post[[i]] <- t(bp_i)
    binom_samples[[i]] <- t(t(bs_i))
  }
  colnames(binom_post[[i]]) <- colnames(binom_samples[[i]]) <- paste(i, 1:binom_n[i], sep = ".")
}
binom_post <- do.call("cbind", binom_post)
binom_samples <- do.call("cbind", binom_samples)

# Sample from Tweedie candidate models
tw_post <- tw_samples <- vector("list", length(tw_n))
for (i in 1:length(tw_n)) {
  mod_i <- as.integer(names(tw_n))[i]
  message(paste0("sampling from tweedie model ", mod_i, "\r"), appendLF = FALSE)
  fit <- predict(tw_models[[mod_i]], newdata = ROMS_full, type = "link", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)
  twp_i <- t(simplify2array(purrr::map2(fit$fit, fit$se.fit, ~exp(rnorm(tw_n[i], .x, .y)))))
  tws_i <- simplify2array(apply(twp_i, 2, \(x) tweedie::rtweedie(length(x), mu = x, power = tw_power[i], phi = tw_scale[i])))
  if (tw_n[i] > 1) {
    tw_post[[i]] <- twp_i
    tw_samples[[i]] <- tws_i 
  } else { 
    tw_post[[i]] <- t(twp_i)
    tw_samples[[i]] <- t(t(tws_i))
  }
  colnames(tw_post[[i]]) <- colnames(tw_samples[[i]]) <- paste(i, 1:tw_n[i], sep = ".")
}
tw_post <- do.call("cbind", tw_post)
tw_samples <- do.call("cbind", tw_samples)

post <- list(binomial = binom_post, tweedie = tw_post)
samples <- list(binomial = binom_samples, tweedie = tw_samples)

saveRDS(post, paste0(save_dir, "posterior_mu.rds"))
saveRDS(samples, paste0(save_dir, "posterior_predictive.rds"))

## Model predictions & SE - ROMS level 2 scale ------------------------------------------

# Average area swept in km2
area_avg <- round(mean(read.csv(here("data", "surveyrep_observed_1982-2022.csv"))$AREA_SWEPT_HA / 100), 5)

dir.create(paste0(save_dir, "level2_projections"))

roms_files <- list.files(here("data", "roms_level2_bc_annual"))

cat("obtaining grid-scale predictions:")

for (i in 1:length(roms_files)) {
  
  cat(paste("   ", gsub(".rds", "", roms_files[i]), "\n"))
  
  roms <- readRDS(here("data", "roms_level2_bc_annual", roms_files[i]))
  
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
  
  saveRDS(roms, paste0(save_dir, "level2_projections/", roms_files[i]))
  
}

# Plot partial effects ------------------------------------------------------------------

binom_smooths <- tw_smooths <- vector("list", length(binom_models))

for (i in 1:length(mod_forms)) {
  
  terms <- attr(terms(mod_forms[[i]]), "term.labels")
  
  bsi <- tsi <- vector("list", length(terms))
  
  for (j in 1:length(terms)) {
    
    plot_var <- all.vars(as.formula(paste("~", terms[j])))
    term_lab <- terms[j]
    
    if (length(plot_var) == 1) {
      
      bsi[[j]] <- plot_smoother(binom_models[[i]], plot_var)$data
      bsi[[j]] <- bsi[[j]] |> mutate(component = "binomial", model = i, predictor = plot_var)
      
      tsi[[j]] <- plot_smoother(tw_models[[i]], plot_var)$data
      tsi[[j]] <- tsi[[j]] |> mutate(component = "tweedie", model = i, predictor = plot_var)
      
    }
    
  }
  
 binom_smooths[[i]] <- do.call("rbind", bsi)
 tw_smooths[[i]] <- do.call("rbind", tsi)
  
}

binom_smooths <- do.call("rbind", binom_smooths)
tw_smooths <- do.call("rbind", tw_smooths)

binom_smooths$weight <- w_binom[binom_smooths$model]
tw_smooths$weight <- w_tw[tw_smooths$model]

smooth_data <- rbind(binom_smooths, tw_smooths)

write.csv(smooth_data, paste0(save_dir, "smoother_fits.csv"), row.names = FALSE)

# Long format of observed predictors for plotting rug
data_long <- pivot_longer(model_data[,unique(binom_smooths$predictor)], cols = unique(binom_smooths$predictor), names_to = "predictor", values_to = "x")
data_long <- data_long |> mutate(predictor = vapply(predictor, \(x) strsplit(x, "_")[[1]][1], character(1)))

# Re-number so that models with higher weights plot on top
smooth_data <- smooth_data |> filter(weight > 0)
smooth_data$model[smooth_data$component == "binomial"] <- rank(w_binom)[smooth_data$model[smooth_data$component == "binomial"]]
smooth_data$model[smooth_data$component == "tweedie"] <- rank(w_tw)[smooth_data$model[smooth_data$component == "tweedie"]]

smooth_plot <- smooth_data |> 
  mutate(
    component = ifelse(component == "binomial", "s(log odds of occurrence)", "s(log CPUE kg/km2)"),
    predictor = vapply(predictor, \(x) strsplit(x, "_")[[1]][1], character(1)), 
  ) |> 
  ggplot(aes(x, fit)) + 
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, fill = weight, group = model), alpha = 0.5) + 
  geom_line(aes(color = weight, group = model)) + 
  geom_rug(aes(y = -Inf), data = data_long, alpha = 0.5) + 
  facet_grid(component~predictor, scales = "free", switch = "both") + 
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "Blues")[-1]) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Blues")[-1]) + 
  guides(
    color = guide_colorbar("model weight", title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "black", ticks.colour = "black"), 
    fill = guide_colorbar("model weight", title.position = "top", title.hjust = 0.5, barwidth = unit(15, "lines"), frame.colour = "black", ticks.colour = "black") 
  ) + 
  theme_bw() + 
  theme(
    strip.placement = "outside", 
    strip.background = element_blank(),
    axis.title = element_blank(), 
    legend.position = "bottom"
  )

ggsave(paste0(save_dir, "model_smoothers.png"), smooth_plot, height = 5, width = 10, units = "in", dpi = 500)

# Average cold pool SVC -----------------------------------------------------------------

# Identify models with cold pool SVCs
svc_models <- which(sapply(mod_forms, \(x) "cold_pool_2C" %in% all.vars(x)))

if (any(w_binom[svc_models] > 0) | any(w_tw[svc_models] > 0)) {
  
  # Create EBS grid to predict on, holding other covariates at their means
  ebs <- get_ebs_shapefile()
  ebs_grid <- ebs |> st_make_grid(10000) |> st_intersection(ebs) |> st_as_sf()
  
  ebs_grid <- ebs_grid |> mutate(
    centroid = st_centroid(x), 
    X = purrr::map_dbl(centroid, ~st_coordinates(.x)[,1]/1000),
    Y = purrr::map_dbl(centroid, ~st_coordinates(.x)[,2]/1000),
    temp_bottom5m = mean(model_data$temp_bottom5m),
    oxygen_bottom5m = mean(model_data$oxygen_bottom5m), 
    pH_bottom5m = mean(model_data$pH_bottom5m), 
    area_swept_km2 = area_avg
  )
  
  # Cold pool extents to predict for
  cold_pool_levels <- c(0, 0.35, 0.7)
  
  svc_data <- vector("list", length(svc_models))
  names(svc_data) <- names(svc_models)
  
  # Loop over extents, store fitted values
  for (i in 1:length(cold_pool_levels)) {
    
    svc_data[[i]] <- ebs_grid |> dplyr::select(x) |> mutate(extent = cold_pool_levels[i])
    
    ebs_grid <- ebs_grid |> mutate(cold_pool_2C = cold_pool_levels[i])
    
    if (any(w_binom[svc_models] > 0)) {
      binom_fit <- simplify2array(lapply(binom_models[svc_models], \(x) predict(x, newdata = ebs_grid, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
      binom_se <- simplify2array(lapply(binom_models[svc_models], \(x) predict(x, newdata = ebs_grid, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
      binom_se <- weighted_se(binom_fit, binom_se, w_binom[svc_models])
      binom_fit <- apply(binom_fit, 1, weighted.mean, w = w_binom[svc_models])
      svc_data[[i]] <- svc_data[[i]] |> mutate(p_occurrence = binom_fit, p_occurrence_se = binom_se)
    }
    
    if (any(w_tw[svc_models] > 0)) {
      tw_fit <- simplify2array(lapply(tw_models[svc_models], \(x) predict(x, newdata = ebs_grid, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE)))
      tw_se <- simplify2array(lapply(tw_models[svc_models], \(x) predict(x, newdata = ebs_grid, type = "response", exclude = "s(year_chr)", newdata.guaranteed = TRUE, se.fit = TRUE)$se.fit))
      tw_se <- weighted_se(tw_fit, tw_se, w_tw[svc_models])
      tw_fit <- apply(tw_fit, 1, weighted.mean, w = w_tw[svc_models]) 
      svc_data[[i]] <- svc_data[[i]] |> mutate(biomass_fit = tw_fit, biomass_se = tw_se)
    }
    
  }
  
  svc_data <- do.call("rbind", svc_data)
  
  saveRDS(svc_data, paste0(save_dir, "cold_pool_SVC.rds"))
  
}

# Annual Summaries ----------------------------------------------------------------------

range_summary <- unique(ROMS_full[,c("year", "sim")])

range_summary[,c(
  "eastings_mean", "eastings_sd", "northings_mean", "northings_sd",
  "eastings_mean_occ", "eastings_sd_occ", "northings_mean_occ", "northings_sd_occ", 
  "area_occupied_mean", "area_occupied_median", "area_occupied_2.5", "area_occupied_25", "area_occupied_75", 
  "area_occupied_97.5", "mahalanobis_binom", "mahalanobis_tweedie"
)] <- NA

# Threshold for predicting presences from probability of occurrence
kappa <- vapply(binom_models, \(x) {
  p <- predict(x, type = "response")
  dismo::threshold(dismo::evaluate(c(p[model_data$present == 1]), c(p[model_data$present == 0])))$kappa
}, numeric(1))
kappa <- c(dismo::threshold(dismo::evaluate(fit_obs$p_occurrence[fit_obs$cpue_kgkm2 > 0], fit_obs$p_occurrence[fit_obs$cpue_kgkm2 == 0]))$kappa, kappa)
names(kappa) <- c( "ensemble", paste0("model", 1:length(binom_models)))
saveRDS(kappa, paste0(save_dir, "kappa.rds"))

for (i in 1:nrow(range_summary)) {
  
  row_sub <- which(ROMS_full$year == range_summary$year[i] & ROMS_full$sim == range_summary$sim[i])
  
  ## Northings and eastings standard deviation
  Xw_sim <- apply(samples$tweedie[row_sub,], MARGIN = 2, \(x) weighted.mean(ROMS_full$X[row_sub], x))
  Yw_sim <- apply(samples$tweedie[row_sub,], MARGIN = 2, \(x) weighted.mean(ROMS_full$Y[row_sub], x))
  
  ## Eastings mean
  range_summary$eastings_mean[i] <- weighted.mean(ROMS_full$X[row_sub], ROMS_fit$biomass_fit[row_sub])
  range_summary$eastings_sd[i] <- sd(Xw_sim)

  ## Northings mean
  range_summary$northings_mean[i] <- weighted.mean(ROMS_full$Y[row_sub], ROMS_fit$biomass_fit[row_sub])
  range_summary$northings_sd[i] <- sd(Yw_sim)
  
  ## Northings and eastings standard deviation, occurrence-weighted
  Xw_sim <- apply(samples$binomial[row_sub,], MARGIN = 2, \(x) weighted.mean(ROMS_full$X[row_sub], x))
  Yw_sim <- apply(samples$binomial[row_sub,], MARGIN = 2, \(x) weighted.mean(ROMS_full$Y[row_sub], x))
  
  ## Eastings mean, occurrence-weighted
  range_summary$eastings_mean_occ[i] <- weighted.mean(ROMS_full$X[row_sub], ROMS_fit$p_occurrence[row_sub])
  range_summary$eastings_sd_occ[i] <- sd(Xw_sim)
  
  ## Northings mean, occurrence-weighted
  range_summary$northings_mean_occ[i] <- weighted.mean(ROMS_full$Y[row_sub], ROMS_fit$p_occurrence[row_sub])
  range_summary$northings_sd_occ[i] <- sd(Yw_sim)
  
  ## Area occupied mean and 95% confidence interval
  range_summary$area_occupied_mean[i] <- sum(ROMS_fit$p_occurrence[row_sub] > kappa["ensemble"])/length(row_sub)
  area_occ_sim <- apply(samples$binomial[row_sub,], 2, \(x) sum(x)/length(x))
  range_summary$area_occupied_2.5[i] <- quantile(area_occ_sim, 0.025)
  range_summary$area_occupied_25[i] <- quantile(area_occ_sim, 0.25)
  range_summary$area_occupied_median[i] <- median(area_occ_sim)
  range_summary$area_occupied_75[i] <- quantile(area_occ_sim, 0.75)
  range_summary$area_occupied_97.5[i] <- quantile(area_occ_sim, 0.975)
  
  ## Mahalanobis distance - binomial component
  range_summary$mahalanobis_binom[i] <- mean(ROMS_fit$mahalanobis_binom[row_sub])
  
  ## Mahalanobis distance - tweedie component
  range_summary$mahalanobis_tweedie[i] <- mean(ROMS_fit$mahalanobis_tweedie[row_sub])
  
}

## Fitted values corresponding to observed data
station_fit <- ROMS_fit |> ungroup() |> 
    mutate(row = 1:nrow(ROMS_fit)) |> 
    dplyr::select(year, station_id, p_occurrence, biomass_fit, row) |> 
    right_join(dplyr::select(model_data, station_id, year, cpue_kgkm2, X, Y), by = c("station_id", "year"))

## Compute northings, eastings, and area occupied for estimates corresponding to observed data
range_summary_obs <- station_fit |> 
  group_by(year) |> 
  summarize(
    northings_survey_mean = weighted.mean(Y, biomass_fit), 
    northings_survey_mean_occ = weighted.mean(Y, p_occurrence),
    eastings_survey_mean  = weighted.mean(X, biomass_fit),
    eastings_survey_mean_occ  = weighted.mean(X, p_occurrence),
    area_occupied_survey_mean = sum(p_occurrence > 0.5)/n(), 
    .groups = "drop"
  )

range_summary_obs$area_occupied_survey_median <- by(samples$binomial[station_fit$row,], station_fit$year, \(x) median(apply(x, 2, \(y) sum(y > 0.5)/length(y))))

## Join to range summary
range_summary <- range_summary |> left_join(range_summary_obs, by = "year")

write.csv(range_summary, paste0(save_dir, "range_summary.csv"), row.names = FALSE)

file.remove(paste0(save_dir, "running"))
file.create(paste0(save_dir, "complete"))
