
#' Download Bering 10K level 2 gridded ROMS output
#'
#' @param var variable name
#' @param type "hindcast", "forecast", or "historical" run
#' @param scenario "SSP126" or "SSP585" for CMIP6 outputs, "RCP45" or "RCP85" for CMIP5 outputs
#' @param start integer start year. NA defaults to all available years.
#' @param end integer end year. NA defaults to all available years.
#' @param version "CMIP6" or "CMIP5"
#' @param earth_model "GFDL", "CESM", or "MIROC"
#' @param crop_ebs if TRUE, crop to EBS/NBS survey boundary
#' @param write_dir path to write downloaded raster. if NA, object is not written to disk.
#'
#' @return A `stars` object
#' @export
get_level2 <- function(var,
                       type = c("hindcast", "forecast", "historical"),
                       scenario = c("SSP126", "SSP585", "RCP45", "RCP85"),
                       start = NA,
                       end = NA,
                       version = c("CMIP6", "CMIP5"),
                       earth_model = c("GFDL", "CESM", "MIROC"),
                       crop_ebs = TRUE,
                       write_dir = ".") {
  
  require("curl")
  require("dplyr")
  require("stringr")
  require("stars")
  
  # Base directory url for Thredds server
  url <- "https://data.pmel.noaa.gov/aclim/thredds"
  
  # Check argument validity
  type <- match.arg(type)
  scenario <- tolower(match.arg(scenario))
  earth_model <- match.arg(earth_model)
  version <- match.arg(version)
  if(version == "CMIP5" & grepl("SSP", scenario)) stop("CMIP5 options are RCP45 or RCP85")
  if(version == "CMIP6" & grepl("RCP", scenario)) stop("CMIP6 options are SSP126 or SSP585")
  if(version == "CMIP6") earth_model <- tolower(earth_model)
  
  # Derived arguments
  roms_model <- switch(version, CMIP5 = "B10K-H16", CMIP6 = "B10K-K20P19")
  var_base <- strsplit(var, "_")[[1]][1]
  
  # Paste together sim directory based on arguments
  sim <- type |> switch( 
    hindcast = paste0(roms_model, "_CORECFS"), 
    historical = paste(roms_model, version, earth_model, "historical", sep = "_"), 
    forecast = paste(roms_model, version, earth_model, scenario, sep = "_")
  )
  
  # List directories on sim page using html scraping
  sim_page <- suppressWarnings(readLines(paste0(url, "/catalog/files/", sim, "/Level2.html")))
  sim_dirs <- sim_page[grepl("/Level2/", sim_page, fixed = TRUE)]
  sim_dirs <- str_match(sim_dirs, "/Level2/(.*?)/catalog.html")[,2]
  
  # Subset directories to download using start and end years
  if(!is.na(start) | !is.na(end)) {
    if (!is.na(start)) if(!(start == as.integer(start))) stop("start should be an integer year")
    if (!is.na(end)) if(!(end == as.integer(end))) stop("end should be an integer year")
    dir_yrs <- lapply(sim_dirs, \(x) as.integer(strsplit(x, "-")[[1]][1]):as.integer(strsplit(x, "-")[[1]][2]))
    start_sub <- vapply(dir_yrs, \(x) any(ifelse(is.na(start), dir_yrs[[1]][1], start) <= x), logical(1))
    end_sub <- vapply(dir_yrs, \(x) any(x <= ifelse(is.na(end), max(dir_yrs[[length(dir_yrs)]]), end)), logical(1))
    sim_dirs <- sim_dirs[start_sub & end_sub]
  }
  
  # Paste urls for all years in sim
  sim_urls <- paste0(url, "/fileServer/", sim, "/Level2/", sim_dirs, "/", sim, "_", sim_dirs, "_average_", var, ".nc")
  
  # Download data
  dir.create(tmp <- tempdir())
  dir.create(save_dir <- paste(tmp, sim, sep = "/"))
  file_names <- paste0(save_dir, "/", var, "_", sim_dirs, ".nc")
  for (i in 1:length(sim_urls)) {
    tryCatch(
      curl::curl_download(sim_urls[i], file_names[i]),
      error = function(x) {
        stop(paste(
          "Download failed: variable likely not avaialable for selected dataset.\n", 
          "Browse available outputs for selected sim by visitng:\n",
          paste0(url, "/catalog/files/", sim, "/Level2.html")
        ))
      }
    )
  }
  
  # Read in rasters, join together
  suppressWarnings(suppressMessages({
    roms <- read_ncdf(file_names[1], var = var_base)
    roms_dates <- st_get_dimension_values(roms, "ocean_time")
    for (i in 2:length(file_names)) {
      roms_i <- read_ncdf(file_names[i], var = var_base)
      roms <- c(roms, roms_i, along = "ocean_time")
      roms_dates <- c(roms_dates, st_get_dimension_values(roms_i, "ocean_time"))
    }
    st_crs(roms) <- "+proj=longlat +datum=WGS84 +no_defs"
  }))
  
  # Fix issue with `stars` where if first raster is single band, all bands receive same date
  roms <- st_set_dimensions(roms, "ocean_time", values = roms_dates)
  
  # Subset time, if applicable
  if(!is.na(start) | !is.na(end)) {
    if (!is.na(start)) {
      start_date <- as.POSIXct(paste0(start, "-01-01 00:00:00"), tz = "UTC")
      dates_in <- which(st_get_dimension_values(roms, "ocean_time") >= start_date)
      roms <- roms |> slice(dates_in, along = "ocean_time")
    }
    if (!is.na(end)) {
      end_date <- as.POSIXct(paste0(end, "-12-31 23:59:00"), tz = "UTC")
      dates_in <- which(st_get_dimension_values(roms, "ocean_time") <= end_date)
      roms <- roms |> slice(dates_in, along = "ocean_time")
    }
  }
  
  # Crop to eastern Bering Sea survey region, if applicable
  if(crop_ebs) {
    suppressWarnings(suppressMessages({
      require("aclim2sdms")
      ebs <- get_ebs_shapefile() |> st_transform("+proj=longlat +datum=WGS84") |> st_shift_longitude()
      roms <- roms |> st_crop(ebs) 
    }))
  }
  
  # Write out stars object, if applicable
  if(!is.na(write_dir)) {
    write_stars(roms, file.path(write_dir, paste0(sim, "_", var, ".tif")), driver = "GTiff")
  }
  
  unlink(tmp, recursive = TRUE, force = TRUE)
  
  roms
  
}


#' Compute week-of-year means for ROMS output
#'
#' @param x A multiband `stars` object returned by `get_level2`
#' @param start integer start year. NA defaults to all available years.
#' @param end integer end year. NA defaults to all available years.
#'
#' @return A `stars` object
#' @export
weight_weeks <- function(x, start = NA, end = NA) {
  
  require("stars")
  require("lubridate")
  
  # Subset time, if applicable
  if(!is.na(start) | !is.na(end)) {
    if (!is.na(start)) {
      if(!(start == as.integer(start))) stop("start should be an integer year")
      start_date <- as.POSIXct(paste0(start, "-01-01 00:00:00"), tz = "UTC")
      dates_in <- which(st_get_dimension_values(x, "ocean_time") >= start_date)
      x <- x |> slice(dates_in, along = "ocean_time")
    }
    if (!is.na(end)) {
      if(!(end == as.integer(end))) stop("end should be an integer year")
      end_date <- as.POSIXct(paste0(end, "-12-31 23:59:00"), tz = "UTC")
      dates_in <- which(st_get_dimension_values(x, "ocean_time") <= end_date)
      x <- x |> slice(dates_in, along = "ocean_time")
    }
  }
  
  # Extract dates from raster
  days <- as.Date(st_get_dimension_values(x, "ocean_time"))
  
  # Matrix to hold distance between raster dates and week-of-year midpoints
  roms_weights <- matrix(NA_real_, nrow = length(days), ncol = 52)
  
  # Loop over raster layers, store distance between raster layer date and week midpoints
  for (i in 1:length(days)) {
    
    # Midpoints for corresponding year
    wk_mdpts <- seq(as.Date(paste0(year(days[i]), "-01-04")), by = "week", length.out = 52)
    
    # Handle special cases - e.g. distance from January dates to week 52 of last year
    if(month(days)[i] == 1) wk_mdpts[52] <- seq(as.Date(paste0(year(days[i]) - 1, "-01-04")), by = "week", length.out = 52)[52]
    if(month(days)[i] == 12) wk_mdpts[1] <- as.Date(paste0(year(days[i]) + 1, "-01-04"))
    
    # Distance between raster layer date and week midpoints (>7 days set to zero)
    roms_weights[i,] <- pmax(1 - abs((days[i] - wk_mdpts) / 7), 0)
    
  }
  
  # Loop over weeks, store weighted mean of applicable raster layers for each week
  for (i in 1:ncol(roms_weights)) {
    
    # Raster layers to compute weighted mean from, and corresponding weights
    weeks <- which(roms_weights[,i] > 0)
    weights <- roms_weights[weeks,i]
    
    if (i == 1) {
      roms_weekly <- x |> slice(weeks, along = "ocean_time") |> 
        st_apply(1:2, weighted.mean, weights = weights, na.rm = TRUE)
    } else {
      roms_week <- x |> slice(weeks, along = "ocean_time") |> 
        st_apply(1:2, weighted.mean, weights = weights, na.rm = TRUE)
      roms_weekly <- c(roms_weekly, roms_week)
    }
    
  }
  
  # Set time dimension to week-of-year
  roms_weekly <- roms_weekly |> merge(name = "week") |> st_set_dimensions("week", values = 1:52)
  
  # Format names and units and return
  names(roms_weekly) <- names(x)
  if ("units" %in% class(x[[names(x)]])) units(roms_weekly[[names(x)]]) <- units(x[[names(x)]])
  roms_weekly
  
}

#' Delta-correct ROMS level 2 outputs
#'
#' @param x ROMS level 2 `stars` object returned by `get_level2`
#' @param hindcast Weekly means for ROMS level 2 hindcast, returned by `weight_weekly`
#' @param historical Weekly means for ROMS level 2 historical run, returned by `weight_weekly`
#' @param lower If applicable, lower threshold for returned variable
#' @param upper If applicable, upper thresdhold for returned variable
#'
#' @return A `stars` object
#' @export
bias_correct <- function(x, hindcast, historical, lower = NA, upper = NA) {
  
  require("stars")
  
  days <- as.Date(st_get_dimension_values(x, "ocean_time"))
  
  for (i in 1:length(days)) {
    
    # Midpoints for corresponding year
    wk_mdpts <- seq(as.Date(paste0(year(days[i]), "-01-04")), by = "week", length.out = 52)
    if(month(days)[i] == 1) wk_mdpts[52] <- seq(as.Date(paste0(year(days[i]) - 1, "-01-04")), by = "week", length.out = 52)[52]
    if(month(days)[i] == 12) wk_mdpts[1] <- as.Date(paste0(year(days[i]) + 1, "-01-04"))
    
    # Week bin corresponding to date
    wk_bin <- which.min(abs(wk_mdpts - days[i]))
    
    # Delta correction
    if (i == 1) {
      bc_rast <- hindcast[,,,wk_bin, drop = TRUE] + (x[,,,i, drop = TRUE] - historical[,,, wk_bin, drop = TRUE])
    } else {
      bc_rast_day <- hindcast[,,,wk_bin, drop = TRUE] + (x[,,,i, drop = TRUE] - historical[,,, wk_bin, drop = TRUE])
      bc_rast <- c(bc_rast, bc_rast_day)
    }
    
  }
  
  # Set dimensions and attributes
  bc_rast <- bc_rast |> merge(name = "ocean_time") |> st_set_dimensions("ocean_time", days)
  names(bc_rast) <- names(x)
  
  # Deal with bounds
  if (!is.na(lower)) bc_rast[bc_rast < lower] <- lower
  if (!is.na(upper)) bc_rast[bc_rast > upper] <- upper
  
  bc_rast
  
}
