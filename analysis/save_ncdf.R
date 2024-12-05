
pkgs <- c("dplyr", "stars", "ncdf4", "here")
sapply(pkgs, require, character.only = TRUE)

dir.create(here("output", "ncdf"))

save_ncdf_utm <- function(obj, path, cellsize = 5e+03) {
  
  # Resample on regular UTM-2 grid
  obj <- obj |> st_transform(crs = "+proj=utm +zone=2 +datum=WGS84 +units=m +no_defs")
  obj <- st_warp(obj, crs = st_crs(obj), cellsize = cellsize)
  
  # Time dimension
  origin <- "1970-01-01"
  dates <- as.double(as.Date(st_get_dimension_values(obj, 3)) - as.Date(origin))
  time_units <- paste("days since", origin, "00:00:00.0 -0:00")
  
  # Define dimensions
  x <- ncdim_def("x", units = "m", longname = "UTM zone 2 eastings", st_get_dimension_values(obj, "x"))
  y <- ncdim_def("y", units = "m", longname = "UTM zone 2 eastings", st_get_dimension_values(obj, "y"))
  time <- ncdim_def("time", units = time_units, dates) 
  
  # Define variables and dimensions
  projname <- "transverse_mercator"
  pocc <- ncvar_def("p_occurrence", units = "", list(x, y, time), missval = NA, longname = "probability of occurrence")
  biom <- ncvar_def("log_biomass", units = "", list(x, y, time), missval = NA, longname = "log of estimated biomass")
  pocc_se <- ncvar_def("p_occurrence_se", units = "", list(x, y, time), missval = NA, longname = "SE(probability of occurrence)")
  biom_se <- ncvar_def("log_biomass_se", units = "", list(x, y, time), missval = NA, longname = "SE(log of estimated biomass)")
  proj_def <- ncvar_def(projname, "1", NULL, NULL, longname = projname, prec="char")
  
  # Create netcdf file
  ncout <- nc_create(path, list(pocc, biom, pocc_se, biom_se, proj_def), force_v4 = TRUE)
  
  # Write variables to netcdf file
  ncvar_put(ncout, pocc, obj$p_occurrence)
  ncvar_put(ncout, biom, obj$biomass_fit)
  ncvar_put(ncout, pocc_se, obj$p_occurrence_se)
  ncvar_put(ncout, biom_se, obj$biomass_se)
  
  # Write dimension and variable attributes to netcdf file
  ncatt_put(ncout,"x", "axis","X")
  ncatt_put(ncout,"x", "standard_name", "projection_x_coordinate")
  ncatt_put(ncout,"x", "_CoordinateAxisType", "GeoX")
  ncatt_put(ncout,"y", "axis", "Y")
  ncatt_put(ncout,"y", "standard_name", "projection_y_coordinate")
  ncatt_put(ncout,"y", "_CoordinateAxisType", "GeoY")
  ncatt_put(ncout, "p_occurrence", "grid_mapping", projname)
  ncatt_put(ncout, "p_occurrence_se", "grid_mapping", projname)
  ncatt_put(ncout, "log_biomass", "grid_mapping", projname)
  ncatt_put(ncout, "log_biomass_se", "grid_mapping", projname)
  
  # Write CRS attributes to netcdf file
  ncatt_put(ncout, projname, "name", projname)
  ncatt_put(ncout, projname,"grid_mapping_name", projname)
  ncatt_put(ncout, projname, "longitude_of_central_meridian", -171)
  ncatt_put(ncout, projname, "false_easting", 500)
  ncatt_put(ncout, projname, "false_northing", 0)
  ncatt_put(ncout, projname, "latitude_of_projection_origin", 0)
  ncatt_put(ncout, projname, "scale_factor_at_central_meridian", 0.9996)
  ncatt_put(ncout, projname, "long_name", "CRS definition")
  ncatt_put(ncout, projname, "longitude_of_prime_meridian", 0)
  ncatt_put(ncout, projname, "semi_major_axis", 6378137)
  ncatt_put(ncout, projname, "inverse_flattening", 298.257223563)
  ncatt_put(ncout, projname, "spatial_ref", "PROJCS[‘unknown’,GEOGCS[‘unknown’,DATUM[‘WGS_1984’,SPHEROID[‘WGS 84’,6378137,298.257223563],AUTHORITY[‘EPSG’,’6326’]],PRIMEM[‘Greenwich’,0,AUTHORITY[‘EPSG’,’8901’]],UNIT[‘degree’,0.0174532925199433]],PROJECTION[‘Transverse_Mercator’],PARAMETER[‘latitude_of_origin’,0],PARAMETER[‘central_meridian’,-171],PARAMETER[‘scale_factor’,0.9996],PARAMETER[‘false_easting’,500],PARAMETER[‘false_northing’,0],UNIT[‘metre’,1,AUTHORITY[‘EPSG’,’9036’]],AXIS[‘Easting’,EAST],AXIS[‘Northing’,NORTH]]")
  ncatt_put(ncout, projname, "crs_wkt", "PROJCS[‘unknown’,GEOGCS[‘unknown’,DATUM[‘WGS_1984’,SPHEROID[‘WGS 84’,6378137,298.257223563],AUTHORITY[‘EPSG’,’6326’]],PRIMEM[‘Greenwich’,0,AUTHORITY[‘EPSG’,’8901’]],UNIT[‘degree’,0.0174532925199433]],PROJECTION[‘Transverse_Mercator’],PARAMETER[‘latitude_of_origin’,0],PARAMETER[‘central_meridian’,-171],PARAMETER[‘scale_factor’,0.9996],PARAMETER[‘false_easting’,500],PARAMETER[‘false_northing’,0],UNIT[‘metre’,1,AUTHORITY[‘EPSG’,’9036’]],AXIS[‘Easting’,EAST],AXIS[‘Northing’,NORTH]]")
  ncatt_put(ncout, projname, "GeoTransform", "-192.6105578478245 5 0 7461.274939415367 0 -5")
  ncatt_put(ncout, projname, "_CoordinateAxisTypes","GeoX GeoY")
  
  # Write global attributes to netcdf file
  ncatt_put(ncout, 0, "title", basename(tools::file_path_sans_ext(path)))
  ncatt_put(ncout, 0, "history", paste("M.C. Goodman,", date()))
  ncatt_put(ncout, 0, "Conventions", "CF-1.5")
  ncatt_put(ncout, 0, "GDAL", "GDAL 3.8.4, released 2024/02/08")
  
  # Close connection to netcdf file
  nc_close(ncout)
  
}

projection_files <- list.files("output", pattern = ".rds", recursive = TRUE)
projection_files <- projection_files[grepl("level2", projection_files)]

for (i in 1:length(projection_files)) {
  
  cat("writing netcdf", i, "/", length(projection_files), "\r")
  projection <- readRDS(here("output", projection_files[i]))
  filename <- here("output", "ncdf", gsub(".rds", ".nc", gsub("/level2_projections/", "-", projection_files[i])))
  suppressMessages(suppressWarnings(save_ncdf_utm(projection, filename)))
  
}
