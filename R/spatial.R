
# Functions for reading in package-associated GIS files -------------------------------------------

#' @title EBS survey shapefiles
#' @param region Either "EBS" for full survey region including NBS, or "SEBS"
#' @param type Whether to return multipolygon corresponding to survey grid ("grid") or region boundary only ("boundary")
#' @source \href{https://github.com/afsc-gap-products/akgfmaps}{QTL Archive}
#'
#' @return An "sf" object
#' @export
get_ebs_shapefile <- function(region = c("EBS", "SEBS"), type = c("boundary", "grid")) {
  
  require("sf")
  
  region <- match.arg(region)
  type <- match.arg(type)
  
  option <- paste(region, type, sep = "_")
  
  option |> switch(
    EBS_grid = sf::st_read(system.file(package = "aclim2sdms", "GIS", "EBS grid"), quiet = TRUE),
    EBS_boundary = sf::st_read(system.file(package = "aclim2sdms", "GIS", "EBS boundary"), quiet = TRUE),
    SEBS_grid = sf::st_read(system.file(package = "aclim2sdms", "GIS", "SEBS grid"), quiet = TRUE), 
    SEBS_boundary = sf::st_read(system.file(package = "aclim2sdms", "GIS", "SEBS boundary"), quiet = TRUE)
  )
  
}


#' @title Alakska coastline shapefile
#' @return An "sf" object
#' @export
get_ak_coast <- function() {
  
  require("sf")
  
  sf::st_read(system.file(package = "aclim2sdms", "GIS", "Alaska Shoreline", "ak_russia.shp"))
  
}


#' @title NOAA 15-arcsecond EBS digitial bedrock elevation model
#' @return A "stars" object
#' @export
get_bathymetry <- function() {
  
  require("stars")
  
  stars::read_stars(system.file(package = "aclim2sdms", "GIS", "etopo_bedrock_15arcsecond.tif"))
  
}

#' @title NOAA 15-arcsecond EBS digitial bedrock elevation model
#' @return A "stars" object
#' @export
get_sediment <- function() {
  
  require("stars")
  
  phi <- stars::read_stars(system.file(package = "aclim2sdms", "GIS", "phi.grd"))
  
  phi |> st_warp(crs = "+proj=longlat +datum=WGS84 +no_defs")
  
}


# Functions for manipulating and plotting spatial data --------------------------------------------

#' @title Function to append UTM columns to data frame with lat / lon columns
#' @description The bulk of this function is duplicated from sdmTMB::add_utm_columns
#' @param data data frame
#' @param ll_cols length-2 character vector; names of longitude and latitude columns (in that order)
#' @param ll_crs Coordinate reference system of lat-long coordinates
#' @param utm_names Names of new UTM columns to be appended to data frame
#' @param utm_crs Coordinate reference system for UTM coordinates
#'
#' @return A copy of the input data frame with UTM columns appended
#' @export
add_utm <- function(data, ll_cols = c("longitude", "latitude"), ll_crs = 4326, utm_names = c("X", "Y"), utm_crs) {
  
  coords <- data |> 
    sf::st_as_sf(crs = ll_crs, coords = ll_cols) |> 
    sf::st_transform(utm_crs) |> 
    sf::st_coordinates() |> 
    as.data.frame() |> 
    dplyr::mutate(X = X / 1000, Y = Y / 1000)
  
  data[[utm_names[1]]] <- coords$X
  data[[utm_names[2]]] <- coords$Y
  
  data
  
}


#' @title Simple function to convert longitude from 0/360 to -180/180 and vice-versa
#' @param x A numeric vector of longitudes
#' @param from Whether to convert from -180/180 or 0/360 
#'
#' @return A numeric vector
#' @export
rotate_lon <- function(x, from = c("-180/180", "0/360")) {
  
  from <- match.arg(from)
  
  from |> switch(
    `-180/180` = (x + 360) %% 360, 
    `0/360` = ((x + 180) %% 360) - 180
  )
  
}
