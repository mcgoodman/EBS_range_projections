
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