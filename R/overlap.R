

#' @title Predator-prey overlap
#' @description Functions to compute spatial predator-prey overlap
#' @rdname overlap
#' @param prey Vector of prey density / biomass values
#' @param pred Vector of predator density / biomass values
#' @param area area of survey corresponding to density / biomass estimates (no units)
#' 
#' @details 
#' Selected functions to compute spatial predator-prey overlap from estimated predator and
#' prey densities. Code is adapted from Carroll et al. (2019), 'A review of methods 
#' for quantifying predator-prey overlap,' Global Ecology and Biogeography. Please see 
#' Carroll et al. for detailed descriptions of provided metrics and their ecological interpretations:
#' - `loc_collocfn`: The local index of collocation
#' - `glob_collocfn`: The global index of collocation
#' - `biomass_overlapfn`: Biomass-weighted overlap
#' - `hurlbert_overlapfn`: Hurlbert's overlap
#' - `AB_overlapfn`: AB ratio
#' - `area_overlapfn`: Area overlap
#' - `range_overlapfn`: Range overlap
#'
#' @return Numeric overlap value
#' @export
loc_collocfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
}

#' @rdname  overlap
#' @export
biomass_overlapfn <- function(prey, pred) {
  sum((prey/max(prey, na.rm = T)) * (pred/max(pred, na.rm = T)), na.rm = T)/sum(prey/max(prey, na.rm = T), na.rm = T)
}

#' @rdname  overlap
#' @export
hurlbert_overlapfn <- function(prey, pred, area) {
  area_occupied <- sum(area[pred > 0 | prey > 0], na.rm = T)
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum((p_pred*p_prey)/(area/area_occupied), na.rm = T)
}

#' @rdname  overlap
#' @export
AB_overlapfn <- function(prey, pred) { 
  mean((pred - mean(pred, na.rm = T)) * (prey - mean(prey, na.rm = T)), na.rm = T)/(mean(pred, na.rm = T) * mean(prey, na.rm = T)) 
}

#' @rdname  overlap
#' @export
bhatta_coeffn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(sqrt(p_prey*p_pred), na.rm = T)
}

#' @rdname  overlap
#' @export
glob_collocfn <- function(prey_x, prey_y, prey, pred_x, pred_y, pred){
  prey_cgx <- sum(prey_x*prey, na.rm = T)/sum(prey, na.rm = T)
  prey_cgy <- sum(prey_y*prey, na.rm = T)/sum(prey, na.rm = T)
  prey_ix <- prey_x - prey_cgx
  prey_iy <- prey_y - prey_cgy
  prey_i <- sqrt(prey_ix^2 + prey_iy^2)
  prey_inert <- sum(prey * (prey_i^2), na.rm = T)/sum(prey, na.rm = T)
  pred_cgx <- sum(pred_x*pred, na.rm = T)/sum(pred, na.rm = T)
  pred_cgy <- sum(pred_y*pred, na.rm = T)/sum(pred, na.rm = T)
  pred_ix <- pred_x - pred_cgx
  pred_iy <- pred_y - pred_cgy
  pred_i <- sqrt(pred_ix^2 + pred_iy^2)
  pred_inert <- sum(pred * (pred_i^2), na.rm = T)/sum(pred, na.rm = T)
  GIC <- (((prey_cgx - pred_cgx)^2+(prey_cgy - pred_cgy)^2)/ (((prey_cgx-pred_cgx)^2+(prey_cgy-pred_cgy)^2)+prey_inert + pred_inert))
  if(!is.na(GIC))
    GIC <- 1-GIC
  else GIC <- 1
  GIC
}

#' @rdname  overlap
#' @export
area_overlapfn <- function(prey, pred, area){
  total_area <- sum(area, na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/total_area
}

#' @rdname  overlap
#' @export
range_overlapfn <- function(prey, pred, area){
  area_prey <- sum(area[prey > 0], na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/area_prey
}
