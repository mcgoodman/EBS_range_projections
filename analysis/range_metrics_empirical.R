
pkgs <- c("here", "dplyr", "tidyr", "purrr", "aclim2sdms")
sapply(pkgs, require, character.only = TRUE)

summarize_range <- function(x) {
  x |> group_by(year) |> 
    summarize(
      centroid_eastings = ifelse(all(cpue_kgkm2 == 0), NA, weighted.mean(X, cpue_kgkm2, na.rm = TRUE)),
      centroid_northings = ifelse(all(cpue_kgkm2 == 0), NA, weighted.mean(Y, cpue_kgkm2, na.rm = TRUE)),
      centroid_eastings_occ = ifelse(all(cpue_kgkm2 == 0), NA, mean(X[present == 1], na.rm = TRUE)),
      centroid_northings_occ = ifelse(all(cpue_kgkm2 == 0), NA, mean(Y[present == 1], na.rm = TRUE)),
      area_occupied = ifelse(all(cpue_kgkm2 == 0), NA, sum(present == 1, na.rm = TRUE)/n())
    )
}

for (i in 1:length(specs$species)) {
  
  species <- specs$species[i]
  survey_file <- specs$survey_file[i]
  threshold <- specs$threshold[i]
  length_bin <- specs$length_bin[i]
  
  if (!is.na(length_bin)) {
    
    source(here("analysis", "load_trawl_size_binned.R"))
    
    save_dir <- paste0(here("output", paste0(gsub(" ", "_", species), paste0("-", length_bin))), "/")
    
    range_empirical <- joined_data |> filter(bin == length_bin) |> summarize_range()
    
  } else {
    
    source(here("analysis", "load_trawl.R"))
    
    save_dir <- paste0(here("output", paste0(gsub(" ", "_", species))), "/")
    
    range_empirical <- summarize_range(joined_data)
  
  }
  
  write.csv(range_empirical, paste0(save_dir, "range_empirical.csv"), row.names = FALSE)
  
}
