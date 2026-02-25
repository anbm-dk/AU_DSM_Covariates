# Fill gaps with flat values

library(terra)
library(magrittr)
library(tools)
library(stringr)
library(dplyr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- root %>%
  paste0(., "/covariates/")

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates_10m/")

cov_files <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

dir_out <- dir_dat %>%
  paste0(., "/new_covariates/") %T>%
  dir.create()

tmpfolder <- paste0(root, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

# Specify covariates

fill_with_zero <- c(
  "mid_slope_position", "normalized_height", "standardized_height", 	
  "micro_aspsd", 	"micro_demmad", 	"micro_edginess", 	"micro_flowsd", 	
  "micro_nmins", 	"micro_ridge_noise", 	"micro_ridginess", 	"micro_saddles", 
  "micro_slopeaspsd", 	"micro_valleyness", 	"soilsuite_baresoilfrequency" 
)

fill_with_half <- c(
  "micro_ridge_slope_index",
  "micro_ridge_valley_index"
)

# Fill half value covariates

cov_fillhalf <- dir_cov %>%
  paste0(., fill_with_half, ".tif") %>%
  rast()

for (i in 1:length(fill_with_half)) {
  ifel(
    test = is.na(cov_fillhalf[[i]]),
    yes = dem*0 + 0.5,
    no = cov_fillhalf[[i]],
    filename = paste0(
      dir_out, fill_with_half[i], ".tif"
    ),
    names = fill_with_half[i],
    datatype = datatype(cov_fillhalf[[i]]),
    overwrite = TRUE,
    gdal = "TILED=YES"
  )
}

# Fill zero value covariates

cov_fillzero <- dir_cov %>%
  paste0(., fill_with_zero, ".tif") %>%
  rast()

for (i in 1:length(fill_with_zero)) {
  ifel(
    test = is.na(cov_fillzero[[i]]),
    yes = dem*0,
    no = cov_fillzero[[i]],
    filename = paste0(
      dir_out, fill_with_zero[i], ".tif"
    ),
    names = fill_with_zero[i],
    datatype = datatype(cov_fillzero[[i]]),
    overwrite = TRUE,
    gdal = "TILED=YES"
  )
}

# END