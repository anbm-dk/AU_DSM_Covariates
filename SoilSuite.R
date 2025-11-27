# Processing the SoilSuite bare soil composite

# Bare soil reflectance (and std):
# - Fill holes in layers
# All layers:
# - Resample to 10 m resolution
# - Mask to DEM extent


# 1: Start up

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

dir_soilsuite <- dir_dat %>%
  paste0(., "/LayersDK/")

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

# List files

files_soilsuite <- dir_soilsuite %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  str_subset(pattern = "confidence", negate = TRUE) %>%
  str_subset(pattern = "mask", negate = TRUE) %>%
  str_subset(pattern = "count", negate = TRUE)

files_bare <- files_soilsuite %>%
  str_subset(pattern = "baresoil", negate = FALSE) %>%
  str_subset(pattern = "frequency", negate = TRUE)

names_bare <- basename(files_bare) %>%
  file_path_sans_ext()

newnames_bare <- names_bare %>%
  str_subset(pattern = "baresoil", negate = FALSE) %>%
  paste0("soilsuite_", .)

files_full <- setdiff(files_soilsuite, files_bare)

names_full <- basename(files_full) %>%
  file_path_sans_ext()

newnames_full <- names_full %>%
  str_subset(pattern = "baresoil", negate = TRUE) %>%
  paste0("soilsuite_MREF_", .) %>%
  c(., 
    names_full %>%
      str_subset(pattern = "baresoil", negate = FALSE) %>%
      paste0("soilsuite_", .)
    )

# Process full extent files

r_full <- files_full %>% rast()

r_full_resampled <- r_full %>%
  terra::project(
    x = .,
    y = dem,
    method = "cubicspline",
    mask = TRUE,
    threads = 10
  )

names(r_full_resampled) <- newnames_full

newfiles_full <- dir_out %>%
  paste0(., newnames_full, ".tif")

for (i in 1:(nlyr(r_full_resampled) - 1)) {
    writeRaster(
      r_full_resampled[[i]],
      filename = newfiles_full[i],
      names = newnames_full[i],
      datatype = "INT2U",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
}

r_full_resampled[[nlyr(r_full_resampled)]] %>%
  round(
    digits = 3,
    filename = newfiles_full[nlyr(r_full_resampled)],
    names = newnames_full[nlyr(r_full_resampled)],
    datatype = "FLT4S",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

# Function to fill gaps

fill_gaps_gauss <- function(
    inrast,
    nsteps,
    include_list = FALSE
) {
  r1 <- rast(ncols = 180, nrows = 180, xmin = 0)
  myfilter1 <- round(
    focalMat(r1, c(1, 2), "Gauss"),
    3
  )
  myfilter2 <- myfilter1
  
  smooth_up_list <- list()
  aggregated_list <- list()
  aggregated_list[[1]] <- c(
    inrast*0 + 1,
    inrast
  )
  names(aggregated_list[[1]]) <- c("count", "mean")
  # Stepwise smoothing and aggregation
  for (i in 2:nsteps) {
    smoothed_down <- terra::focal(
      aggregated_list[[i - 1]],
      w = myfilter1,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
    aggregated_list[[i]] <- terra::aggregate(
      smoothed_down,  
      fun = "mean",
      na.rm = TRUE
    )
  }
  # Stepwise disaggregation, merging and smoothing
  smooth_up_list[[nsteps]] <- aggregated_list[[nsteps]]
  for (i in (nsteps - 1):1) {
    # Disaggregate by 2
    splitted <- terra::project(
      x = smooth_up_list[[i + 1]],
      y = aggregated_list[[i]],
      method = "near"
    )
    # Merge with aggregated layers
    merged <- terra::merge(
      x = aggregated_list[[i]],
      y = splitted
    )
    # Smoothing
    smooth_up_list[[i]] <- terra::focal(
      merged,
      w = myfilter2,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
  }
  # Divide mean values by the number of cells, to get a weighted mean
  final_lyr <- smooth_up_list[[1]][[2]] / smooth_up_list[[1]][[1]]
  out <- list()
  # Merge with input layer
  out$final <- terra::merge(
    inrast,
    final_lyr,
    wopt = list(datatype = datatype(inrast))
  )
  out$aggregated_list <- aggregated_list
  out$smooth_up_list <- smooth_up_list
  return(out)
}

# Process bare soil files

r_bare <- files_bare %>% rast()


# Fill gaps

for (j in 1:nlyr(r_bare)) {
  r <- r_bare[[i]]

  r2 <- fill_gaps_gauss(
    r,
    nsteps = 9
  )

  r3 <- mask(
    r2$final,
    dem,
    filename = paste0(dir_out, newnames_bare[i], ".tif"),
    names = newnames_bare[i],
    datatype = "INT2U",
    gdal = "TILED=YES"
  )
  
  rm(r, r2)
  tmpFiles(remove = TRUE)
}

# END