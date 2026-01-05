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

terraOptions(
  tempdir = tmpfolder,
  memfrac = 0.05
)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

source("Fill_raster_gaps.R")
source("Focal_density.R")

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
  c(
    .,
    names_full %>%
      str_subset(pattern = "baresoil", negate = FALSE) %>%
      paste0("soilsuite_", .)
  )

# # Process bare soil files [ok]
#
# r_bare <- files_bare %>%
#   rast()
#
#
# # Fill gaps
#
# for (i in 1:nlyr(r_bare)) {
#   r_filled <- r_bare[[i]] %>%
#     terra::clamp(
#       lower = 0,
#       upper = 10000,
#       values = FALSE,
#     ) %>%
#     fill_gaps_gauss(
#       nsteps = 9,
#       weighted = TRUE
#     )  %>%
#     terra::project(
#       x = .,
#       y = dem,
#       method = "cubicspline",
#       mask = TRUE,
#       threads = 10
#     ) %>%
#     mask(
#       mask = dem
#     )
#
#   names(r_filled) <- newnames_bare[i]
#   varnames(r_filled) <- newnames_bare[i]
#
#   r_filled %>%
#     terra::clamp(
#       lower = 0,
#       upper = 10000,
#       values = TRUE,
#       filename = paste0(dir_out, newnames_bare[i], ".tif"),
#       names = newnames_bare[i],
#       datatype = "INT2U",
#       gdal = "TILED=YES"
#     )
#
#   tmpFiles(remove = TRUE)
# }

# # Process full extent files [ok]
#
# r_full <- files_full %>% rast()
#
# newfiles_full <- dir_out %>%
#   paste0(., newnames_full, ".tif")
#
# for (i in 1:(nlyr(r_full) - 1)) {
# # for (i in 12:(nlyr(r_full) - 1)) {
#   r_filled_i <- r_full[[i]] %>%
#     terra::clamp(
#       lower = 0,
#       upper = 10000,
#       values = FALSE,
#     ) %>%
#     fill_gaps_gauss(
#       nsteps = 3,
#       weighted = TRUE
#     )
#
#   r_resampled_i <- r_filled_i %>%
#     terra::project(
#       x = .,
#       y = dem,
#       method = "cubicspline",
#       mask = TRUE,
#       threads = 10
#     ) %>%
#     mask(
#       mask = dem
#     )
#
#   names(r_resampled_i) <- newnames_full[i]
#   varnames(r_resampled_i) <- newnames_full[i]
#
#   r_resampled_i %>%
#     terra::clamp(
#       lower = 0,
#       upper = 10000,
#       values = TRUE,
#       filename = newfiles_full[i],
#       names = newnames_full[i],
#       datatype = "INT2U",
#       overwrite = TRUE,
#       gdal = "TILED=YES"
#     )
#
#   tmpFiles(remove = TRUE)
# }

# # Bare soil frequency [done]
#
# r_full_freq_resampled <- r_full[[nlyr(r_full)]] %>%
#   terra::project(
#     x = .,
#     y = dem,
#     method = "cubicspline",
#     mask = TRUE,
#     threads = 10
#   ) %>%
#   mask(
#     mask = dem
#   )
#
# names(r_full_freq_resampled) <- newnames_full[nlyr(r_full)]
# varnames(r_full_freq_resampled) <- newnames_full[nlyr(r_full)]
#
# r_full_freq_resampled %>%
#   terra::math(
#     "round",
#     digits = 3,
#     filename = newfiles_full[nlyr(r_full)],
#     names = newnames_full[nlyr(r_full)],
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )


# # Focal density filter for bare soil extent [done]
#
# r_bare_mask <- paste0(
#   dir_soilsuite,
#   "/mask.tif"
# ) %>%
#   rast() %>%
#   subst(
#     from = c(2:3),
#     to = c(NA, NA)
#   )
#
# newname_bare_mask <- "soilsuite_bare_mask_density"
# newfile_bare_mask <- dir_out %>%
#   paste0(., newname_bare_mask, ".tif")
#
# r_bare_mask_dens <- r_bare_mask %>%
#   focal_density(
#     nsteps = 9
#   )
#
# r_bare_mask_dens_resampled <- r_bare_mask_dens$final %>%
#   terra::project(
#     x = .,
#     y = dem,
#     method = "cubicspline",
#     mask = TRUE,
#     threads = 10
#   ) %>%
#   mask(
#     mask = dem
#   )
#
# names(r_bare_mask_dens_resampled) <- newname_bare_mask
# varnames(r_bare_mask_dens_resampled) <- newname_bare_mask
#
# r_bare_mask_dens_resampled %>%
#   terra::math(
#     "round",
#     digits = 3,
#     filename = newfile_bare_mask,
#     names = newname_bare_mask,
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )

# END
