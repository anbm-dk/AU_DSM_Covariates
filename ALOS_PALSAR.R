# ALOS/PALSAR

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

dir_alos <- dir_dat %>%
  paste0(., "/AlosPalsar/")

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

files_alos <- dir_alos %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

r_alos <- files_alos %>%
  rast()

# New names

newnames <- files_alos %>%
  basename() %>%
  file_path_sans_ext()

newfiles_string <- tmpfolder %>%
  paste0("/", newnames, ".tif")


# Resample and mask

for (i in 1:nlyr(r_alos)) {
  r_alos[[i]] %>%
    resample(
      y = dem,
      method = "cubicspline",
      threads = 18
    ) %>%
    mask(
      mask = dem,
      filename = newfiles_string[i],
      names = newnames[i],
      datatype = "INT2U",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
}

# Old code
# # HH
# 
# r_hh <- files_alos %>%
#   str_subset(pattern = "HH", negate = FALSE) %>%
#   rast()
# 
# hh_mean <- mean(r_hh, na.rm = TRUE) %>%
#   resample(y = dem) %>%
#   mask(
#     mask = dem,
#     filename = newfiles_string[1],
#     names = newnames[1],
#     datatype = "INT2U",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# 
# hh_sd <- stdev(r_hh, na.rm = TRUE) %>%
#   resample(y = dem) %>%
#   mask(
#     mask = dem,
#     filename = newfiles_string[2],
#     names = newnames[2],
#     datatype = "INT2U",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# 
# # HV
# 
# r_hv <- files_alos %>%
#   str_subset(pattern = "HV", negate = FALSE) %>%
#   rast()
# 
# hv_mean <- mean(r_hv, na.rm = TRUE) %>%
#   resample(y = dem) %>%
#   mask(
#     mask = dem,
#     filename = newfiles_string[3],
#     names = newnames[3],
#     datatype = "INT2U",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# 
# hv_sd <- stdev(r_hv, na.rm = TRUE) %>%
#   resample(y = dem) %>%
#   mask(
#     mask = dem,
#     filename = newfiles_string[4],
#     names = newnames[4],
#     datatype = "INT2U",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )

# END
