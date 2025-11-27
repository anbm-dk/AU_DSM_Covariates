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
  paste0(., "/AlosPalsar/raw")

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

# New names

newnames <- paste0(
  "alos_palsar_2017_2021_sar_", 
  c("HH", "HH", "HV", "HV"),
  c("_mean", "_sd")
  )

newfiles_string <- tmpfolder %>%
  paste0("/", newnames, ".tif")

# HH

r_hh <- files_alos %>%
  str_subset(pattern = "HH", negate = FALSE) %>%
  rast()

hh_mean <- mean(r_hh, na.rm = TRUE) %>%
  resample(y = dem) %>%
  mask(
    mask = dem,
    filename = newfiles_string[1],
    names = newnames[1],
    datatype = "INT2U",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

hh_sd <- stdev(r_hh, na.rm = TRUE) %>%
  resample(y = dem) %>%
  mask(
    mask = dem,
    filename = newfiles_string[2],
    names = newnames[2],
    datatype = "INT2U",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

# HV

r_hv <- files_alos %>%
  str_subset(pattern = "HV", negate = FALSE) %>%
  rast()

hv_mean <- mean(r_hv, na.rm = TRUE) %>%
  resample(y = dem) %>%
  mask(
    mask = dem,
    filename = newfiles_string[3],
    names = newnames[3],
    datatype = "INT2U",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

hv_sd <- stdev(r_hv, na.rm = TRUE) %>%
  resample(y = dem) %>%
  mask(
    mask = dem,
    filename = newfiles_string[4],
    names = newnames[4],
    datatype = "INT2U",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

# END