# Remame layers to match filenames

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

dir_out <- dir_dat %>%
  paste0(., "new_covariates/")

cov_files <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

tmpfolder <- paste0(dir_dat, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

# Check names

cov_filesnames_base_sans_ext <- cov_files %>%
  basename() %>%
  file_path_sans_ext()

cov <- cov_files %>% rast()

cov_layernames <- names(cov)

cbind(cov_filesnames_base_sans_ext, cov_layernames)

names_diff <- cov_filesnames_base_sans_ext != cov_layernames

ind_change <- c(1:length(cov_layernames))[names_diff]

cbind(cov_filesnames_base_sans_ext, cov_layernames)[ind_change, ]

# Rewrite covariate rasters with new layer names

for (i in 1:length(ind_change)) {
  ind_i <- ind_change[i]
  
  rast_i <- rast(cov_files[ind_i])
  
  dtyp_i <- datatype(rast_i)
  
  filename_new_i <- paste0(
    dir_out,
    cov_filesnames_base_sans_ext[ind_i],
    ".tif"
  )
  
  writeRaster(
    rast_i,
    filename = filename_new_i,
    names = cov_filesnames_base_sans_ext[ind_i],
    datatype = dtyp_i,
    overwrite = TRUE,
    gdal = "TILED=YES"
  )
}

# Check names for removed covariates

dir_removed <- dir_dat %>%
  paste0("covariates_removed/")

cov_removed_files <- dir_removed %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

removed_filesnames_base_sans_ext <- cov_removed_files %>%
  basename() %>%
  file_path_sans_ext()

cov_removed <- cov_removed_files %>% rast()

cov_removed_layernames <- names(cov_removed)

cbind(removed_filesnames_base_sans_ext, cov_removed_layernames)

names_removed_diff <- removed_filesnames_base_sans_ext != cov_removed_layernames

ind_change_removed <- c(1:length(cov_removed_layernames))[names_removed_diff]

cbind(
  removed_filesnames_base_sans_ext, 
  cov_removed_layernames)[ind_change_removed, ]


# Rewrite removed rasters with new layer names

for (i in 1:length(ind_change_removed)) {
  ind_i <- ind_change_removed[i]
  
  rast_i <- rast(cov_removed_files[ind_i])
  
  dtyp_i <- datatype(rast_i)
  
  filename_new_i <- paste0(
    dir_out,
    removed_filesnames_base_sans_ext[ind_i],
    ".tif"
  )
  
  writeRaster(
    rast_i,
    filename = filename_new_i,
    names = removed_filesnames_base_sans_ext[ind_i],
    datatype = dtyp_i,
    overwrite = TRUE,
    gdal = "TILED=YES"
  )
}

# END