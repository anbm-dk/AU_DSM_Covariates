# Rename indicators


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

dir_input <- dir_dat %>%
  paste0(., "/New_covariates_input/")

tmpfolder <- paste0(dir_dat, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

codes_lu <- data.frame(
  code = c(1, 2, 3, 4, 5),
  name = c("urban", "agriculture", "nature_dry", "nature_wet", "water")
)

codes_lu

codes_landscape <- data.frame(
  code = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
  name = c(
    "aeolian", "postglacial_marine", "reclaimed", "lateglacial_marine",
    "weichsel_terminal_moraine", "weichsel_subglacial_moraine",
    "weichsel_kettled_moraine", "weichsel_tunnel_valley",
    "weichsel_outwash", "saalian_moraine", "marsh", "bedrock"
  )
)

codes_landscape

codes_geology <- data.frame(
  code = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  name = c(
    "aeolian_sand", "freshwater_clay", "freshwater_sand", "freshwater_peat",
    "marine_clay", "marine_sand", "till_clay", "till_sand",
    "glaciofluvial_clay", "glaciofluvial_sand", "other_unknown"
  )
)

codes_geology

codes_georegion <- data.frame(
  code = c(1, 2, 3, 4, 6, 8, 10),
  name = c(
    "nw_jutland", "ne_jutland", "w_jutland", "s_cent_jutland",
    "n_cent_jutland", "e_denmark", "s_denmark"
  )
)

codes_georegion

geology_ind <- grepl(
  "geology",
  cov_files
)

geology_files <- cov_files[geology_ind]

geology_files

landscape_ind <- grepl(
  "landscape",
  cov_files
)

landscape_files <- cov_files[landscape_ind]

georeg_ind <- grepl(
  "georeg_",
  cov_files
)

georeg_files <- cov_files[georeg_ind]

lu_ind <- grepl(
  "lu_",
  cov_files
)

lu_files <- cov_files[lu_ind]

## Rename crisp indicators

# Rename geology

# Rename georegions

# Rename landscapes

# Rename lu

# Rename basemap

## Rename fuzzy indicators

# END
