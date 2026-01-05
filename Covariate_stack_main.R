### Collection and standardization of 10 m covariates for digital soil mapping
### Main script for overview and planning.

# To do:

# First steps:
#  1: Fix NA values in fuzzy landscape elements and geology. [ok]
#  2: Add layers from SoilSuite. [ok]
#  3: Add ALOS/PALSAR. [ok]
#  4: ADK wetlands as covariate (fuzzy), based on raw version without landscape
#     elements [ok].

# Second steps:
#  5: Peat probabilities and point probabilities based on Jupiter and the Ochre
#     database. [!]
#  6: Microtopography. [!]
#  7: 100 m DEM derivatives?
#  8: Addtional DEM layers with sea surfaces filled in?

# Wait:
#  9: Add new wetland extent as a covariate (fuzzy?). [wait]
# 10: Alpha Earth layers? (Ask Julian/Sebastian) [wait]
# 11: Peat types from Giri's work. (fuzzy) (Amélie?) [wait]

# Standard processing:
# 12: Standardize names.
# 13: Round off values.
# 14: Fill gaps (especially bare soil products).
# 15: Redo tiles.

# Startup

library(terra)
library(tidyverse)
library(magrittr)
library(tidyr)
library(dplyr)
library(magrittr)
library(tools)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- root %>%
  paste0(., "/covariates/")

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates_10m/")

dir_out <- dir_dat %>%
  paste0(., "/new_covariates/") %T>%
  dir.create()

files_cov <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

r_cov <- files_cov %>%
  rast()


# Check covariate names

names_cov <- names(r_cov)

basenames_cov <- files_cov %>%
  basename() %>%
  file_path_sans_ext()

sum(basenames_cov != names_cov)

names_needfix <- names_cov[basenames_cov != names_cov]
names_fixed <- basenames_cov[basenames_cov != names_cov]

names_needfix
names_fixed

# for (i in 1:length(names_needfix)) {
#   name_i <- names_needfix[i]
#
#   r_i <- terra::subset(r_cov, names_needfix[i])
#
#   names(r_i) <- names_fixed[i]
#
#   datatype(r_i)
#
#   writeRaster(
#     r_i,
#     filename = paste0(dir_out, names_fixed[i], ".tif"),
#     datatype = datatype(r_i),
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# Add new covariates to overview table

covariates_overview_old <- root %>%
  paste0(., "/cov_categories_20240304.csv") %>%
  read.table(header = TRUE, sep = ";")

setdiff(names_cov, covariates_overview_old$name)

setdiff(covariates_overview_old$name, names_cov)

# Style for code

library(styler)

style_dir(getwd())

### END ###
