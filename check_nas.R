# Check for NAs in layers

library(terra)
library(magrittr)
library(dplyr)
library(stringr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- root %>%
  paste0(., "/covariates/")

tmpfolder <- paste0(root, "/Temp/") %T>%
  dir.create()

terraOptions(
  tempdir = tmpfolder
)

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates_10m/")

cov_files <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

cov <- rast(cov_files)

cov_nas <- app(
  cov,
  function(x) {
    if(is.na(sum(x, na.rm = TRUE))) {
      out <- NA
    } else {
      out <- sum(is.na(x))
    }
    return(out)
  },
  filename = paste0(tmpfolder, "/cov_nas.tif")
)

# END