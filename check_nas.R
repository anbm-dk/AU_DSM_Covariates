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

# cov_nas <- app(
#   cov,
#   function(x) {
#     if(is.na(sum(x, na.rm = TRUE))) {
#       out <- NA
#     } else {
#       out <- sum(is.na(x))
#     }
#     return(out)
#   },
#   filename = paste0(tmpfolder, "/cov_nas.tif")
# )

cov_nas <- rast(paste0(tmpfolder, "/cov_nas.tif"))

# as.data.frame(cov_nas)


# Some layer needs to be masked, as it covers areas outside the dem. Which one?
# A random sample doesn't give us much of a hint.

# spatSample(cov, 100) %>%
#   lapply(function(x) {sum(is.na(x))}) %>%
#   unlist() %>%
#   sort()

# Mask using the dem

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

# cov_nas_masked <- terra::mask(
#   cov_nas,
#   mask = dem,
#   filename = paste0(tmpfolder, "/cov_nas_masked.tif")
# )

cov_nas_masked <- rast(paste0(tmpfolder, "/cov_nas_masked.tif"))

# hist(cov_nas_masked)

# cov_nas_masked2 <- ifel(cov_nas_masked == 0, NA, cov_nas_masked)

# cov_nas_masked2

# writeRaster(
#   cov_nas_masked2,
#   filename = paste0(tmpfolder, "/cov_nas_masked2.tif")
# )

cov_nas_masked2 <- rast(paste0(tmpfolder, "/cov_nas_masked2.tif"))

# hist(cov_nas_masked2)

# as.data.frame(cov_nas_masked2)

# cov_masked <- mask(
#   cov,
#   cov_nas_masked2
# )
# 
# cov_notna_summary <- terra::global(cov_masked, "notNA")
# 
# cov_notna_summary
# 
# saveRDS(cov_notna_summary, paste0(dir_dat, "cov_notna_summary.Rds"))

# Inverse mask using dem, to find areas outside dem.

# cov_nas_masked_inverse <- terra::mask(
#     cov_nas,
#     mask = dem,
#     inverse = TRUE,
#     filename = paste0(tmpfolder, "/cov_nas_masked_inverse.tif")
#   )
# 
# spatSample(cov_nas_masked_inverse, 10)

# Apparently, there are no nonNA values outside the dem.

# END