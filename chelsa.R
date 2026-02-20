# Process chelsa data

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

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

dir_out <- dir_dat %>%
  paste0(., "/new_covariates/") %T>%
  dir.create()

tmpfolder <- paste0(root, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

# Load chelsa data

dir_chelsa <- dir_dat %>%
  paste0(., "/chelsa/")

chelsa_files <- dir_chelsa %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

chelsa_rast <- rast(chelsa_files)

# get extent of denmark

ext_chelsa_dk <- ext(cov) %>%
  extend(20000) %>%
  project(crs(cov), crs(chelsa_rast))

chelsa_crop <- crop(
  chelsa_rast,
  ext_chelsa_dk
)

plot(chelsa_crop[[1]])

as.data.frame(
  chelsa_crop
) %>%
  sample_n(10)

global(
  chelsa_crop,
  "min"
)

global(
  chelsa_crop,
  "max"
)

chelsa_decimals <- c(
  2, 1, 3, 
  1, 2, 2,
  1, 2, 2,
  2, 2, 0,
  1, 1, 1,
  1, 1, 1,
  1
)

chelsa_extraname <- c(
  "MAT", 	"MDTR", 	"ISOT", 	"TS", 	"TWARMM", 	"TCOLDM", "AMTR", 	"TWETQ", 	
  "TDRYQ", 	"TWARMQ", 	"TCOLDQ", 	"AP", 	"PWETM", 	"PDRYM", 	"PS", 	"PWETQ", 	
  "PDRYQ", 	"PWARMQ", 	"PCOLDQ"
)

chelsa_newnames <- paste0(
  "CHELSA_bio", str_pad(1:19, 2, pad = "0"), "_1981_2010_",
  chelsa_extraname
) %>% tolower()

r_na <- rast(ncols = 180, nrows = 180, xmin = 0)
mygaussmat <- focalMat(r_na, c(1, 1), "Gauss")
mygaussmat2 <- focalMat(r_na, c(1, 2), "Gauss")



# xy ratio of cells
918 / 514
# [1] 1.785992

myfilter <- terra::focalMat(r_na, c(10, 22), "Gauss") %>%
  terra::rast() %>%
  terra::aggregate(fact = c(5, 9), fun = "sum") %>%
  terra::values() %>%
  as.vector() %>%
  matrix(ncol = 9) %>%
  round(3)

chelsa_crop[[1]] %>%
  focal(myfilter) %>%
  plot()

for (i in 1:length(chelsa_newnames)) {
  chelsa_crop[[i]] %>%
    focal(myfilter) %>%
    project(
      dem,
      method = "cubicspline",
      mask = TRUE
    ) %>%
    mask(
      dem
    ) %>%
    math(
      "round",
      digits = chelsa_decimals[i],
      filename = paste0(dir_out, chelsa_newnames[i], ".tif"),
      names = chelsa_newnames[i],
      datatype = "FLT4S",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
}



# END