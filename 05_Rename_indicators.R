# Rename indicators
# Pad zeroes
# Append names
# Basemap doesn't require renaming

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

tmpfolder <- paste0(dir_dat, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs


# Codes and class names for each source map

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

codes_georeg <- data.frame(
  code = c(1, 2, 3, 4, 6, 8, 10),
  name = c(
    "nw_jutland", "ne_jutland", "w_jutland", "s_cent_jutland",
    "n_cent_jutland", "e_denmark", "s_denmark"
  )
)

codes_georeg

geology_ind <- grepl(
  "geology",
  cov_files
)

fuzzy_geology_files <- cov_files[geology_ind]

fuzzy_geology_files

landscape_ind <- grepl(
  "landscape",
  cov_files
)

fuzzy_landscape_files <- cov_files[landscape_ind]

georeg_ind <- grepl(
  "georeg_",
  cov_files
)

fuzzy_georeg_files <- cov_files[georeg_ind]

lu_ind <- grepl(
  "lu_",
  cov_files
)

fuzzy_lu_files <- cov_files[lu_ind]



# Functions for renaming

get_indicator_name <- function(files) {
  out <- files %>%
    basename() %>%
    file_path_sans_ext() %>%
    strsplit("_") %>%
    do.call(rbind, .) %>%
    magrittr::extract(1 , 1:2) %>%
    as.vector() %>%
    paste(collapse = "_")
  return(out)
}

get_indicator_numbers <- function(files) {
  out <- files %>%
    basename() %>%
    file_path_sans_ext() %>%
    strsplit("_") %>%
    do.call(rbind, .) %>%
    magrittr::extract(, 3) %>%
    as.numeric()
  return(out)
}

rename_indicator_files <- function(
    files,
    codes,
    numwidth
) {
  dir_files <- files %>%
    dirname() %>%
    magrittr::extract(1) %>%
    paste0(., "/")
  
  catname <- files %>%
    get_indicator_name()
  
  file_numbers <- files %>%
    get_indicator_numbers()
  
  names_sorted <- codes$name[rank(file_numbers)]
  
  numbers_pad <- file_numbers %>%
    str_pad(width = numwidth, pad = "0")
  
  newfiles <- paste0(
    dir_files, catname, "_", numbers_pad, "_", names_sorted, 
    ".tif"
  )
  
  return(newfiles)
}

## Rename fuzzy indicators ##

# Rename fuzzy geology

fuzzy_geology_newfiles <- rename_indicator_files(
  fuzzy_geology_files,
  codes = codes_geology,
  numwidth = 2
)

fuzzy_geology_newfiles

file.rename(
  fuzzy_geology_files,
  fuzzy_geology_newfiles
)

# Rename fuzzy georegions

fuzzy_georeg_newfiles <- rename_indicator_files(
  fuzzy_georeg_files,
  codes = codes_georeg,
  numwidth = 2
)

fuzzy_georeg_newfiles

file.rename(
  fuzzy_georeg_files,
  fuzzy_georeg_newfiles
)

# Rename fuzzy landscapes

fuzzy_landscape_newfiles <- rename_indicator_files(
  fuzzy_landscape_files,
  codes = codes_landscape,
  numwidth = 2
)

fuzzy_landscape_newfiles

file.rename(
  fuzzy_landscape_files,
  fuzzy_landscape_newfiles
)

# Rename fuzzy lu

fuzzy_lu_newfiles <- rename_indicator_files(
  fuzzy_lu_files,
  codes = codes_lu,
  numwidth = 1
)

fuzzy_lu_newfiles

file.rename(
  fuzzy_lu_files,
  fuzzy_lu_newfiles
)

## Rename crisp indicators ##

dir_crisp <- dir_dat %>%
  paste0(., "covariates_removed/")

crisp_files_all <- dir_crisp %>%
  list_files_with_exts(exts = "tif") %>%
  str_subset(pattern = "crisp")

# Rename crisp geology

crisp_geology_files <- crisp_files_all %>%
  str_subset(pattern = "geology")

crisp_geology_newfiles <- rename_indicator_files(
  crisp_geology_files,
  codes = codes_geology,
  numwidth = 2
)

crisp_geology_newfiles

file.rename(
  crisp_geology_files,
  crisp_geology_newfiles
)

# Rename crisp georegions

crisp_georeg_files <- crisp_files_all %>%
  str_subset(pattern = "georeg")

crisp_georeg_newfiles <- rename_indicator_files(
  crisp_georeg_files,
  codes = codes_georeg,
  numwidth = 2
)

crisp_georeg_newfiles

file.rename(
  crisp_georeg_files,
  crisp_georeg_newfiles
)

# Rename crisp landscapes

crisp_landscape_files <- crisp_files_all %>%
  str_subset(pattern = "landscape")

crisp_landscape_newfiles <- rename_indicator_files(
  crisp_landscape_files,
  codes = codes_landscape,
  numwidth = 2
)

crisp_landscape_newfiles

file.rename(
  crisp_landscape_files,
  crisp_landscape_newfiles
)

# Rename crisp lu

crisp_lu_files <- crisp_files_all %>%
  str_subset(pattern = "_lu_")

crisp_lu_newfiles <- rename_indicator_files(
  crisp_lu_files,
  codes = codes_lu,
  numwidth = 1
)

crisp_lu_newfiles

file.rename(
  crisp_lu_files,
  crisp_lu_newfiles
)

## Fill out NA values in the crisp rasters using the DEM ##

crisp_files_all <- dir_crisp %>%
  list_files_with_exts(exts = "tif") %>%
  str_subset(pattern = "crisp")

dir_out <- dir_dat %>%
  paste0(., "new_covariates/")

dem_zero <- dem*0

i <- 1

for (i in 1:length(crisp_files_all)) {
  name_i <- basename(crisp_files_all[i]) %>%
    file_path_sans_ext()
  
  rast_i <- rast(crisp_files_all[i])
  
  dtyp_i <- datatype(rast_i)
  
  rast_i %>%
    terra::cover(y = dem_zero) %>%
    terra::mask(
      mask = dem_zero,
      filename = paste0(dir_out, "/", name_i, ".tif"),
      names = name_i,
      datatype = dtyp_i,
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
}

# END