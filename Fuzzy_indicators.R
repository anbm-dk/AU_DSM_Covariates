# Indicator rasters with fuzzy boundaries for categorial covariates

# Process these layers:
# crisp_adk_wetlands.tif - 1:20,000   - smallest units about 15 m across [ok]
# geology_         - 1:25,000   - smallest units about 25 m across [ok]
# landscape_       - 1:100,000  - smallest units about 100 m across [ok]
# georeg_          - 1:100,000  - the uncertainty seems to reach 500 m in some cases [ok]
# lu_              - 10 m - (Corine LU has a scale of 1:00,000, but the basemap has 10 m resolution) [ok]
# (Use a sigma of less than 1 for lu.)
# imk             - 10 m, use half sigma for fuzzification [ok]
# cwl_10m_  # Already processed in ArcGIS (original resolution 20 m) [removed from stack]
# nature types:   - 10 m - (the basemap has 10 m resolution) [  ]

# Also rename crisp indicators to highlight differences.

# lu and imk will only need new names for the crisp layers. They will not
# need to be processed again

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

dir_input <- dir_dat %>%
  paste0(., "/New_covariates_input/")

tmpfolder <- paste0(root, "/Temp/") %T>%
  dir.create()

terraOptions(tempdir = tmpfolder)

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs

my_focal_weights <- focalMat(
  dem,
  c(10, 20),
  type = c("Gauss")
)

# Process (ADK) wetlands layer [ok]

# wetlands_crisp <- dir_input %>%
#   paste0(., "/crisp_adk_wetlands.tif") %>%
#   rast()
#
# names(wetlands_crisp) <- "crisp_adk_wetlands"
#
# drylands_crisp <- 1 - wetlands_crisp
#
# wl_twolayers_crisp <- c(wetlands_crisp, drylands_crisp)
#
#
# wl_twolayers_fuzzy <- focal(
#   wl_twolayers_crisp,
#   w = my_focal_weights,
#   na.policy = "omit",
#   na.rm = TRUE
# )
#
# wl_twolayers_fuzzy_sum <- sum(wl_twolayers_fuzzy)
# wl_fuzzy_norm <- wl_twolayers_fuzzy[[1]] / wl_twolayers_fuzzy_sum
# wl_fuzzy_norm_round <- round(wl_fuzzy_norm, digits = 2)
# names(wl_fuzzy_norm_round) <- "fuzzy_adk_wetlands"
#
# writeRaster(
#   wl_fuzzy_norm_round,
#   filename = paste0(tmpfolder, "/fuzzy_adk_wetlands.tif"),
#   datatype = "FLT4S",
#   overwrite = TRUE,
#   gdal = "TILED=YES"
#   )


# Process geological map [ok]
# 2025-11-24: I only need to rename previous layers [ok]

# geology_ind <- grepl(
#   "geology",
#   cov_files
# )
#
# geology_files <- cov_files[geology_ind] %>%
#   # str_subset(pattern = "crisp", negate = TRUE) %>%
#   str_subset(pattern = "fuzzy", negate = TRUE)
#
# geology_files
#
# geology_crisp <- geology_files %>%
#   rast()
#
# geology_crisp
# names(geology_crisp)
#
# geology_names <- geology_files %>%
#   basename() %>%
#   file_path_sans_ext()
#
# names(geology_crisp) <- geology_names
#
# geology_files_tmp <- paste0(tmpfolder, geology_names, ".tif")
#
# geology_files_tmp
#
# for(i in 1:nlyr(geology_crisp)) {
#   writeRaster(
#     geology_crisp[[i]],
#     filename = geology_files_tmp[i],
#     names = geology_names[i],
#     datatype = "INT2U",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }
#
# # Renaming layers
# # geology_newnames <- geology_files %>%
# #   basename() %>%
# #   file_path_sans_ext() %>%
# #   paste0("crisp_", .)
# #
# # geology_newnames
# #
# # geology_newfiles <- paste0(dir_cov, geology_newnames, ".tif")
# #
# # geology_files
# # geology_newfiles
# #
# # file.rename(
# #   from = geology_files,
# #   to = geology_newfiles
# # )
#
# geology_sum <- sum(geology_crisp)
# geology_res <- 1 - geology_sum
# names(geology_res) <- "geology_res"
# geology_crisp_full <- c(geology_crisp, geology_res)
#
# geology_fuzzy <- focal(
#   geology_crisp_full,
#   w = my_focal_weights,
#   na.policy = "omit",
#   na.rm = TRUE
# )
#
#
# geology_fuzzy_sum <- sum(geology_fuzzy)
# geology_fuzzy_norm <- geology_fuzzy / geology_fuzzy_sum
# geology_fuzzy_norm_round <- round(geology_fuzzy_norm, digits = 2)
# geology_names <- names(geology_fuzzy_norm_round)
# geology_names_fuzzy <- paste0("fuzzy_", geology_names)
# geology_files_fuzzy <- paste0(tmpfolder, geology_names_fuzzy, ".tif")
# names(geology_fuzzy_norm_round) <- geology_names_fuzzy
#
# for (i in 1:nlyr(geology_fuzzy_norm_round)) {
#   writeRaster(
#     geology_fuzzy_norm_round[[i]],
#     filename = geology_files_fuzzy[[i]],
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# Function for remaining layers

fuzzify_indicators <- function(
  x = NULL,
  aggregation_factor = 1,
  residual_layer = TRUE,
  local_filter = NULL,
  n_digits = 3,
  outfolder = NULL
) {
  x_sum <- sum(x)

  if (residual_layer) {
    x_res <- 1 - x_sum
    names(x_res) <- "x_res"
    x_crisp_full <- c(x, x_res)
  } else {
    x_crisp_full <- x
  }

  if (aggregation_factor > 1) {
    x_crisp_full <- terra::aggregate(
      x_crisp_full,
      fact = aggregation_factor,
      fun = "sum",
      na.rm = TRUE # NB!
    )

    x_fuzzy <- focal(
      x_crisp_full,
      w = local_filter,
      na.policy = "all",
      na.rm = TRUE
    )

    x_fuzzy <- terra::resample(
      x = x_fuzzy,
      y = x,
      method = "cubicspline"
    )

    x_fuzzy <- mask(
      x_fuzzy,
      mask = x[[1]]
    )
  } else {
    x_fuzzy <- focal(
      x_crisp_full,
      w = local_filter,
      na.policy = "omit",
      na.rm = TRUE
    )
  }

  x_fuzzy_sum <- sum(x_fuzzy)
  x_fuzzy_norm <- x_fuzzy / x_fuzzy_sum
  x_fuzzy_norm_round <- signif(x_fuzzy_norm, digits = n_digits)
  x_names <- names(x_fuzzy_norm_round)
  # x_names_fuzzy <- paste0("fuzzy_", x_names)
  x_names_fuzzy <- str_replace(
    x_names,
    pattern = "crisp",
    replacement = "fuzzy"
  )
  x_files_fuzzy <- paste0(outfolder, x_names_fuzzy, ".tif")
  names(x_fuzzy_norm_round) <- x_names_fuzzy

  for (i in 1:nlyr(x_fuzzy_norm_round)) {
    writeRaster(
      x_fuzzy_norm_round[[i]],
      filename = x_files_fuzzy[[i]],
      datatype = "FLT4S",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
  }

  invisible(NULL)
}

# Function to rename rasters with "crisp" and write new layers

rename_crisp <- function(
  filenames,
  outfolder = NULL
) {
  x <- filenames %>% rast()

  x_names_crisp <- filenames %>%
    basename() %>%
    file_path_sans_ext() %>%
    paste0("crisp_", .)

  new_files <- paste0(outfolder, x_names_crisp, ".tif")

  for (i in 1:nlyr(x)) {
    writeRaster(
      x[[i]],
      filename = new_files[i],
      names = x_names_crisp[i],
      datatype = "INT2U",
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
  }
  out <- new_files %>% rast()

  return(out)
}

# # Rename landscape elements [ok]
#
# landscape_ind <- grepl(
#   "landscape",
#   cov_files
# )
#
# landscape_files <- cov_files[landscape_ind] %>%
#   str_subset(pattern = "crisp", negate = TRUE) %>%
#   str_subset(pattern = "fuzzy", negate = TRUE)
#
# landscape_crisp <- rename_crisp(
#   landscape_files,
#   outfolder = tmpfolder
#   )
#
# # Process landscape elements [ok]
#
# fuzzify_indicators(
#   landscape_crisp,
#   aggregation_factor = 5,
#   local_filter = my_focal_weights,
#   n_digits = 2,
#   outfolder = tmpfolder
# )

# # Rename georegions [ok]
#
# georeg_ind <- grepl(
#   "georeg_",
#   cov_files
# )
#
# georeg_files <- cov_files[georeg_ind] %>%
#   str_subset(pattern = "crisp", negate = TRUE) %>%
#   str_subset(pattern = "fuzzy", negate = TRUE)
#
# georeg_crisp <- rename_crisp(
#   georeg_files,
#   outfolder = tmpfolder
#   )
#
# # Process georegions [ok]
#
# fuzzify_indicators(
#   georeg_crisp,
#   aggregation_factor = 5,
#   local_filter = my_focal_weights,
#   n_digits = 2,
#   outfolder = tmpfolder
# )


# Sigma experiment

r1 <- matrix(
  sample(c(0, 1), 400, replace = TRUE),
  nrow = 20
) %>%
  rast()

# plot(r1)
#
# r2 <- focal(r1, w = my_focal_weights, na.rm = TRUE)
#
# plot(r2)

# halfsigma <- focalMat(r1, d = c(0.56, 2), type = 'Gauss')

halfsigma <- focalMat(r1, d = c(0.5, 1), type = "Gauss")

halfsigma

# r3 <- focal(r1, w = halfsigma, na.rm = TRUE)
#
# library(viridis)
#
# plot(r1, col = cividis(100))
# plot(r3, col = cividis(100))
#
# hist(r3, breaks = 100)
#
# plot(r3 - r1)
#
# r3 - r1

# Rename LU [ok]

# lu_ind <- grepl(
#   "lu_",
#   cov_files
# )
#
# lu_files <- cov_files[lu_ind] %>%
#   str_subset(pattern = "crisp", negate = TRUE) %>%
#   str_subset(pattern = "fuzzy", negate = TRUE)
#
# lu_crisp <- rename_crisp(
#   lu_files,
#   outfolder = tmpfolder
#   )

# Process LU [ok]

# lu_ind <- grepl(
#   "lu_",
#   cov_files
# )
#
# lu_crisp <- cov_files[lu_ind] %>% rast()
#
# fuzzify_indicators(
#   lu_crisp,
#   local_filter = halfsigma,
#   n_digits = 2,
#   outfolder = tmpfolder
# )

# # Rename IMK [ok]
#
# imk_ind <- grepl(
#   "imk_",
#   cov_files
# )
#
# imk_files <- cov_files[imk_ind] %>%
#   str_subset(pattern = "crisp", negate = TRUE) %>%
#   str_subset(pattern = "fuzzy", negate = TRUE)
#
# imk_crisp <- rename_crisp(
#   imk_files,
#   outfolder = tmpfolder
# )

# Process imk [ok]

# imk_files <- grep(
#   "imk_",
#   cov_files,
#   value = TRUE
# )
#
# for (i in 1:length(imk_files)) {
#   r <- rast(imk_files[i])
#
#   newname <- names(r) %>% paste0("fuzzy_", .)
#
#   r_fuzzy <- x_fuzzy <- focal(
#     r,
#     w = halfsigma,
#     na.policy = "omit",
#     na.rm = TRUE
#   )
#
#   r_fuzzy_round <- round(
#     r_fuzzy,
#     digits = 2
#   )
#
#   names(r_fuzzy_round) <- newname
#
#   writeRaster(
#     r_fuzzy_round,
#     filename = paste0(tmpfolder, "/", newname, ".tif"),
#     datatype = "FLT4S",
#     overwrite = TRUE,
#     gdal = "TILED=YES"
#   )
# }

# Rename_nature types from basemap and write new rasters [ok]
# Clip crisp and fuzzy rasters afterwards using dem extent

dir_nature <- dir_dat %>%
  paste0(., "/BasemapForMapping/")

dir_nature_crisp <- dir_dat %>%
  paste0(., "/basemap_crisp/") %T>%
  dir.create()

dir_nature_fuzzy <- dir_dat %>%
  paste0(., "/basemap_fuzzy/") %T>%
  dir.create()

dir_nature_masked <- dir_dat %>%
  paste0(., "/basemap_both_masked/") %T>%
  dir.create()

nature_names_original <- dir_nature %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

nature_newnames_crisp <- nature_names_original %>%
  str_replace_all(" ", "_") %>%
  basename() %>%
  file_path_sans_ext() %>%
  tolower() %>%
  paste0("crisp_basemap_", .)

nature_newnames_fuzzy <- nature_names_original %>%
  str_replace_all(" ", "_") %>%
  basename() %>%
  file_path_sans_ext() %>%
  tolower() %>%
  paste0("fuzzy_basemap_", .)

nature_newfiles_crisp1 <- nature_newnames_crisp %>%
  paste0(dir_nature_crisp, ., ".tif")

# Load crisp files, rename, rewrite [ok]

nature_crisp <- nature_names_original %>%
  rast()

names(nature_crisp) <- nature_newnames_crisp
varnames(nature_crisp) <- nature_newnames_crisp

# for (i in 1:nlyr(nature_crisp)) {
#     writeRaster(
#       nature_crisp[[i]],
#       filename = nature_newfiles_crisp1[i],
#       datatype = "INT2U",
#       overwrite = TRUE,
#       gdal = "TILED=YES"
#     )
# }

# # Process fuzzy nature types [ok]

nature_crisp1 <- dir_nature_crisp %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  rast()

# fuzzify_indicators(
#   nature_crisp1,
#   local_filter = halfsigma,
#   n_digits = 2,
#   outfolder = dir_nature_fuzzy
# )

# Mask crisp layers [ok]

for (i in 1:nlyr(nature_crisp1)) {
  outname_i <- nature_newfiles_crisp1[i] %>%
    basename() %>%
    file_path_sans_ext()

  masked_i <- terra::mask(
    nature_crisp1[[i]],
    mask = dem
  )

  names(masked_i) <- outname_i
  varnames(masked_i) <- outname_i

  writeRaster(
    masked_i,
    filename = paste0(
      dir_nature_masked, outname_i, ".tif"
    ),
    datatype = "INT2U",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

  tmpFiles(remove = TRUE)
}

# Mask fuzzy layers

nature_fuzzy <- dir_nature_fuzzy %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  rast()
# %>%
#   subset(-nlyr(.))

for (i in 1:nlyr(nature_fuzzy)) {
  outname_i <- nature_newnames_fuzzy[i]

  masked_i <- terra::mask(
    nature_fuzzy[[i]],
    mask = dem
  )

  names(masked_i) <- outname_i
  varnames(masked_i) <- outname_i

  writeRaster(
    masked_i,
    filename = paste0(
      dir_nature_masked, outname_i, ".tif"
    ),
    datatype = "FLT4S",
    overwrite = TRUE,
    gdal = "TILED=YES"
  )
}

# END
