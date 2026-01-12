# Layers for historical peat delineation

library(terra)
library(magrittr)
library(dplyr)
library(stringr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/covariates/")

tmpfolder <- paste0(dir_dat, "/Temp/") %T>% dir.create()

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates_10m/")

cov_files <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

terraOptions(tempdir = tmpfolder)

# Load points

Jupiter_pts <- paste0(
  root, 
  "/Historical_peat_delineation_2025/Input_points/Jupiter_pts_processed.rds"
  ) %>%
  readRDS() %>%
  vect(
    geom = c("x", "y"), 
    crs = mycrs, 
    keepgeom = TRUE
  )

ochre_pts <- paste0(
  root, 
  "/Historical_peat_delineation_2025/Input_points/ochre_pts_processed.rds"
) %>%
  readRDS() %>%
  vect(
    geom = c("x", "y"), 
    crs = mycrs, 
    keepgeom = TRUE
  )

LU_1700_pts <- paste0(
  dir_dat,
  "/LU_18thcentury_points/LU_points_18thCentury.shp"
  ) %>% vect()

plot(Jupiter_pts, "is_peat")
plot(ochre_pts, "is_peat")
plot(LU_1700_pts, "LU_txt")

# Load covariates

cov_all <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  rast()

names(cov_all)

covnames_topo_clim <- c(
  "chelsa_bio01_1981_2010_10m",
  "chelsa_bio12_1981_2010_10m",
  "convergence_index",
  "cos_aspect_radians",
  "cross_sectional_curvature",
  "detrended_3_mean",
  "dhm2015_terraen_10m",
  "flooded_depth_10m_mean",
  "flow_accumulation",
  "hillyness",
  "longitudinal_curvature",
  "maximal_curvature",
  "mid_slope_positon",
  "minimal_curvature",
  "normalized_height",
  "positive_openness",
  "profile_curvature",
  "profile_curvature2", 
  "rvb_bios",
  "rvb_fot",
  "saga_wetness_index",
  "sin_aspect_radians",
  "slope",
  "slope_height",
  "standardized_height",
  "tangential_curvature",
  "total_curvature",
  "valley_depth",
  "vdtochn",                             
  "vdtochngt0"
)

cov_names_selected <- c(
  cov_all %>%
    names() %>%
    str_subset(pattern = "georeg", negate = FALSE),
  cov_all %>%
    names() %>%
    str_subset(pattern = "geology", negate = FALSE),
  cov_all %>%
    names() %>%
    str_subset(pattern = "landscape", negate = FALSE),
  cov_all %>%
    names() %>%
    str_subset(pattern = "ogc_", negate = FALSE),
  covnames_topo_clim
)

# Add more covariates for LU?


# Check number of background samples for JUpiter/ochre DB



# END