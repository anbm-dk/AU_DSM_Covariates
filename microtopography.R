# Microtography

# Notes
# I am using 0.8 m resolution instead of the original 0.4 resolution, as the the
# highest resolution seems to contain linear artifacts. Maybe related to
# overlaps between flights?
# Also, resampling the DEM to 0.8  m increases the speed of computation.

# I should first aggregate the intermediate results to 5 m resolution.
# The merge the tiles, smooth the rasters and aggregate to 10 m.

library(terra)
library(magrittr)

dir_code <- getwd()
root <- dirname(dir_code)
dir_dat <- paste0(root, "/covariates/")

tmpfolder <- paste0(dir_dat, "/Temp/") %T>% dir.create()
terraOptions(
  memfrac = 0.02,
  tempdir = tmpfolder
)

# Folder for output tiles

dir_microdem_tiles <- dir_dat %>%
  paste0(., "/microdem_tiles/") %T>%
  dir.create()

dir_tiles_nmins <- dir_microdem_tiles %>%
  paste0(., "/nmins/") %T>%
  dir.create()

dir_tiles_demmad <- dir_microdem_tiles %>%
  paste0(., "/demmad/") %T>%
  dir.create()

dir_tiles_aspsd <- dir_microdem_tiles %>%
  paste0(., "/aspsd/") %T>%
  dir.create()

dir_tiles_flowsd <- dir_microdem_tiles %>%
  paste0(., "/flowsd/") %T>%
  dir.create()

dir_tiles_slopeaspsd <- dir_microdem_tiles %>%
  paste0(., "/slopeaspsd/") %T>%
  dir.create()

dir_tiles_valleyness <- dir_microdem_tiles %>%
  paste0(., "/valleyness/") %T>%
  dir.create()

dir_tiles_ridginess <- dir_microdem_tiles %>%
  paste0(., "/ridginess/") %T>%
  dir.create()

dir_tiles_ridge_noise <- dir_microdem_tiles %>%
  paste0(., "/ridge_noise/") %T>%
  dir.create()

dir_tiles_ridge_slope_index <- dir_microdem_tiles %>%
  paste0(., "/ridge_slope_index/") %T>%
  dir.create()

dir_tiles_ridge_valley_index <- dir_microdem_tiles %>%
  paste0(., "/ridge_valley_index/") %T>%
  dir.create()

dir_tiles_saddles <- dir_microdem_tiles %>%
  paste0(., "/saddles/") %T>%
  dir.create()

dir_tiles_edginess <- dir_microdem_tiles %>%
  paste0(., "/edginess/") %T>%
  dir.create()

dir_microdem_merged <- dir_dat %>%
  paste0(., "/microdem_merged/") %T>%
  dir.create()

# Function to update mosaic

update_running_mosaic <- function(zip_mosaic_path, national_path, layer_name) {
  if (!file.exists(zip_mosaic_path)) {
    return(invisible(NULL))
  }

  if (!file.exists(national_path)) {
    # First chunk: just copy/rename by writing once via terra
    r <- terra::rast(zip_mosaic_path)
    terra::writeRaster(r, national_path, overwrite = TRUE, names = layer_name)
  } else {
    # Merge existing national with new chunk and overwrite national
    r_nat <- terra::rast(national_path)
    r_zip <- terra::rast(zip_mosaic_path)
    r_new <- terra::merge(r_nat, r_zip)
    terra::writeRaster(r_new, national_path, overwrite = TRUE, names = layer_name)
  }

  # Delete the intermediate zip mosaic
  unlink(zip_mosaic_path, force = TRUE)
  invisible(TRUE)
}


# Focal matrix

myfocalmat <- matrix(c(1, 1, 1, 1, NA, 1, 1, 1, 1), nrow = 3)

r_na <- rast(ncols = 180, nrows = 180, xmin = 0)
mygaussmat <- focalMat(r_na, c(1, 1), "Gauss")
mygaussmat2 <- focalMat(r_na, c(1, 2), "Gauss")
mygaussmat5 <- focalMat(r_na, c(2.5, 5), "Gauss")


# Function to aggregate to 5 m resolution
# Include "fun" as an argument (mean/sum/sd)

agg_5m <- function(
  x,
  fun_focal1,
  fun_focal2,
  fun_agg,
  decimals
) {
  out <- x %>%
    focal(mygaussmat, fun = fun_focal1, na.rm = TRUE) %>%
    aggregate(fact = 2, fun = fun_agg) %>%
    focal(mygaussmat, fun = fun_focal2, na.rm = TRUE) %>%
    aggregate(fact = 2, fun = fun_agg) %>%
    focal(mygaussmat, fun = fun_focal2, na.rm = TRUE) %>%
    resample(r5, "bilinear") %>%
    round(decimals)
  return(out)
}

# Function to aggregate from 5 m to 10 m resolution.

agg_5to10m <- function(
  x,
  fun_focal,
  fun_agg,
  decimals
) {
  out <- x %>%
    focal(mygaussmat, fun = fun_focal, na.rm = TRUE) %>%
    aggregate(fact = 2, fun = fun_agg) %>%
    round(decimals)
  return(out)
}


# Function to calculate differences in aspect

differ_aspect <- function(x) {
  library(circhelp)

  asp <- terrain(x, "aspect", unit = "radians")

  dem_focal <- focal(x, mygaussmat, fun = "mean", na.rm = TRUE)

  asp_focal <- dem_focal %>%
    terrain("aspect", unit = "radians")

  asp_diff2 <- angle_diff_rad(asp, asp_focal)

  asp_diff_agg <- asp_diff2 %>%
    agg_5m(
      fun_focal1 = "sd",
      fun_focal2 = "mean",
      fun_agg = "mean",
      decimals = 3
    )

  return(asp_diff_agg)
}

# Build Gabor kernels

library(OpenImageR)
library(magrittr)

init_gb <- GaborFeatureExtract$new()

gb_f <- init_gb$gabor_filter_bank(
  scales = 1, orientations = 8,
  gabor_rows = 3, gabor_columns = 3,
  plot_data = TRUE
)

# Normalize kernels like you do (abs -> mean center -> L1 norm)
gabor_kernels <- lapply(1:8, function(i) {
  k <- gb_f$gabor_real[[i]] %>%
    abs() %>%
    subtract(mean(.)) %>%
    divide_by(sum(abs(.)))
  k
})

# Cluster initialization

library(future)
library(future.apply)
library(parallel)

n_workers <- min(6, max(1, parallel::detectCores() - 1))
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
future::plan(future::cluster, workers = cl)

parallel::clusterExport(
  cl = cl,
  list("tmpfolder"),
  envir = environment()
)

parallel::clusterEvalQ(cl, {
  library(terra)
  library(magrittr)
  library(Rcpp)
  library(circhelp)
  library(flowpackage) # for Rcpp functions
})

# Path for DEMs

dem_zips <- list.files(
  "O:/AUIT_Geodata/Denmark/Digital_elevation_models/Lidar/DHM_2023_working",
  pattern = ".zip",
  full.names = TRUE
)

for (i in 1:length(dem_zips)) {
# for (i in 300:301) {
  unlink(
    list.files(tmpfolder, full.names = TRUE),
    recursive = TRUE,
    force = TRUE
  )

  unzip(
    dem_zips[i],
    exdir = tmpfolder
  )

  rasters <- list.files(
    tmpfolder,
    pattern = "\\.tif$",
    full.names = TRUE
  )

  # Parallel workflow

  # Function to process tiles
  process_tile <- function(j, rasters, mygaussmat, myfocalmat, gabor_kernels) {
    demtile0 <- terra::rast(rasters[j])
    r5 <- terra::rast(terra::ext(demtile0), resolution = 5)

    agg_5m <- function(x, fun_focal1, fun_focal2, fun_agg, decimals) {
      x %>%
        terra::focal(mygaussmat, fun = fun_focal1, na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = fun_agg) %>%
        terra::focal(mygaussmat, fun = fun_focal2, na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = fun_agg) %>%
        terra::focal(mygaussmat, fun = fun_focal2, na.rm = TRUE) %>%
        terra::resample(r5, "bilinear") %>%
        round(decimals)
    }

    differ_aspect <- function(x) {
      asp <- terra::terrain(x, "aspect", unit = "radians")
      dem_focal <- terra::focal(x, mygaussmat, fun = "mean", na.rm = TRUE)
      asp_focal <- terra::terrain(dem_focal, "aspect", unit = "radians")
      asp_diff2 <- circhelp::angle_diff_rad(asp, asp_focal)

      asp_diff2 %>%
        agg_5m(
          fun_focal1 = "sd",
          fun_focal2 = "mean",
          fun_agg = "mean",
          decimals = 3
        )
    }

    demtile <- demtile0 %>%
      terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      terra::aggregate(fact = 2, fun = "mean")

    # Local depressions
    focalmin <- terra::focal(demtile, myfocalmat, fun = "min", na.rm = TRUE)
    lessthanmin <- demtile < focalmin
    mins <- lessthanmin %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "sum",
        decimals = 1
      )

    # Roughness
    focalmeans <- terra::focal(demtile, mygaussmat, fun = "mean", na.rm = TRUE)
    demdiffs <- abs(focalmeans - demtile)
    demmad <- demdiffs %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "median",
        decimals = 3
      )

    # Aspect differences
    aspsd <- differ_aspect(demtile)

    # Local flow accumulation tendency
    tileflow <- terra::focalCpp(
      demtile,
      w = 3,
      fun = localflowCpp,
      fillvalue = NA_real_
    )


    flowsd <- (tileflow^2) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 4
      ) %>%
      sqrt() %>%
      round(2)

    # Slope-aspect SD
    tileslope_sin <- terra::terrain(demtile, "slope", unit = "radians") %>%
      sin()

    tileasp <- terra::terrain(demtile, "aspect", unit = "radians")
    tileaspcos <- cos(tileasp)
    tileaspsin <- sin(tileasp)

    tileasp_focal <- terra::terrain(focalmeans, "aspect", unit = "radians")
    tileaspsin_focal <- sin(tileasp_focal)
    tileaspcos_focal <- cos(tileasp_focal)

    tileslope_smooth <- terra::terrain(focalmeans, "slope", unit = "radians")
    tileslope_sin_smooth <- sin(tileslope_smooth)

    mystack3 <- c(
      tileaspcos, tileaspsin, tileslope_sin,
      tileaspcos_focal, tileaspsin_focal, tileslope_sin_smooth
    )

    angletest5 <- terra::app(mystack3, chordAngleCpp)

    slopeaspsd <- (angletest5^2) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      ) %>%
      sqrt() %>%
      round(3)

    # -------------------
    # Gabor products (integrated)
    # -------------------

    # Apply each Gabor kernel via focal()
    outrasters <- lapply(gabor_kernels, function(k) {
      terra::focal(demtile, w = k)
    })

    outras_stack <- terra::rast(outrasters)

    maxgab <- max(outras_stack)
    mingab <- min(outras_stack) * -1

    # SD across the 8 orientation responses
    sdgab <- terra::app(outras_stack, fun = sd)

    # Smooth/aggregate helper used repeatedly in your gabor block
    smooth_to_r5 <- function(x) {
      x %>%
        terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = "mean") %>%
        terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = "mean") %>%
        terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
        terra::resample(r5, "bilinear")
    }

    valleyness <- smooth_to_r5(maxgab)
    ridginess <- smooth_to_r5(mingab)

    ridge_noise <- outras_stack %>%
      abs() %>%
      mean() %>%
      smooth_to_r5()

    make_negative <- function(x) {
      out <- x * -1
      return(out)
    }

    ridge_slope_index <- (sdgab * log(tileslope_smooth)) %>%
      terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      terra::subst(from = -Inf, to = 0) %>%
      terra::aggregate(fact = 2, fun = "mean") %>%
      terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      terra::aggregate(fact = 2, fun = "mean") %>%
      terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      terra::resample(r5, "bilinear") %>%
      make_negative()

    ridge_valley_index <- (maxgab / (maxgab + mingab)) %>%
      terra::subst(from = NaN, to = 0.5) %>%
      terra::subst(from = Inf, to = 0.5) %>%
      terra::subst(from = -Inf, to = 0.5) %>%
      smooth_to_r5()
    
    saddles <- min(
        max(mingab, 0),
        max(maxgab, 0)
      ) %>%
      smooth_to_r5()
    
    edginess <- outras_stack %>%
      abs() %>%
      prod() %>%
      raise_to_power(1 / 8) %>%
      smooth_to_r5()

    out_paths <- list(
      mins = file.path(tmpfolder, sprintf("tile_%05d_mins.tif", j)),
      demmad = file.path(tmpfolder, sprintf("tile_%05d_demmad.tif", j)),
      aspsd = file.path(tmpfolder, sprintf("tile_%05d_aspsd.tif", j)),
      flowsd = file.path(tmpfolder, sprintf("tile_%05d_flowsd.tif", j)),
      slopeaspsd = file.path(tmpfolder, sprintf("tile_%05d_slopeaspsd.tif", j)),
      valleyness = file.path(tmpfolder, sprintf("tile_%05d_valleyness.tif", j)),
      ridginess = file.path(tmpfolder, sprintf("tile_%05d_ridginess.tif", j)),
      ridge_noise = file.path(tmpfolder, sprintf("tile_%05d_ridge_noise.tif", j)),
      ridge_slope_index = file.path(tmpfolder, sprintf("tile_%05d_ridge_slope_index.tif", j)),
      ridge_valley_index = file.path(tmpfolder, sprintf("tile_%05d_ridge_valley_index.tif", j)),
      saddles = file.path(tmpfolder, sprintf("tile_%05d_saddles.tif", j)),
      edginess = file.path(tmpfolder, sprintf("tile_%05d_edginess.tif", j))
    )

    terra::writeRaster(mins, out_paths$mins, overwrite = TRUE)
    terra::writeRaster(demmad, out_paths$demmad, overwrite = TRUE)
    terra::writeRaster(aspsd, out_paths$aspsd, overwrite = TRUE)
    terra::writeRaster(flowsd, out_paths$flowsd, overwrite = TRUE)
    terra::writeRaster(slopeaspsd, out_paths$slopeaspsd, overwrite = TRUE)
    terra::writeRaster(valleyness, out_paths$valleyness, overwrite = TRUE)
    terra::writeRaster(ridginess, out_paths$ridginess, overwrite = TRUE)
    terra::writeRaster(ridge_noise, out_paths$ridge_noise, overwrite = TRUE)
    terra::writeRaster(ridge_slope_index, out_paths$ridge_slope_index, overwrite = TRUE)
    terra::writeRaster(ridge_valley_index, out_paths$ridge_valley_index, overwrite = TRUE)
    terra::writeRaster(saddles, out_paths$saddles, overwrite = TRUE)
    terra::writeRaster(edginess, out_paths$edginess, overwrite = TRUE)

    # Return only filenames
    return(out_paths)
  }

  res_list <- future_lapply(
    X = seq_along(rasters),
    FUN = process_tile,
    rasters = rasters,
    mygaussmat = mygaussmat,
    myfocalmat = myfocalmat,
    gabor_kernels = gabor_kernels,
    future.seed = TRUE
  )

  tile_files <- res_list

  # Helper: merge a set of tile rasters (filenames) into one zip mosaic file
  merge_tiles_to_zip <- function(paths, out_file, layer_name) {
    paths <- unlist(paths, use.names = FALSE)
    paths <- paths[file.exists(paths)]
    if (length(paths) == 0) {
      return(invisible(NULL))
    }

    # Build SpatRasterCollection from filenames (lazy on disk)
    s <- terra::sprc(paths)

    # Merge and write to disk
    terra::merge(
      s,
      filename = out_file,
      overwrite = TRUE
    ) |> setNames(layer_name)

    invisible(out_file)
  }

  # Build output filenames (one per product per zip)
  zip_nmins_path <- file.path(
    dir_tiles_nmins, paste0("nmins_", sprintf("%03d", i), ".tif"))
  zip_demmad_path <- file.path(
    dir_tiles_demmad, paste0("demmad_", sprintf("%03d", i), ".tif"))
  zip_aspsd_path <- file.path(
    dir_tiles_aspsd, paste0("aspsd_", sprintf("%03d", i), ".tif"))
  zip_flowsd_path <- file.path(
    dir_tiles_flowsd, paste0("flowsd_", sprintf("%03d", i), ".tif"))
  zip_slopeaspsd_path <- file.path(
    dir_tiles_slopeaspsd, paste0("slopeaspsd_", sprintf("%03d", i), ".tif"))
  zip_valleyness_path <- file.path(
    dir_tiles_valleyness, paste0("valleyness_", sprintf("%03d", i), ".tif"))
  zip_ridginess_path <- file.path(
    dir_tiles_ridginess, paste0("ridginess_", sprintf("%03d", i), ".tif"))
  zip_ridge_noise_path <- file.path(
    dir_tiles_ridge_noise, paste0("ridge_noise_", sprintf("%03d", i), ".tif"))
  zip_ridge_slope_idx_path <- file.path(
    dir_tiles_ridge_slope_index,
    paste0("ridge_slope_index_", sprintf("%03d", i), ".tif"))
  zip_ridge_valley_idx_path <- file.path(
    dir_tiles_ridge_valley_index, 
    paste0("ridge_valley_index_", sprintf("%03d", i), ".tif"))
  zip_saddles_path <- file.path(
    dir_tiles_saddles, 
    paste0("saddles_", sprintf("%03d", i), ".tif"))
  zip_edginess_path <- file.path(
    dir_tiles_edginess, 
    paste0("edginess_", sprintf("%03d", i), ".tif"))

  # Merge each product across tiles for this zip
  merge_tiles_to_zip(lapply(tile_files, `[[`, "mins"), 
                     zip_nmins_path, "nmins")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "demmad"), 
                     zip_demmad_path, "demmad")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "aspsd"), 
                     zip_aspsd_path, "aspsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "flowsd"), 
                     zip_flowsd_path, "flowsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "slopeaspsd"), 
                     zip_slopeaspsd_path, "slopeaspsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "valleyness"), 
                     zip_valleyness_path, "valleyness")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridginess"), 
                     zip_ridginess_path, "ridginess")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_noise"), 
                     zip_ridge_noise_path, "ridge_noise")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_slope_index"), 
                     zip_ridge_slope_idx_path, "ridge_slope_index")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_valley_index"), 
                     zip_ridge_valley_idx_path, "ridge_valley_index")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "saddles"), 
                     zip_saddles_path, "saddles")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "edginess"), 
                     zip_edginess_path, "edginess")

  # # Optional: after merging, you can delete the tile-level files to keep temp small
  # # (Only do this if you're sure the merge succeeded.)
  # all_tile_paths <- unique(unlist(tile_files, use.names = FALSE))
  # unlink(all_tile_paths, force = TRUE)

  saveRDS(i, file = paste0(dir_dat, "/microdem_loop_current_i.rds"))
  
  gc()
}

parallel::stopCluster(cl)

# Update mosaic (too timeconsuming do this afterwards)

# Example paths for the running national mosaics (unmasked; mask at the end)
nat_nmins_path <- file.path(
  dir_microdem_merged, "micro_nmins_unmasked.tif"
)
nat_demmad_path <- file.path(
  dir_microdem_merged, "micro_demmad_unmasked.tif"
)
nat_aspsd_path <- file.path(
  dir_microdem_merged, "micro_aspsd_unmasked.tif"
)
nat_flowsd_path <- file.path(
  dir_microdem_merged, "micro_flowsd_unmasked.tif"
)
nat_slopeaspsd_path <- file.path(
  dir_microdem_merged, "micro_slopeaspsd_unmasked.tif"
)
nat_valleyness_path <- file.path(
  dir_microdem_merged, "micro_valleyness_unmasked.tif"
)
nat_ridginess_path <- file.path(
  dir_microdem_merged, "micro_ridginess_unmasked.tif"
)
nat_ridge_noise_path <- file.path(
  dir_microdem_merged, "micro_ridge_noise_unmasked.tif"
)
nat_ridge_slope_idx_path <- file.path(
  dir_microdem_merged, "micro_ridge_slope_index_unmasked.tif"
)
nat_ridge_valley_idx_path <- file.path(
  dir_microdem_merged, "micro_ridge_valley_index_unmasked.tif"
)
nat_saddles_path <- file.path(
  dir_microdem_merged, "micro_saddles_unmasked.tif"
)
nat_edginess_path <- file.path(
  dir_microdem_merged, "micro_edginess_unmasked.tif"
)

# Merge tiles, aggregate to 10 m and mask

dir_tiles_aspsd %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  sprc() %>%
  merge() %>%
  agg_5to10m(
    fun_focal = "mean",
    fun_agg = "mean",
    decimals = 3
  ) %>%
  terra::mask(
    mask = dem,
    filename = file.path(
      dir_microdem_merged,
      "micro_aspsd.tif"
    ),
    names = "micro_aspsd",
    overwrite = TRUE
  )

dir_tiles_nmins %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  sprc() %>%
  merge() %>%
  agg_5to10m(
    fun_focal = "mean",
    fun_agg = "sum",
    decimals = 1
  ) %>%
  terra::mask(
    mask = dem,
    filename = file.path(
      dir_microdem_merged,
      "micro_nmins.tif"
    ),
    names = "micro_nmins",
    overwrite = TRUE
  )

dir_tiles_demmad %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  sprc() %>%
  merge() %>%
  agg_5to10m(
    fun_focal = "mean",
    fun_agg = "median",
    decimals = 3
  ) %>%
  terra::mask(
    mask = dem,
    filename = file.path(
      dir_microdem_merged,
      "micro_demmad.tif"
    ),
    names = "micro_demmad",
    overwrite = TRUE
  )

dir_tiles_flowsd %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  ) %>%
  sprc() %>%
  merge() %>%
  agg_5to10m(
    fun_focal = "mean",
    fun_agg = "mean",
    decimals = 4
  ) %>%
  terra::mask(
    mask = dem,
    filename = file.path(
      dir_microdem_merged,
      "micro_flowsd.tif"
    ),
    names = "micro_flowsd",
    overwrite = TRUE
  )

# Mask and Write results to files


mask_and_write(
  nat_slopeaspsd_path,
  file.path(
    dir_microdem_merged,
    "micro_slopeaspsd.tif"
  ),
  "micro_slopeaspsd"
)

mask_and_write(
  nat_valleyness_path,
  file.path(
    dir_microdem_merged,
    "micro_valleyness.tif"
  ),
  "micro_valleyness"
)

mask_and_write(
  nat_ridginess_path,
  file.path(
    dir_microdem_merged,
    "micro_ridginess.tif"
  ),
  "micro_ridginess"
)

mask_and_write(
  nat_ridge_noise_path,
  file.path(
    dir_microdem_merged,
    "micro_ridge_noise.tif"
  ),
  "micro_ridge_noise"
)

mask_and_write(
  nat_ridge_slope_idx_path,
  file.path(
    dir_microdem_merged,
    "micro_ridge_slope_index.tif"
  ),
  "micro_ridge_slope_index"
)

mask_and_write(
  nat_ridge_valley_idx_path,
  file.path(
    dir_microdem_merged,
    "micro_ridge_valley_index.tif"
  ),
  "micro_ridge_valley_index"
)

mask_and_write(
  nat_saddles_path,
  file.path(
    dir_microdem_merged,
    "micro_saddles.tif"
  ),
  "micro_saddles"
)

mask_and_write(
  nat_edginess_path,
  file.path(
    dir_microdem_merged,
    "micro_edginess.tif"
  ),
  "micro_edginess"
)

# Inspect results from first loop

plot(demtile)
plot(mins)
plot(demmad)
plot(aspsd)
plot(flowsd)

# Plot tiles for second loop

# alltiles_dem <- rasters %>%
#   sprc() %>%
#   merge() %>%
#   aggregate(fact = 25, fun = "mean")
# plot(alltiles_dem)
plot(alltiles_nmins) # 1 decimal
plot(alltiles_demmad) # 3 decimals
plot(alltiles_aspsd) # 3 decimals
plot(alltiles_flowsd) # 2 decimals
plot(alltiles_slopeaspsd) # 3 decimals

# Plot merged results

plot(micro_nmins)
plot(micro_demmad)
plot(micro_aspsd)
plot(micro_flowsd)
plot(micro_flowsd)


# END
