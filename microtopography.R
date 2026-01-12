# Microtography

# Notes
# I am using 1.25 m resolution instead of the original 0.4 resolution, as the the
# highest resolution seems to contain linear artifacts. Maybe related to
# overlaps between flights?
# Also, resampling the DEM to 1.25  m increases the speed of computation.

# I should first aggregate the intermediate results to 5 m resolution.
# The merge the tiles, smooth the rasters and aggregate to 10 m.

library(terra)
library(magrittr)
library(flowpackage)

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
# for (i in 14) {
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
    r1p25 <- terra::rast(terra::ext(demtile0), resolution = 1.25)

    agg_5m <- function(x, fun_focal1, fun_focal2, fun_agg, decimals) {
      x %>%
        terra::focal(mygaussmat, fun = fun_focal1, na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = fun_agg) %>%
        terra::focal(mygaussmat, fun = fun_focal2, na.rm = TRUE) %>%
        terra::aggregate(fact = 2, fun = fun_agg) %>%
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
      terra::aggregate(fact = 2, fun = "mean") %>%
      terra::focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
      terra::resample(r1p25, "bilinear")

    # Local depressions
    focalmin <- terra::focal(demtile, myfocalmat, fun = "min", na.rm = TRUE)
    focaln <- terra::focal(
      demtile*0 + 1,
      3,
      fun = "sum",
      na.rm = FALSE
      )
    lessthanmin <- (demtile < focalmin)*(1-is.na(focaln))
    agg_n <- (!is.na(focaln)*1) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "sum",
        decimals = 3
      )
    agg_mins <- lessthanmin %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "sum",
        decimals = 3
      )
    mins <- (agg_mins * (5^2/1.25^2)) / agg_n

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
    tileslope <- terra::terrain(demtile, "slope", unit = "radians")
    tileslope_sin <- tileslope %>%
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

    valleyness <- maxgab %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      )
    ridginess <- mingab %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      )

    ridge_noise <- outras_stack %>%
      abs() %>%
      mean() %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      )

    make_negative <- function(x) {
      out <- x * -1
      return(out)
    }
    
    # The multiplication and subtraction use arbitrary parameters
    param_mult1 <- 2
    param_add <- -2
    param_mult2 <- 0.5
    
    mylogodds <- (
      (
        terra::clamp(log(sdgab), lower = -15) - 
          terra::clamp(
            log(tan(tileslope_smooth)), lower = -10
          )*param_mult1 + param_add
        )*param_mult2
    ) %>%
      terra::subst(from = NaN, to = 0) %>%
      terra::subst(from = Inf, to = 10) %>%
      terra::subst(from = -Inf, to = -10) 

    ridge_slope_index <- (exp(mylogodds)/(1 + exp(mylogodds))) %>%
      mask(mask = tileslope) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 2
      )

    ridge_valley_index <- (maxgab / (maxgab + mingab)) %>%
      terra::subst(from = NaN, to = 0.5) %>%
      terra::subst(from = Inf, to = 0.5) %>%
      terra::subst(from = -Inf, to = 0.5) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 3
      )
    
    saddles_abs_min <- min(
      max(mingab, 0),
      max(maxgab, 0)
    )
    saddles_abs_max <- max(
      max(mingab, 0),
      max(maxgab, 0)
    )
  
    saddles <- (1-((saddles_abs_max - saddles_abs_min)/saddles_abs_max)) %>%
      terra::subst(from = NaN, to = 0.5) %>%
      terra::subst(from = Inf, to = 0.5) %>%
      terra::subst(from = -Inf, to = 0.5) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      )
    
    edginess <- outras_stack %>%
      abs() %>%
      prod() %>%
      raise_to_power(1 / 8) %>%
      agg_5m(
        fun_focal1 = "mean",
        fun_focal2 = "mean",
        fun_agg = "mean",
        decimals = 5
      )

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

# Load DEM

mycrs <- "EPSG:25832"

dir_cov <- dir_dat %>%
  paste0(., "/covariates_10m/")

cov_files <- dir_cov %>%
  list.files(
    pattern = "\\.tif$",
    full.names = TRUE
  )

dem_ind <- grepl(
  "dhm",
  cov_files
)

dem <- cov_files[dem_ind] %>% rast()

crs(dem) <- mycrs


# Merge tiles, aggregate to 10 m and mask

merge_agg_mask <- function(
    tiledir,
    fun_focal,
    fun_agg,
    decimals,
    mask,
    outname
) {
  tiledir %>%
    list.files(
      pattern = "\\.tif$",
      full.names = TRUE
    ) %>%
    sprc() %>%
    merge() %>%
    agg_5to10m(
      fun_focal = fun_focal,
      fun_agg = fun_agg,
      decimals = decimals
    ) %>%
    terra::mask(
      mask = mask,
      filename = file.path(
        dir_microdem_merged,
        paste0(outname, ".tif")
      ),
      names = outname,
      overwrite = TRUE,
      gdal = "TILED=YES"
    )
}

merge_agg_mask(
  tiledir   = dir_tiles_aspsd,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 3,
  mask      = dem,
  outname   = "micro_aspsd"
)

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
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

merge_agg_mask(
  tiledir   = dir_tiles_nmins,
  fun_focal = "mean",
  fun_agg   = "sum",
  decimals  = 1,
  mask      = dem,
  outname   = "micro_nmins"
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
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

merge_agg_mask(
  tiledir   = dir_tiles_demmad,
  fun_focal = "mean",
  fun_agg   = "median",
  decimals  = 3,
  mask      = dem,
  outname   = "micro_demmad"
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
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

merge_agg_mask(
  tiledir   = dir_tiles_flowsd,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 4,
  mask      = dem,
  outname   = "micro_flowsd"
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
    overwrite = TRUE,
    gdal = "TILED=YES"
  )

merge_agg_mask(
  tiledir   = dir_tiles_slopeaspsd,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_slopeaspsd"
)

merge_agg_mask(
  tiledir   = dir_tiles_valleyness,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_valleyness"
)

merge_agg_mask(
  tiledir   = dir_tiles_ridginess,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_ridginess"
)

merge_agg_mask(
  tiledir   = dir_tiles_ridge_noise,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_ridge_noise"
)

merge_agg_mask(
  tiledir   = dir_tiles_ridge_slope_index,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 2,
  mask      = dem,
  outname   = "micro_ridge_slope_index"
)

merge_agg_mask(
  tiledir   = dir_tiles_ridge_valley_index,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 3,
  mask      = dem,
  outname   = "micro_ridge_valley_index"
)

merge_agg_mask(
  tiledir   = dir_tiles_saddles,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_saddles"
)

merge_agg_mask(
  tiledir   = dir_tiles_edginess,
  fun_focal = "mean",
  fun_agg   = "mean",
  decimals  = 5,
  mask      = dem,
  outname   = "micro_edginess"
)


# END
