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
  # memfrac = 0.02,
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

dir_microdem_merged <- dir_dat %>%
  paste0(., "/microdem_merged/") %T>%
  dir.create()

# Function to update mosaic

update_running_mosaic <- function(zip_mosaic_path, national_path, layer_name) {
  if (!file.exists(zip_mosaic_path)) return(invisible(NULL))
  
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

myfocalmat <- matrix(c(1,1,1,1,NA,1,1,1,1), nrow = 3)

r_na <- rast(ncols = 180, nrows = 180, xmin=0)
mygaussmat <- focalMat(r_na, c(1, 1), "Gauss")
mygaussmat2 <- focalMat(r_na, c(1, 2), "Gauss")
mygaussmat5 <- focalMat(r_na, c(2.5, 5), "Gauss")

# Local (pseudo) flow accumulation

localflow <- function(x) {
  x_center <- x[ceiling(length(x)/2)]
  out <- sum((x > x_center), na.rm = TRUE) - sum((x < x_center), na.rm = TRUE)
  return(out)
}

# C++ implementation of the above

library(Rcpp)

cppFunction( 
  "#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector localflowCpp(NumericVector x, std::size_t ni, std::size_t nw) {
  NumericVector out(ni);
  if (nw == 0) return out;            // nothing to do

  std::size_t start = 0;
  const std::size_t center_offset = (nw - 1) / 2; // window size is odd in focal

  for (std::size_t i = 0; i < ni; ++i) {
    std::size_t end = start + nw;

    double center = x[start + center_offset];
    int greater_count = 0;
    int smaller_count = 0;

    if (!NumericVector::is_na(center)) {
      for (std::size_t j = start; j < end; ++j) {
        double v = x[j];
        if (NumericVector::is_na(v)) continue;
        if (v > center)      ++greater_count;
        else if (v < center) ++smaller_count;
      }
    }
    // if center is NA, both counts stay 0 — same as R's na.rm=TRUE behavior
    out[i] = static_cast<double>(greater_count - smaller_count);
    start = end;
  }
  return out;
}"
)

# # Function to calculate the angle based on the chord (original version)
# 
# cppFunction( 
#   '#include <Rcpp.h>
# #include <cmath>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# double chordAngleCpp(NumericVector x) {
#   if (x.size() < 6) {
#     stop("Input vector must have at least 6 elements.");
#   }
# 
#   double x1 = x[0], x2 = x[1], x3 = x[2],
#          x4 = x[3], x5 = x[4], x6 = x[5];
# 
#   // Propagate NA like R would
#   if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
#       NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
#       NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
#     return NA_REAL;
#   }
# 
#   // crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
#   double t1 = x1 * x3 - x4 * x6;
#   double t2 = x2 * x3 - x5 * x6;
#   double crd = std::sqrt(t1 * t1 + t2 * t2);
# 
#   // out <- asin(crd/2) * 2  (in radians)
#   double arg = crd / 2.0;
#   // Clamp to avoid tiny numerical excursions beyond [-1, 1]
#   if (arg > 1.0) arg = 1.0;
#   if (arg < -1.0) arg = -1.0;
# 
#   return 2.0 * std::asin(arg);
# }'
# )


# # Function to calculate the angle based on the chord (modified to discount differences in slope)
# 
# cppFunction( 
#   '#include <Rcpp.h>
# #include <cmath>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# double chordAngleCpp(NumericVector x) {
#   if (x.size() < 6) {
#     stop("Input vector must have at least 6 elements.");
#   }
# 
#   double x1 = x[0], x2 = x[1], x3 = x[2],
#          x4 = x[3], x5 = x[4], x6 = x[5];
# 
#   // Propagate NA like R would
#   if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
#       NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
#       NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
#     return NA_REAL;
#   }
# 
#   // Use the minimum of x3 and x6
#   double xm = std::min(x3, x6);
# 
#   // crd <- sqrt((x[1]*xm - x[4]*xm)^2 + (x[2]*xm - x[5]*xm)^2)
#   double t1 = x1 * xm - x4 * xm;
#   double t2 = x2 * xm - x5 * xm;
#   double crd = std::sqrt(t1 * t1 + t2 * t2);
# 
#   // out <- asin(crd/2) * 2  (in radians)
#   double arg = crd / 2.0;
#   // Clamp to avoid tiny numerical excursions beyond [-1, 1]
#   if (arg > 1.0) arg = 1.0;
#   if (arg < -1.0) arg = -1.0;
# 
#   return 2.0 * std::asin(arg);
# }'
# )


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
  
  # x <- demtile
  
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

parallel::clusterEvalQ(cl, {
  library(terra)
  library(magrittr)
  library(Rcpp)
  library(circhelp)
}
)

#  Compile the C++ functions ON EACH WORKER
#    (cppFunction creates compiled code in the current R session, so each worker needs it.)
parallel::clusterEvalQ(
  cl, 
  {
    Rcpp::cppFunction(
      code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double chordAngleCpp(NumericVector x) {
  if (x.size() < 6) stop("Input vector must have at least 6 elements.");

  double x1 = x[0], x2 = x[1], x3 = x[2],
         x4 = x[3], x5 = x[4], x6 = x[5];

  if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
      NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
      NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
    return NA_REAL;
  }

  double t1 = x1 * x3 - x4 * x6;
  double t2 = x2 * x3 - x5 * x6;
  double crd = std::sqrt(t1 * t1 + t2 * t2);

  double arg = crd / 2.0;
  if (arg > 1.0) arg = 1.0;
  if (arg < -1.0) arg = -1.0;

  return 2.0 * std::asin(arg);
}
'
    )
  }
)

parallel::clusterEvalQ(cl, {
  Rcpp::cppFunction(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector localflowCpp_worker(NumericVector x, int ni, int nw) {
  NumericVector out(ni);
  if (nw <= 0) return out;

  int start = 0;
  const int center_offset = (nw - 1) / 2;

  for (int i = 0; i < ni; ++i) {
    int end = start + nw;

    double center = x[start + center_offset];
    int greater_count = 0;
    int smaller_count = 0;

    if (!NumericVector::is_na(center)) {
      for (int j = start; j < end; ++j) {
        double v = x[j];
        if (NumericVector::is_na(v)) continue;
        if (v > center)      ++greater_count;
        else if (v < center) ++smaller_count;
      }
    }

    out[i] = (double)(greater_count - smaller_count);
    start = end;
  }

  return out;
}
')
})


# Path for DEMs

dem_zips <- list.files(
  "O:/AUIT_Geodata/Denmark/Digital_elevation_models/Lidar/DHM_2023_working",
  pattern = ".zip",
  full.names = TRUE
)

# for (i in 1:length(dem_zips)) {
  for (i in 300) {
  
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
    pattern = ".tif",
    full.names = TRUE
  )
  
  # Parallel workflow
  
  # Function to process tiles
  process_tile <- function(j, rasters, mygaussmat, myfocalmat, gabor_kernels) {
    
    library(Rcpp)
    
    if (!exists("localflowCpp_worker", mode = "function")) {
      Rcpp::cppFunction(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector localflowCpp_worker(NumericVector x, int ni, int nw) {
  NumericVector out(ni);
  if (nw <= 0) return out;

  int start = 0;
  const int center_offset = (nw - 1) / 2;

  for (int i = 0; i < ni; ++i) {
    int end = start + nw;

    double center = x[start + center_offset];
    int greater_count = 0;
    int smaller_count = 0;

    if (!NumericVector::is_na(center)) {
      for (int j = start; j < end; ++j) {
        double v = x[j];
        if (NumericVector::is_na(v)) continue;
        if (v > center)      ++greater_count;
        else if (v < center) ++smaller_count;
      }
    }

    out[i] = (double)(greater_count - smaller_count);
    start = end;
  }

  return out;
}
')
    }
  
  if (!exists("chordAngleCpp", mode = "function")) {
    Rcpp::cppFunction(
      code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double chordAngleCpp(NumericVector x) {
  if (x.size() < 6) stop("Input vector must have at least 6 elements.");

  double x1 = x[0], x2 = x[1], x3 = x[2],
         x4 = x[3], x5 = x[4], x6 = x[5];

  if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
      NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
      NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
    return NA_REAL;
  }

  double t1 = x1 * x3 - x4 * x6;
  double t2 = x2 * x3 - x5 * x6;
  double crd = std::sqrt(t1 * t1 + t2 * t2);

  double arg = crd / 2.0;
  if (arg > 1.0) arg = 1.0;
  if (arg < -1.0) arg = -1.0;

  return 2.0 * std::asin(arg);
}
'
    )
  }
  
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
      fun = localflowCpp_worker,
      fillvalue = NA_real_
    )
    
    
    flowsd <- (tileflow ^ 2) %>%
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
    
    slopeaspsd <- (angletest5 ^ 2) %>%
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
    ridginess  <- smooth_to_r5(mingab)
    
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
      terra::subst(from = NaN,  to = 0) %>%
      terra::subst(from = Inf,  to = 0) %>%
      terra::subst(from = -Inf, to = 0) %>%
      smooth_to_r5()
    
    out_paths <- list(
      mins              = file.path(tmpfolder, sprintf("tile_%05d_mins.tif", j)),
      demmad           = file.path(tmpfolder, sprintf("tile_%05d_demmad.tif", j)),
      aspsd             = file.path(tmpfolder, sprintf("tile_%05d_aspsd.tif", j)),
      flowsd            = file.path(tmpfolder, sprintf("tile_%05d_flowsd.tif", j)),
      slopeaspsd        = file.path(tmpfolder, sprintf("tile_%05d_slopeaspsd.tif", j)),
      valleyness        = file.path(tmpfolder, sprintf("tile_%05d_valleyness.tif", j)),
      ridginess         = file.path(tmpfolder, sprintf("tile_%05d_ridginess.tif", j)),
      ridge_noise       = file.path(tmpfolder, sprintf("tile_%05d_ridge_noise.tif", j)),
      ridge_slope_index = file.path(tmpfolder, sprintf("tile_%05d_ridge_slope_index.tif", j)),
      ridge_valley_index= file.path(tmpfolder, sprintf("tile_%05d_ridge_valley_index.tif", j))
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
    if (length(paths) == 0) return(invisible(NULL))
    
    # Build SpatRasterCollection from filenames (lazy on disk)
    s <- terra::sprc(paths)
    
    # Merge and write to disk
    terra::merge(
      s,
      filename = out_file,
      overwrite = TRUE
    ) |> terra::setNames(layer_name)
    
    invisible(out_file)
  }
  
  # Build output filenames (one per product per zip)
  zip_nmins_path            <- file.path(dir_tiles_nmins,             paste0("nmins_",             sprintf("%03d", i), ".tif"))
  zip_demmad_path           <- file.path(dir_tiles_demmad,            paste0("demmad_",            sprintf("%03d", i), ".tif"))
  zip_aspsd_path            <- file.path(dir_tiles_aspsd,             paste0("aspsd_",             sprintf("%03d", i), ".tif"))
  zip_flowsd_path           <- file.path(dir_tiles_flowsd,            paste0("flowsd_",            sprintf("%03d", i), ".tif"))
  zip_slopeaspsd_path       <- file.path(dir_tiles_slopeaspsd,        paste0("slopeaspsd_",        sprintf("%03d", i), ".tif"))
  zip_valleyness_path       <- file.path(dir_tiles_valleyness,        paste0("valleyness_",        sprintf("%03d", i), ".tif"))
  zip_ridginess_path        <- file.path(dir_tiles_ridginess,         paste0("ridginess_",         sprintf("%03d", i), ".tif"))
  zip_ridge_noise_path      <- file.path(dir_tiles_ridge_noise,       paste0("ridge_noise_",       sprintf("%03d", i), ".tif"))
  zip_ridge_slope_idx_path  <- file.path(dir_tiles_ridge_slope_index, paste0("ridge_slope_index_", sprintf("%03d", i), ".tif"))
  zip_ridge_valley_idx_path <- file.path(dir_tiles_ridge_valley_index,paste0("ridge_valley_index_",sprintf("%03d", i), ".tif"))
  
  # Merge each product across tiles for this zip
  merge_tiles_to_zip(lapply(tile_files, `[[`, "mins"),              zip_nmins_path,            "nmins")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "demmad"),           zip_demmad_path,           "demmad")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "aspsd"),             zip_aspsd_path,            "aspsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "flowsd"),            zip_flowsd_path,           "flowsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "slopeaspsd"),        zip_slopeaspsd_path,       "slopeaspsd")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "valleyness"),        zip_valleyness_path,       "valleyness")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridginess"),         zip_ridginess_path,        "ridginess")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_noise"),       zip_ridge_noise_path,      "ridge_noise")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_slope_index"), zip_ridge_slope_idx_path,  "ridge_slope_index")
  merge_tiles_to_zip(lapply(tile_files, `[[`, "ridge_valley_index"),zip_ridge_valley_idx_path, "ridge_valley_index")
  
  # Optional: after merging, you can delete the tile-level files to keep temp small
  # (Only do this if you're sure the merge succeeded.)
  all_tile_paths <- unique(unlist(tile_files, use.names = FALSE))
  unlink(all_tile_paths, force = TRUE)
  
  gc()
  
  
  # Serial workflow
  # rasterlist_mins <- list()
  # rasterlist_demmad <- list()
  # rasterlist_aspsd <- list()
  # rasterlist_flowsd <- list()
  # rasterlist_slopeaspsd <- list()
  # 
  # for (j in 1:length(rasters)) {
  #   # for (j in 1) {
  #   
  #   j <- 1  # NB
  #   
  #   demtile0 <- rast(rasters[j])
  #   
  #   r2 <- rast(ext(demtile0), resolution = 1)
  #   
  #   r5 <- rast(ext(demtile0), resolution = 5)
  #   
  #   r10 <- rast(ext(demtile0), resolution = 10)
  #   
  #   # demtile <- resample(demtile, r2, "average")
  #   
  #   # Aggregate to 0.8 m resolution
  #   
  #   demtile <- demtile0 %>%
  #     focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
  #     aggregate(fact = 2, fun = "mean")
  #   
  #   # Local depressions
  #   
  #   focalmin <- focal(demtile, myfocalmat, fun = "min", na.rm = TRUE)
  #   
  #   lessthanmin <- demtile < focalmin
  #   
  #   mins <- lessthanmin %>%
  #     agg_5m(
  #       fun_focal1 = "mean", 
  #       fun_focal2 = "mean",
  #       fun_agg = "sum",
  #       decimals = 1
  #     )
  #   
  #   # mins <- aggregate(lessthanmin, fact = 10, fun = "sum")
  #   
  #   rasterlist_mins[[j]] <- mins
  #   
  #   # Roughness
  #   
  #   focalmeans <- focal(demtile, mygaussmat, fun = "mean", na.rm = TRUE)
  #   
  #   demdiffs <- abs(focalmeans - demtile)
  #   
  #   demmad <- demdiffs %>%
  #     agg_5m(
  #       fun_focal1 = "mean", 
  #       fun_focal2 = "mean",
  #       fun_agg = "median",
  #       decimals = 3
  #     )
  #   
  #   rasterlist_demmad[[j]] <- demmad
  #   
  #   # Differences in aspect
  #   
  #   aspsd <- differ_aspect(demtile)
  #   
  #   rasterlist_aspsd[[j]] <- aspsd
  #   
  #   # Tendency towards local flow accumulation
  #   
  #   tileflow <- focalCpp(
  #     demtile,
  #     w = 3,
  #     fun = localflowCpp,
  #     fillvalue = NA_real_
  #     )
  #   
  #   flowsd <- tileflow %>%
  #     raise_to_power(2) %>%
  #     agg_5m(
  #       fun_focal1 = "mean", 
  #       fun_focal2 = "mean",
  #       fun_agg = "mean",
  #       decimals = 4
  #     ) %>%
  #     sqrt() %>%
  #     round(2)
  #   
  #   # Note: It is also possible to do channels/ridges and the mean value.
  #   # However, this will have to wait.
  #   
  #   rasterlist_flowsd[[j]] <- flowsd
  #   
  #   # Standard deviation of slope aspect (giving less weight to flat areas)
  #   
  #   tileslope_sin <- terrain(demtile, "slope", unit="radians") %>% sin()
  # 
  #   tileasp <- terrain(demtile, "aspect", unit="radians")
  #   tileaspcos <- tileasp %>% cos()
  #   tileaspsin <- tileasp %>% sin()
  #   
  #   tileasp_focal <- terrain(focalmeans, "aspect", unit="radians")
  #   tileaspsin_focal <- tileasp_focal %>% sin()
  #   tileaspcos_focal <- tileasp_focal %>% cos()
  #   
  #   tileslope_smooth <- terrain(focalmeans, "slope", unit="radians") 
  #   tileslope_sin_smooth <- sin(tileslope_smooth)
  # 
  #   mystack3 <- c(
  #     tileaspcos, tileaspsin, tileslope_sin,
  #     tileaspcos_focal, tileaspsin_focal, tileslope_sin_smooth
  #   )
  # 
  #   angletest5 <- app(mystack3, chordAngleCpp)
  #   
  #   slopeaspsd <- angletest5 %>%
  #     raise_to_power(2) %>%
  #     agg_5m(
  #       fun_focal1 = "mean", 
  #       fun_focal2 = "mean",
  #       fun_agg = "mean",
  #       decimals = 5
  #     ) %>%
  #     sqrt() %>%
  #     round(3)
  #   
  #   rasterlist_slopeaspsd[[j]] <- slopeaspsd
  #   
  #   # Gabor filter products
  #   
  #   
  #   
  #   # Flow direction products
  #   
  #   
  #   print(paste0("i = ", i, "; j = ", j))
  # }
  
  
  # Update mosaic
  
  # Example paths for the running national mosaics (unmasked; mask at the end)
  nat_nmins_path            <- file.path(
    dir_microdem_merged, "micro_nmins_unmasked.tif")
  nat_demmad_path           <- file.path(
    dir_microdem_merged, "micro_demmad_unmasked.tif")
  nat_aspsd_path            <- file.path(
    dir_microdem_merged, "micro_aspsd_unmasked.tif")
  nat_flowsd_path           <- file.path(
    dir_microdem_merged, "micro_flowsd_unmasked.tif")
  nat_slopeaspsd_path       <- file.path(
    dir_microdem_merged, "micro_slopeaspsd_unmasked.tif")
  nat_valleyness_path       <- file.path(
    dir_microdem_merged, "micro_valleyness_unmasked.tif")
  nat_ridginess_path        <- file.path(
    dir_microdem_merged, "micro_ridginess_unmasked.tif")
  nat_ridge_noise_path      <- file.path(
    dir_microdem_merged, "micro_ridge_noise_unmasked.tif")
  nat_ridge_slope_idx_path  <- file.path(
    dir_microdem_merged, "micro_ridge_slope_index_unmasked.tif")
  nat_ridge_valley_idx_path <- file.path(
    dir_microdem_merged, "micro_ridge_valley_index_unmasked.tif")
  
  # Update each running mosaic; these calls also delete the per-zip file they ingest
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_nmins, 
                                paste0("nmins_", sprintf("%03d", i), ".tif")),
    national_path   = nat_nmins_path,
    layer_name      = "micro_nmins"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_demmad, 
                                paste0("demmad_", sprintf("%03d", i), ".tif")),
    national_path   = nat_demmad_path,
    layer_name      = "micro_demmad"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_aspsd, 
                                paste0("aspsd_", sprintf("%03d", i), ".tif")),
    national_path   = nat_aspsd_path,
    layer_name      = "micro_aspsd"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_flowsd, 
                                paste0("flowsd_", sprintf("%03d", i), ".tif")),
    national_path   = nat_flowsd_path,
    layer_name      = "micro_flowsd"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_slopeaspsd, 
                                paste0("slopeaspsd_", sprintf("%03d", i), ".tif")),
    national_path   = nat_slopeaspsd_path,
    layer_name      = "micro_slopeaspsd"
  )
  
  # Gabor products
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_valleyness, 
                                paste0("valleyness_", sprintf("%03d", i), ".tif")),
    national_path   = nat_valleyness_path,
    layer_name      = "micro_valleyness"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_ridginess, 
                                paste0("ridginess_", sprintf("%03d", i), ".tif")),
    national_path   = nat_ridginess_path,
    layer_name      = "micro_ridginess"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_ridge_noise, 
                                paste0("ridge_noise_", sprintf("%03d", i), ".tif")),
    national_path   = nat_ridge_noise_path,
    layer_name      = "micro_ridge_noise"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_ridge_slope_index, 
                                paste0("ridge_slope_index_", sprintf("%03d", i), ".tif")),
    national_path   = nat_ridge_slope_idx_path,
    layer_name      = "micro_ridge_slope_index"
  )
  
  update_running_mosaic(
    zip_mosaic_path = file.path(dir_tiles_ridge_valley_index, 
                                paste0("ridge_valley_index_", sprintf("%03d", i), ".tif")),
    national_path   = nat_ridge_valley_idx_path,
    layer_name      = "micro_ridge_valley_index"
  )
  
  gc()
}

parallel::stopCluster(cl)

# Mask and Write results to files

dem <- terra::rast(file.path(dir_dat, "Sampling_Input", "dhm2015_terraen_10m.tif"))

mask_and_write <- function(unmasked_path, final_path, layer_name) {
  r <- terra::rast(unmasked_path)
  terra::mask(
    r, dem,
    filename = final_path,
    names = layer_name,
    overwrite = TRUE
  )
}

mask_and_write(
  nat_nmins_path,            
  file.path(
    dir_microdem_merged, 
    "micro_nmins.tif"
  ),           
  "micro_nmins"
)

mask_and_write(
  nat_demmad_path,          
  file.path(
    dir_microdem_merged, 
    "micro_demmad.tif"
  ),         
  "micro_demmad"
)

mask_and_write(
  nat_aspsd_path,           
  file.path(
    dir_microdem_merged, 
    "micro_aspsd.tif"
  ),           
  "micro_aspsd"
)

mask_and_write(
  nat_flowsd_path,          
  file.path(
    dir_microdem_merged, 
    "micro_flowsd.tif"
  ),          
  "micro_flowsd"
)

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


# dem <- dir_dat %>%
#   paste0(., "/Sampling_Input/dhm2015_terraen_10m.tif") %>%
#   rast()
# 
# micro_nmins <- dir_tiles_nmins %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_nmins.tif"
#     ),
#     names = "micro_nmins",
#     overwrite = TRUE
#   )
# 
# micro_demmad <- dir_tiles_demmad %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_demmad.tif"
#     ),
#     names = "micro_demmad",
#     overwrite = TRUE
#   )
#     
# 
# micro_aspsd <- dir_tiles_aspsd %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_aspsd.tif"
#     ),
#     names = "micro_aspsd",
#     overwrite = TRUE
#   )
# 
# 
# micro_flowsd <- dir_tiles_flowsd %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_flowsd.tif"
#     ),
#     names = "micro_flowsd",
#     overwrite = TRUE
#   )
# 
# micro_slopeaspsd <- dir_tiles_slopeaspsd %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_slopeaspsd.tif"
#     ),
#     names = "micro_slopeaspsd",
#     overwrite = TRUE
#   )
# 
# micro_valleyness <- dir_tiles_valleyness %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_valleyness.tif"
#     ),
#     names = "micro_valleyness",
#     overwrite = TRUE
#   )
# 
# micro_ridginess <- dir_tiles_ridginess %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_ridginess.tif"
#     ),
#     names = "micro_ridginess",
#     overwrite = TRUE
#   )
# 
# micro_ridge_noise <- dir_tiles_ridge_noise %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_ridge_noise.tif"
#     ),
#     names = "micro_ridge_noise",
#     overwrite = TRUE
#   )
# 
# micro_ridge_slope_index <- dir_tiles_ridge_slope_index %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_ridge_slope_index.tif"
#     ),
#     names = "micro_ridge_slope_index",
#     overwrite = TRUE
#   )
# 
# micro_ridge_valley_index <- dir_tiles_ridge_valley_index %>%
#   list.files(".tif", full.names = TRUE) %>%
#   sprc() %>%
#   merge() %>%
#   mask(
#     mask = dem,
#     filename = paste0(
#       dir_microdem_merged,
#       "micro_ridge_valley_index.tif"
#     ),
#     names = "micro_ridge_valley_index",
#     overwrite = TRUE
#   )


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
plot(alltiles_nmins)   # 1 decimal
plot(alltiles_demmad)   # 3 decimals
plot(alltiles_aspsd)  # 3 decimals
plot(alltiles_flowsd)  # 2 decimals
plot(alltiles_slopeaspsd) # 3 decimals

# Plot merged results

plot(micro_nmins)
plot(micro_demmad)
plot(micro_aspsd)
plot(micro_flowsd)
plot(micro_flowsd)

# Old stuff

# # Test Gabor filters
# 
# library(OpenImageR)
# 
# init_gb <- GaborFeatureExtract$new()
# 
# gb_f <- init_gb$gabor_filter_bank(
#   scales = 1, orientations = 8, gabor_rows = 3,
#   gabor_columns = 3, plot_data = TRUE
# )
# 
# 
# outrasters <- list()
# 
# for (i in 1:8) {
#   mygabfilter <- gb_f$gabor_real [[i]] %>%
#     abs() %>%
#     subtract(mean(.)) %>%
#     divide_by(sum(abs(.)))
#   # divide_by(sd(.))
#   
#   r_gabor <- focal(demtile, w = mygabfilter)
#   
#   outrasters[[i]] <- r_gabor
# }
# 
# plot(r_gabor)
# 
# outrasters %>% rast() %>% plot()
# 
# mean_gab <- outrasters %>% rast() %>% mean() 
# 
# plot(mean_gab)
# 
# maxgab <- outrasters %>% rast() %>% max()
# mingab <- outrasters %>% rast() %>% min() %>% multiply_by(-1)
# 
# plot(maxgab)  # Valleys
# plot(mingab)  # Ridges
# 
# plot(mingab - maxgab)
# 
# sdgab <- outrasters %>% rast() %>% app(function(x) {sd(x)})
# 
# plot(sdgab)  # Rillyness
# 
# plot(mean_gab + sdgab)
# plot(mean_gab - sdgab)
# 
# maxgab2 <- ifel(mean_gab > 0, maxgab, 0)
# mingab2 <- ifel(mean_gab < 0, mingab, 0)
# 
# plot(maxgab2)
# plot(mingab2)
# 
# extremegab <- mingab2 - maxgab2
# 
# plot(extremegab)
# 
# plot(extremegab + mean_gab)
# 
# extremegab2 <- ifel(mean_gab > 0, -sdgab, sdgab)
# 
# plot(extremegab2)
# 
# plot(extremegab2 + mean_gab)
# 
# plot(mingab)
# 
# extremegab3 <- ifel(maxgab > (-1*mingab), maxgab, mingab)
# 
# plot(extremegab3)
# 
# maxgab %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Interesting  - "valleyness"
# 
# mingab %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Very similar to maxgab  - "ridginess"
# 
# extremegab %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Possibly useful  - "enhanced ridge index"
# # However, the alternative below is better, so do not use.
# 
# extremegab2 %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Possibly useful, more variable than the previous one
# # Also "enhanced ridge index", still too flat, do not use
# 
# mean_gab %>%
#   multiply_by(-1) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Possibly useful  # "ridge index", a bit too flat, so do not use
# 
# 
# mean_gab %>%
#   multiply_by(-1) %>%
#   is_greater_than(0) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Possibly useful  # "ridge probability"
# # Probably not as useful as using the individual rasters, so skip for now.
# # See example with individual rasters below.
# 
# sdgab %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Interesting - similar to mingab and maxgab  "ridge noise"
# # However, there is a bettwe alternative below, so do not use
# 
# plot(maxgab - mingab)
# 
# (maxgab - mingab) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # could also be useful - skip for now, do not use
# 
# outrasters %>%
#   rast() %>%
#   median()  %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # could also be useful - skip for now, do not use
# 
# outrasters %>%
#   rast() %>%
#   is_greater_than(0)  %>%
#   mean() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Looks interesting - could also be useful
# # "Valley probability", but skip for now, do not use
# 
# outrasters %>%
#   rast() %>%
#   abs()  %>%
#   mean() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()
# # Better alternativ to sdgab - "ridge noise" - use this one, use
# 
# mean_gab %>%
#   is_less_than(sdgab) %>% 
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()
# # Probably not useful, do not use
# 
# mean_gab %>%
#   multiply_by(-1) %>%
#   is_less_than(sdgab) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()
# # Probably not useful, do not use
# 
# mean_gab %>%
#   abs() %>%
#   is_less_than(sdgab) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()
# # Probably not useful, do not use
# 
# (sdgab*log(tileslope_smooth)) %>% 
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>% 
#   subst(from = -Inf, to = 0) %>%
#   aggregate(fact = 2, fun = "mean") %>% 
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   multiply_by(-1) %>%
#   plot()  # Ridge slope index, use
# 
# which.lyr(rast(outrasters) == extremegab3) %>%
#   subtract(1) %>%
#   multiply_by(2*pi/8) %>%
#   cos() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Ridge direction, interesting, but do not use
# 
# plot((mean_gab - mingab)/(maxgab + mingab))
# 
# ((mean_gab - mingab)/(maxgab + mingab)) %>%
#   subst(from = NaN, to = 0) %>%
#   subst(from = Inf, to = 0) %>%
#   subst(from = -Inf, to = 0) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Local Ridge/valley index, the next one is better, so do not use
# 
# ((maxgab)/(maxgab + mingab)) %>%
#   subst(from = NaN, to = 0) %>%
#   subst(from = Inf, to = 0) %>%
#   subst(from = -Inf, to = 0) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Local Ridge/valley index, use


# # Test flow direction 
# 
# tileflowdir <- terrain(demtile, "flowdir")
# 
# flowmat <- c(32, 64, 128, 16, NA, 1, 8, 4, 2) %>%
#   rev() %>%
#   matrix(., nrow = 3) %>%
#   t()
# 
# flowmat_flat <- c(32, 64, 128, 16, NA, 1, 8, 4, 2) %>%
#   rev()
# 
# focalflow <- focal(tileflowdir, 3, function(x) {sum(x == flowmat_flat, na.rm = TRUE)})
# 
# plot(focalflow)
# 
# focalflow %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Degree of convergent flow, do not use
# 
# tileslope_tan <- terrain(demtile, "slope", unit = "radians") %>% tan()
# 
# focalflow %>%
#   add(1) %>%
#   log() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Degree of convergent flow, log scale, use
# 
# focalflow %>%
#   add(1) %>%
#   divide_by(tileslope_tan) %>%
#   subst(from = Inf, to = 1) %>%
#   log() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Local TWI, use
# 
# focalflow %>%
#   add(1) %>%
#   log() %>%
#   subtract(1) %>%
#   abs() %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r5, "bilinear") %>%
#   plot()  # Overall flow convergence/divergence, use


# Test asp diff with slope
# plot(demtile)
# 
# tileslope <- terrain(demtile, "slope", unit="radians") 
# 
# tileslope_cotan <- tan(pi/2 - tileslope)
# 
# tileslope_sin <- sin(tileslope)
# 
# plot(tileslope_cotan)
# 
# tileasp <- terrain(demtile, "aspect", unit="radians")
# 
# tileaspcos <- tileasp %>% cos()
# 
# plot(tileaspcos)
# 
# tileaspsin <- tileasp %>% sin()
# 
# plot(tileaspsin)
# 
# tileasp_focal <- terrain(focalmeans, "aspect", unit="radians")
# 
# tileaspsin_focal <- tileasp_focal %>% sin()
# 
# plot(tileaspsin_focal)
# 
# tileaspcos_focal <- tileasp_focal %>% cos()
# 
# plot(tileaspcos_focal)
# 
# # tileslope_smooth <- focal(tileslope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# tileslope_smooth <- terrain(focalmeans, "slope", unit="radians") 
# 
# # tileslope_cotan_smooth <- focal(tileslope_cotan, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# tileslope_cotan_smooth <- tan(pi/2 - tileslope_smooth)
# 
# tileslope_sin_smooth <- sin(tileslope_smooth)
# 
# plot(tileaspsin - tileaspsin_focal)
# plot(tileaspcos - tileaspcos_focal)
# 
# plot(tileaspsin*tileslope_sin_smooth)
# 
# plot(tileaspsin*tileslope_sin - tileaspsin_focal*tileslope_sin_smooth)
# 
# plot(tileslope_cotan_smooth)
# 
# mystack3 <- c(
#   tileaspcos, tileaspsin, tileslope_sin,
#   tileaspcos_focal, tileaspsin_focal, tileslope_sin_smooth
# )
# 
# angletest4 <- app(mystack3, function(x) {
#   a <- c(x[1], x[2], x[3])
#   b <- c(x[4], x[5], x[6])
#   dot.prod <- a%*%b 
#   norm.a <- norm(a,type = "2")
#   norm.b <- norm(b,type = "2")
#   theta <- acos(dot.prod / (norm.a * norm.b))
#   out <- as.numeric(theta)
#   return(out)
# })
# 
# plot(angletest4)
# 
# angletest5 <- app(mystack3, function(x) {
#   crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
#   out <- asin(crd/2)*2
#   return(out)
# }
# )
# 
# plot(angletest5)
# 
# aggregate(angletest5, fact = 10, na.rm = TRUE) %>% plot()
# 
# angletest5 <- app(mystack3, function(x) {
#   z <- (x[3] + x[6])/2
#   a <- c(x[1], x[2], z)
#   b <- c(x[4], x[5], z)
#   dot.prod <- a%*%b 
#   norm.a <- norm(a,type = "2")
#   norm.b <- norm(b,type = "2")
#   theta <- acos(dot.prod / (norm.a * norm.b))
#   out <- as.numeric(theta)
#   return(out)
# })
# 
# angletest5 <- app(mystack3, chordAngleCpp)
# 
# plot(angletest5)
# 
# cppFunction( 
#   '#include <Rcpp.h>
# #include <cmath>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# double chordAngleCpp(NumericVector x) {
#   if (x.size() < 6) {
#     stop("Input vector must have at least 6 elements.");
#   }
# 
#   double x1 = x[0], x2 = x[1], x3 = x[2],
#          x4 = x[3], x5 = x[4], x6 = x[5];
# 
#   // Propagate NA like R would
#   if (NumericVector::is_na(x1) || NumericVector::is_na(x2) ||
#       NumericVector::is_na(x3) || NumericVector::is_na(x4) ||
#       NumericVector::is_na(x5) || NumericVector::is_na(x6)) {
#     return NA_REAL;
#   }
# 
#   // crd <- sqrt((x[1]*x[3]-x[4]*x[6])^2 + (x[2]*x[3]-x[5]*x[6])^2)
#   double t1 = x1 * x3 - x4 * x6;
#   double t2 = x2 * x3 - x5 * x6;
#   double crd = std::sqrt(t1 * t1 + t2 * t2);
# 
#   // out <- asin(crd/2) * 2  (in radians)
#   double arg = crd / 2.0;
#   // Clamp to avoid tiny numerical excursions beyond [-1, 1]
#   if (arg > 1.0) arg = 1.0;
#   if (arg < -1.0) arg = -1.0;
# 
#   return 2.0 * std::asin(arg);
# }'
# )
# 
# angletest5 <- app(mystack3, chordAngleCpp)
# 
# plot(angletest5)
# 
# aggregate(angletest5, fact = 10, na.rm = TRUE) %>% plot()
# 
# angletest5 %>%
#   raise_to_power(2) %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   aggregate(fact = 2, fun = "mean") %>%
#   focal(mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   resample(r10, "bilinear") %>%
#   sqrt() %>%
#   round(3) %>% plot()
# 
# angle <- function(x,y){
#   dot.prod <- x%*%y 
#   norm.x <- norm(x,type="2")
#   norm.y <- norm(y,type="2")
#   theta <- acos(dot.prod / (norm.x * norm.y))
#   as.numeric(theta)
# }


# Tests of different approaches

# plot(terrain(demtile, "slope"))
# 
# terrain(demtile, "slope") %>%
#   aggregate( fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# plot(terrain(aggregate(demtile, fact = 25, fun = "mean", na.rm = TRUE), "slope"))
# 
# terrain(demtile, "TRIrmsd") %>%
#   aggregate(fact = 25, fun = "median", na.rm = TRUE) %>%
#   plot()
# 
# terrain(demtile, "slope") %>%
#   aggregate(fact = 25, fun = "median", na.rm = TRUE) %>%
#   plot()
# 
# 
# 
# tileasp_diff_agg <- differ_aspect(demtile)
# 
# plot(tileasp_diff_agg)
# 
# detrend_asp_diff <- differ_aspect(demdiffs)
# 
# plot(detrend_asp_diff)
# 
# # Aspect
# 
# tileslope <- terrain(demtile, "slope", unit = "radians") %>%
#   sin()
# 
# tileasp <- terrain(demtile, "aspect", unit = "radians")
# 
# plot(tileasp)
# 
# # sine
# 
# tileaspsin <- tileasp %>% sin()
# 
# plot(tileaspsin)
# 
# tileaspsin_focal <- focal(tileaspsin, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspsin_focal)
# 
# tileaspsin_diff <- tileaspsin - tileaspsin_focal
# 
# plot(tileaspsin_diff)
# 
# tileaspsin_agg <- aggregate(tileaspsin_diff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileaspcos_agg)
# 
# # cosine
# 
# tileaspcos <- tileasp %>% cos()
# 
# plot(tileaspcos)
# 
# tileaspcos_focal <- focal(tileaspcos, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspcos_focal)
# 
# tileaspcos_diff <- tileaspcos - tileaspcos_focal
# 
# plot(tileaspcos_diff)
# 
# tileaspcos_agg <- aggregate(tileaspcos_diff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileaspcos_agg)
# 
# # Combine aggregated sine and cosine
# 
# tileaspagg_diff <- atan2(tileaspsin_agg, tileaspcos_agg)
# 
# plot(tileaspagg_diff)
# 
# # Combine raw sine and cosine
# 
# library(circhelp)
# 
# tileasp2 <- atan2(tileaspsin, tileaspcos)
# 
# plot(tileasp2)
# 
# tileasp_focal <- atan2(tileaspsin_focal, tileaspcos_focal)
# 
# plot(tileasp_focal)
# 
# tileasp_diff2 <- angle_diff_rad(tileasp2, tileasp_focal) %>% abs()
# 
# plot(tileasp_diff2)
# 
# aggregate(
#   tileasp_diff2^2, fact = 25, fun = "mean", na.rm = TRUE) %>%
#   sqrt() %>%
#   plot()
# 
# tileasp_diff <- atan2(tileaspsin_diff, tileaspcos_diff)
# 
# plot(tileasp)
# 
# plot(tileasp_diff)
# 
# plot(abs(tileasp_diff))
# 
# plot(abs(tileasp_diff) > pi*0.5)
# 
# aggregate(
#   abs(tileasp_diff) > pi*0.5, fact = 25, fun = "sum", na.rm = TRUE) %>%
#   plot()
# 
# aggregate(
#   abs(tileasp_diff), fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# 
# tileasp_diff_agg <- aggregate(
#   tileasp_diff^2, fact = 25, fun = "mean", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileasp_diff_agg)  # This one has promise
# 
# tileasp_diff_sd <- aggregate(
#   tileasp_diff^2, fact = 25, fun = "sd", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(tileasp_diff_sd)
# 
# # Slope
# 
# plot(tileslope)
# 
# focalslope <- focal(tileslope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# slopediff <- tileslope - focalslope
# 
# plot(slopediff)
# 
# slopediff_agg <- aggregate(
#   slopediff^2, fact = 25, fun = "median", na.rm = TRUE) %>%
#   sqrt()
# 
# plot(slopediff_agg)
# 
# # Different focal filters
# 
# # Laplacian filter
# 
# dem_lap <- focal(demtile, matrix(c(0,1,0,1,-4,1,0,1,0), nrow=3), fun = "sum")
# 
# plot(dem_lap)
# 
# dem_lap_agg <- aggregate(
#   dem_lap, fact = 25, fun = "mean", na.rm = TRUE)
# 
# plot(dem_lap_agg)
# 
# # Sobel filter
# 
# fx <- matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3)
# 
# fy <- matrix(c(1,0,-1,2,0,-2,1,0,-1), nrow=3)
# 
# dem_sob <- focal(demtile, fx, fun = "sum")
# 
# plot(dem_sob)
# 
# # Asp also using slope
# 
# plot(tileasp_diff*tileslope)
# 
# tileaspslope_diff_agg <- aggregate(tileasp_diff*tileslope, fact = 25, fun = "median", na.rm = TRUE)
# 
# plot(tileaspslope_diff_agg)
# 
# tileaspsin_slope <- tileasp %>% sin() %>% multiply_by(tileslope)
# 
# tileaspsin_slope_focal <- focal(tileaspsin_slope, mygaussmat, fun = "mean", na.rm = TRUE)
# 
# plot(tileaspsin_slope - tileaspsin_slope_focal)
# 
# aggregate((tileaspsin_slope - tileaspsin_slope_focal)^2, fact = 25, fun = "median", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# terrain(demtile, "aspect", unit = "radians") %>%
#   cos() %>%
#   aggregate(fact = 25, fun = "mean", na.rm = TRUE) %>%
#   plot()
# 
# demtile %>%
#   aggregate(fact = 25, fun = "mean", na.rm = TRUE) %>%
#   terrain("aspect", unit = "radians") %>%
#   cos() %>%
#   plot()
# 
# # using smooth dem as baseline
# 
# smoothdem_aspsin <- focal(demtile, mygaussmat, fun = "mean", na.rm = TRUE) %>%
#   terrain("aspect", unit = "radians") %>%
#   sin()
#   
# plot(smoothdem_aspsin)
# 
# plot(tileaspsin - smoothdem_aspsin)
# 
# aggregate((tileaspsin - smoothdem_aspsin)^2, fact = 25, fun = "median", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# # Local flow accumulation
# 
# localflow <- function(x) {
#   xmat <- matrix(x, nrow = sqrt(length(x))) %>% rast()
#   flow <- terrain(xmat, "flowdir") %>% flowAccumulation()
#   out <- flow %>% values() %>% unlist() %>% magrittr::extract(5)
#   return(out)
# }
# 
# f <- system.file("ex/elev.tif", package="terra")
# r <- rast(f)
# 
# plot(r)
# 
# r_focalized <- focal(r, w = 3, localflow)
# 
# plot(r_focalized)
# 
# mynumbers <- c(1:9)
# 
# localflow(mynumbers)
# 
# localflow_R <- function(x) {
#   x_center <- x[ceiling(length(x)/2)]
#   out <- sum((x > x_center), na.rm = TRUE) - sum((x < x_center), na.rm = TRUE)
#   return(out)
# }
# 
# r_focalized <- focal(demtile, w = 3, localflow)
# 
# plot(r_focalized)
# 
# aggregate(r_focalized^2, fact = 25, fun = "mean", na.rm = TRUE) %>% sqrt() %>% plot()
# 
# aggregate(abs(r_focalized), fact = 25, fun = "mean", na.rm = TRUE) %>% plot()
# 
# aggregate(r_focalized, fact = 25, fun = "sd", na.rm = TRUE) %>% plot()
# 
# raw <- focalCpp(demtile, w=3, fun=localflowCpp, fillvalue=NA_real_)

# END