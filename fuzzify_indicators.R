# Function to fuzzify indicators


fuzzify_indicators <- function(
  x = NULL,
  aggregation_factor = 1,
  residual_layer = TRUE,
  local_filter = NULL,
  n_digits = 3,
  outfolder = NULL,
  mask = NULL
) {
  if (!is.null(mask)) {
    x <- terra::cover(
      x = x,
      y = mask * 0
    )
  }

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

  if (!is.null(mask)) {
    x_fuzzy <- terra::mask(
      x = x_fuzzy,
      mask = mask
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

# END
