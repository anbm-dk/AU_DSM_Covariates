# Focal density

focal_density <- function(
    inrast,
    nsteps,
    include_list = FALSE
) {
  r1 <- rast(ncols = 180, nrows = 180, xmin = 0)
  myfilter1 <- round(
    focalMat(r1, c(1, 2), "Gauss"),
    3
  )
  myfilter2 <- myfilter1
  
  smooth_up_list <- list()
  aggregated_list <- list()
  # aggregated_list[[1]] <- c(
  #   inrast*0 + 1,    # Number of cells
  #   inrast           # Cell values
  # )
  aggregated_list[[1]] <- inrast*0 + 1
  # names(aggregated_list[[1]]) <- c("count", "mean")
  # Stepwise smoothing and aggregation
  for (i in 2:nsteps) {
    smoothed_down <- terra::focal(
      aggregated_list[[i - 1]],
      w = myfilter1,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE,
      expand = TRUE,
      fillvalue = 0
    )
    aggregated_list[[i]] <- terra::aggregate(
      smoothed_down,  
      fun = "mean",
      na.rm = TRUE
    )
  }
  # Stepwise disaggregation, merging and smoothing
  smooth_up_list[[nsteps]] <- aggregated_list[[nsteps]]
  for (i in (nsteps - 1):1) {
    # Disaggregate by 2
    splitted <- terra::resample(
      x = smooth_up_list[[i + 1]],
      y = aggregated_list[[i]],
      method = "near"
    )
    # Merge with aggregated layers
    # merged <- terra::merge(
    #   x = aggregated_list[[i]],
    #   y = splitted
    # )
    merged <- max(c(
      aggregated_list[[i]],
      splitted
    ), na.rm = TRUE)
    # Smoothing
    smooth_up_list[[i]] <- terra::focal(
      merged,
      w = myfilter2,
      fun = "sum",
      na.policy = "all",
      na.rm = TRUE
    )
  }
  # Divide mean values by the number of cells, to get a weighted mean
  # final_lyr <- smooth_up_list[[1]][[2]] / smooth_up_list[[1]][[1]]
  final_lyr <- smooth_up_list[[1]]
  out <- list()
  # Merge with input layer
  dtyp_inrast <- datatype(inrast)
  if(dtyp_inrast != "") {
    my_wopt <- list(datatype = datatype(inrast))
  } else {
    my_wopt <- list()
  }
  out$final <- terra::merge(
    aggregated_list[[1]],
    final_lyr,
    wopt = my_wopt
  )
  if (include_list) {
    out$aggregated_list <- aggregated_list
    out$smooth_up_list <- smooth_up_list
  }
  return(out)
}

f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
plot(r)


r_dens <- focal_density(r, 6)
# r_dens <- fill_gaps_gauss(r * 0 + 1, 6)


plot(r_dens$final)

# END