# Function to fill gaps in raster

fill_gaps_gauss <- function(
    inrast,
    nsteps,
    include_list = FALSE,
    weighted
) {
  r1 <- rast(ncols = 180, nrows = 180, xmin = 0)
  # myfilter <- round(
  #   focalMat(r1, c(1, 2), "Gauss"),
  #   3
  # )
  myfilter <- terra::focalMat(r1, c(9, 22), "Gauss") %>%
    terra::rast() %>%
    terra::aggregate(fact = 9, fun = "sum") %>%
    terra::values() %>%
    as.vector() %>%
    matrix(ncol = 5) %>%
    round(3)
  
  smooth_up_list <- list()
  aggregated_list <- list()
  aggregated_list[[1]] <- c(
    inrast*0 + 1,    # Number of cells
    inrast           # Cell values
  )
  names(aggregated_list[[1]]) <- c("count", "mean")
  # Stepwise smoothing and aggregation
  for (i in 2:nsteps) {
    if (weighted) {
      product_i <- aggregated_list[[i - 1]][[1]] * aggregated_list[[i - 1]][[2]]
      aggregated_i <- c(
        aggregated_list[[i - 1]][[1]],
        product_i
      ) %>%
        terra::focal(
        w = myfilter,
        fun = "mean",
        na.policy = "omit",
        na.rm = TRUE
      ) %>%
        terra::aggregate(
          fun = "sum",
          na.rm = TRUE
        )
      aggregated_list[[i]] <- c(
        aggregated_i[[1]],
        aggregated_i[[2]] / aggregated_i[[1]]
        # aggregated_i[[2]]
      )
    } else {
      smoothed_down <- terra::focal(
        aggregated_list[[i - 1]],
        w = myfilter,
        fun = "mean",
        na.policy = "all",
        na.rm = TRUE
      )
      aggregated_list[[i]] <- terra::aggregate(
        smoothed_down,  
        fun = "mean",
        na.rm = TRUE
      )
    }
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
    if (weighted) {
      splitted[[1]] <- splitted[[1]] / 4
      # Merge with aggregated layers
      merged <- terra::merge(
        x = aggregated_list[[i]],
        y = splitted
      )
      product_i <- merged[[1]] * merged[[2]]
      # Smoothing
      smoothed_i <- c(
        merged[[1]],
        product_i
      ) %>%
        terra::focal(
          w = myfilter,
          fun = "mean",
          na.policy = "omit",
          na.rm = TRUE
        )
      
      smooth_up_list[[i]] <- c(
        smoothed_i[[1]],
        smoothed_i[[2]] / smoothed_i[[1]]
      )
      
    } else {
      # Merge with aggregated layers
      merged <- terra::merge(
        x = aggregated_list[[i]],
        y = splitted
      )
      # Smoothing
      smooth_up_list[[i]] <- terra::focal(
        merged,
        w = myfilter,
        fun = "mean",
        na.policy = "omit",
        na.rm = TRUE
      )
    }
    
  }
  # Divide mean values by the number of cells, to get a weighted mean
  if (!weighted) {
    final_lyr <- smooth_up_list[[1]][[2]] / smooth_up_list[[1]][[1]]
  } else {
    final_lyr <- smooth_up_list[[1]][[2]]
    # final_lyr <- smooth_up_list[[1]][[2]] / smooth_up_list[[1]][[1]]
  }
  dtyp_inrast <- datatype(inrast)
  if(dtyp_inrast != "") {
    my_wopt <- list(datatype = datatype(inrast))
  } else {
    my_wopt <- list()
  }
  if (include_list) {
    # Merge with input layer
    out <- list()
    out$final <- terra::merge(
      inrast,
      final_lyr,
      wopt = my_wopt
    )
    out$aggregated_list <- aggregated_list
    out$smooth_up_list <- smooth_up_list
  } else {
    out <- terra::merge(
      inrast,
      final_lyr,
      wopt = my_wopt
    )
  }
  return(out)
}

f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
plot(r)

r[1200] <- 0

plot(r)

r_fill <- fill_gaps_gauss(r, 7, weighted = FALSE)
r_fill_w <- fill_gaps_gauss(r, 7, weighted = TRUE)
r_fill_w_list <- fill_gaps_gauss(r, 7, weighted = TRUE, include_list = TRUE)

plot(r_fill)
plot(r_fill_w)

plot(r_fill_w_list$aggregated_list[[3]])
plot(r_fill_w_list$smooth_up_list[[3]])

plot(r_fill_w_list$smooth_up_list[[1]])

# END