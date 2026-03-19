# Function for cropping covariates

cropstack <- function(
  x, # list of files
  y, # extent, polygon or raster
  folder, # target folder
  mask = FALSE # Use y for masking?
) {
  for (i in 1:length(x)) {
    require(terra)
    require(magrittr)
    require(dplyr)
    require(tools)
    r <- x[i] %>%
      terra::rast(.)
    dtype <- terra::datatype(r)
    outname_base <- x[i] %>%
      base::basename(.) %>%
      tools::file_path_sans_ext(.) %>%
      base::make.names(.)
    outfile <- outname_base %>%
      base::paste0(folder, "/", ., ".tif")
    terra::crop(
      x = r,
      y = y,
      mask = mask,
      names = outname_base,
      datatype = dtype,
      filename = outfile,
      overwrite = TRUE
    )
  }
}

# END
