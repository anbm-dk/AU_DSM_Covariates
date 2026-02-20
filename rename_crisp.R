# Function to rename crisp indicator layers

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


# END