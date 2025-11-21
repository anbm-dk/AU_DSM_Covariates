# Fix extent for 30 m covariates

library(terra)
library(tidyverse)
library(magrittr)
library(dplyr)

cov_files_30m <- list.files(
  "O:/Tech_AGRO/Jord/DSM/Covariates/national_30m/Covariates_v1",
  pattern = "\\.tif$",
  full.names = TRUE
)

cov_names_30m <- cov_list_30m %>%
  lapply(names) %>%
  unlist()

cov_list_30m <- cov_files_30m %>%
  lapply(
    function(x) {
      out <- rast(x)
      crs(out) <- "EPSG:25832"
      return(out)
      }
    )

cov_list_30m

cov_ext_30m <- cov_list_30m %>%
  lapply(
    function(x) {
      out <- ext(x) %>%
        as.vector()
      return(out)
    }
  ) %>%
  bind_rows()

cov_ext_median <- cov_ext_30m %>%
  apply(2, median)

iscorrext_ext <- cov_ext_30m %>%
  as.matrix() %>%
  t() %>%
  equals(cov_ext_median) %>%
  t() %>%
  apply(
    1,
    sum
  ) %>%
  equals(4)
  
hascorrext_ext1 <- cov_list_30m[[which(iscorrext_ext)[1]]] %>% rast()

cov_list_30m %>%
  .[which(!iscorrext_ext)]

cov_names_30m %>%
  .[which(!iscorrext_ext)]
# "bluespot" "rootzone" "wetland"

decimals_dtyp_corrected <- data.frame(
  name = cov_names_30m %>%
    .[which(!iscorrext_ext)],
  decimals = c(2, 1, 0),
  dtyp = c("FLT4S", "FLT4S", "INT2U")
  
)

for (i in which(!iscorrext_ext)) {
  rast_name <- names(cov_list_30m[[i]])
  
  decimals_i <- decimals_dtyp_corrected %>%
    filter(name == rast_name) %>%
    select(decimals) %>%
    unlist() %>%
    unname()
  
  rast_resamp <- resample(
    x = cov_list_30m[[i]],
    y = hascorrext_ext1
    ) %>%
    round(digits = decimals_i)
  
  names(rast_resamp) <- rast_name
  varnames(rast_resamp) <- rast_name
  
  cov_list_30m[[i]] <- rast_resamp
}

cov_list_30m

cov_stack_30m <- cov_list_30m %>% rast()

dtypes <- datatype(cov_stack_30m)

dir_newstack <- getwd() %>%
  dirname() %>%
  paste0(., "/covariates_30m_20251121/") %T>%
  dir.create()

dir_newstack

for(i in 1:length(cov_list_30m)) {
  
  rast_name <- names(cov_list_30m[[i]])
  
  dtyp_i <- dtypes[i]
  
  if (dtyp_i == "") {
    dtyp_i <- decimals_dtyp_corrected %>%
      filter(name == rast_name) %>%
      select(dtyp) %>%
      unlist() %>%
      unname()
  }
  
  writeRaster(
    x = cov_list_30m[[i]],
    filename = dir_newstack %>%
      paste0(., rast_name, ".tif"),
    overwrite = TRUE,
    names = rast_name,
    datatype = dtyp_i,
    gdal = "TILED=YES"
  )
}



# END