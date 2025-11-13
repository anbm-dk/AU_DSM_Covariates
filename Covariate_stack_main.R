### Collection and standardization of 10 m covariates for digital soil mapping
### Main script for overview and planning.

# To do:
#  1: Fix NA values in fuzzy landscape elements and geology. [!]
#  2: Add new wetland extent as a covariate (fuzzy?).
#  3: ADK wetlands as covariate (fuzzy), based on raw version without landscape
#     elements.
#  4: Add layers from SoilSuite. [!]
#  5: Add ALOS/PALSAR. [!]
#  6: Alpha Earth layers? (Ask Julian/Sebastian)
#  7: Peat types from Giri's work. (fuzzy) (Amélie?)
#  8: Peat probabilities and point probabilities based on Jupiter and the Ochre
#     database. [!]
#  9: Microtopography. [!]
# 10: Addtional DEM layers with sea surfaces filled in?
# 11: 100 m DEM derivatives?
# 12: Standardize names.
# 13: Round off values.
# 14: Fill gaps (especially bare soil products).
# 15: Redo tiles.

# Startup

library(terra)
library(tidyverse)
library(magrittr)
library(tidyr)
library(dplyr)


### END ###