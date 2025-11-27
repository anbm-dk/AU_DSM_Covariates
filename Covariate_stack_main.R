### Collection and standardization of 10 m covariates for digital soil mapping
### Main script for overview and planning.

# To do:

# First steps:
#  1: Fix NA values in fuzzy landscape elements and geology. [ok]
#  2: Add layers from SoilSuite. [!]
#  3: Add ALOS/PALSAR. [!]
#  4: ADK wetlands as covariate (fuzzy), based on raw version without landscape
#     elements [ok].

# Second steps:
#  5: Peat probabilities and point probabilities based on Jupiter and the Ochre
#     database. [!]
#  6: Microtopography. [!]
#  7: 100 m DEM derivatives?
#  8: Addtional DEM layers with sea surfaces filled in?

# Wait:
#  9: Add new wetland extent as a covariate (fuzzy?). [wait]
# 10: Alpha Earth layers? (Ask Julian/Sebastian) [wait]
# 11: Peat types from Giri's work. (fuzzy) (Amélie?) [wait]

# Standard processing:
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