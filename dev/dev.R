# Scripts use to run pkg development


# Documentation -----------------------------------------------------------
# Package-level documentation
use_package_doc(open = rlang::is_interactive())

# Function documentation
devtools::document()

# Vignettes
usethis::use_vignette("EFHSDM")
devtools::build_vignettes()


# Dependencies ------------------------------------------------------------
# ***Need to add everything that gets called w library() or require() here.
# akgfmaps dependencies
usethis::use_package("ggplot2", type = "Imports", min_version = NULL)
usethis::use_package("viridis", type = "Imports", min_version = NULL) #used by multiple fns
usethis::use_package("stars", type = "Imports", min_version = NULL)
usethis::use_package("sf", type = "Imports", min_version = NULL)
usethis::use_package("gridExtra", type = "Imports", min_version = NULL)
usethis::use_package("patchwork", type = "Imports", min_version = NULL)
usethis::use_package("MASS", type = "Imports", min_version = NULL)
usethis::use_package("scales", type = "Imports", min_version = NULL)
usethis::use_package("labeling", type = "Imports", min_version = NULL)

#GamModel dependencies
usethis::use_package("mgcv", type = "Imports", min_version = NULL) #used by multiple fns
usethis::use_package("raster", type = "Imports", min_version = NULL) #used by multiple fns
usethis::use_package("PresenceAbsence", type = "Imports", min_version = NULL)

#LoadMap dependencies
usethis::use_package("rgdal", type = "Imports", min_version = NULL)
usethis::use_package("sp", type = "Imports", min_version = NULL)
usethis::use_package("gstat", type = "Imports", min_version = NULL)
#usethis::use_package("viridis", type = "Imports", min_version = NULL)
#usethis::use_package("mgcv", type = "Imports", min_version = NULL)
#usethis::use_package("raster", type = "Imports", min_version = NULL)


# Packages in development
usethis::use_dev_package("akgfmaps", type = "Imports", remote = "git::https://github.com/sean-rohan-NOAA/akgfmaps.git@master")
#devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)

# Datasets ----------------------------------------------------------------

use_data() #FILL IN with example datasets that need to be stored with the pkg
