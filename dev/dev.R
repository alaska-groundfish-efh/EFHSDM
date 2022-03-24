# Scripts use to run pkg development


# Documentation -----------------------------------------------------------
# Package-level documentation
use_package_doc(open = rlang::is_interactive())
usethis::use_vignette("EFHSDM") # package-level vignette

# Function documentation
usethis::use_vignette("BasicsOfFiveYearReview")

# Build documentation
devtools::document()
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
usethis::use_package("magrittr", type = "Imports", min_version = NULL)

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

#Maxent dependencies
#usethis::use_package("raster", type = "Imports", min_version = NULL)
usethis::use_package("PresenceAbsence", type = "Imports", min_version = NULL)
usethis::use_package("maxnet", type = "Imports", min_version = NULL)
usethis::use_package("ENMeval", type = "Imports", min_version = NULL)

# Xtable dependencies
usethis::use_package("xtable", type = "Imports", min_version = NULL)
usethis::use_package("XML", type = "Imports", min_version = NULL)


# Packages in development
usethis::use_dev_package("akgfmaps", type = "Imports", remote = "git::https://github.com/sean-rohan-NOAA/akgfmaps.git@master")
#devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)


# Ignore meatgrinder file -------------------------------------------------
usethis::use_build_ignore(files = "R/Meatgrinder5.R")


# Datasets ----------------------------------------------------------------
# Set up sysdata.rda, which includes data that should be automatically loaded with the package and will not be avail to users outside the package

# These are all in the current example; more will probably need to be added
region_data_all <- read.csv("Y:/RACE_EFH_Variables/Trawl_Models/GOA/all_GOA_data_2021.csv")
GOA_bathy <- raster::raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Bathy")
GOA_btemp <- raster::raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Btemp")
GOA_slope <- raster::raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Slope")
GOA_sponge <- raster::raster("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Spongefactor")

GOA_lat <- raster::init(GOA_bathy, v = "y")
GOA_lat <- raster::mask(lat, GOA_bathy, overwrite = F)
GOA_lon <- raster::init(GOA_bathy, v = "x")
GOA_lon <- raster::mask(lon, GOA_bathy, overwrite = F)

usethis::use_data(region_data_all)
usethis::use_data(GOA_bathy)
usethis::use_data(GOA_btemp)
usethis::use_data(GOA_slope)
usethis::use_data(GOA_sponge)
usethis::use_data(GOA_lat)
usethis::use_data(GOA_lon)
#save(region_data_all, GOA_bathy, GOA_btemp, GOA_slope, GOA_sponge, GOA_lat, GOA_lon, file = here::here("R","sysdata.rda"))


# Build package -----------------------------------------------------------

devtools::build()
# Build pkg including vignettes
