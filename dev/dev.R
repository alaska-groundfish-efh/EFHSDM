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
usethis::use_package("terra", type = "Imports", min_version = NULL)

#GamModel dependencies
usethis::use_package("mgcv", type = "Imports", min_version = NULL) #used by multiple fns
usethis::use_package("raster", type = "Imports", min_version = NULL) #used by multiple fns
usethis::use_package("PresenceAbsence", type = "Imports", min_version = NULL)

#LoadMap dependencies
#usethis::use_package("rgdal", type = "Imports", min_version = NULL)
#usethis::use_package("sp", type = "Imports", min_version = NULL)
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
usethis::use_dev_package("akgfmaps", type = "Imports", remote = "git::https://github.com/afsc-gap-products/akgfmaps.git@master")
#devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)


# Ignore meatgrinder file -------------------------------------------------
usethis::use_build_ignore(files = "R/Meatgrinder5.R")


# Datasets ----------------------------------------------------------------
# Advice from Sean:
# There are a couple of ways to make .rda data set load into the environment
# One is to put everything in data and set LazyData: true in the description
# The other is ensure that all of the data sets that should be exported have @export in the .R file for your data
# In coldpool I have LazyData: false and export the individual data sets: https://github.com/afsc-gap-products/coldpool/blob/2019492cb8db3bc71e086075406217ebba073419/R/data.R#L1

# In FishStatsUtils, Jim has LazyData: true and puts all of the lazy load .rdas in in /data/
#   https://github.com/James-Thorson-NOAA/FishStatsUtils/tree/main/data
# and
# https://github.com/James-Thorson-NOAA/FishStatsUtils/blob/51c7713682ecc8a30ddb3d1f949b7c241451e75b/DESCRIPTION#L57



# Set up sysdata.rda, which includes data that should be automatically loaded with the package and will not be avail to users outside the package

# These are all in the current example; more will probably need to be added
region_data_all <- read.csv("Y:/RACE_EFH_Variables/Trawl_Models/GOA/all_GOA_data_2021.csv")
GOA_bathy <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Bathy")
GOA_btemp <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Btemp")
GOA_slope <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Slope")
GOA_sponge <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Spongefactor")

GOA_lat <- terra::init(GOA_bathy, v = "y")
GOA_lat <- terra::mask(lat, GOA_bathy, overwrite = F)
GOA_lon <- terra::init(GOA_bathy, v = "x")
GOA_lon <- terra::mask(lon, GOA_bathy, overwrite = F)

# region.data <- subset(region_data_all, year >= 2012)
# region.data$sponge <- as.integer(region.data$sponge > 0)
# region.data$logarea <- log(region.data$area)

# this step formulates lon and lat so they can easily be used as variables, and stacks them in one object

#Multiply raster by one to get "free-floating" raster w no filenames attached (good because we don't want paths in there)
raster_stack <- 1*terra::rast(list(GOA_lon, GOA_lat, GOA_bathy, GOA_btemp, GOA_slope, GOA_sponge))

names(raster_stack) <- c("lon", "lat", "bdepth", "btemp", "slope", "sponge")

usethis::use_data(region_data_all)
usethis::use_data(GOA_bathy)
usethis::use_data(GOA_btemp)
usethis::use_data(GOA_slope)
usethis::use_data(GOA_sponge)
usethis::use_data(GOA_lat)
usethis::use_data(GOA_lon)
usethis::use_data(raster_stack)
#save(region_data_all, GOA_bathy, GOA_btemp, GOA_slope, GOA_sponge, GOA_lat, GOA_lon, file = here::here("R","sysdata.rda"))


# Build package -----------------------------------------------------------

devtools::build()
# Build pkg including vignettes
