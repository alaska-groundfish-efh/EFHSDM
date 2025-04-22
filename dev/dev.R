# Scripts use to run pkg development
library(usethis)
library(devtools)


# Set up package skeleton -------------------------------------------------

usethis::create_package()
usethis::use_git()

# Documentation -----------------------------------------------------------

usethis::use_package_doc(open = rlang::is_interactive()) # Package-level documentation
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
usethis::use_package("viridis", type = "Imports", min_version = NULL) # used by multiple fns
usethis::use_package("stars", type = "Imports", min_version = NULL)
usethis::use_package("sf", type = "Imports", min_version = NULL)
usethis::use_package("gridExtra", type = "Imports", min_version = NULL)
usethis::use_package("scales", type = "Imports", min_version = NULL)
usethis::use_package("magrittr", type = "Imports", min_version = NULL)
usethis::use_package("terra", type = "Imports", min_version = NULL)

# GamModel dependencies
usethis::use_package("mgcv", type = "Imports", min_version = NULL) # used by multiple fns
usethis::use_package("PresenceAbsence", type = "Imports", min_version = NULL)

# LoadMap dependencies
usethis::use_package("viridis", type = "Imports", min_version = NULL)


# Maxent dependencies
usethis::use_package("PresenceAbsence", type = "Imports", min_version = NULL)
usethis::use_package("maxnet", type = "Imports", min_version = NULL)
# usethis::use_package("ENMeval", type = "Imports", min_version = NULL)

# Xtable dependencies
usethis::use_package("xtable", type = "Imports", min_version = NULL)
# usethis::use_package("XML", type = "Imports", min_version = NULL)


# Packages in development
usethis::use_dev_package("akgfmaps", type = "Imports", remote = "git::https://github.com/afsc-gap-products/akgfmaps.git@master")
# devtools::install_github("sean-rohan-noaa/akgfmaps", build_vignettes = TRUE)

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
region_data_goa <- read.csv("Y:/RACE_EFH_Variables/Trawl_Models/GOA/all_GOA_data_2021.csv")
GOA_bathy <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Bathy.grd")
GOA_btemp <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Btemp.grd")
GOA_slope <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Slope.grd")
GOA_sponge <- terra::rast("Y:/RACE_EFH_variables/Variables/Variables_GOA_1km/Spongefactor.grd")

GOA_lat0 <- terra::init(GOA_bathy, fun = "y")
GOA_lat <- terra::mask(GOA_lat0, GOA_bathy, overwrite = F)
GOA_lon0 <- terra::init(GOA_bathy, fun = "x")
GOA_lon <- terra::mask(GOA_lon0, GOA_bathy, overwrite = F)

# region.data <- subset(region_data_goa, year >= 2012)
# region.data$sponge <- as.integer(region.data$sponge > 0)
# region.data$logarea <- log(region.data$area)

# this step formulates lon and lat so they can easily be used as variables, and stacks them in one object

# Multiply raster by one to get "free-floating" raster w no filenames attached (good because we don't want paths in there)
raster_stack <- 1 * terra::rast(list(GOA_lon, GOA_lat, GOA_bathy, GOA_btemp, GOA_slope, GOA_sponge))

names(raster_stack) <- c("lon", "lat", "bdepth", "btemp", "slope", "sponge")

# Multiply them all by 1
GOA_bathy <- terra::wrap(GOA_bathy)
GOA_lon <- terra::wrap(GOA_lon)
GOA_lat <- terra::wrap(GOA_lat)
GOA_btemp <- terra::wrap(GOA_btemp)
GOA_slope <- terra::wrap(GOA_slope)
GOA_sponge <- terra::wrap(GOA_sponge)

raster_stack <- terra::wrap(raster_stack) # Makes PackedSpatRaster object, need to unwrap to use

# Save new raster objects
# save(GOA_bathy, file = here::here("data","GOA_bathy.rda"))
# save(GOA_lon, file = here::here("data","GOA_lon.rda"))
# save(GOA_lat, file = here::here("data","GOA_lat.rda"))
# save(GOA_btemp, file = here::here("data","GOA_btemp.rda"))
# save(GOA_slope, file = here::here("data","GOA_slope.rda"))
# save(GOA_sponge, file = here::here("data","GOA_sponge.rda"))
# save(raster_stack, file = here::here("data","raster_stack.rda"))


usethis::use_data(region_data_goa, overwrite = TRUE)
usethis::use_data(region_data_ebs, overwrite = TRUE)
usethis::use_data(region_data_ai, overwrite = TRUE)

# Only do this if you want to save each raster layer separately
usethis::use_data(GOA_bathy, overwrite = TRUE)
usethis::use_data(GOA_btemp, overwrite = TRUE)
usethis::use_data(GOA_slope, overwrite = TRUE)
usethis::use_data(GOA_sponge, overwrite = TRUE)
usethis::use_data(GOA_lat, overwrite = TRUE)
usethis::use_data(GOA_lon, overwrite = TRUE)

# Do this if you want to save the full raster stack as one stacked object
usethis::use_data(raster_stack, overwrite = TRUE)

save(region_data_goa, GOA_bathy, GOA_btemp, GOA_slope, GOA_sponge, GOA_lat, GOA_lon, file = here::here("R", "sysdata.rda"))


# Document datasets -------------------------------------------------------
devtools::document(roclets = c("rd", "collate", "namespace"))



# Setup testing -----------------------------------------------------------
usethis::use_testthat()

# Example script for drafting a test
# usethis::use_test("name")
usethis::use_test("name")

# Build package -----------------------------------------------------------

devtools::build()
# Build pkg including vignettes
