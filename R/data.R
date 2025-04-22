#' Bathymetry data from the GOA
#'
#' A raster layer containing bathymetry covariates for SDMs.
#'
#' @format A raster layer with 12 slots
#' \describe{
#'   \item{name}{character string with path to raster layer}
#'   \item{datanotation}{Raster file native format}
#' }
#' @source Compiled by Jeremy Harris using 1-km resolution data from Mark Zimmerman
"GOA_bathy"

#' GOA region data for abundance for EFH species and lifestages included in the 2023 Review
#'
#' The variables are as follows:
#' @format a data frame with 8,507 rows and 184 variables (each species/lifestage gets a column):
#' \describe{
#'   \item{hauljoin}{A numeric value for joining catch data between catch and haul tables.}
#'   \item{region}{EFH area (GOA)}
#'   \item{subregion}{Subregion.}
#'   \item{vessel}{Vessel ID from the survey; matches RACEBASE, FOSS, others}
#'   \item{cruise}{Cruise ID number}
#'   \item{haul}{Haul number}
#' }
#' @usage data(region_data_goa)
"region_data_goa"

#' EBS region data for abundance for EFH species and lifestages included in the 2023 Review
#'
#' The variables are as follows:
#' @format a data frame with 15,592 rows and 184 variables (each species/lifestage gets a column):
#' \describe{
#'   \item{hauljoin}{A numeric value for joining catch data between catch and haul tables.}
#'   \item{region}{EFH area (EBS)}
#'   \item{subregion}{Subregion.}
#'   \item{vessel}{Vessel ID from the survey; matches RACEBASE, FOSS, others}
#'   \item{cruise}{Cruise ID number}
#'   \item{haul}{Haul number}
#' }
#' @usage data(region_data_ebs)
"region_data_ebs"

#' AI region data for abundance for EFH species and lifestages included in the 2023 Review
#'
#' The variables are as follows:
#' @format a data frame with 5,347 rows and 184 variables (each species/lifestage gets a column):
#' \describe{
#'   \item{hauljoin}{A numeric value for joining catch data between catch and haul tables.}
#'   \item{region}{EFH area (AI)}
#'   \item{subregion}{Subregion.}
#'   \item{vessel}{Vessel ID from the survey; matches RACEBASE, FOSS, others}
#'   \item{cruise}{Cruise ID number}
#'   \item{haul}{Haul number}
#' }
#' @usage data(region_data_ai)
"region_data_ai"

#' GOA covariate data
#'
#' A terra SpatRaster containing environmental covariates for the GOA,
#' as used in the 2023 EFH review.
#'
#' @format A terra PackedSpatRaster with 6 slots
#' \describe{
#'   \item{lon}{Longitude}
#'   \item{lat}{Latitude}
#'   \item{bdepth}{Bottom depth in meters}
#'   \item{btemp}{Bottom temperature in deg C}
#'   \item{slope}{The rate of change in bathymetry over a defined area, measured as the first derivative of the bathymetric surface, in degrees of incline (Horn 1981, Dolan and Lucieer 2014.)}
#'   \item{sponge}{CPUE of sponge from bottom trawl surveys. Usually, this is transformed into sponge presence-absence when fitting EFH SDMs.}
#' }
#' @source Compiled by Jeremy Harris using 1-km resolution data from Mark Zimmerman
#' @docType data
#' @keywords datasets
#' @name raster_stack
#' @usage data(raster_stack)
"raster_stack"


