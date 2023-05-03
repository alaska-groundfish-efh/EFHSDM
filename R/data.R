#' Bathymetry data from the GOA
#'
#' A raster layer containing bathymetry covariates for SDMs.
#'
#' @format A raster layer with 12 slots
#' \describe{
#'   \item{name}{character string with path to raster layer}
#'   \item{datanotation}{Raster file native format}
#'   ...
#' }
#' @source Compiled by Jeremy Harris using 1-km resolution data from Mark Zimmerman
"GOA_bathy"

#' Region data for abundance for all the EFH species included in the 2022 Review
#'
#' The variables are as follows:
#'
#' \itemize{
#'   \item hauljoin.
#'   \item region. EFH area (EBS, GOA, or AI)
#'   \item subregion.
#'   \item vessel. Vessel ID from the survey; matches RACEBASE, FOSS, others
#'   \item cruise. Cruise ID number.
#'   \item haul. Haul number.
#'   ...
#' }
#'
#' @docType data
#' @keywords datasets
#' @name region_data_all
#' @usage data(region_data_all)
#' @format A data frame with 15592 rows and 184 variables
NULL
