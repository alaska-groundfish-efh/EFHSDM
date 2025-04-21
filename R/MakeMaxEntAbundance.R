#' Make abundance/distribution raster from MaxEnt model
#'
#' @description Make an abundance/distribution raster from any dismo or maxnet model. Depending on the settings, it will apply various masks and the cloglog transformation.
#' @param model a fitted maxnet model
#' @param maxent.stack a raster stack containing all covariates
#' @param scale.fac numeric; scale factor that will multiply the model predictions
#' @param land raster; a mask to be applied to the output usually landmasses
#' @param mask raster; an additional mask to apply to the output, if needed
#' @param type character; use "maxnet" for prob suitable habitat or "cloglog" for approximate abundance
#' @param clamp Logical; just leave this as F; should the covariates be restricted to the range observe in fitting the model?
#' @param filename a filename to save results; "" writes to active memory
#'
#' @return a raster map with the desired prediction
#' @export
#' @importFrom terra values
#' @importFrom terra setValues
#'
#' @examples

MakeMaxEntAbundance <- function(model,
                                maxent.stack,
                                scale.fac = 1,
                                land = NULL,
                                mask = NULL,
                                type = "cloglog",
                                clamp = F,
                                filename = "") {
  # correct a common mistake
  if (is.null(filename) || is.na(filename)) {
    filename <- ""
  }

  # Identify the type and make the main prediction
  type <- tolower(type)
  # Somewhat counterintuitive, but this is the type if using the cloglog link to make an abundance estimate
  if (type == "cloglog") {
    # since ENMeval 2.0, they got rid of the useful function and I need to make the predictions the long way
    dat <- terra::values(maxent.stack)

    # using predict with maxnet will quietly remove the NAs, so need to track them manually
    na.spots <- which(apply(X = dat, MARGIN = 1, FUN = function(x) {
      return(any(is.na(x)))
    }))
    dat.spots <- which(seq(1:nrow(dat)) %in% na.spots == F)

    preds <- stats::predict(model, dat[dat.spots, ], type = "link")
    preds2 <- exp(preds + model$ent) * scale.fac
    new.vals <- vector(length = nrow(dat))
    new.vals[na.spots] <- NA
    new.vals[dat.spots] <- preds2
    habitat.prediction <- terra::setValues(x = terra::rast(maxent.stack[[1]]), values = new.vals) # Creates a SpatRaster object. terra needs to call just one raster in the stack bc it's just setting dimensions.
  }
  # this makes a habitat suitability map from a maxnet model
  if (type == "maxnet") {
    dat <- terra::values(maxent.stack)

    # using predict with maxnet will quietly remove the NAs, so need to track them manually
    na.spots <- which(apply(X = dat, MARGIN = 1, FUN = function(x) {
      return(any(is.na(x)))
    }))
    dat.spots <- which(seq(1:nrow(dat)) %in% na.spots == F)

    preds <- stats::predict(model, dat[dat.spots, ], type = "cloglog")
    new.vals <- vector(length = nrow(dat))
    new.vals[na.spots] <- NA
    new.vals[dat.spots] <- preds
    habitat.prediction <- terra::setValues(x = terra::rast(maxent.stack[[1]]), values = new.vals)
  }
  # need to add a check to see about the strange problems with the EBS
  if (is.null(land) == F) {
    habitat.prediction <- raster::extend(x = habitat.prediction, y = land)
  }

  # For some reason, the crs info isn't always carrying over
  terra::crs(habitat.prediction) <- terra::crs(maxent.stack)
  if (filename != "") {
    terra::writeRaster(x = habitat.prediction, filename = filename, overwrite = TRUE)
  }

  # Apply additional masks if necessary
  if (is.null(land) == F) {
    habitat.prediction <- terra::mask(habitat.prediction, land,
      inverse = T, overwrite = TRUE,
      filename = filename
    )
  }
  # Apply additional masks if necessary
  if (is.null(mask) == F) {
    habitat.prediction <- terra::mask(habitat.prediction, mask,
      overwrite = TRUE,
      filename = filename
    )
  }
  return(habitat.prediction)
}
