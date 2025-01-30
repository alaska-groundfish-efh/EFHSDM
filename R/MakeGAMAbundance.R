#' Make GAM abundance
#'
#' @description This function makes abundance rasters for the GAMs. When an offset is present in the model, generates a map for the mean value of the offset
#' @param model a mgcv GAM model object
#' @param r.stack  raster stack with all the necessary covariates
#' @param scale.factor numeric; scale factor that will multiply the model predictions
#' @param filename filename for the output
#' @param land raster; a mask to be applied to the output usually landmasses
#' @param mask raster; an additional mask to apply to the output, if needed
#'
#' @return raster; usually representing abundance
#' @export
#'
#' @examples
MakeGAMAbundance <- function(model,
                             r.stack,
                             scale.factor = 1,
                             filename = "",
                             land = NULL,
                             mask = NULL) {

  # correct a common mistake
  if (is.na(filename) | is.null(filename)) {
    filename <- ""
  }

  model.terms <- AutodetectGAMTerms(model)
  if (model$family$family == "ziplss") {
    dterms <- model.terms[[1]]
    pterms <- model.terms[[2]]
    model.terms <- unique(rbind(dterms, pterms))
  }

  # now, need to add a check for the offset, and make a dummy raster if there is one
  if ("offset" %in% model.terms$type) {
    off.name <- model.terms$term[which(model.terms$type == "offset")]
    off.val <- ifelse(is.list(model$offset), mean(model$offset[[1]]), mean(model$offset))
    off.raster <- terra::rast(ext = terra::ext(r.stack), crs = terra::crs(r.stack), nrow = terra::nrow(r.stack), ncol = terra::ncol(r.stack),
                              vals = off.val
    )
    names(off.raster) <- off.name
    r.stack <- c(r.stack, off.raster)
  }

  # will also need to detect factors and format them into a list
  gam.factors <- model.terms$term[model.terms$type == "factor"]
  gam.factors2 <- list()
  if (length(gam.factors) > 0) {
    for (t in 1:length(gam.factors)) {
      range <- terra::subset(x = r.stack, subset = which(names(r.stack) == gam.factors[t]))
      gam.factors2[[t]] <- sort(unique(stats::na.omit(terra::values(range))))
    }
    names(gam.factors2) <- gam.factors
  } else {
    gam.factors2 <- NULL
  }

  ## Detect the link function and make the predictions
  pred.type <- ifelse(model$family$link == "cloglog", "link", "response")[1]

  # predict raster function malfunctions if ziplss model has no factors, so need a special case
  if(model$family$family == "ziplss" & is.null(gam.factors2)){
    r.vals<-as.data.frame(terra::values(r.stack))
    if ("offset" %in% model.terms$type) {
      r.vals<-cbind(r.vals,data.frame(off.val))
      names(r.vals[ncol(r.vals)])<-off.name
    }
    pred.vals<-mgcv::predict.gam(model,newdata=r.vals,type="response")
    predict.raster<-terra::setValues(r.stack, pred.vals)
  }else{
    predict.raster <- terra::predict(r.stack, model, factors = gam.factors2, progress = "text",
                                     overwrite = TRUE, type = pred.type)
  }

  # Detects and apply the cloglog approximation if appropriate
  if (model$family$link[1] == "cloglog") {
    predict.raster <- exp(predict.raster)
  }

  predict.raster <- predict.raster * scale.factor

  if (filename != "") {
    terra::writeRaster(x = predict.raster, filename = filename, overwrite = TRUE)
  }

  # now apply masks if appropriate
  if (is.null(land) == F) {
    predict.raster <- terra::mask(predict.raster, land,
                                  inverse = TRUE, overwrite = TRUE,
                                  filename = filename
    )
  }
  if (is.null(mask) == F) {
    predict.raster <- terra::mask(predict.raster, mask,
                                  overwrite = TRUE,
                                  filename = filename
    )
  }
  return(predict.raster)
}
