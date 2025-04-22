#' Make ensemble abundance
#' @description Make the ensemble into a raster.
#' @details The abundance rasters should be a list, and must match the ordering of the weights. After this, use the FindEFHbreaks function and others to generate an EFH map as normal.
#' @param model.weights a vector of model weights
#' @param abund.list list of abundance rasters corresponding to model weights
#' @param filename character; a file name for the abundance raster
#'
#' @return raster formed from weighted average of abund.list
#' @export
#'
#' @examples
MakeEnsembleAbundance<-function(model.weights,
                                abund.list,
                                filename=""){

  # intialize a raster with appropriate properties
  good.abund<-which(is.na(abund.list)==F)
  new.abund<-terra::rast(abund.list[[good.abund[1]]])
  new.abund<-terra::setValues(new.abund,values = 0)

  #To make the raster, we'll need another loop
  for(i in 1:length(model.weights)){
    if(model.weights[i]>0){
      new.abund<-new.abund+abund.list[[i]]*model.weights[i]
    }
  }
  if(filename!=""){
    terra::writeRaster(x = new.abund,filename = filename,overwrite=T)
  }
  return(new.abund)
}
