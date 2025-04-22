#' Get ensemble variance
#'
#' @description This function computes the variance in abundance based on a formula from Burnham and Anderson
#' @param model.weights a vector of model weights
#' @param variance.list list of rasters of the variance in model predictions
#' @param abund.list list of rasters of the predicted abundance from each model
#' @param ensemble.abund raster of ensemble predicted abundance
#'
#' @return raster of ensemble predicted variance in predictions
#' @export
#'
#' @examples
GetEnsembleVariance<-function(model.weights,variance.list,abund.list,ensemble.abund){

  keepers<-which(model.weights>0)
  weights2<-model.weights[keepers]
  variance.list2<-list()
  abund.list2<-list()

  for(k in 1:length(keepers)){
    variance.list2[[k]]<-variance.list[[keepers[k]]]
    abund.list2[[k]]<-abund.list[[keepers[k]]]
  }

  data.spots<-which(is.na(terra::values(abund.list2[[1]]))==F)
  if(length(abund.list2)>1){
    for(r in 2:length(abund.list2)){
      spots<-which(is.na(terra::values(abund.list2[[r]]))==F)
      data.spots<-data.spots[data.spots%in%spots]
    }
  }
  e.abund<-as.numeric(terra::extract(ensemble.abund,data.spots, raw=TRUE))

  dat<-matrix(nrow=length(data.spots),ncol=length(keepers))

  for(m in 1:length(keepers)){
    a.dat<-terra::extract(abund.list2[[m]],data.spots, raw=TRUE)
    v.dat<-terra::extract(variance.list2[[m]],data.spots, raw=TRUE)

    dat[,m]<-weights2[m]*sqrt(v.dat+(a.dat-e.abund)^2)
  }

  std.error<-apply(dat,MARGIN = 1,FUN = sum)

  out.raster<-abund.list2[[1]]
  val.vec<-terra::values(out.raster)
  val.vec[data.spots]<-std.error

  out.raster<-terra::setValues(out.raster,val.vec)
  return(out.raster)
}

