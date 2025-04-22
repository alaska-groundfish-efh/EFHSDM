#' Make variance rasters
#' @description This function makes a non-parametric estimate of spatial variance based on the cross validation models. In order to do this, it needs to hold a lot of data in memory at once, so this can take awhile.
#' @param model.list list of models produced by the CV folds
#' @param raster.stack raster stack of the covariates used for the model
#' @param model.type character; the type of model ("maxnet","cloglog","hgam","gam")
#' @param scale.factor numeric; a scale factor to be applied
#' @param efh.break numeric; the EFH breakpoint for the full model, optionally creates an extra map
#' @importFrom terra extract
#' @importFrom stats density
#'
#' @return raster of estimate non-parametric variance in model predictions
#' @export
#'
#' @examples
MakeVarianceRasters<-function(model.list,            # a list of models for each cv fold
                              raster.stack,          # a stack of covariates for the model
                              model.type,            # the type of model
                              scale.factor = 1,      # should the output be scaled
                              efh.break=NA){         # the efh breakpoint from the full model


  # to speed things up and preserve memory, we're going to separate out the NAs
  data.spots<-which(is.na(terra::values(raster.stack[[1]]))==F)
  for(r in 2:terra::nlyr(raster.stack)){
    spots<-which(is.na(terra::values(raster.stack[[r]]))==F)
    data.spots<-data.spots[data.spots%in%spots]
  }
  data<-terra::extract(raster.stack,data.spots)

  #check for empty models
  list.index<-1
  model.list2<-list()
  for(m in 1:length(model.list)){
    if(is.list(model.list[[m]])){
      model.list2[[list.index]]<-model.list[[m]]
      list.index<-list.index+1
    }
  }

  out.data<-matrix(nrow=nrow(data),ncol = length(model.list2))

  # kind of complicated method of detecting an offset
  if(model.type!="maxnet"){
    terms<-AutodetectGAMTerms(model.list[[1]],hgam = "d")
    has.offset<-"offset"%in%terms$type
    facs<-terms$term[terms$type=="factor"]
    offset.name<-terms$term[terms$type=="offset"]

    # check for an offset, kind of kludgy solution for now
    if(model.type=="hgam" & has.offset){
      data<-data.frame(data,mean(model.list2[[1]]$offset[[1]]))
      pterms<-AutodetectGAMTerms(model.list[[1]],hgam = "p")
      facs<-unique(c(facs,pterms$term[pterms$type=="factor"]))
    }
    if(model.type%in%c("cloglog","gam")){
      data<-data.frame(data,mean(model.list2[[1]]$offset))
    }
    names(data)[ncol(data)]<-offset.name

    for(f in 1:length(facs)){
      data[,facs[f]]<-as.factor(data[,facs[f]])
    }
  }


  #progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(model.list2), style = 3)

  # loop through and get predictions for each model
  for(m in 1:length(model.list2)){
    if(model.type=="maxnet"){
      out.data[,m]<-exp(stats::predict(model.list2[[m]],newdata=data,type="link")+model.list2[[m]]$entropy)
    }
    if(model.type=="cloglog"){
      out.data[,m]<-exp(mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="link"))
    }
    if(model.type=="hgam"){
      out.data[,m]<-mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="response")
    }
    if(model.type=="gam"){
      out.data[,m]<-mgcv::predict.gam(object = model.list2[[m]],newdata=data,type="response")
    }
    utils::setTxtProgressBar(pb, m)
  }
  close(pb)
  # now tally things up

  raster.template<-terra::rast(raster.stack[[1]]) #isolate only one of the layers in the raster stack for the template
  #raster.template<-terra::rast(raster.stack)
  names(raster.template) <- "layer"

  variances<-apply(X = out.data*scale.factor,MARGIN = 1,FUN = stats::var)
  var.vec<-rep(NA,times=terra::ncell(raster.stack))
  var.vec[data.spots]<-variances
  var.raster<-terra::setValues(raster.template,values = var.vec)

  if(is.na(efh.break)==F){
    percents<-apply(X = out.data>efh.break,MARGIN = 1,FUN = sum)/ncol(out.data)
    per.vec<-rep(NA,times=terra::ncell(raster.stack))
    per.vec[data.spots]<-percents
    per.raster<-terra::setValues(raster.template,values = per.vec)
    return(list(var.raster,per.raster))
  }else{
    return(var.raster)
  }
}


