#' Make model ensemble
#' @description Make an ensemble model vector with weights for each model
#' @details This function calculates the weights for an ensemble model using a very simple rmse method. Future work may allow for more complex methods. One can supply a list of dataset (old way), or one can more sensibly just supply a named vector of RMSE values (better way)
#' @param ensemble.list list of data frames from which rmse can be calculated, not recommended
#' @param rmse named vector of rmse values to calculate weights
#' @param names.vec vector of names corresponding to the list or rmse vector
#' @param minimum numeric; drops any model given less weight than this value
#' @param model.types vector of names of models, only one model of a given type will advance
#'
#' @return a named vector with weights for each model
#' @export
#'
#' @examples
MakeEnsemble<-function(ensemble.list=NA,
                       rmse=NA,
                       names.vec=NA,
                       minimum=NA,
                       model.types=NULL ){


  # if the rmse vector isn't supplied, compute it from the ensemble.list
  if(any(is.na(rmse))){
    if(is.na(names.vec)){
      names.vec<-names(ensemble.list)
    }

    rmse<-vector(length=length(ensemble.list))

    # find the rmse for each set
    for(m in 1:length(ensemble.list)){
      ensemble.dat<-stats::na.omit(ensemble.list[[m]])
      if(nrow(ensemble.dat)>0){
        rmse[m]<-sqrt(sum((ensemble.dat$abund-ensemble.dat$pred)^2)/nrow(ensemble.dat))
      }else{
        rmse[m]<-NA
      }
    }
  }else{
    if(is.na(names.vec)){
      names.vec<-names(rmse)
    }
  }

  # these days we are using the square of RMSE to weight things
  rmse2<-rmse^2

  # weight it, by the inverse
  weights<-(1/(rmse2))/sum(1/(rmse2),na.rm=T)
  names(weights)<-names.vec
  weights[is.na(weights)]<-0

  # apply any adjustments
  if(length(model.types)>0){
    for(m in unique(model.types)){
      index<-which(model.types==m)
      weights2<-weights[index]
      keep<-index[which.max(weights2)]
      drop<-index[index!=keep]
      weights[drop]<-0
    }
    weights<-weights/sum(weights)
  }

  # drop anything that's below the cutoff
  if(is.na(minimum)==F){
    included<-which(weights>=minimum)
    set.to.zero<-which(weights<minimum)
    weights2<-weights[included]
    weights2<-weights2/sum(weights2)
    weights[set.to.zero]<-0
    weights[included]<-weights2
  }
  return(weights)
}
