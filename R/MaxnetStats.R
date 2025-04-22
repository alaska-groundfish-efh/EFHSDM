#' Get MaxEnt stats
#' @description This function institutes a jackknife estimate of the importance of each covariate for a maxnet model. Can sometimes run slowly if there are convergence issues.
#' @param model a fitted maxnet model
#' @param maxnet2d a list of vectors of variables to be treated as 2D; used for consistency with the GAMs
#' @param regmult numeric; a regularization penalty multiplier
#' @param data data frame; typically the same data used to fit the model
#' @param species character; the name of the column for the dependent variable in the model
#'
#' @return a named vector with estimates of the relative deviance explained by each term
#' @export
#'
#' @examples
MaxnetStats<-function(model,                 # a maxnet model
                      maxnet2d=NULL,         # a list of terms that are to be considered jointly
                      regmult=1,             # a multiplier for regularization
                      data,                  # the data used to fit that model
                      species){              # the species or column of data for the dependent variable

  # This first section that is to find how many coefficients each covariate uses
  covar.names<-names(model$samplemeans)

  if(length(maxnet2d)>0){
    coef.vec2<-vector(length=length(maxnet2d))
    vars1d<-covar.names[covar.names%in%unlist(maxnet2d)==F]
    maxnet2d.name<-unlist(lapply(X = maxnet2d,FUN = function(x){paste(x,collapse="*")}))
  }else{
    vars1d<-covar.names
  }

  #detect factors, and remove them from the vars1d vector
  min.names<-names(model$varmin)
  facs<-vars1d[vars1d%in%min.names==F]
  vars1d<-vars1d[vars1d%in%min.names]

  pb <- utils::txtProgressBar(min = 0, max = length(c(maxnet2d,vars1d,facs)), style = 3)
  pb.i<-1
  # This part makes the estimate of deviance explained for each covariate
  # each category of covariate needs to be handled separately
  if(length(maxnet2d)>0){
    dev.vec2d<-rep(NA,length=length(maxnet2d))
    for(i in 1:length(maxnet2d)){
      #print(paste0("testing deviance for ",maxnet2d.name[i]))
      d.covars<-covar.names[covar.names%in%maxnet2d[[i]]==F]

      try(test.model<-FitMaxnet(data=data,species = species,vars = d.covars,facs = facs,regmult = regmult))
      if(exists("test.model")){
        dev.vec2d[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
        rm(test.model)
        utils::setTxtProgressBar(pb, pb.i)
        pb.i<-pb.i+1
      }else{
        close(pb)
        break
      }
    }
    names(dev.vec2d)<-maxnet2d.name
  }else{
    dev.vec2d<-NULL
    maxnet2d.name<-NULL
  }

  # 1D covariates
  dev.vec<-rep(NA,length=length(vars1d))
  for(i in 1:length(vars1d)){
    #print(paste0("testing deviance for ",vars1d[i]))
    try(test.model<-FitMaxnet(data=data,species = species,vars = c(unlist(maxnet2d),vars1d[-i]),facs = facs,regmult = regmult))
    if(exists("test.model")){
      dev.vec[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i<-pb.i+1
    }else{
      close(pb)
      break
    }
  }
  names(dev.vec)<-vars1d

  # factors
  fac.vec<-rep(NA,length=length(facs))
  for(i in 1:length(facs)){
    #print(paste0("testing deviance for ",facs[i]))
    try(test.model<-FitMaxnet(data=data,species = species,vars = c(unlist(maxnet2d),vars1d),facs = facs[-i],regmult = regmult))
    if(exists("test.model")){
      fac.vec[i]<-test.model$dev.ratio[length(test.model$dev.ratio)]
      rm(test.model)
      utils::setTxtProgressBar(pb, pb.i)
      pb.i<-pb.i+1
    }else{
      close(pb)
      break
    }
  }
  close(pb)

  names(fac.vec)<-facs

  out.dev.vec<-c(dev.vec2d,dev.vec,fac.vec)
  if(sum(is.na(out.dev.vec))==0){
    dev.lost<-1-out.dev.vec/model$dev.ratio[length(model$dev.ratio)]
    dev.exp<-dev.lost/sum(dev.lost,na.rm=T)*100

    #sometimes you end up with negative deviance, so correct that
    if(min(dev.exp)<0){
      dev.exp<-dev.exp-min(dev.exp)
      dev.exp<-dev.exp/sum(dev.exp)*100
    }
  }else{
    dev.exp<-rep(NA,times=length(c(maxnet2d,vars1d,facs)))
  }
  names(dev.exp)<-c(maxnet2d.name,vars1d,facs)
  return(dev.exp)
}
