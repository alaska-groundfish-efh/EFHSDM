#' Get MaxEnt coefficients
#'
#' @description Get the number of coefficients for each variable in a maxnet model
#' @details If you want things to match the gams, supply a list for maxnet2d
#' @param model a fitted maxnet model
#' @param maxnet2d a list of vectors of variables to be treated as 2D; used for consistency with the GAMs
#'
#' @return a named vector listing the number of coefficients that apply to each term
#' @export
#'
#' @examples
MaxnetCoefs<-function(model,
                      maxnet2d=NULL){
  # This first section that is to find how many coefficients each covariate uses
  covar.names<-names(model$samplemeans)
  beta.names<-strsplit(names(model$betas),split=c("\\(|\\:|\\^|\\)"))

  if(length(maxnet2d)>0){
    coef.vec2<-vector(length=length(maxnet2d))
    vars1d<-covar.names[covar.names%in%unlist(maxnet2d)==F]
    maxnet2d.name<-unlist(lapply(X = maxnet2d,FUN = function(x){paste(x,collapse="*")}))

    for(i in 1:length(maxnet2d)){
      coef2d1<-unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=maxnet2d[[i]][1]))
      coef2d2<-unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=maxnet2d[[i]][2]))

      coef.vec2[i]<-sum((coef2d1+coef2d2)>0)
      names(coef.vec2)<-maxnet2d.name
    }
  }else{
    coef.vec2<-NULL
    vars1d<-covar.names
  }

  coef.vec<-vector(length=length(vars1d))
  for(i in 1:length(vars1d)){
    coef.vec[i]<-sum(unlist(lapply(X=beta.names,FUN=function(x,y){y%in%x},y=vars1d[i])))
  }
  names(coef.vec)<-vars1d
  return(c(coef.vec2,coef.vec))
}

