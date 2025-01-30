#' Calculate PDE
#'
#' @description A quick function to calculate the PDE
#' @param pred vector of predictions
#' @param obs vector of observations
#'
#' @return returns the estimate percent deviance explained, assuming a Poisson distribution
#' @export
#'
#' @examples
PDE<-function(obs,pred){
  term1<-obs*log(obs/mean(obs))
  term1[is.nan(term1)]<-0
  term2<-obs-mean(obs)

  nulldev<-2*sum(term1-term2)

  pred[pred<.00001]<-.00001
  term1<-obs*log(obs/pred)
  term1[is.nan(term1)]<-0
  term2<-obs-pred

  pdev<-2*sum(term1-term2)
  return(1-(pdev/nulldev))
}
