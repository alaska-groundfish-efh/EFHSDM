#' Calculate RMSE
#'
#' @description A quick function to calculate the RMSE.
#' @param pred vector of predictions
#' @param obs vector of observations
#'
#' @return returns RMSE of obs and preds
#' @export
#'
#' @examples
RMSE<-function(pred,obs){
  keep<-which(is.na(pred)==F & is.na(obs)==F)
  return(sqrt(sum((pred[keep]-obs[keep])^2)/length(pred[keep])))
}


