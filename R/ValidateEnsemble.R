#' Validate ensemble
#' @description This function makes residual plots for the ensemble model and returns a helpful data frame of ensemble predictions
#' @param pred.list list of data frames containing observations and predictions from constituent models
#' @param model.weights named vector of weights for each model
#' @param method character; method for the cor() function
#' @param output logical; should the data frame be returned, or just plots
#' @param make.plots logical; should the plots be output, or just the data frame
#' @param key character; a column name present in pred.list that works as a case identifier; usually "hauljoin"
#' @param latlon logical; should lat and lon be included in the output table, must be present in pred.list
#' @param group character; name of a column to be included in the output table representing CV folds
#' @param histogram logical; should histograms be plotted (they often don't look nice)
#'
#' @return data frame of ensemble observations vs predictions (pred and cvpred are the same for the ensemble)
#' @export
#'
#' @examples
ValidateEnsemble<-function(pred.list,
                           model.weights,
                           method="pearson",
                           output=T,
                           make.plots=T,
                           key=NA,
                           latlon=T,
                           group=NA,
                           histogram=F){

  ensemble.preds<-data.frame(abund=pred.list[[1]]$abund,pred=0)

  if(latlon){
    ensemble.preds<-cbind(pred.list[[1]][,c("lon","lat")],ensemble.preds)
  }
  if(is.na(group)==F){
    ensemble.preds<-cbind(pred.list[[1]][,group],ensemble.preds)
    names(ensemble.preds)[1]<-group
  }
  if(is.na(key)==F){
    ensemble.preds<-cbind(pred.list[[1]][,key],ensemble.preds)
    names(ensemble.preds)[1]<-key
  }

  # Now use the weights to find out the predictions
  for(m in 1:length(pred.list)){
    if(model.weights[m]>0){
      ensemble.preds$pred<-ensemble.preds$pred+pred.list[[m]]$pred*model.weights[m]
    }
  }
  ensemble.preds$prob<-1-stats::dpois(0,ensemble.preds$pred)
  ensemble.preds$cvprob<-ensemble.preds$prob

  keepers<-which(is.na(ensemble.preds$pred)==F & is.infinite(ensemble.preds$pred)==F)

  # compute the summary stats for the output
  ensemble.preds2<-ensemble.preds
  if(method=="spearman"){
    ensemble.preds2$abund<-rank(ensemble.preds$abund)
    ensemble.preds2$pred<-rank(ensemble.preds$pred)
  }

  ensemble.rmse<-sqrt(sum((ensemble.preds$abund-ensemble.preds$pred)^2,na.rm=T)/nrow(ensemble.preds))
  regr <- stats::lm(ensemble.preds2$pred[keepers]~ensemble.preds2$abund[keepers])
  rsqr <- summary(regr)$r.squared

  # ploting portion
  if(make.plots==T){
    old.par<-graphics::par()[c("mfcol","family","mar","xaxs","yaxs")]
    graphics::par(mfcol = c(2+as.integer(histogram),1), family = "sans", mar = c(4,4,3,1))

    stats::qqnorm((ensemble.preds2$abund - ensemble.preds2$pred))
    stats::qqline((ensemble.preds2$abund - ensemble.preds2$pred))

    if(histogram){graphics::hist((ensemble.preds2$abund - ensemble.preds2$pred), xlab = "Residuals", main = "")}

    pred.max <- ifelse(method=="pearson",stats::quantile(ensemble.preds2$pred[keepers],probs=.99,na.rm=T),nrow(ensemble.preds2))
    abund.max <- ifelse(method=="pearson",stats::quantile(ensemble.preds2$pred[keepers],probs=.99,na.rm=T),nrow(ensemble.preds2))
    plot.max<-base::max(pred.max,abund.max)*1.1

    plot(y=ensemble.preds2$pred[keepers], x=ensemble.preds2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
         ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
         xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
    graphics::abline(coef = c(0,1), lty = 2)
    graphics::abline(regr,col=2)
    graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(rsqr,2)), pos = 4)

    suppressWarnings(graphics::par(old.par))
  }

  print(paste("Ensemble",method,"Rsq =",round(rsqr,3)))
  print(paste("Ensemble RMSE =",round(ensemble.rmse,2)))

  ensemble.preds$error<-ensemble.preds$pred-ensemble.preds$abund

  if(output==T){return(ensemble.preds)}
}

