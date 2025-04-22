#' Make cross-validation plots
#'
#' @description Use the model and the cross validated errors to make some residual plots
#' @details This is low priority, but could be much better.
#' @param error.data data frame containing observation, predictions, and CV predictions
#' @param method character; a method to be passed to the cor function
#' @param make.hist should the histograms be plotted, they sometimes fail and may need to be turned off
#'
#' @return nothing, but creates some plots
#' @export
#'
#' @examples
MakeCrossValidationPlots<-function(error.data,           # a data frame, typically from the CrossvalidateModel function
                                   method="pearson",     # a method for the residuals, accepts pearson and spearman
                                   make.hist=T){         # should the histograms be plotted, they sometimes fail and may need to be turned off

  #remove any bad data
  keepers<-unique(which(is.infinite(error.data$pred)==F & is.infinite(error.data$cvpred)==F))
  if(length(keepers)<nrow(error.data)){
    warning(paste0(nrow(error.data)-length(keepers)," Infinite values detected and removed; estimates of error may be too low"))
  }

  error.data2<-error.data
  # compute the summary stats for the output
  if(method=="spearman"){
    error.data2$abund<-rank(error.data2$abund)
    error.data2$pred<-rank(error.data2$pred)
    error.data2$cvpred<-rank(error.data2$cvpred)
  }

  main.regr<-stats::lm(error.data2$pred[keepers]~error.data2$abund[keepers])
  main.r2<-summary(main.regr)$r.squared
  main.rmse<-sqrt(sum((stats::na.omit(error.data$abund[keepers]-error.data$pred[keepers]))^2)/nrow(stats::na.omit(error.data[keepers,])))

  #now need to do the cv tests, which should already be in a nice format from the CV function
  cv.regr<-stats::lm(error.data2$cvpred[keepers]~error.data2$abund[keepers])
  cv.r2<-summary(cv.regr)$r.squared
  cv.rmse<-sqrt(sum((stats::na.omit(error.data$abund[keepers]-error.data$cvpred[keepers]))^2)/nrow(stats::na.omit(error.data[keepers,])))

  print(paste("Full model",method,"Rsq =",round(main.r2,2)))
  print(paste("CV",method,"Rsq =",round(cv.r2,2)))
  print(paste("Full model RMSE =",round(main.rmse,3)))
  print(paste("CV RMSE =",round(cv.rmse,3)))

  # make all the plots
  old.par<-graphics::par()[c("mfcol","family","mar","xaxs","yaxs")]
  graphics::par(mfcol = c(ifelse(make.hist==T,3,2),2), family = "sans", mar = c(4,4,3,1))

  stats::qqnorm((error.data$pred[keepers] - error.data$abund[keepers]), main = "Model Predictions")
  stats::qqline((error.data$pred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){graphics::hist((error.data$pred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}
  pred.max <- ifelse(method=="pearson",stats::quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",stats::quantile(error.data2$pred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  plot(y=error.data2$pred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  graphics::abline(coef = c(0,1), lty = 2)
  graphics::abline(main.regr,col=2)
  graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(main.r2,2)), pos = 4)


  #Plots for test/CV data
  pred.max <- ifelse(method=="pearson",stats::quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  abund.max <- ifelse(method=="pearson",stats::quantile(error.data2$cvpred[keepers],probs=.99,na.rm=T),nrow(error.data2))
  plot.max<-max(pred.max,abund.max)*1.1

  stats::qqnorm((error.data$cvpred[keepers] - error.data$abund[keepers]), main = "Test data")
  stats::qqline((error.data$cvpred[keepers] - error.data$abund[keepers]))
  if(make.hist==T){graphics::hist((error.data$cvpred[keepers] - error.data$abund[keepers]), xlab = "Residuals", main = "")}

  plot(y=error.data2$cvpred[keepers], x=error.data2$abund[keepers], ylim = c(0,plot.max), xlim = c(0,plot.max),
       ylab = ifelse(method=="pearson","Predicted","Predicted Ranks"),
       xlab = ifelse(method=="pearson","Observed","Observed Ranks"),main = "", pch = 20)
  graphics::abline(coef = c(0,1), lty = 2)
  graphics::abline(cv.regr,col=2)
  graphics::text(1, plot.max*.9, paste(method,"R-squared = ", signif(cv.r2,2)), pos = 4)
  suppressWarnings(graphics::par(old.par))
}

