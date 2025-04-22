#' Make deviance explained table (for ensemble)
#'
#' @description Make an HTML table of deviance explained and fit metrics for a model ensemble.
#' @details For now, you need to include each of the elements for the ensemble table to work, so no shortcuts. It is kind of cumbersome and unlikely to be useful for others. May soon be deprecated.
#' @param model.names vector of names for the models
#' @param ensemble logical; is an ensemble included
#' @param weights vector of weights for each model other than the ensemble
#' @param converge.vec vector of T/F for whether models converged
#' @param abund.check.vec vector of T/F for whether models passed the abundance check
#' @param cor.method character; a method for calculating the r squared values
#' @param preds.table data frame with abundance observations and predictions for each model and the ensemble
#' @param scale.facs vector of scale factors for each model other than the ensemble
#' @param efh.breaks vector of efh breaks for each model including the ensemble
#' @param areas vector of areas for each model, including the ensemble
#' @param filename character; filename to write the html table to
#'
#' @return does not return anything, but instead writes a html table at filename
#' @export
#'
#' @examples
MakeEnsembleXtable<-function(model.names=c("maxnet","cloglog","hpoisson","poisson","negbin"),
                             ensemble=T,
                             weights=NULL,
                             converge.vec=NULL,
                             abund.check.vec=NULL,
                             cor.method="pearson",       #
                             preds.table=NULL,
                             scale.facs=NULL,
                             efh.breaks=NULL,
                             areas=NULL,
                             filename){

  # this part detects what data are available and formats things appropriately, sturdy to missing data
  dim.names<-"Model"

  if(sum(c(is.null(preds.table)==F,is.na(preds.table)==F))>1){
    N.col<-2
    dim.names<-c(dim.names,"N")
  }
  if(is.logical(converge.vec)){
    converge.col<-length(dim.names)+1
    dim.names<-c(dim.names,"Converged")
  }
  if(is.numeric(scale.facs)){
    sc.col<-length(dim.names)+1
    dim.names<-c(dim.names,"Scale Factor")
  }
  if(is.logical(abund.check.vec)){
    abund.check.col<-length(dim.names)+1
    dim.names<-c(dim.names,"Abund check")
  }
  if(is.numeric(weights)){
    w.col<-length(dim.names)+1
    dim.names<-c(dim.names,"weight")
  }
  if(sum(c(is.null(preds.table)==F,is.na(preds.table)==F))>1){
    pred.cols<-(length(dim.names)+1):(length(dim.names)+4)
    dim.names<-c(dim.names,"RMSE","Correlation","AUC","PDE")
  }
  if(is.numeric(efh.breaks)){
    break.col<-length(dim.names)+1
    dim.names<-c(dim.names,"EFH break")
  }

  if(is.numeric(areas)){
    area.col<-length(dim.names)+1
    dim.names<-c(dim.names,"Area")
  }

  n.models<-length(model.names)
  model.rows<-2:(n.models+1)

  etable<-array("",dim=c(3+n.models,length(dim.names)))

  # now go back and fill the rows in
  etable[,1]<-c("Model",model.names,".","ensemble")
  etable[1,]<-dim.names

  if(sum(c(is.null(preds.table)==F,is.na(preds.table)==F))>1){
    N.vec<-rep(NA,length(model.names))
    cor.vec<-rep(NA,length(model.names))
    rmse.vec<-rep(NA,length(model.names))
    auc.vec<-rep(NA,length(model.names))
    pde.vec<-rep(NA,length(model.names))
    for(m in 1:length(model.names)){
      if(model.names[m]%in%preds.table$Model){
        pred.dat<-stats::na.omit(subset(preds.table,Model==model.names[m]))
        pred.dat$pred[which(is.infinite(pred.dat$pred))]<-10^10
        pred.cvpred<-ifelse(RMSE(pred = pred.dat$pred,obs = pred.dat$abund)<
                              RMSE(pred = pred.dat$cvpred,obs = pred.dat$abund),"cvpred","pred")

        prob.cvprob<-ifelse(pred.cvpred=="cvpred","cvprob","prob")

        N.vec[m]<-sum(pred.dat$abund>0,na.rm=T)
        rmse.vec[m]<-RMSE(pred = pred.dat[,pred.cvpred],obs = pred.dat$abund)
        cor.vec[m]<-stats::cor(pred.dat$abund,round(pred.dat[,pred.cvpred],2),method=cor.method)
        auc.vec[m]<-PresenceAbsence::auc(data.frame(1:nrow(pred.dat),pred.dat$abund,pred.dat[,prob.cvprob]))[[1]]
        pde.vec[m]<-PDE(obs = pred.dat$abund,pred = pred.dat[,pred.cvpred])
      }
    }

    digs<-ifelse(min(rmse.vec,na.rm=T)>1000,0,ifelse(min(rmse.vec,na.rm = T)<10,2,1))

    etable[model.rows,N.col]<-N.vec
    etable[model.rows,pred.cols[1]]<-format(round(rmse.vec,digs),nsmall=digs,big.mark = ",")
    etable[model.rows,pred.cols[2]]<-format(round(cor.vec, 2),nsmall=2)
    etable[model.rows,pred.cols[3]]<-format(round(auc.vec, 2),nsmall=2)
    etable[model.rows,pred.cols[4]]<-format(round(pde.vec, 2),nsmall=2)

    if(ensemble){

      ensemble.dat<-stats::na.omit(subset(preds.table,Model=="ensemble"))
      etable[n.models+3,N.col]<-sum(ensemble.dat$abund>0,na.rm=T)
      etable[n.models+3,pred.cols[1]]<-format(round(RMSE(obs=ensemble.dat$abund,pred=ensemble.dat$pred),digs),nsmall=digs,big.mark = ",")
      etable[n.models+3,pred.cols[2]]<-format(round(stats::cor(ensemble.preds$abund,round(ensemble.preds$pred,2),method=cor.method), 2),nsmall=2)
      etable[n.models+3,pred.cols[3]]<-format(round(PresenceAbsence::auc(data.frame(1:nrow(ensemble.dat),ensemble.preds$abund,ensemble.preds$prob))[[1]], 2),nsmall=2)
      etable[n.models+3,pred.cols[4]]<-format(round(PDE(obs = ensemble.preds$abund,pred = ensemble.preds$pred), 2),nsmall=2)
    }
  }
  if(is.logical(converge.vec)){
    converge.vec2<-rep(F,length(model.names))
    for(m in 1:length(model.names)){
      converge.vec2[m]<-converge.vec[which(names(converge.vec)==model.names[m])]
    }
    etable[model.rows,converge.col]<-converge.vec2
    if(ensemble){etable[n.models+3,converge.col]<-TRUE}
  }
  if(is.logical(abund.check.vec)){
    abund.check.vec2<-rep(F,length(model.names))
    for(m in 1:length(model.names)){
      abund.check.vec2[m]<-abund.check.vec[which(names(abund.check.vec)==model.names[m])]
    }
    etable[model.rows,abund.check.col]<-abund.check.vec2
    if(ensemble){etable[n.models+3,abund.check.col]<-TRUE}
  }
  if(is.numeric(scale.facs)){
    scale.fac.vec<-rep(NA,length(model.names))
    for(m in 1:length(model.names)){
      if(model.names[m]%in%names(scale.facs)){scale.fac.vec[m]<-scale.facs[which(names(scale.facs)==model.names[m])]}
    }
    etable[model.rows,sc.col]<-format(round(scale.fac.vec,2),nsmall=2)
    if(ensemble){etable[n.models+3,sc.col]<-1}
  }
  if(is.numeric(weights)){
    weight.vec<-rep(0,length(model.names))
    for(m in 1:length(model.names)){
      if(model.names[m]%in%names(weights)){weight.vec[m]<-weights[which(names(weights)==model.names[m])]}
    }
    etable[model.rows,w.col]<-format(round(weight.vec,2),nsmall=2)
    if(ensemble){etable[n.models+3,w.col]<-1}
  }
  if(is.numeric(efh.breaks)){
    breaks.vec<-rep(NA,length(model.names))
    for(m in 1:length(model.names)){
      if(model.names[m]%in%names(efh.breaks)){breaks.vec[m]<-efh.breaks[which(names(efh.breaks)==model.names[m])]}
    }
    etable[model.rows,break.col]<-format(round(breaks.vec,2),nsmall=2)
    if(ensemble){etable[n.models+3,break.col]<-format(round(efh.breaks[which(names(efh.breaks)=="ensemble")],2),nsmall=2)}
  }
  if(is.numeric(areas)){
    area.vec<-rep(NA,length(model.names))
    for(m in 1:length(model.names)){
      if(model.names[m]%in%names(areas)){area.vec[m]<-areas[which(names(areas)==model.names[m])]}
    }
    etable[model.rows,area.col]<-format(round(area.vec,-2),big.mark = ",")
    if(ensemble){etable[n.models+3,area.col]<-format(round(areas[which(names(areas)=="ensemble")],-2),big.mark = ",")}
  }

  etable[trimws(etable)=="NA"]<-"--"

  ensemble.table <- xtable::xtable(etable)
  xtable::print.xtable(ensemble.table, type = "html", file = filename,
                       include.rownames = getOption("xtable.include.rownames", FALSE), html.table.attributes = 2,
                       include.colnames = getOption("xtable.include.colnames", FALSE), hline.after = getOption("xtable.hline.after",
                                                                                                               c(-1, 1, nrow(ensemble.table))))
}
