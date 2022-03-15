# This script contains some functions for making the tables for summary stats

rpackages <- c("xtable", "XML")

which_not_installed <- which(rpackages %in% rownames(installed.packages()) == FALSE)

if(length(which_not_installed) > 1){
  install.packages(rpackages[which_not_installed], dep = TRUE)
}
rm(rpackages,which_not_installed)

require(xtable)
require(XML)

#' Make a table with constituent model statistics
#'
#' @description Make a table of deviance explained for a model object.
#' @details This is a slightly over-complicated function, but it automatically adjusts the table to whatever data is supplied to it. Basically, first it makes a lot of notes as to where things should go, then it takes a second pass construct the table and put everything in the right place. Known issue: does not calculate the degrees of freedom for factors.
#' @param model a fitted model object produced by mgcv::gam() or maxnet::maxnet()
#' @param rsq.method character; a method for the cor() function
#' @param model.type character; the type of model used, does not autodetect options: maxnet, cloglog, hgam, gam
#' @param regmult numeric; a regularization constant for maxnet models, ignored otherwise
#' @param scale.factor numeric; a scale factor for abundance predictions
#' @param efh.break numeric; a breakpoint or threshold for EFH area
#' @param nice.names data frame; table linking abbreviated covariate names to nicer versions
#' @param filename character; a file name for the html table produced
#' @param abund vector of abundance observations
#' @param train vector of predictions from the model trained ont he full data set
#' @param test vector of out-of-bag predictions
#' @param group vector of group labels
#' @param forecast vector of predictions from a third group (not used much anymore)
#' @param ncoefs vector of coefficients used in maxnet models
#' @param devs vector of deviance explained estimates
#' @param area numeric; total EFH area
#'
#' @return no return value; writes an html table at filename
#' @export
#'
#' @examples
MakeXtable<-function(model,             # a model object
                     rsq.method="pearson", # a method for calculating the r squared values
                     model.type,        # the type of model, basically just checks if its "maxnet"
                     regmult=1,         # if the model is a maxnet, then record the regularization multiplier
                     scale.factor=1,    # a sclae factor used to multiple the predictions
                     efh.break=NULL,    # a break point to be used for the accuracy table
                     nice.names=NULL,   # a table with variable names matched to nicer names for printing
                     filename,          # a filename to save the table
                     abund=NULL,        # a vector of abundance
                     train=NULL,        # a vector of training set predictions
                     test=NULL,         # a vector of testing or crossvalidation predictions
                     group=NULL,        # a vector of groups that apply to the other data columns
                     forecast=NULL,     # a set of forecasting results
                     ncoefs=NULL,       # a named vector of coefficients from the MaxnetStats function
                     devs=NULL,         # a named vector of deviance %s from one of the ___Stats functions
                     area=NULL){        # a single value for the area of the EFH

  dim.names1<-""
  dim.names2<-"Model term"
  if(length(devs)>0){if(is.na(devs)){devs<-NULL}}

  if(tolower(model.type)=="maxnet"){
    if(is.null(ncoefs)==F){
      dim.names1<-c(dim.names1,"")
      dim.names2<-c(dim.names2,"Coefficients")
      term.names<-names(ncoefs)
      nterms<-ncoefs
    }
    if(is.null(devs)==F){
      term.names<-names(devs)
      dev.col=length(dim.names1)+1
      dim.names1<-c(dim.names1,"")
      dim.names2<-c(dim.names2,"% Explained Deviance")
      term.order<-order(devs,decreasing = T)
      edf.order<-term.order
    }else{
      term.order<-1:length(model$samplemeans)
      edf.order<-term.order
    }
    if(is.null(ncoefs)==F & is.null(devs)==F){
      term.order<-order(devs,ncoefs,decreasing = T)
      edf.order<-term.order
    }

    dim.names1<-c(dim.names1,"")
    dim.names2<-c(dim.names2,"Total Deviance Explained")
    tot.dev<-format(signif(max(model$dev.ratio)*100,3),nsmall=1)
    totdev.col<-length(dim.names1)
  }else{
    # it'll be helpful to keep the names attached
    term.table<-AutodetectGAMTerms(model,hgam="d")
    term.table<-subset(term.table,type!="offset")

    dim.names1<-c(dim.names1,"")
    dim.names2<-c(dim.names2,"EDF")

    # a little work to make sure all the names and stats match up
    nterms <- round(summary(model)$s.table[1:nrow(term.table[term.table$type!="factor",]),1],1)
    edf.names<-vector(length=length(nterms))
    for(n in 1:length(nterms)){
      name0<-strsplit(names(nterms)[n],split="[()]")[[1]]
      name1<-strsplit(name0[2],split=",")[[1]]
      edf.names[n]<-paste(name1,collapse="*")
    }

    if(is.null(devs)==F){
      dev.col=length(dim.names1)+1
      dim.names1<-c(dim.names1,"")
      dim.names2<-c(dim.names2,"% Explained Deviance")
      term.order<-order(devs,decreasing = T)
      term.names<-names(devs)
    }else{
      term.order<-1:length(model$samplemeans)
      term.names<-term.table$term
    }
    dim.names1<-c(dim.names1,"")
    dim.names2<-c(dim.names2,"Total Deviance Explained")
    totdev.col<-length(dim.names1)

    tot.dev<-format(signif(summary(model)$dev.expl*100, 3),nsmall=1)

    #check that edf is in the same order as the devs
    edf.order<-match(edf.names,names(devs))
  }
  # set up the confusion matrix now, put it into a table later
  train.table<-NULL
  if(is.null(train)==F){
    train.col<-c(length(dim.names1)+1,length(dim.names1)+2)
    dim.names1<-c(dim.names1,"Training data","Training data")
    dim.names2<-c(dim.names2,"rsq","rmse")

    if(is.null(efh.break)==F){
      obs.efh<-as.integer(abund>efh.break)
      pred.efh<-as.integer(train>efh.break)

      n.nas<-sum(is.na(pred.efh))
      denom<-length(obs.efh)-n.nas

      t.pos<-round(sum(as.integer(obs.efh==1 & pred.efh==1),na.rm=T)/denom*100,1)
      t.neg<-round(sum(as.integer(obs.efh==0 & pred.efh==0),na.rm=T)/denom*100,1)
      f.pos<-round(sum(as.integer(obs.efh==0 & pred.efh==1),na.rm=T)/denom*100,1)
      f.neg<-round(sum(as.integer(obs.efh==1 & pred.efh==0),na.rm=T)/denom*100,1)

      train.table<-array(dim=c(7,3),data="")
      train.table[,1]<-c("Training Data","","Obs. EFH","Obs. non-EFH","","# Samples","Total Accuracy %")
      train.table[,2]<-c("","Pred. EFH",t.pos,f.pos,"",length(train),t.pos+t.neg)
      train.table[,3]<-c("","Pred. non-EFH",f.neg,t.neg,"","","")
    }
  }
  test.table<-NULL
  if(is.null(test)==F){
    test.col<-c(length(dim.names1)+1,length(dim.names1)+2)
    dim.names1<-c(dim.names1,"Test data","Test data")
    dim.names2<-c(dim.names2,"rsq","rmse")

    if(is.null(efh.break)==F){
      obs.efh<-as.integer(abund>efh.break)
      pred.efh<-as.integer(test>efh.break)

      n.nas<-sum(is.na(pred.efh))
      denom<-length(obs.efh)-n.nas

      t.pos<-round(sum(as.integer(obs.efh==1 & pred.efh==1),na.rm=T)/denom*100,1)
      t.neg<-round(sum(as.integer(obs.efh==0 & pred.efh==0),na.rm=T)/denom*100,1)
      f.pos<-round(sum(as.integer(obs.efh==0 & pred.efh==1),na.rm=T)/denom*100,1)
      f.neg<-round(sum(as.integer(obs.efh==1 & pred.efh==0),na.rm=T)/denom*100,1)

      test.table<-array(dim=c(7,3),data="")
      test.table[,1]<-c("CV Testing Data","","Obs. EFH","Obs. non-EFH","","# Samples","Total Accuracy %")
      test.table[,2]<-c("","Pred. EFH",t.pos,f.pos,"",length(test),t.pos+t.neg)
      test.table[,3]<-c("","Pred. non-EFH",f.neg,t.neg,"","","")
    }
  }
  if(is.null(forecast)==F){
    forecast.col<-c(length(dim.names1)+1,length(dim.names1)+2)
    dim.names1<-c(dim.names1,"Forecast data","Forecast data")
    dim.names2<-c(dim.names2,"rsq","rmse")
  }
  if(is.null(area)==F){
    area.col=length(dim.names1)+1
    dim.names1<-c(dim.names1,"")
    dim.names2<-c(dim.names2,"Area")
  }
  # Covers a weird edge case if the table is very small
  if(length(dim.names1)<7 & is.na(group)==F & is.null(group)==F){
    diff<-7-length(dim.names1)
    dim.names1<-c(dim.names1,rep("",diff))
    dim.names2<-c(dim.names2,rep("",diff))
  }
  # Fix the names
  if(is.null(nice.names)==F){
    term.names<-nice.names[match(term.names,nice.names[,1]),2]
  }

  # now we have a set of everything we need, so make the table and start filling it in
  table1 <- array("", dim = c(length(term.names)+2,length(dim.names1)))
  table1[1,] <- dim.names1
  table1[2,] <- dim.names2
  table1[3:(length(term.names)+2),1]<-term.names[term.order]
  table1[3,totdev.col]<-tot.dev

  if(is.null(ncoefs)==F | model.type!="maxnet"){
    table1[3:(length(term.names)+2),2]<-nterms[term.order]
    table1[is.na(table1[,2]),2]<-1
  }
  if(is.null(devs)==F){table1[3:(length(term.names)+2),dev.col]<-format(round(devs[term.order],1),nsmall = 1)}
  if(is.null(train)==F){
    train.dat<-na.omit(train)
    table1[3,train.col[1]]<-round(cor(abund,train,method=rsq.method)^2,3)
    table1[3,train.col[2]]<-format(round(RMSE(pred = train,obs = abund),1),nsmall=1,big.mark = ",")

  }
  if(is.null(test)==F){
    test.dat<-na.omit(test)
    table1[3,test.col[1]]<-round(cor(abund,test,method=rsq.method)^2,3)
    table1[3,test.col[2]]<-format(round(RMSE(pred = test,obs = abund),1),nsmall=1,big.mark = ",")
  }
  if(is.null(forecast)==F){
    forecast.dat<-na.omit(forecast)
    table1[3,forecast.col[1]]<-round(cor(abund,forecast,method=rsq.method)^2,3)
    table1[3,forecast.col[2]]<-format(round(RMSE(pred = forecast,obs = abund),1),nsmall=1,big.mark = ",")
  }
  if(is.null(area)==F){table1[3,area.col]<-format(round(area,-2),big.mark = ",")}

  if(is.null(model$family$family)==F && model$family$family=="quasipoisson"){
    theta.table<-array("",dim=c(2,length(dim.names1)))
    theta.table[1,1]<-"."
    theta.table[2,1]<-"Theta"
    theta.table[2,2]<-format(round(model$scale,1),nsmall=1)
    table1<-rbind(table1,theta.table)
  }
  if(is.null(model$family$family)==F && tolower(strsplit(model$family$family,split="[()]")[[1]][1])=="negative binomial"){
    theta.table<-array("",dim=c(2,length(dim.names1)))
    theta.table[1,1]<-"."
    theta.table[2,1]<-"Theta"
    theta.table[2,2]<-round(as.numeric(strsplit(model$family$family,split="[()]")[[1]][2]),3)
    table1<-rbind(table1,theta.table)
  }

  if(model.type=="maxnet" & is.numeric(regmult)){
    reg.table<-array("",dim=c(2,length(dim.names1)))
    reg.table[1,1]<-"."
    reg.table[2,1:2]<-c("Regularization Constant",format(round(regmult,1),nsmall=1))
    table1<-rbind(table1,reg.table)
  }

  if(is.numeric(scale.factor) & scale.factor!=1){
    scale.table<-array("",dim=c(2,length(dim.names1)))
    scale.table[1,1]<-"."
    scale.table[2,1:2]<-c("Scale Factor",format(round(scale.factor,2),nsmall=1))
    table1<-rbind(table1,scale.table)
  }

  # make a spacer so the table is formatted nice and readable, if necessary
  spacer<-NULL
  if(is.na(group)==F & is.null(group)==F | is.null(train)==F |is.null(test)==F){
    spacer<-array("",dim=c(2,length(dim.names1)))
    spacer[,1]<-"."
  }


  # now go and add in a table for the confusion matrices
  acc.table<-NULL
  if(is.null(train)==F |is.null(test)==F & is.null(efh.break)==F){
    acc.table<-array(data="",dim=c(8,length(dim.names1)))
    acc.table[1,1]<-"Confusion Matrix"
    if(is.null(train)==F){
      acc.table[2:8,1:3]<-train.table
      if(is.null(test)==F){
        acc.table[2:8,5:7]<-test.table
      }
    }else{
      acc.table[2:8,1:3]<-test.table
    }
  }

  cv.table<-NULL
  #now, add in a table for the cross-validation info
  if(is.na(group)==F & is.null(group)==F){

    cv.rows<-length(unique(group))+3
    cv.table<-array("",dim = c(cv.rows,length(dim.names1)))
    cv.table[1,1]<-"Cross Validation Results"
    cv.table[2,1:6]<-c("Fold","N","Pres","rsq","mean error","rmse")
    cv.table[3:cv.rows,1]<-c(unique(group),"Totals/Averages")
    cv.table[3:cv.rows,2]<-c(table(group),sum(table(group)))
    cv.table[3:cv.rows,3]<-c(aggregate(x = abund,by=list(group),FUN=function(x){sum(x>0)})[,2],
                             sum(abund>0))
    cv.table[3:cv.rows,5]<-format(signif(c(aggregate(x=test-abund,by=list(group),FUN="mean",na.rm=T)$x,mean(test-abund,na.rm=T)),3),nsmall=2)
    cv.table[cv.rows,4]<-format(round(cor(abund,test,method=rsq.method)^2, 3),nsmall=3)
    cv.table[cv.rows,6]<-format(round(RMSE(pred = test,obs = abund),2),nsmall=1,big.mark = "'")
    # now need a loop
    for(g in 1:length(unique(group))){
      gabund<-abund[group==unique(group)[g]]
      gtest<-test[group==unique(group)[g]]
      if(length(gtest)>0){
        cv.table[2+g,4]<-format(round(cor(gabund,gtest,method=rsq.method)^2, 3),nsmall=3)
        cv.table[2+g,6]<-format(round(RMSE(pred = gtest,obs = gabund),2),nsmall=1,big.mark = ",")
      }else{
        cv.table[2+g,4]<-"NA"
        cv.table[2+g,6]<-"NA"
      }
    }
  }
  model.table <- xtable::xtable(rbind(table1,spacer,acc.table,spacer,cv.table))

  xtable::print.xtable(model.table, type = "html", file = filename,
                       include.rownames = getOption("xtable.include.rownames", FALSE), html.table.attributes = 2,
                       include.colnames = getOption("xtable.include.colnames", FALSE), hline.after = getOption("xtable.hline.after",
                                                                                                               c(-1, 1, nrow(table1))))
}



#' Make a table of summary statistics (for ensemble)
#'
#' @description Make an HTML table of summary stats and fit metrics for a model ensemble.
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
        pred.dat<-na.omit(subset(preds.table,Model==model.names[m]))
        pred.dat$pred[which(is.infinite(pred.dat$pred))]<-10^10
        pred.cvpred<-ifelse(RMSE(pred = pred.dat$pred,obs = pred.dat$abund)<
                              RMSE(pred = pred.dat$cvpred,obs = pred.dat$abund),"cvpred","pred")

        prob.cvprob<-ifelse(pred.cvpred=="cvpred","cvprob","prob")

        N.vec[m]<-sum(pred.dat$abund>0,na.rm=T)
        rmse.vec[m]<-RMSE(pred = pred.dat[,pred.cvpred],obs = pred.dat$abund)
        cor.vec[m]<-cor(pred.dat$abund,round(pred.dat[,pred.cvpred],2),method=cor.method)
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
      ensemble.dat<-na.omit(subset(preds.table,Model=="ensemble"))

      etable[n.models+3,N.col]<-sum(ensemble.dat$abund>0,na.rm=T)
      etable[n.models+3,pred.cols[1]]<-format(round(RMSE(obs=ensemble.dat$abund,pred=ensemble.dat$pred),digs),nsmall=digs,big.mark = ",")
      etable[n.models+3,pred.cols[2]]<-format(round(cor(ensemble.preds$abund,ensemble.preds$pred,method=cor.method), 2),nsmall=2)
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




#' Make a deviance table
#'
#' @description Make an HTML table of deviance explained for a model ensemble.
#' @details Makes a simple html table of deviance explained values for models and the ensemble.
#' @param model.names character vector with names for the models, should match order of dev.list
#' @param model.types character vector with types of model, should match order of dev.list
#' @param nice.names dataframe matching abbreviated names to nicer versions
#' @param dev.list a list of named vectors giving the deviance explained values for each model
#' @param model.weights numeric vector with weights for each model, should match order of dev.list
#' @param filename a filename to use for saving the output
#'
#' @return does not return anything, but writes two tables at based on the supplied filename
#' @export
#'
#' @examples
MakeDevianceTable<-function(model.names=NULL,
                            model.types=NULL,
                            nice.names=NULL,
                            dev.list,
                            model.weights=NULL,
                            filename){



  if(is.null(model.weights)||is.na(model.weights)){
    good.models<-1:length(dev.list)
  }else{
    good.models<-which(model.weights>0)
  }
  if(is.null(model.names)||is.na(model.names)){
    model.names<-paste0("Model ",1:length(dev.list))
  }
  good.names<-model.names[good.models]

  columns<-"Model"
  if(is.null(model.types)==F&&is.na(model.types)==F){
    t.col<-2
    columns<-c(columns,"Type")
  }
  if(is.null(model.weights)==F &&is.na(model.weights)==F){
    w.col<-length(columns)+1
    columns<-c(columns,"Weight")
  }

  dev.names<-unique(names(unlist(dev.list)))
  dev.names<-dev.names[dev.names!=""]
  n.terms<-length(dev.names)

  # if names aren't supplied, use the terms from the model
  if(is.null(nice.names)){
    nice.names<-data.frame(var=dev.names,name=dev.names)
  }

  model.rows<-2:(length(good.models)+1)
  term.col<-(length(columns)+1):(length(columns)+n.terms)

  columns<-c(columns,nice.names$name[match(dev.names,nice.names$var)])

  dtable<-array(dim=c(length(good.names)+2,length(columns)))

  dtable[1,]<-c(columns)

  for(d in 1:length(good.models)){
    dtable[1+d,1]<-good.names[d]

    if("Type"%in%columns){
      dtable[1+d,2]<-model.types[good.models[d]]
    }
    if("Weight"%in%columns){
      dtable[1+d,w.col]<-round(model.weights[good.models[d]],3)
    }
    model.dev<-dev.list[[good.models[d]]]
    for(t in 1:n.terms){
      term<-dev.names[t]
      d.term<-which(names(model.dev)==term)
      if(length(d.term>0)){
        dtable[(1+d),term.col[t]]<-round(model.dev[d.term],1)
      }else{
        dtable[(1+d),term.col[t]]<-"--"
      }
      if(is.na(dtable[(1+d),term.col[t]])){dtable[(1+d),term.col[t]]<-"--"}
    }
  }

  dtable[nrow(dtable),1]<-"Mean"

  mean.matrix<-matrix(nrow=length(model.rows),ncol=length(term.col),data = as.numeric(dtable[model.rows,term.col]),byrow = F)

  term.means<-round(apply(mean.matrix,MARGIN = 2,FUN = mean,na.rm=T),1)

  dtable[nrow(dtable),term.col]<-term.means
  term.order<-order(term.means,decreasing = T)

  dtable2<-dtable
  dtable2[,term.col]<-dtable[,term.col[term.order]]

  if("Weight"%in%columns){
    w.mean<-array(dim=c(1,ncol(dtable2)))

    term.means2<-round(apply(mean.matrix*model.weights[good.models],MARGIN = 2,FUN = sum,na.rm=T),1)

    term.order<-order(term.means2,decreasing = T)

    dtable2<-dtable
    dtable2[,term.col]<-dtable[,term.col[term.order]]

    w.mean[1,1]<-"Weighted Mean"
    w.mean[,term.col]<-term.means2[term.order]

    dtable2<-rbind(dtable2,w.mean)
  }

  deviance.table <- xtable::xtable(dtable2)
  xtable::print.xtable(deviance.table, type = "html", file = filename,
                       include.rownames = getOption("xtable.include.rownames", FALSE), html.table.attributes = 2,
                       include.colnames = getOption("xtable.include.colnames", FALSE), hline.after = getOption("xtable.hline.after",
                                                                                                               c(-1, 1, nrow(deviance.table))))

  # Make a second table

  rescale.col2<-as.numeric(dtable2[nrow(dtable2),4:ncol(dtable2)])
  rescale.col3<-cumsum(dtable2[nrow(dtable2),4:ncol(dtable2)])

  dtable3<-array(dim=c(1+length(term.means),3))
  dtable3[,1]<-c("Covariate",dtable2[1,4:ncol(dtable2)])
  dtable3[,2]<-c("% Contribution",format(round(100*rescale.col2/sum(rescale.col2),1),nsmall=1))
  dtable3[,3]<-c("Cumulative %",format(round(100*rescale.col3/sum(rescale.col2),1),nsmall=1))

  deviance.table2 <- xtable::xtable(dtable3)
  xtable::print.xtable(deviance.table2, type = "html", file = paste0(strsplit(filename,split = ".html")[[1]],"2.html"),
                       include.rownames = getOption("xtable.include.rownames", FALSE), html.table.attributes = 2,
                       include.colnames = getOption("xtable.include.colnames", FALSE), hline.after = getOption("xtable.hline.after",
                                                                                                               c(-1, 1, nrow(deviance.table2))))
}
