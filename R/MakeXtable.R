#' Make deviance explained table with constituent model statistics
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
#' @param train vector of predictions from the model trained on the full data set
#' @param test vector of out-of-bag predictions
#' @param group vector of group labels
#' @param forecast vector of predictions from a third group (not used much anymore)
#' @param ncoefs vector of coefficients used in maxnet models
#' @param devs vector of deviance explained estimates
#' @param area numeric; total EFH area
#'
#' @return no return value; writes an html table at filename
#' @export
#' @importFrom stats aggregate
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

  # fix user input
  if(length(devs)==1 && is.na(devs[1])){devs<-NULL}
  if(length(group)==1 && is.na(group)){group<-NULL}

  dim.names1<-""
  dim.names2<-"Model term"

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
  if(length(dim.names1)<7 & is.null(group)==F){
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
    train.dat<-stats::na.omit(train)
    table1[3,train.col[1]]<-round(stats::cor(abund,train,method=rsq.method)^2,3)
    table1[3,train.col[2]]<-format(round(RMSE(pred = train,obs = abund),1),nsmall=1,big.mark = ",")

  }
  if(is.null(test)==F){
    test.dat<-stats::na.omit(test)
    table1[3,test.col[1]]<-round(stats::cor(abund,test,method=rsq.method)^2,3)
    table1[3,test.col[2]]<-format(round(RMSE(pred = test,obs = abund),1),nsmall=1,big.mark = ",")
  }
  if(is.null(forecast)==F){
    forecast.dat<-stats::na.omit(forecast)
    table1[3,forecast.col[1]]<-round(stats::cor(abund,forecast,method=rsq.method)^2,3)
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
  if(is.null(group)==F | is.null(train)==F |is.null(test)==F){
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
  if(is.null(group)==F){

    cv.rows<-length(unique(group))+3
    cv.table<-array("",dim = c(cv.rows,length(dim.names1)))
    cv.table[1,1]<-"Cross Validation Results"
    cv.table[2,1:6]<-c("Fold","N","Pres","rsq","mean error","rmse")
    cv.table[3:cv.rows,1]<-c(unique(group),"Totals/Averages")
    cv.table[3:cv.rows,2]<-c(table(group),sum(table(group)))
    cv.table[3:cv.rows,3]<-c(stats::aggregate(x = abund,by=list(group),FUN=function(x){sum(x>0)})[,2],
                             sum(abund>0))
    cv.table[3:cv.rows,5]<-format(signif(c(aggregate(x=test-abund,by=list(group),FUN="mean",na.rm=T)$x,mean(test-abund,na.rm=T)),3),nsmall=2)
    cv.table[cv.rows,4]<-format(round(stats::cor(abund,test,method=rsq.method)^2, 3),nsmall=3)
    cv.table[cv.rows,6]<-format(round(RMSE(pred = test,obs = abund),2),nsmall=1,big.mark = "'")
    # now need a loop
    for(g in 1:length(unique(group))){
      gabund<-abund[group==unique(group)[g]]
      gtest<-test[group==unique(group)[g]]
      if(length(gtest)>0){
        cv.table[2+g,4]<-format(round(stats::cor(gabund,gtest,method=rsq.method)^2, 3),nsmall=3)
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
