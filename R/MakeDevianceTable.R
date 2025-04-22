#' Make deviance table
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
