# # This is going to be a new, complete version that runs all models, including any recent fixes such as the maxnet tuning routine.
# # Will also incorporate the publication figures.
#
# source("G:/Harris/EFH_copy/Functions_LoadMap.R")
# source("G:/Harris/EFH_copy/Functions_GamModel.R")
# source("G:/Harris/EFH_copy/Functions_Maxent.R")
# source("G:/Harris/EFH_copy/Functions_Xtable.R")
# source("G:/Harris/EFH_copy/Functions_Ensemble.R")
#
# computer.name<-"JHlaptop"
# region<-F
#
# # Use this to specify a particular species or subset of species and some other control parameters
# # set species.vec and region.vec to NA in order to run everything
# species.vec<-NA                              # use the column abbreviations to select specific species
# region.vec<-NA
# stop.early<-F                                                  # useful if you want to stop and check results
# update.table<-T
# rerun<-F
#
# EFH.path<-"Y:/RACE_EFH_variables"
#
#
# # set of nicer names to be used in figure labels
# # Use name for most figures, but use name2 for deviance tables
# nice.names<-data.frame(stringsAsFactors = F,
#                        var=c("lon","lat","bdepth","slope","aspectE","aspectN","curve","btemp","speed","tmax","BPI","phi","vitalrate",
#                              "color","sponge","coral","pen","area","bcurrentU","bcurrentV","bcurrentUSD","bcurrentVSD","lon*lat",
#                              "bcurrentU*bcurrentV","bcurrentUSD*bcurrentVSD","aspectE*aspectN","rocky"),
#                        name=c("Longitude","Latitude","Depth (m)","Slope (degrees)","Slope Eastness","Slope Northness",
#                               "Terrain Curvature","Temperature (C)","Current speed (m/s)",
#                               "Tidal Maximum (cm/s)","BPI","Sediment Grain size (phi)","Growth Potential (g/g/day)",
#                               "Ocean color","Sponge presence","Coral presence","Pennatulacean presence","Area Swept",
#                               "Current Velocity East (m/s)","Current Velocity North (m/s)","Current Velocity East SD (m/s)",
#                               "Current Velocity North SD (m/s)","Position (lon,lat)","Current Velocity (east,north)",
#                               "Current Velocity SD","Slope Aspect (east,north)","Bottom Rockiness (%)"))
#
# nice.names2<-data.frame(stringsAsFactors = F,
#                         var=c("lon","lat","bdepth","slope","aspectE","aspectN","curve","btemp","speed","tmax","BPI","phi","vitalrate",
#                               "color","sponge","coral","pen","area","bcurrentU","bcurrentV","bcurrentUSD","bcurrentVSD","lon*lat",
#                               "bcurrentU*bcurrentV","bcurrentUSD*bcurrentVSD","aspectE*aspectN","rocky"),
#                         name=c("longitude","latitude","bottom depth","slope","aspect east","aspect north","curvature",
#                                "bottom temperature","current speed",
#                                "tidal maximum","BPI","phi","growth potential",
#                                "color","sponge presence","coral presence","pennatulacean presence","area swept",
#                                "current east","current north","current east SD",
#                                "current north SD","position","current",
#                                "current SD","slope aspect","rockiness"))
#
# #Do you want to group by spatial zone or random fold
# group.var<-"Folds"
#
# #make ensemble versions
# make.ensemble<-T
#
# done<-F
#
# while(done==F){
#   # now, going to figure out which species to do next
#   progress.table<-read.csv(file=paste0(EFH.path,"/ModelProgress_offset_reruns.csv"),stringsAsFactors = F)
#
#   if(is.na(species.vec)){species.vec<-unique(progress.table$Abbreviation)}
#   if(is.na(region.vec)){region.vec<-c("AI","EBS","GOA")}
#
#   unclaimed<-which(is.na(progress.table$Claimed) & progress.table$Region%in%region.vec & progress.table$Abbreviation%in%species.vec)
#
#   if(length(unclaimed)==0 & rerun==F){
#     print("No unclaimed models; ending routine!!!")
#     done<-T
#   }else{
#
#     # alternately, you can just set a value for i here to run a single species
#     i<-unclaimed[1]
#     if(rerun){
#       if(length(species.vec)>1 | length(region.vec)>1){stop("Rerun option can only accept one species at a time")}
#       i=which(progress.table$Abbreviation==species.vec & progress.table$Region==region.vec)
#       done<-T
#     }
#     progress.table$Claimed[i]<-as.character(Sys.time())
#     progress.table$Computer[i]<-computer.name
#     if(update.table){
#       try(write.csv(x = progress.table,file =paste0(EFH.path,"/ModelProgress_offset_reruns.csv"),row.names = F))
#       try(write.csv(x = progress.table,file =paste0(EFH.path,"/ProgressChecker.csv"),row.names = F))
#     }
#     region0<-progress.table$Region[i]
#     s<-progress.table$Abbreviation[i]
#     figure.name<-progress.table$Figure_Name[i]
#     styear<-progress.table$Start_Year[i]
#
#     #behold the main loop
#     print(paste0("Starting species: ",figure.name,"; Number ",i," of ",nrow(progress.table)))
#     #Check data and load if necessary
#     if(region0!=region){
#       region<-region0
#
#       LoadEFHData(region,raster.path=paste0(EFH.path,"/Variables"))
#       LoadMap(region,parameter.file=paste0(EFH.path,"/Map_Settings.csv"),covariate.raster=raster.stack,
#               coast.file=paste0(EFH.path,"/shapefiles"),
#               GOA.mask = paste0(EFH.path,"/shapefiles"))
#
#       if(region=="GOA"){raster.stack<-mask(raster.stack,GOA.mask)}
#
#       region.data<-read.csv(paste0(EFH.path,"/Trawl_Models2/",region,"/all_",region,"_data_2021.csv"))
#
#       # these are sizes for png files for various figures
#       png.width<-6
#       png.height<-ifelse(region=="EBS",6,5)
#       legend.pos<-ifelse(region=="EBS","topleft",ifelse(region=="GOA","bottom","top"))
#       legend.horiz<-ifelse(region=="EBS",F,T)
#
#       # small adjustment for the REBS complex
#       region.data$j_rebs<-region.data$j_rebs+region.data$j_rough+region.data$j_bspot
#       region.data$a_rebs<-region.data$a_rebs+region.data$a_rough+region.data$a_bspot
#     }
#
#     species.data0<-subset(region.data,year>=styear)
#
#     if(is.na(progress.table$Lifestage[i])){
#       species.path<-paste0(EFH.path,"/Trawl_Models2/",region,"/",progress.table$Species[i])
#     }else{
#       species.path<-paste0(EFH.path,"/Trawl_Models2/",region,"/",progress.table$Lifestage[i],"_",progress.table$Species[i])
#     }
#
#     # Assign the folds from here
#     seed<-progress.table$Species_Code[i]
#     set.seed(seed)
#
#     random.folds<-rep(LETTERS[1:10],length.out=nrow(species.data0))
#     species.data0$Folds<-sample(random.folds,size = nrow(species.data0),replace = F)
#
#     n.pres<-sum(species.data0[,s]>0)
#
#     if(n.pres<50){
#       print(paste("WARNING: Insufficient samples for",region,figure.name,". Proceeding to next species/stage",sep = " "))
#       species.data<-species.data0
#       maxnet.abund.check<-F
#       cloglog.abund.check<-F
#       hpoisson.abund.check<-F
#       poisson.abund.check<-F
#       negbin.abund.check<-F
#       ensemble.success<-F
#       start.time<-Sys.time()
#     }
#     if(n.pres>=50){
#       # if there are at least 50 total hauls with presence, we will at least run some kind of model
#       pres.table<-table(species.data0[species.data0[,s]>0,group.var])
#       smallest<-min(n.pres-pres.table)
#
#       #New reallocation loop, basically, if one fold has too many of the presences, then the largest fold
#       # donates a few cases to the smallest fold, loop continues until model always has at least 40 to train with
#       if(smallest<40){
#         species.data01<-species.data0
#         while(smallest<40){
#           least.zone<-names(pres.table)[which.min(pres.table)]
#           most.zone<-names(pres.table)[which.max(pres.table)]
#
#           n.donate<-40-smallest
#           donors<-which(species.data01[,group.var]==most.zone & species.data01[,s]>0)
#           winners<-sample(x = donors,size = n.donate,replace = F)
#           species.data01[winners,group.var]<-least.zone
#
#           pres.table<-table(species.data0[species.data0[,s]>0,group.var])
#           smallest<-min(n.pres-pres.table)
#         }
#         species.data<-species.data01
#       }else{
#         species.data<-species.data0
#       }
#
#       # FOr now, we will leave the weights of the SFIs in the data set, but we are using it as a binary variable for the analysis
#       species.data$sponge<-as.factor(as.integer(species.data$sponge>0))
#       species.data$coral<-as.factor(as.integer(species.data$coral>0))
#       species.data$pen<-as.factor(as.integer(species.data$pen>0))
#
#       # also remember to log the area
#       species.data$logarea<-log(species.data$area)
#
#       covars2d<-list(c("lon","lat"),c("bcurrentU","bcurrentV"),c("bcurrentUSD","bcurrentVSD"))
#       if(region=="EBS"){
#         covars<-covars.vec[covars.vec%in%c("bdepth","slope","aspectE","aspectN","curve","btemp","tmax","phi","BPI")]
#         maxnet.covars<-c("bcurrentU","bcurrentV","bcurrentUSD","bcurrentVSD","bdepth",
#                          "slope","aspectE","aspectN","curve","btemp","tmax","phi","BPI")
#       }else{
#         covars<-covars.vec[covars.vec%in%c("bdepth","slope","aspectE","aspectN","curve","btemp","tmax","rocky","BPI")]
#         maxnet.covars<-c("bcurrentU","bcurrentV","bcurrentUSD","bcurrentVSD","bdepth",
#                          "slope","aspectE","aspectN","curve","btemp","tmax","rocky","BPI")
#       }
#       cofactors<-c("sponge","coral","pen")
#       maxnet2d<-list(c("bcurrentU","bcurrentV"),c("bcurrentUSD","bcurrentVSD"))
#
#       # Now can start running the models
#       start.time<-Sys.time()
#       dir.create(species.path)
#
#       png(paste0(species.path,"/dotplot.png"),width = png.width,height = png.height,res=300,units = "in")
#       plotDots(train.data = species.data,species = s,legPosition = legend.pos,horiz = legend.horiz,
#                figure.name = figure.name)
#       dev.off()
#
#       # this looks cumbersome, but automates the gam formulas
#       basic.gam.table<-data.frame(type=c(rep("smooth",length(covars2d)+length(covars)),rep("factor",length(cofactors)),"offset"),
#                                   dims=c(rep(2,length(covars2d)),rep(1,length(covars)+length(cofactors)+1)),
#                                   term=c(unlist(lapply(covars2d,FUN=function(x){return(x[1])})),covars,cofactors,"logarea"),
#                                   term2=c(unlist(lapply(covars2d,FUN=function(x){return(x[2])})),rep(NA,length(covars)+length(cofactors)+1)),
#                                   bs=c(rep("ds",length(covars2d)),rep("tp",length(covars)),rep(NA,length(cofactors)+1)),
#                                   k=c(rep(10,length(covars2d)),rep(4,length(covars)),rep(NA,length(cofactors)+1)),
#                                   m=c(rep(1,length(covars2d)+length(covars)),rep(NA,length(cofactors)+1)),
#                                   m2=c(rep(.5,length(covars2d)),rep(NA,length(covars)+length(cofactors)+1)))
#
#       alt.gam.table<-basic.gam.table
#       alt.gam.table$bs[alt.gam.table$bs=="tp"]<-"cr"
#
#       basic.gam.formula<-AssembleGAMFormula(yvar = s, gam.table = basic.gam.table)
#       alt.gam.formula<-AssembleGAMFormula(yvar = s, gam.table = alt.gam.table)
#       basic.hgam.formula<-AssembleGAMFormula(yvar = s, gam.table = basic.gam.table,hgam=T)
#       alt.hgam.formula<-AssembleGAMFormula(yvar = s, gam.table = alt.gam.table,hgam=T)
#
#       ##############################
#       #Maxnet
#       print(paste(Sys.time(),"Starting the Maxnet model for",region,figure.name,sep=" "))
#       maxnet.converge0<-rep(F,6)
#       maxnet.abund.check0<-rep(F,6)
#       maxnet.model.list<-list()
#       maxnet.abund.list<-list()
#       maxnet.cv.model.list<-list()
#       maxnet.error.list<-list()
#       maxnet.pred.list<-list()
#       maxnet.scale.vec<-rep(NA,6)
#       maxnet.rmse0<-rep(NA,6)
#
#       # loop through and check different multiplication constants
#       for(r in 1:6){
#         r.mult<-c(.5,1,1.5,2,2.5,3)[r]
#
#         # check out a prospective maxnet model
#         try(maxnet.model0<-FitMaxnet(data = species.data,species = s,vars = maxnet.covars,facs=cofactors,
#                                      regmult = r.mult,reduce = T))
#
#         if(exists("maxnet.model0")){
#           maxnet.converge0[r]<-T
#
#           maxnet.scale.vec[r]<-mean(species.data[,s])/
#             mean(exp(predict(maxnet.model0,newdata=species.data,type="link")+maxnet.model0$ent))
#           maxnet.abund.list[[r]]<-MakeMaxEntAbundance(model = maxnet.model0,maxent.stack = raster.stack,
#                                                       scale.fac = maxnet.scale.vec[r],type = "cloglog",
#                                                       land = ak.raster,filename = "")
#
#           maxnet.abund.check0[r]<-cellStats(maxnet.abund.list[[r]],max)<(max(species.data[,s])*10)
#         }
#
#         # if it converges, do the checks and crossvalidation
#         if(maxnet.converge0[r] & maxnet.abund.check0[r]){
#           maxnet.model.list[[r]]<-maxnet.model0
#
#           maxnet.cv<-CrossValidateModel(model = maxnet.model0,data = species.data,species = s,group = group.var,
#                                         model.type = "maxnet",key="hauljoin",scale.preds = T,regmult = r.mult)
#           maxnet.error.list[[r]]<-maxnet.cv[[1]]
#           maxnet.cv.model.list[[r]]<-maxnet.cv[[2]]
#
#           maxnet.rmse0[r]<-max(c(RMSE(pred = maxnet.error.list[[r]]$pred,obs = maxnet.error.list[[r]]$abund),
#                                  RMSE(pred = maxnet.error.list[[r]]$cvpred,obs = maxnet.error.list[[r]]$abund)))
#
#           # will also discard a model any of the cv folds fails, as it becomes impossible to calculate RMSE
#           if(any(is.na(maxnet.cv[[1]]$cvpred))){
#             maxnet.model.list[[r]]<-NA
#             maxnet.cv.model.list[[r]]<-NA
#             maxnet.error.list[[r]]<-NA
#             maxnet.pred.list[[r]]<-NA
#             maxnet.rmse0[r]<-NA
#           }
#         }else{
#           maxnet.model.list[[r]]<-NA
#           maxnet.cv.model.list[[r]]<-NA
#           maxnet.error.list[[r]]<-NA
#           maxnet.pred.list[[r]]<-NA
#           maxnet.rmse0[r]<-NA
#         }
#         rm(maxnet.model0)
#       }
#
#       # as long as some of the constants managed to converge, carry something forward
#       if(sum(is.na(maxnet.rmse0))!=6){
#
#         maxnet.model<-maxnet.model.list[[which.min(maxnet.rmse0)]]
#         maxnet.abund<-maxnet.abund.list[[which.min(maxnet.rmse0)]]
#
#         writeRaster(maxnet.abund,filename = paste0(species.path,"/maxnet_abundance"),overwrite=T)
#
#         maxnet.cv.models<-maxnet.cv.model.list[[which.min(maxnet.rmse0)]]
#         maxnet.errors<-maxnet.error.list[[which.min(maxnet.rmse0)]]
#         maxnet.scale<-maxnet.scale.vec[which.min(maxnet.rmse0)]
#         opt.mult<-c(.5,1,1.5,2,2.5,3)[which.min(maxnet.rmse0)]
#
#         maxnet.rmse<-min(maxnet.rmse0,na.rm = T)
#
#         maxnet.effects<-GetMaxnetEffects(model = maxnet.model,cv.models = maxnet.cv.models,
#                                          vars = "all",add.entropy = T,scale = "log",data=species.data,
#                                          maxnet2d = maxnet2d,scale.factor = maxnet.scale)
#         maxnet.converge<-T
#         maxnet.abund.check<-T
#
#         saveRDS(maxnet.model,file=paste0(species.path,"/maxnet_model.rds"))
#         saveRDS(maxnet.effects,file=paste0(species.path,"/maxnet_effects.rds"))
#       }else{
#         print(paste("Maxnet model for",s,"unsuccessful -- moving on",sep=" "))
#         maxnet.converge<-sum(maxnet.converge0)>0
#         maxnet.abund.check<-sum(maxnet.abund.check0)>0
#       }
#
#       ##############################
#       #cloglog PAGAM
#       print(paste(Sys.time(),"Starting the Cloglog GAM for",region,figure.name,sep=" "))
#       try(cloglog.model<-FitGAM(gam.formula = basic.gam.formula,data = species.data,
#                                 family.gam = "binomial",link.fx = "cloglog",reduce = T,select = T,
#                                 verbose=F))
#
#       cloglog.converge<-exists("cloglog.model") & any(is.infinite(predict(cloglog.model,type="response")))==F &
#         any(is.na(predict(cloglog.model,type="response")))==F
#       if(cloglog.converge){
#         cloglog.converge<-T
#         cloglog.scale<-mean(species.data[,s])/mean(exp(predict(cloglog.model,type="link")))
#         cloglog.abund<-MakeGAMAbundance(model = cloglog.model,r.stack = raster.stack,scale.factor = cloglog.scale,
#                                         land = ak.raster,filename = "")
#
#         cloglog.abund.check<-cellStats(cloglog.abund,max)<(max(species.data[,s])*10)
#       }else{
#         cloglog.abund.check<-F
#       }
#
#       if(cloglog.abund.check==F){
#         rm(cloglog.abund.check,cloglog.abund,cloglog.scale)
#
#         print("TPS cloglog model failed abundance test; Trying alternate version")
#         try(cloglog.model<-FitGAM(gam.formula = alt.gam.formula,data = species.data,species = s,
#                                   family.gam = "binomial",link.fx = "cloglog",reduce = T,
#                                   select = T,verbose=F))
#         cloglog.converge<-exists("cloglog.model") & any(is.infinite(predict(cloglog.model,type="response")))==F &
#           any(is.na(predict(cloglog.model,type="response")))==F
#         if(cloglog.converge){
#           cloglog.scale<-mean(species.data[,s])/mean(exp(predict(cloglog.model,type="link")))
#           cloglog.abund<-MakeGAMAbundance(model = cloglog.model,r.stack = raster.stack,scale.factor = cloglog.scale,
#                                           land = ak.raster,filename = "")
#           cloglog.abund.check<-cellStats(cloglog.abund,max)<(max(species.data[,s])*10)
#
#         }else{
#           cloglog.converge<-F
#         }
#       }
#
#       if(cloglog.abund.check==F){
#         print("both versions of cloglog model unacceptable; moving to hgam")
#         rm(cloglog.model,cloglog.abund,cloglog.scale)
#       }
#
#       if(cloglog.converge & cloglog.abund.check){
#         print(paste0(Sys.time(),"Cloglog model successful for ",figure.name," starting evaluations"))
#         cloglog.cv<-CrossValidateModel(model = cloglog.model,model.type = "cloglog",data = species.data,scale.preds = T,
#                                        species = s,group = group.var,key="hauljoin")
#         cloglog.errors<-cloglog.cv[[1]]
#         cloglog.cv.models<-cloglog.cv[[2]]
#
#         writeRaster(cloglog.abund,paste0(species.path,"/cloglog_abundance"),overwrite=T)
#
#         cloglog.effects<-GetGAMEffects(model = cloglog.model,data = species.data,vars = "all",scale = "log",
#                                        cv.model.list = cloglog.cv.models,scale.factor = cloglog.scale)
#
#         cloglog.rmse<-max(c(RMSE(pred = cloglog.errors$pred,obs = cloglog.errors$abund),
#                             RMSE(pred = cloglog.errors$cvpred,obs = cloglog.errors$abund)))
#
#         cloglog.success<-ifelse(any(is.na(cloglog.errors$cvpred)),F,T)
#
#         saveRDS(cloglog.model,file=paste0(species.path,"/cloglog_model.rds"))
#         saveRDS(cloglog.effects,file=paste0(species.path,"/cloglog_effects.rds"))
#       }else{
#         print(paste("Cloglog gam for",s,"unsuccessful -- moving on",sep=" "))
#       }
#
#       ##############################
#       #HurdleGAM time
#       print(paste(Sys.time(),"Starting the Hurdle model for",region,figure.name,sep=" "))
#       try(hpoisson.model<-FitHurdleGAM(density.formula = basic.hgam.formula[[1]],prob.formula = basic.hgam.formula[[2]],
#                                        data = species.data,verbose = F,select = T,reduce = T))
#
#       hpoisson.converge<-exists("hpoisson.model") & any(is.infinite(predict(hpoisson.model,type="response")))==F &
#         any(is.na(predict(hpoisson.model,type="response")))==F
#
#       if(hpoisson.converge){
#         hpoisson.scale<-mean(species.data[,s])/mean(predict(hpoisson.model,type="response"))
#         hpoisson.abund<-MakeGAMAbundance(model = hpoisson.model,r.stack = raster.stack,scale.factor = hpoisson.scale,
#                                          land = ak.raster,filename = "")
#
#         hpoisson.abund.check<-cellStats(hpoisson.abund,max)<(max(species.data[,s])*10)
#       }else{
#         hpoisson.abund.check<-F
#       }
#
#       if(hpoisson.abund.check==F){
#         rm(hpoisson.model,hpoisson.abund,hpoisson.scale)
#
#         print("TPS hurdle model failed abundance test; Trying alternate version")
#         try(hpoisson.model<-FitHurdleGAM(density.formula = alt.hgam.formula[[1]],prob.formula = alt.hgam.formula[[2]],
#                                          data = species.data,verbose = F,select = T,reduce = T))
#
#         hpoisson.converge<-exists("hpoisson.model") & any(is.infinite(predict(hpoisson.model,type="response")))==F &
#           any(is.na(predict(hpoisson.model,type="response")))==F
#         if(hpoisson.converge){
#           hpoisson.scale<-mean(species.data[,s])/mean(predict(hpoisson.model,type="response"))
#           hpoisson.abund<-MakeGAMAbundance(model = hpoisson.model,r.stack = raster.stack,scale.factor = hpoisson.scale,
#                                            land = ak.raster,filename = "")
#           hpoisson.abund.check<-cellStats(hpoisson.abund,max)<(max(species.data[,s])*10)
#
#         }
#       }
#       if(hpoisson.abund.check==F){
#         print("both versions of hpoisson model unacceptable; moving to gams")
#         rm(hpoisson.model,hpoisson.abund,hpoisson.scale)
#       }
#
#       if(hpoisson.converge & hpoisson.abund.check){
#         print(paste0(Sys.time(),"hpoisson model successful for ",figure.name," starting evaluations"))
#         hpoisson.cv<-CrossValidateModel(model = hpoisson.model,model.type = "hgam",data = species.data,
#                                         species = s,group = group.var,key="hauljoin",scale.preds = T)
#
#         hpoisson.errors<-hpoisson.cv[[1]]
#         hpoisson.cv.models<-hpoisson.cv[[2]]
#
#         writeRaster(hpoisson.abund,paste0(species.path,"/hpoisson_abundance"),overwrite=T)
#
#         hpoisson.effects<-GetGAMEffects(model = hpoisson.model,data=species.data,cv.model.list = hpoisson.cv.models,
#                                         scale.factor = hpoisson.scale,vars = "all",scale = "log")
#
#         hpoisson.rmse<-max(c(RMSE(pred = hpoisson.errors$pred,obs = hpoisson.errors$abund),
#                              RMSE(pred = hpoisson.errors$cvpred,obs = hpoisson.errors$abund)))
#
#         hpoisson.success<-ifelse(any(is.na(hpoisson.errors$cvpred)),F,T)
#
#         saveRDS(hpoisson.model,file=paste0(species.path,"/hpoisson_model.rds"))
#         saveRDS(hpoisson.effects,file=paste0(species.path,"/hpoisson_effects.rds"))
#       }else{
#         print(paste("Hurdle gam for",s,"unsuccessful -- moving on",sep=" "))
#       }
#
#
#       ###############
#       # Normal GAM
#       print(paste(Sys.time(),"Starting the poisson model for",region,figure.name,sep=" "))
#       try(poisson.model<-FitGAM(gam.formula = basic.gam.formula,data = species.data,verbose = F,
#                                 reduce = T,select = T,family.gam = "poisson"))
#
#       poisson.converge<-exists("poisson.model") & any(is.infinite(predict(poisson.model,type="response")))==F &
#         any(is.na(predict(poisson.model,type="response")))==F
#
#       if(poisson.converge){
#         poisson.scale<-mean(species.data[,s])/mean(predict(poisson.model,type="response"))
#         poisson.abund<-MakeGAMAbundance(model = poisson.model,r.stack = raster.stack,scale.factor = poisson.scale,
#                                         land = ak.raster,filename = "")
#
#         poisson.abund.check<-cellStats(poisson.abund,max)<(max(species.data[,s])*10)
#       }else{
#         poisson.abund.check<-F
#       }
#
#       if(poisson.abund.check==F){
#         rm(poisson.model,poisson.abund,poisson.scale)
#
#         print("TPS poisson model failed abundance test; Trying alternate version")
#         try(poisson.model<-FitGAM(gam.formula = alt.gam.formula,data = species.data,verbose = F,
#                                   select = T,reduce = T,family.gam = "poisson"))
#
#         poisson.converge<-exists("poisson.model") & any(is.infinite(predict(poisson.model,type="response")))==F &
#           any(is.na(predict(poisson.model,type="response")))==F
#         if(poisson.converge){
#           poisson.scale<-mean(species.data[,s])/mean(predict(poisson.model,type="response"))
#           poisson.abund<-MakeGAMAbundance(model = poisson.model,r.stack = raster.stack,scale.factor = poisson.scale,
#                                           land = ak.raster,filename = "")
#           poisson.abund.check<-cellStats(poisson.abund,max)<(max(species.data[,s])*10)
#
#         }
#       }
#       if(poisson.abund.check==F){
#         print("both versions of poisson model unacceptable; moving to gams")
#         rm(poisson.model,poisson.abund,poisson.scale)
#       }
#
#       if(poisson.converge & poisson.abund.check){
#         print(paste0(Sys.time(),"GAM model successful for ",figure.name," starting evaluations"))
#         poisson.cv<-CrossValidateModel(model = poisson.model,model.type = "gam",species = s,data = species.data,
#                                        group = group.var,key="hauljoin",scale.preds = T)
#         poisson.errors<-poisson.cv[[1]]
#         poisson.cv.models<-poisson.cv[[2]]
#
#         writeRaster(poisson.abund,paste0(species.path,"/poisson_abundance"),overwrite=T)
#
#         poisson.effects<-GetGAMEffects(model = poisson.model,data = species.data,vars = "all",scale = "log",
#                                        cv.model.list = poisson.cv.models,scale.factor = poisson.scale)
#
#         poisson.rmse<-max(c(RMSE(pred = poisson.errors$pred,obs = poisson.errors$abund),
#                             RMSE(pred = poisson.errors$cvpred,obs = poisson.errors$abund)))
#
#         poisson.success<-ifelse(any(is.na(poisson.errors$cvpred)),F,T)
#
#         saveRDS(poisson.model,file=paste0(species.path,"/poisson_model.rds"))
#         saveRDS(poisson.effects,file=paste0(species.path,"/poisson_effects.rds"))
#       }else{
#         print(paste("Gam for",s,"unsuccessful -- moving on",sep=" "))
#       }
#
#       ###############
#       # Let's try a negative binomial poisson GAM
#       print(paste(Sys.time(),"Starting the Negative Binomial GAM model for",region,figure.name,sep=" "))
#       try(negbin.model<-FitGAM(gam.formula = basic.gam.formula,data = species.data,verbose = F,select = T,
#                                reduce = T,family.gam = "nb"))
#
#       negbin.converge<-exists("negbin.model") & any(is.infinite(predict(negbin.model,type="response")))==F &
#         any(is.na(predict(negbin.model,type="response")))==F
#       if(negbin.converge){
#         negbin.converge<-T
#         negbin.scale<-mean(species.data[,s])/mean(predict(negbin.model,type="response"))
#         negbin.abund<-MakeGAMAbundance(model = negbin.model,r.stack = raster.stack,scale.factor = negbin.scale,
#                                        land = ak.raster,filename = "")
#
#         negbin.abund.check<-cellStats(negbin.abund,max)<(max(species.data[,s])*10)
#       }else{
#         negbin.abund.check<-F
#       }
#
#       if(negbin.abund.check==F){
#         rm(negbin.model,negbin.abund,negbin.scale)
#
#         print("TPS negbin model failed abundance test; Trying alternate version")
#         try(negbin.model<-FitGAM(gam.formula = alt.gam.formula,data = species.data,verbose = F,select = T,
#                                  reduce = T,family.gam = "nb"))
#         negbin.converge<-exists("negbin.model") & any(is.infinite(predict(negbin.model,type="response")))==F &
#           any(is.na(predict(negbin.model,type="response")))==F
#         if(negbin.converge){
#           negbin.scale<-mean(species.data[,s])/mean(predict(negbin.model,type="response"))
#           negbin.abund<-MakeGAMAbundance(model = negbin.model,r.stack = raster.stack,scale.factor = negbin.scale,
#                                          land = ak.raster,filename = "")
#           negbin.abund.check<-cellStats(negbin.abund,max)<(max(species.data[,s])*10)
#
#         }
#       }
#
#       if(negbin.abund.check==F){
#         print("both versions of negbin model unacceptable; moving to gams")
#         rm(negbin.model,negbin.abund,negbin.scale)
#       }
#
#       if(negbin.converge & negbin.abund.check){
#         print(paste0(Sys.time(),"NB-GAM model successful for ",figure.name," starting evaluations"))
#         negbin.cv<-CrossValidateModel(model = negbin.model,model.type = "gam",species = s,data = species.data,
#                                       group = group.var,key="hauljoin",scale.preds = T)
#
#         negbin.errors<-negbin.cv[[1]]
#         negbin.cv.models<-negbin.cv[[2]]
#
#         writeRaster(negbin.abund,paste0(species.path,"/negbin_abundance"),overwrite=T)
#
#         negbin.effects<-GetGAMEffects(model = negbin.model,data = species.data,vars = "all",scale = "log",
#                                       cv.model.list = negbin.cv.models,scale.factor = negbin.scale)
#
#         negbin.rmse<-max(c(RMSE(pred = negbin.errors$pred,obs = negbin.errors$abund),
#                            RMSE(pred = negbin.errors$cvpred,obs = negbin.errors$abund)))
#
#         negbin.success<-ifelse(any(is.na(negbin.errors$cvpred)),F,T)
#
#         saveRDS(negbin.model,file=paste0(species.path,"/negbin_model.rds"))
#         saveRDS(negbin.effects,file=paste0(species.path,"/negbin_effects.rds"))
#       }else{
#         print(paste("Negative Binomial Gam for",s,"unsuccessful -- moving on",sep=" "))
#       }
#
#       #############################################################
#       # Special rules for exceptional cases
#       if(region=="GOA" & s=="a_atka"){maxnet.abund.check<-F}
#
#       #############################################################
#       #Now figure out which models were run and make summary tables and plots
#       models<-c("maxnet","cloglog","hpoisson","poisson","negbin")
#
#       model.vec<-NULL
#       model.errors<-list()
#       model.scales<-vector()
#       effects.list<-list()
#       list.index<-1
#       model.types<-NULL
#       rmse.vec<-NULL
#
#       model.converge.vec<-c(maxnet=maxnet.converge,cloglog=cloglog.converge,hpoisson=hpoisson.converge,
#                             poisson=poisson.converge,negbin=negbin.converge)
#       model.abund.check.vec<-c(maxnet=maxnet.abund.check,cloglog=cloglog.abund.check,hpoisson=hpoisson.abund.check,
#                                poisson=poisson.abund.check,negbin=negbin.abund.check)
#
#       # Loop to pull info for each successful model
#       for(m in models){
#         if(get(paste0(m,".abund.check"))){
#           model.vec<-c(model.vec,m)
#           model.types<-c(model.types,ifelse(m%in%c("poisson","negbin"),"gam",ifelse(m=="hpoisson","hgam",m)))
#           model.errors[[list.index]]<-get(paste0(m,".errors"))
#           rmse.vec<-c(rmse.vec,get(paste0(m,".rmse")))
#           model.scales<-c(model.scales,get(paste0(m,".scale")))
#           effects.list[[list.index]]<-get(paste0(m,".effects"))
#           list.index<-list.index+1
#         }
#       }
#       names(rmse.vec)<-model.vec
#       names(model.errors)<-model.vec
#       names(model.scales)<-model.vec
#
#       # now construct the weights to see which models we need to run additional steps on
#       print(paste(Sys.time(),"Making ensemble for",region,figure.name,sep = " "))
#       model.weights<-MakeEnsemble(rmse = rmse.vec,minimum = .1,model.types = model.types)
#
#       # now loop back to work with only the good models
#       model.list<-list()
#       dev.list<-list()
#       cv.model.list<-list()
#       area.vec<-vector(length=length(model.vec))
#       pred.data<-data.frame(x=NA)
#       model.breaks<-vector(length=length(model.vec))
#       abund.list<-list()
#       efh.list<-list()
#       var.list<-list()
#
#
#       # now loop through a run the extra tests for each one
#       for(m in 1:length(model.vec)){
#         model.name<-model.vec[m]
#         if(model.weights[m]>0){
#
#           print(paste0(Sys.time(),": Running additional evaluations for ",figure.name," for ",model.vec[m]," model"))
#
#           model.list[[m]]<-get(paste0(model.name,".model"))
#           cv.model.list[[m]]<-get(paste0(model.name,".cv.models"))
#           abund.list[[m]]<-get(paste0(model.name,".abund"))
#
#
#           if(model.types[m]=="maxnet"){
#             dev.list[[m]]<-MaxnetStats(model=model.list[[m]],data = species.data,species = s,regmult = opt.mult,
#                                        maxnet2d=maxnet2d)
#             maxnet.coefs<-MaxnetCoefs(model.list[[m]],maxnet2d=maxnet2d)
#           }
#           if(model.types[m]%in%c("cloglog","hgam","gam")){
#             dev.list[[m]]<-GAMStats(model=model.list[[m]],data=species.data)
#           }
#           if(is.na(dev.list[[m]][1])==F){
#             good.terms<-names(sort(dev.list[[m]],decreasing = T)[1:min(9,length(dev.list[[m]]))])
#           }else{
#             if(model.types[m]=="maxnet"){
#               good.terms<-names(model.list[[m]]$samplemeans)
#               for(t in 1:length(maxnet2d)){
#                 good2d<-which(good.terms%in%maxnet2d[[t]])
#                 if(length(good2d)>0){
#                   good.terms[good2d[1]]<-paste(maxnet2d[[t]],collapse = "*")
#                 }
#                 if(length(good2d)>1){
#                   good.terms<-good.terms[-good2d[2]]
#                 }
#               }
#               good.terms<-good.terms[1:min(length(good.terms),9)]
#             }else{
#               good.terms<-names(effects.list[[m]])[1:min(length(effects.list[[m]]),9)]
#             }
#           }
#
#
#           m.breaks<-FindEFHbreaks(abund.raster = abund.list[[m]],method = "percentile",threshold=.0513,
#                                   quantiles = c(.05,.25,.5,.75))
#           model.breaks[m]<-m.breaks[2]
#           efh.list[[m]] <- cut(abund.list[[m]], breaks = m.breaks, overwrite = TRUE,
#                                filename = paste0(species.path,"/",model.name,"_efh"))
#
#           area.vec[m]<-sum(getValues(efh.list[[m]])>1,na.rm=T)
#           names(area.vec)[m]<-model.name
#
#           var.list[[m]]<-MakeVarianceRasters(model.list = cv.model.list[[m]],raster.stack = raster.stack,
#                                              model.type = model.types[m],scale.factor = model.scales[m])
#           writeRaster(x = var.list[[m]],filename = paste0(species.path,"/",model.name,"_abund_variance"),overwrite=T)
#
#
#           png(filename = paste0(species.path,"/",model.name,"_abund_stdev.png"),width = png.width,height = png.height,res=600,units="in")
#           plotAbundance(map=var.list[[m]],back.col = NA,legend.text = .8,label.size=.8,
#                         outline=T,outline.lwd=.75,legend.name = "Standard Deviation of Predicted Abundance")
#           AddGrid(legPosition = NA,horiz=legend.horiz,grid=F,depth=T,land.col = "grey30",
#                   axistext.size=.8)
#           dev.off()
#
#
#           png(filename = paste0(species.path,"/",model.name,"_residuals.png"),width = 6,height = 6,res=600,units="in")
#           MakeCrossValidationPlots(error.data = model.errors[[m]],method = "spearman",make.hist = F)
#           dev.off()
#
#           #More plots
#           png(filename = paste0(species.path,"/",model.name,"_abundance.png"),width = png.width,height = png.height,res=600,units="in")
#           plotAbundance(map=abund.list[[m]],back.col = NA,legend.text = .8,label.size=.8,
#                         outline=T,outline.lwd=.75)
#           AddGrid(legPosition = NA,horiz=legend.horiz,grid=F,depth=T,land.col = "grey30",
#                   axistext.size=.8)
#           dev.off()
#
#           png(filename = paste0(species.path,"/",model.name,"_efh.png"),width = png.width,height = png.height,res=600,units="in")
#           plotEFH(map = efh.list[[m]],label.size = .8,outline=T,outline.lwd=.75)
#           AddGrid(legPosition = legend.pos,horiz=legend.horiz,grid=F,grid.col = "grey70",
#                   depth=T,land.col = "grey30",dlabel.size = .3,dlabels = T,axistext.size = .8,
#                   legend.size = .8)
#           dev.off()
#
#
#           png(filename = paste0(species.path,"/",model.name,"_effects.png"),width = 6,height = 6,units = "in",res = 600)
#           plotEffects(effects = effects.list[[m]],nice.names = nice.names,vars = good.terms,land = ak.coast)
#           dev.off()
#
#           if(model.vec[m]=="maxnet"){
#             MakeXtable(model = maxnet.model,model.type = "maxnet",abund = model.errors[[m]]$abund,rsq.method = "spearman",
#                        train = model.errors[[m]]$pred,test=model.errors[[m]]$cvpred,scale.factor = model.scales[m],
#                        ncoefs = maxnet.coefs,devs=dev.list[[m]],area=area.vec[m],group = model.errors[[m]][,group.var],
#                        nice.names = nice.names,efh.break = model.breaks[m],regmult = opt.mult,
#                        filename = paste0(species.path,"/maxnet_table.html"))
#           }else{
#             MakeXtable(model = model.list[[m]],model.type = model.types[m],abund=model.errors[[m]]$abund,
#                        train = model.errors[[m]]$pred,test=model.errors[[m]]$cvpred,scale.factor = model.scales[m],
#                        nice.names = nice.names,rsq.method = "spearman",
#                        devs = dev.list[[m]],area = area.vec[m],group = model.errors[[m]][,group.var],
#                        efh.break = model.breaks[m],filename=paste0(species.path,"/",model.name,"_table.html"))
#           }
#         }else{
#           model.list[[m]]<-get(paste0(model.name,".model"))
#           dev.list[[m]]<-NA
#           abund.list[[m]]<-NA
#           model.breaks[m]<-NA
#           efh.list[[m]]<-NA
#           # save a few things anyway, for potential use with some versions of the effects plot
#           cv.model.list[[m]]<-get(paste0(model.name,".cv.models"))
#           area.vec[m]<-NA
#           var.list[[m]]<-NA
#         }
#       }
#       print("Constructing the ensemble")
#       # now make the actual ensemble
#       ensemble.abund<-MakeEnsembleAbundance(model.weights=model.weights,abund.list = abund.list,
#                                             filename = paste0(species.path,"/ensemble_abundance"))
#       ensemble.breaks<-FindEFHbreaks(abund.raster = ensemble.abund,method = "percentile",threshold = .0513,
#                                      data = species.data,quantiles = c(.05,.25,.5,.75))
#       ensemble.efh<-cut(ensemble.abund, breaks = ensemble.breaks, overwrite = TRUE,
#                         filename = paste0(species.path,"/ensemble_efh"))
#       area.vec<-c(area.vec,ensemble=sum(getValues(ensemble.efh)>1,na.rm=T))
#
#       model.breaks<-c(model.breaks,ensemble=ensemble.breaks[2])
#
#       #figure out the weighted average deviance explained
#       dev.names=NULL
#       for(d in 1:length(dev.list)){
#         dev.names<-unique(c(dev.names,names(dev.list[[d]])))
#       }
#       dev.dat<-data.frame(name=dev.names,as.data.frame(matrix(nrow=length(dev.names),ncol=sum(model.weights>0))),
#                           stringsAsFactors = F)
#       dev.index<-2
#       for(d in 1:length(dev.list)){
#         if(model.weights[d]>0){
#           if(is.na(dev.list)[d]==F){
#             match.vec<-match(names(dev.list[[d]]),dev.dat$name)
#             dev.dat[match.vec,dev.index]<-dev.list[[d]]
#             dev.dat[is.na(dev.dat[dev.index]),dev.index]<-0
#             dev.index<-dev.index+1
#
#           }else{
#             dev.dat[,dev.index]<-NA
#             dev.index<-dev.index+1
#           }
#         }
#       }
#       if(ncol(dev.dat)>2){
#         ave.devs<-apply(dev.dat[,2:ncol(dev.dat)],MARGIN = 1,FUN = "weighted.mean",w=model.weights[model.weights>0],na.rm=T)
#       }else{
#         ave.devs<-dev.dat[,2]
#       }
#       good.terms<-dev.dat$name[order(ave.devs,decreasing = T)[1:min(9,length(ave.devs))]]
#
#       en.eff<-GetEnsembleEffects(effects.list = effects.list,model.weights = model.weights,
#                                  vars = "all",scale = "log")
#       png(filename = paste0(species.path,"/ensemble_effects.png"),
#           width = 8,height = 10,units = "in",res = 600)
#       plotEffects(effects = en.eff,nice.names = nice.names,vars = good.terms,land = ak.coast)
#       dev.off()
#
#       saveRDS(en.eff,file=paste0(species.path,"/effects_data.rds"))
#
#       png(filename = paste0(species.path,"/ensemble_residuals.png"),width = 6,height = 6,res=600,units="in")
#       ensemble.preds<-ValidateEnsemble(pred.list = model.errors,model.weights = model.weights,output = T,make.plots = T,
#                                        key = "hauljoin",latlon = T,group = "Folds",method="spearman")
#       dev.off()
#
#       #Calculate remaining statistics for the ensemble
#       ensemble.acc<-EFHAccuracy(data=ensemble.preds,efh.break = ensemble.breaks)
#
#       en.auc.data<-data.frame(1:nrow(ensemble.preds),ensemble.preds$abund,ensemble.preds$prob)
#       ensemble.auc<-PresenceAbsence::auc(en.auc.data)[1,1]
#
#       ensemble.DevExpP<-PDE(obs = ensemble.preds$abund,pred = ensemble.preds$pred)
#
#       png(filename = paste0(species.path,"/ensemble_abundance.png"),width = png.width,height = png.height,res=600,units="in")
#       plotAbundance(map=ensemble.abund,back.col = NA,legend.text = .8,label.size=.8,
#                     outline=T,outline.lwd=.75)
#       AddGrid(legPosition = NA,horiz=legend.horiz,grid=F,depth=T,land.col = "grey30",
#               axistext.size=.8)
#       dev.off()
#
#       png(filename = paste0(species.path,"/ensemble_efh.png"),width = png.width,height = png.height,res=600,units="in")
#       plotEFH(map = ensemble.efh,label.size = .8,outline=T,outline.lwd=.75)
#       AddGrid(legPosition = legend.pos,horiz=legend.horiz,grid=F,grid.col = "grey70",
#               depth=T,land.col = "grey30",dlabel.size = .3,dlabels = T,axistext.size = .8,
#               legend.size = .8)
#       dev.off()
#
#       #ensemble variance
#       ensemble.var<-GetEnsembleVariance(model.weights = model.weights,variance.list = var.list,abund.list = abund.list,
#                                         ensemble.abund = ensemble.abund)
#       writeRaster(x = ensemble.var,filename = paste0(species.path,"/ensemble_variance"),overwrite=T)
#
#       png(filename = paste0(species.path,"/ensemble_abund_stdev.png"),width = png.width,height = png.height,res=600,units="in")
#       plotAbundance(map=ensemble.var,back.col = NA,legend.text = .8,label.size=.8,
#                     outline=T,outline.lwd=.75,legend.name = "Standard Deviation of Predicted Abundance")
#       AddGrid(legPosition = NA,horiz=legend.horiz,grid=F,depth=T,land.col = "grey30",
#               axistext.size=.8)
#       dev.off()
#
#       # Ensemble deviance table
#       MakeDevianceTable(model.names = model.vec,model.types = model.types,dev.list = dev.list,model.weights = model.weights,
#                         filename = paste0(species.path,"/deviance_table.html"),nice.names=nice.names2)
#
#       # set up some stuff for the species tables
#       ensemble.preds$cvpred<-ensemble.preds$pred
#       ensemble.preds$cverror<-ensemble.preds$error
#
#       model.errors[[list.index]]<-ensemble.preds
#       names(model.errors)[list.index]<-"ensemble"
#
#       # save the cv errors and predictions, in case you want to analyze them later
#       model.vec<-c(model.vec,"ensemble")
#
#       error.table<-data.frame(Model=model.vec[1],model.errors[[1]],stringsAsFactors = F)
#       for(f in 2:length(model.errors)){
#         error.table<-rbind(error.table,data.frame(Model=model.vec[f],model.errors[[f]],stringsAsFactors = F))
#       }
#       write.csv(x = error.table,file = paste0(species.path,"/cv_error_data.csv"),row.names = F)
#
#       MakeEnsembleXtable(weights = model.weights,preds.table = error.table,converge.vec = model.converge.vec,
#                          abund.check.vec=model.abund.check.vec,scale.facs = model.scales,areas = area.vec,
#                          efh.breaks=model.breaks,cor.method="spearman",
#                          filename = paste0(species.path,"/ensemble_table.html"))
#
#       ensemble.success<-T
#     }
#
#     # update the progress table
#     progress.table<-read.csv(file=paste0(EFH.path,"/ModelProgress_offset_reruns.csv"),stringsAsFactors = F)
#
#     progress.table$N[i]<-n.pres
#
#     progress.table$Maxnet[i]<-maxnet.abund.check
#     progress.table$Cloglog[i]<-cloglog.abund.check
#     progress.table$Hpoisson[i]<-hpoisson.abund.check
#     progress.table$Poisson[i]<-poisson.abund.check
#     progress.table$Negbin[i]<-negbin.abund.check
#     progress.table$Ensemble[i]<-ensemble.success
#
#     progress.table$Completed[i]<-as.character(Sys.time())
#     if(update.table){
#       try(write.csv(progress.table,file=paste0(EFH.path,"/ModelProgress_offset_reruns.csv"),row.names = F))
#       try(write.csv(x = progress.table,file =paste0(EFH.path,"/ProgressChecker.csv"),row.names = F))
#     }
#     if(stop.early){
#       stop("Stopping Early due to user selected request!!")
#     }
#     #remove anything you don't want to carry over to the next species/life-stage
#     rm(start.time,n.pres,
#        maxnet.model,cloglog.model,hpoisson.model,poisson.model,negbin.model,model.list,
#        maxnet.abund,cloglog.abund,hpoisson.abund,poisson.abund,negbin.abund,ensemble.abund,
#        maxnet.breaks,cloglog.breaks,hpoisson.breaks,poisson.breaks,negbin.breaks,ensemble.breaks,model.breaks,
#        maxnet.efh,cloglog.efh,hpoisson.efh,poisson.efh,negbin.efh,ensemble.efh,efh.list,
#        maxnet.cv,cloglog.cv,hpoisson.cv,poisson.cv,negbin.cv,
#        maxnet.errors,cloglog.errors,hpoisson.errors,poisson.errors,negbin.errors,
#        maxnet.rmse,cloglog.rmse,hpoisson.rmse,poisson.rmse,negbin.rmse,
#        maxnet.cv.models,cloglog.cv.models,hpoisson.cv.models,poisson.cv.models,negbin.cv.models,
#        maxnet.var,cloglog.var,hpoisson.var,poisson.var,negbin.var,ensemble.var,
#        maxnet.variance,cloglog.variance,hpoisson.variance,poisson.variance,negbin.variance,ensemble.variance,
#        maxnet.error.raster,cloglog.error.raster,hpoisson.error.raster,poisson.error.raster,negbin.error.raster,ensemble.error.raster,
#        maxnet.preds,cloglog.preds,hpoisson.preds,poisson.preds,negbin.preds,ensemble.preds,
#        max.stats,cloglog.dev,hpoisson.dev,poisson.dev,negbin.dev,
#        maxnet.area,cloglog.area,hpoisson.area,poisson.area,negbin.area,ensemble.area,
#        maxnet.effects,cloglog.effects,hpoisson.effects,poisson.effects,negbin.effects,
#        maxnet.converge,cloglog.converge,hpoisson.converge,poisson.converge,negbin.converge,ensemble.success,model.converge.vec,
#        maxnet.abund.check,cloglog.abund.check,hpoisson.abund.check,poisson.abund.check,negbin.abund.check,
#        maxnet.abund.check0,maxnet.converge0,model.abund.check.vec,maxnet.scale.vec,
#        maxnet.scale,cloglog.scale,hpoisson.scale,poisson.scale,negbin.scale,
#        cv.model.list,effects.list,maxnet.abund.list,maxnet.cv.model.list,maxnet.model.list,maxnet.pred.list,
#        dev.dat,dev.list,maxstats,maxnet.coefs,maxnet2d,var.list,rmse.vec,
#        model.types,model.vec,train.list,test.list,pred.data,area.vec,abund.list,list.index,auc.vec,model.weights,
#        weights,abund.list,drop.terms,good.terms,good.terms0,en.eff,en.eff2,error.table,model.errors,model.scales)
#     print(paste0(Sys.time()," finished models for species: ",s))
#   }
# }
#
#
