## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##  
##  Centre of Marine Sciences [ ccmar.ualg.pt ]
##  Faro, Portugal
##
## ---------------------------------------------------------------------
##
##  Machine Learning Species Distribution Modelling [ Ver.301 ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

for( scenario in scenariosToPredict) {
  
  # ---------------------
  # Environmental layers
  
  if( ! dir.exists(paste0(resultsDirectory,"/Predictions/",scenario,"/")) ) { dir.create( paste0(resultsDirectory,"/Predictions/",scenario,"/") , recursive = TRUE) }
  
  rasterLayers <- list.files(gsub("Baseline",scenario,climateLayersDirectory),pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
  rasterLayers <- stack(rasterLayers[as.vector(sapply(tolower(dataLayers),function(x) { which( grepl(x,tolower(rasterLayers))) } ))])
  
  if( scenario == scenariosToPredict[1] ) {
    shape <- subset(rasterLayers,1); shape[!is.na(shape)] <- 1; names(shape) <- "layer"; shape[ shape == 1] <- 0
  }

  occurrenceRecords <- modelData@coords[modelData@pa == 1,]
  rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer, minDepth=ifelse( is.null(minDepth), "NULL" , ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer)) , maxDepth=ifelse( is.null(maxDepth), "NULL" , maxDepth+maxDepthBuffer),intertidal)
  names(rasterLayers) <- dataLayersName
  
  # ---------------------

  modelDataList <- unlist(sapply(algorithms, function(x) { list.files(paste0(resultsDirectory,"/Models/"), pattern = paste0(x,".RData"), full.names = TRUE) } ))
  if( reduceModelcomplexity ) { modelDataList <- modelDataList[grepl("reducedModel",modelDataList)] } 
  if(! reduceModelcomplexity ) { modelDataList <- modelDataList[!grepl("reducedModel",modelDataList)] } 
  ensemblePrediction <- ensembleModels(modelDataList,rasterLayers)
  ensemble <- ensemblePrediction$prediction
  ensembleSD <- ensemblePrediction$prediciton.sd

  # ---------------------

  if( scenario == "Baseline" ) {
    accur <- predictedPerformance(ensemble,modelData,reclassificationIndex)
    reclassificationThreshold <- accur$threshold
    save( accur , file=paste0(resultsDirectory,"/SummaryModels/ensemblePerformance.RData"))
  }

  ensembleReclass <- reclassifyPredicted(ensemble,reclassificationThreshold)

  # ---------------------
  
  save( ensemble , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensemble.RData"), compress=TRUE, compression_level=6)
  save( ensembleSD , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleSD.RData"), compress=TRUE, compression_level=6)
  save( ensembleReclass , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleReclass.RData"), compress=TRUE, compression_level=6)
  
  # ---------------------
  # Reachability
  
  if(scenario == "Baseline") {
    ensembleReclassReach <- reachableRange(ensembleReclass,dispersalFactor=2)
  }
  
  if(scenario != "Baseline") {
    ensembleReclassReach <- reachableRange(ensembleReclass,dispersalFactor=potentialDispersalFactor,nonReachableCells)
  }

  save( ensembleReclassReach , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleReclassReachable.RData"), compress=TRUE, compression_level=6)
  
  if(scenario == "Baseline") {
    accur <- predictedPerformance(ensemble*ensembleReclassReach,modelData,reclassificationIndex)
    save( accur , file=paste0(resultsDirectory,"/SummaryModels/ensemblePerformanceReachable.RData"))
    nonReachableCells <- Which(ensembleReclassReach == 0 & ensembleReclass == 1, cells=TRUE)
  }
  
  # ---------------------------------------------------------------------
  # ---------------------------------------------------------------------
  # Project to global extent

  ensembleGlobal <- ensembleReclassGlobal <- ensembleReclassReachGlobal <- shape
  
  ensembleDF <- as.data.frame(ensemble, xy=T, na.rm=T)
  ensembleReclassDF <- as.data.frame(ensembleReclass, xy=T, na.rm=T)
  ensembleReclassReachDF <- as.data.frame(ensembleReclassReach, xy=T, na.rm=T)
  
  ensembleGlobal[ cellFromXY( ensembleGlobal , ensembleDF[,c("x","y")] ) ] <- ensembleDF[,"layer"]
  ensembleReclassGlobal[ cellFromXY( ensembleReclassGlobal , ensembleReclassDF[,c("x","y")] ) ] <- ensembleReclassDF[,"layer"]
  ensembleReclassReachGlobal[ cellFromXY( ensembleReclassReachGlobal , ensembleReclassReachDF[,c("x","y")] ) ] <- ensembleReclassReachDF[,"layer"]
  
  save( ensembleGlobal , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleGlobal.RData"), compress=TRUE, compression_level=6)
  save( ensembleReclassGlobal , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleReclassGlobal.RData"), compress=TRUE, compression_level=6)
  save( ensembleReclassReachGlobal , file=paste0(resultsDirectory,"/Predictions/",scenario,"/ensembleReclassReachableGlobal.RData"), compress=TRUE, compression_level=6)
  
  # ---------------------
  # ensemble tipping points
  
  if(scenario == "Baseline") {
    
    ensembleReclass.i <- ensembleReclassReach
    ensembleReclass.i[ ensembleReclass.i == 0] <- NA
    
    ensembleTP <- data.frame()
    
    for(i in 1:length(names(rasterLayers))) {
      
      predictor.i <- names(subset(rasterLayers,i))
      res.i <- getValues(mask(subset(rasterLayers,i) , ensembleReclass.i))
      res.i <- res.i[!is.na(res.i)]
      
      if( monotonicity[predictor.i] == 1 ) { res.q <- quantile(res.i, probs=c(0.05)); res.i <- min(res.i) }
      if( monotonicity[predictor.i] == 0 ) { res.q <- quantile(res.i, probs=c(0.5)); res.i <- mean(res.i) }
      if( monotonicity[predictor.i] == -1 ) { res.q <- quantile(res.i, probs=c(0.95)); res.i <- max(res.i) }
      
      ensembleTP <- rbind(ensembleTP,data.frame(predicor = predictor.i , 
                                                tp = res.i,
                                                tp.quantile = res.q,
                                                row.names = NULL
      ))
      
    }
    
    save( ensembleTP , file=paste0(resultsDirectory,"/SummaryModels/EnsembleTippingPoints.RData"))
    
  }
  
  
  # ---------------------
  
  if(scenario == "Baseline") {
    
    depthData <- loadRData(file=paste0(resultsDirectory,"/SummaryModels/","depthData.RData"))
    depthDataPredict <- data.frame(variable=c("Min","Mean","Max","q95"),matrix(NA,nrow=4,ncol=length(scenariosToPredict)+2))
    colnames(depthDataPredict) <- c("variable","Known","Observed",paste0("Pred",scenariosToPredict))
    depthDataPredict[1,2] <- depthData$minKnownDepth
    depthDataPredict[3,2] <- depthData$maxKnownDepth
    depthDataPredict[1,3] <- depthData$minObservedDepth
    depthDataPredict[3,3] <- depthData$maxObservedDepth
    depthDataPredict[4,3] <- depthData$q95ObservedDepth
    
  }

  if( length(Which(ensembleReclassGlobal == 1,cells=T)) > 0) {
    
    occurrenceRecordsDepths <- abs( raster::extract(raster(bathymetryDataLayer),xyFromCell(ensembleReclassGlobal,Which(ensembleReclassGlobal == 1,cells=T))) )
    
  }
  
  if( length(Which(ensembleReclassGlobal == 1,cells=T)) == 0) {
    
    occurrenceRecordsDepths <- NA
    
  }
  
  depthDataPredict[1,paste0("Pred",scenario)] <- min(occurrenceRecordsDepths,na.rm=T)
  depthDataPredict[2,paste0("Pred",scenario)] <- mean(occurrenceRecordsDepths,na.rm=T)
  depthDataPredict[3,paste0("Pred",scenario)] <- max(occurrenceRecordsDepths,na.rm=T)
  depthDataPredict[4,paste0("Pred",scenario)] <- quantile(occurrenceRecordsDepths,probs=0.95,na.rm=T)
  depthDataPredict[depthDataPredict == Inf | depthDataPredict == -Inf | depthDataPredict == "NaN" ] <- NA
  
  save(depthDataPredict,file=paste0(resultsDirectory,"/SummaryModels/","depthDataPredicted.RData"))
  
}

# End of Code
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------