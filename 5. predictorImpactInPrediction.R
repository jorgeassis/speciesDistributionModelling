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



REVIEW




# Read data

modelDataList <- list.files(paste0(resultsDirectory,"/RData/"), pattern="reducedModel", full.names = TRUE)

# -----------------------------

rasterLayersBaseline <- list.files(dataLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
rasterLayersBaseline <- stack(rasterLayersBaseline[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayersBaseline)) } ))])
rasterLayersBaseline <- processLayers(rasterLayersBaseline,occurrenceRecords,regionBuffer, minDepth=ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer) , maxDepth=maxDepth+maxDepthBuffer,intertidal)
names(rasterLayersBaseline) <- dataLayersName

# -----------------------------

for( scenario in scenariosToPredict[scenariosToPredict != "Baseline"] ) {
  
  # ---------------------
  
  rasterLayers <- list.files(gsub("Baseline",scenario,dataLayersDirectory),pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
  rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
  
  rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer, minDepth=ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer) , maxDepth=maxDepth+maxDepthBuffer,intertidal)
  names(rasterLayers) <- dataLayersName
  
  # ---------------------
  
  ensemblePrediction <- list()
  
  for( algorithmType in algorithms ) {
    
    load( modelDataList[grepl(algorithmType,modelDataList)] )
    
    if( algorithmType == "BRT" ) { ensemblePrediction <- c(ensemblePrediction,list( predictDistribution(rasterLayers,modelFull$model , reclassToOne=FALSE) )) }
    if( algorithmType == "MBOOST" ) { ensemblePrediction <- c(ensemblePrediction,list( predictDistribution(rasterLayers,modelFull$model , reclassToOne=FALSE))) }
    if( algorithmType == "MPNN") { stop("!") }
    
  }
  
  beginCluster(nCores)
  finalEnsemble <- clusterR(stack(unlist(ensemblePrediction)), calc, args=list(mean))
  endCluster()
  
  # ---------------------
  
  predictorImpactGain <- list()
  predictorImpactLoss <- list()
  
  for(l in 1:nlayers(rasterLayers)) {
    
    rasterLayers.l <- rasterLayers
    rasterLayers.l <- dropLayer(rasterLayers.l,l)
    rasterLayers.l <- addLayer(rasterLayers.l,subset(rasterLayersBaseline,l))
    
    ensemblePrediction.l <- list()
    
    for( algorithmType in algorithms ) {
      
      load( modelDataList[grepl(algorithmType,modelDataList)] )
      
      if( algorithmType == "BRT" ) { ensemblePrediction.l <- c(ensemblePrediction.l,list( predictDistribution(rasterLayers.l,modelFull$model , reclassToOne=FALSE) )) }
      if( algorithmType == "MBOOST" ) { ensemblePrediction.l <- c(ensemblePrediction.l,list( predictDistribution(rasterLayers.l,modelFull$model , reclassToOne=FALSE))) }
      if( algorithmType == "MPNN") { stop("!") }
      
    }
    
    beginCluster(nCores)
    finalEnsemble.l <- clusterR(stack(unlist(ensemblePrediction.l)), calc, args=list(mean))
    endCluster()
    
    finalEnsemble.l.gain <- finalEnsemble.l - finalEnsemble
    finalEnsemble.l.gain[finalEnsemble.l.gain < 0] <- 0
    
    finalEnsemble.l.loss <- finalEnsemble.l - finalEnsemble
    finalEnsemble.l.loss[finalEnsemble.l.loss > 0] <- 0
    
    predictorImpactGain <- c(predictorImpact,list( finalEnsemble.l.gain ))
    predictorImpactLoss <- c(predictorImpact,list( finalEnsemble.l.loss ))
    
  }
  
  # ---------------------
  
  predictorImpactLossSum <- calc(stack(predictorImpactLoss),sum)
  predictorImpactLossStand <- list()
  
  for( l in 1:length(predictorImpact)) {
    
    predictorImpactLoss.l <- predictorImpactLoss[[l]] / predictorImpactLossSum
    predictorImpactLoss.l[predictorImpactLoss.l > 1] <- 1
    predictorImpactLossStand <- c(predictorImpactLossStand,list(predictorImpactLoss.l))
    
  }
  
  predictorImpactLossStand <- stack(predictorImpactLossStand)
  predictorImpactLossStand <- predictorImpactLossStand * finalEnsemble
  names(predictorImpactLossStand) <- names(rasterLayers)
  
  # ---------------------
  
  predictorImpactGainSum <- calc(stack(predictorImpactGain),sum)
  predictorImpactGainStand <- list()
  
  for( l in 1:length(predictorImpact)) {
    
    predictorImpactGain.l <- predictorImpactGain[[l]] / predictorImpactGainSum
    predictorImpactGain.l[predictorImpactGain.l > 1] <- 1
    predictorImpactGainStand <- c(predictorImpactGainStand,list(predictorImpactGain.l))
    
  }
  
  predictorImpactGainStand <- stack(predictorImpactGainStand)
  predictorImpactGainStand <- predictorImpactGainStand * finalEnsemble
  names(predictorImpactGainStand) <- names(rasterLayers)
  
  # ---------------------
  
  save( finalEnsemble , file=paste0(resultsDirectory,"/RData/reducedModelEnsemble.RData"))

}
  

