## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------

simplifyType <- "negativeRelativeImportance"
keepPredictorsThreshold <- 0.05

## -----------------------

if( simplifyType == "sortRelativeImportance" ) {
  
    # Version by sorting per variable and testing difference in performance
  
    if( exists("modelFull") ) { modelFullBk <- modelFull }
    
    allButOne <- FALSE
    
    relativeImportance <- relativeContribution$dataFrame
    relativeImportance <- relativeImportance[relativeImportance$Predictor %in% names(rasterLayers) ,]
    
    if( simplifyType == "deviance" ) {
      
      predictionModel <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
      predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
      observed <- speciesData$PA
      observed <- observed[!is.na(predicted)]
      predicted <- predicted[!is.na(predicted)]
      modelFullDE <- Dsquared(glm(predicted~observed,family=binomial()))
      
    }
    
    if( simplifyType == "auc" ) {
      
      modelFullDE <- modelFull$auc
      
    }
    
    lossFunction <- data.frame(predictorLoss=rep(NA,length(names(rasterLayers))-1),loss=NA,loss0.01=NA,loss0.99=NA)
    sortLossfunction <- as.numeric(sapply(as.character(relativeImportance[sort(relativeImportance[,2],decreasing = FALSE, index.return=TRUE)$ix,1]),function(x){ which(names(rasterLayers) == x )}))
    
    for(v in 1:(ifelse( length(names(rasterLayers)) >= 4 , length(names(rasterLayers)) - 1 , length(names(rasterLayers))-1 )) ) {
      
      rasterLayersBk <- rasterLayers
      rasterLayers <- subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] )
      
      source("Dependencies/modelMe.R")
      model <- modelFull
      
      if( simplifyType == "deviance" ) {
        
        lab <- "(Deviance Explained)"
        
        if( length(names(rasterLayers)) == 1 ) { 
          dummy.raster <- rasterLayers
          dummy.raster[1:length(rep(1:2,length(dummy.raster) / 2))] <- rep(1:2,length(dummy.raster) / 2)
          names(dummy.raster) <- "layer"
          rasterLayers <- stack(rasterLayers,dummy.raster) 
          }
        
        predictionModel <- predictDistribution(rasterLayers,model$model, reclassToOne=FALSE)
        predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
        observed <- speciesData$PA
        observed <- observed[!is.na(predicted)]
        predicted <- predicted[!is.na(predicted)]
        partialModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
        partialModelDEConf <- 0 
        partialModelDEConf$t <- c(0,0,0,0,0)
        
      }
      
      if( simplifyType == "auc" ) {
        
        lab <- "(Area Under the Curve)"
        
        partialModelDE <- mean(model$cv)
        partialModelDEConf <- list(t=c(model$cv,model$cv))
        
      }
      
      rasterLayers <- rasterLayersBk
      
      lossFunction[v,1] <- v
      lossFunction[v,2] <- modelFullDE - partialModelDE
      lossFunction[v,3] <- modelFullDE - partialModelDE - (sd(partialModelDEConf$t) / sqrt(length(partialModelDEConf$t) ))
      lossFunction[v,4] <- modelFullDE - partialModelDE + (sd(partialModelDEConf$t) / sqrt(length(partialModelDEConf$t) ))
      
    }
    
    keepPredictors <- which( lossFunction[,4] >= keepPredictorsThreshold)
    
    if( length(keepPredictors) == 0) { keepPredictors <- 1 }
    if( ! 1 %in% keepPredictors) { removePredictors <- 1:(min(keepPredictors)-1) }
    if(   1 %in% keepPredictors) { removePredictors <- numeric(0) }
    
    colorsPredictors <- rep("Black",nrow(lossFunction))
    colorsPredictors[removePredictors] <- "Red"
    
    lossFunctionPlot <- ggplot(lossFunction, aes(x=predictorLoss, y=loss)) + 
      geom_line(linetype = "dashed",color="Gray") + #  ylim(-0.5,0.5) +
      geom_point(color=colorsPredictors, size=2.5) +
      geom_hline(aes(yintercept=0  ),color="Black", size=0.2) + 
      xlab("Variables removed (number)") + ylab(paste0("Loss function ",lab)) +
      geom_ribbon(aes(ymin = lossFunction[,3],ymax = lossFunction[,4]), alpha = 0.3)
    
    if( length(removePredictors) > 0 ) {
      
      rasterLayersBk <- rasterLayers
      dataLayersMonotonocityBk <- dataLayersMonotonocity
      
      if(allButOne) { if(length(removePredictors) > 1) { removePredictors <- removePredictors[-length(removePredictors)]} }
      v <- ifelse( length(keepPredictors) > 0 , min(keepPredictors) - 1 , max(removePredictors) ) 
      rasterLayers <- subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] )
      dataLayersMonotonocity <- data.frame(t( sapply( names(rasterLayers.v) ,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )))
      
      source("Dependencies/modelMe.R")
      
      keptVariables <- names(subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] ))
      removedVariables <- names(subset(rasterLayers, (1:length(names(rasterLayers)))[sortLossfunction[1:v]] ))
      
      rasterLayers <- rasterLayersBk
      dataLayersMonotonocity <- dataLayersMonotonocityBk
        
    }
    
    if( length(removePredictors) == 0 ) { 
      
      modelFull <- modelFullBk
      keptVariables <- names(rasterLayers)
      removedVariables <- character(0)
      
    }
    
}

## -----------------------
## -----------------------

if( simplifyType == "negativeRelativeImportance" ) {
  
  # Simple version by removing variables with negative importance
  
  lossFunctionPlot <- NULL
  rasterLayersBk <- rasterLayers
  keptVariables <- names(subset(rasterLayers, which(names(rasterLayers) %in% relativeContributionDF[which(relativeContributionDF$relImportance > 0),"Predictor"]) ))
  removedVariables <- names(subset(rasterLayers, which(names(rasterLayers) %in% relativeContributionDF[which(relativeContributionDF$relImportance < 0),"Predictor"]) ))
  
  rasterLayers <- subset(rasterLayers,which(names(rasterLayers) %in% relativeContributionDF[which(relativeContributionDF$relImportance > 0),"Predictor"]))

  source("Dependencies/modelMe.R")
  
}

## -----------------------

modelReduction <- list(model=modelFull,plot=lossFunctionPlot,keptVariables=keptVariables,removedVariables=removedVariables)

## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------
