## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------

getContributionType <- "modelBased" # testBased

if( getContributionType == "modelBased" ) {
  
  if( modelType == "XGBOOST" ) {
    
    relativeImportance = as.data.frame(xgb.importance( modelFull$model$feature_names, model = modelFull$model)[,1:2])
    rownames(relativeImportance) <- NULL
    colnames(relativeImportance) <- c("Predictor","relImportance")
    relativeImportance$relImportance <- ( (relativeImportance$relImportance ) / sum(relativeImportance$relImportance) ) * 100
    
    missingPredictors <- modelFull$model$feature_names[!modelFull$model$feature_names %in% relativeImportance$Predictor]
    
    if(length(missingPredictors) > 0) {
      relativeImportance <- rbind(relativeImportance,data.frame(Predictor=missingPredictors,relImportance=0))
    }
    relativeImportanceStats <- relativeImportance
    colnames(relativeImportanceStats) <- c("Predictor","relImportanceSignif")
    relativeImportanceStats$relImportanceSignif <- FALSE
    
  }
  
  if( modelType == "BRT" ) {
    
    relativeImportance <- summary(modelFull$model)
    rownames(relativeImportance) <- NULL
    colnames(relativeImportance) <- c("Predictor","relImportance")
    relativeImportanceStats <- relativeImportance
    colnames(relativeImportanceStats) <- c("Predictor","relImportanceSignif")
    relativeImportanceStats$relImportanceSignif <- FALSE
    
  }
  
  if( modelType == "MBOOST" ) {
    
    relativeImportance <- data.frame( Predictor = names(rasterLayers), relImportance=NA)
    
    for( p.i in 1:length(names(rasterLayers))) {
      relativeImportance[p.i,2] <- as.numeric(varimp(modelFull$model)[p.i])
    }
    
    rownames(relativeImportance) <- NULL
    relativeImportance$relImportance <- ((relativeImportance$relImportance - min(relativeImportance$relImportance) ) / sum(relativeImportance$relImportance) ) * 100
    relativeImportanceStats <- relativeImportance
    colnames(relativeImportanceStats) <- c("Predictor","relImportanceSignif")
    relativeImportanceStats$relImportanceSignif <- FALSE
    
  }
  
}

## ----------

if( getContributionType == "testBased" ) {

        varImpType <- "deviance" # boyce deviance auc tss // cvIndex
        
        if( exists("modelFull") ) { modelFullBk <- modelFull }
        
        accurModelFull <- accuracyPredicted( prediction ,speciesData ,cvIndex)
        modelFullDE <- accurModelFull[,varImpType]
        
        relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=NA)
        relativeImportanceStats <- data.frame(Predictor=names(rasterLayers),relImportanceSignif=NA)
        
        for(v in 1:length(names(rasterLayers)) ) {
          
          rasterLayersBk <- rasterLayers
          rasterLayers <- subset(rasterLayers, (1:length(names(rasterLayers)))[-v] )
          
          source("Dependencies/modelMe.R")
          
          if( !is.null(modelFull$model)) { 
            predictionModelPartial <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE) 
            accurModelPartial <- accuracyPredicted( predictionModelPartial ,speciesData,cvIndex)
            partialModelDE <- accurModelPartial[,varImpType]
          }
          
          if( is.null(modelFull$model)) { 
            partialModelDE <- 0 
          }
          
          relativeImportance[v,2] <- modelFullDE - partialModelDE
          
          if( length(modelFull$cv) > 0 ) { 
            relativeImportanceStats[v,2] <- wilcox.test(modelFullBk$cv,modelFull$cv,alternative = c("greater"), paired = ifelse(length(modelFull$cv) == length(model$c),TRUE,FALSE))$p.value < 0.05
          }
          
          rasterLayers <- rasterLayersBk
          
        }
        
        relativeImportance[,2] <- (relativeImportance[,2] / sum( abs(relativeImportance[,2]) )) * 100

        
}

## ----------

relativeImportancePlot <- ggplot(relativeImportance[sort(relativeImportance[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
  geom_bar( aes(x= reorder(Predictor, relImportance) , y=relImportance), stat="identity", fill="black", alpha=0.5) +
  coord_flip() + theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)

relativeImportance <- data.frame( relativeImportance , significative= sapply(relativeImportance$Predictor , function(x) { relativeImportanceStats[relativeImportanceStats$Predictor == x,2] } ))
row.names(relativeImportance) <- NULL

relativeContribution <- list(dataFrame=relativeImportance,plot=relativeImportancePlot)

if(exists("modelFullBk")) { modelFull <- modelFullBk }

## -----------------------

closeAllConnections()
gc(reset=TRUE)