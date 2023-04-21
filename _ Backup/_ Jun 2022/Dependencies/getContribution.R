## -----------------------
closeAllConnections()
gc(reset=TRUE)
## -----------------------

if(exists("modelFull")) { modelFullBk <- modelFull }

if( varImpType == "deviance" ) {
  
  predictionModel <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
  predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
  observed <- speciesData$PA
  observed <- observed[!is.na(predicted)]
  predicted <- predicted[!is.na(predicted)]
  modelFullDE <- Dsquared(glm(predicted~observed,family=binomial()))
  
}
  
if( varImpType == "auc" ) {
  
  modelFullDE <- modelFull$auc
  
}

relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=NA)
relativeImportanceStats <- data.frame(Predictor=names(rasterLayers),relImportanceSignif=NA)

for(v in 1:length(names(rasterLayers)) ) {
  
  rasterLayersBk <- rasterLayers
  rasterLayers <- subset(rasterLayers, (1:length(names(rasterLayers)))[-v] )
  
  source("Dependencies/modelMe.R")
  
  model <- modelFull
  rasterLayers <- rasterLayersBk

  if( varImpType == "deviance" ) {
    
    predictionModel <- predictDistribution(rasterLayers,model$model, reclassToOne=FALSE)
    predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
    observed <- speciesData$PA
    observed <- observed[!is.na(predicted)]
    predicted <- predicted[!is.na(predicted)]
    partialModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
    
  }
  
  if( varImpType == "auc" ) {
    
    partialModelDE <- model$auc
    
  }
  
  relativeImportance[v,2] <- modelFullDE - partialModelDE
  relativeImportanceStats[v,2] <- wilcox.test(modelFull$cv,model$cv,alternative = c("greater"), paired = ifelse(length(modelFull$cv) == length(model$c),TRUE,FALSE))$p.value < 0.05
  
}

relativeImportance[,2] <- (relativeImportance[,2] / sum( abs(relativeImportance[,2]) )) * 100

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
