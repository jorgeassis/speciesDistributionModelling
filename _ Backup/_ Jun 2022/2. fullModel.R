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

# Read data

load(paste0(resultsDirectory,"/RData/","corrDistances.RData"))
load(paste0(resultsDirectory,"/RData/","speciesData.RData"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Test for correlations / collinearity

pdf(file = paste0(resultsDirectory,"/pairsPlot.pdf"), width=14, height=14 )
pairs(rasterLayers)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

crossValidation <- dataPartitioning(speciesData,rasterLayers,cvType,cvRemoveEdges,meanCorrDistance * cvKFolds,cvKFolds)
crossValidationCentroids <- data.frame(as.data.frame(gCentroid(crossValidationPlot,byid=TRUE)),ID=crossValidationPlot$ID)

crossValidationPlot <- ggplot(data = worldMap) + geom_sf(color = NA, fill = "Gray") +
  geom_polygon(data = crossValidationPlot , fill = NA, colour = "Black" , size = 0.2,  aes(long, lat, group = group)) +
  geom_text(data=crossValidationCentroids, aes(label = ID, x = x, y = y), size = 3) +
  coord_sf(xlim = c(min(speciesData$Lon) - 5 , max(speciesData$Lon) + 5), ylim = c(min(speciesData$Lat) - 5, max(speciesData$Lat) + 5), expand = FALSE)

pdf(file = paste0(resultsDirectory,"/crossValidationPlot.pdf"), width=14 )
print(crossValidationPlot)  
dev.off()

# -----------------------

for( modelType in algorithms ) {

  # Perform full model
  
  source("Dependencies/modelMe.R")
  
  prediction <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
  accur <- accuracyPredicted( prediction ,speciesData, cvIndex )
  accur <- data.frame( accur ,aucCV=mean(modelFull$cv),aucSD=sd(modelFull$cv))
  threshold <- accur$threshold
  
  save( modelFull , file=paste0( resultsDirectory, "/RData/fullModel",modelType,".RData"))
  save( accur , file=paste0(resultsDirectory,"/summaryModel/",modelType,"Performance.RData"))
  
  # Variable contribution
  
  varImpType <- "auc" # deviance
  source("Dependencies/getContribution.R")
  
  pdf(file = paste0(resultsDirectory,"/summaryModel/variableContributionFull",modelType,".pdf"), onefile=FALSE, width=10, height=7 )
  print(relativeContribution$plot)
  dev.off()
  
  relativeContributionDF <- relativeContribution$dataFrame
  save( relativeContributionDF , file=paste0(resultsDirectory,"/summaryModel/",modelType,"Contrib.RData"))

  # Partial Plots
  
  summaryModel <- relativeContribution$dataFrame
  summaryModel <- summaryModel[sort(summaryModel$relImportance,decreasing = TRUE,index.return=TRUE)$ix,]
  varsToPlot <- as.character(summaryModel[summaryModel$relImportance > 0,1])
  varsToPlotTippingPoint <- numeric(0)
  
  pdf(file = paste0(resultsDirectory,"/summaryModel/partialPlots",modelType,".pdf"), onefile=FALSE, width=8, height=10 )
  par(mfrow=c(round((length(names(rasterLayers)) + 1) / 2),2))
  for( var.n in sapply(summaryModel[,1],function(x){ which(names(rasterLayers) == x)}) ) {
    contribution <- summaryModel[summaryModel$Predictor == names(rasterLayers)[var.n],2]
    modelPlot(modelFull$model,rasterLayers,distribution="binomial",variable.to.plot=var.n,occurrenceRecords,print.limiting=FALSE,auto.limiting=TRUE,distribution.threshold=threshold,distribution.predicted=prediction,val.limiting=NULL,export.data.to.plot=FALSE,plot.x.lim=NULL,plot.y.lim=NULL)
    varsToPlotTippingPoint <- c(varsToPlotTippingPoint,tippingPoint)
  }
  dev.off()
  par(mfrow=c(1,1))
  
  # Export tipping points
  
  TippingPoints <- data.frame(Predictor=as.character(summaryModel[,1]),tippingPoint=varsToPlotTippingPoint)
  save(TippingPoints,file=paste0(resultsDirectory,"/summaryModel/",modelType,"TippingPoints.RData"))
  
  # Simplify model
  
  source("Dependencies/simplifyModel.R")
  
  modelFull <- modelReduction$model
  keptVariables <- modelReduction$keptVariables
  removedVariables <- modelReduction$removedVariables
  
  prediction <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
  accur <- accuracyPredicted( prediction ,speciesData, cvIndex ) # 
  accur <- data.frame( accur ,aucCV=mean(modelFull$cv),aucSD=sd(modelFull$cv))
  threshold <- accur$threshold
  
  save( modelFull , file=paste0( resultsDirectory, "/RData/reducedModel",modelType,".RData"))
  save( accur , file=paste0(resultsDirectory,"/summaryModel/",modelType,"ReducedPerformance.RData"))
  save( keptVariables,removedVariables , file=paste0(resultsDirectory,"/summaryModel/",modelType,"Reduction.RData"))
  
}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# End of Code