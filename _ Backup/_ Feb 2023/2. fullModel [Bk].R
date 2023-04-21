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

load(paste0(resultsDirectory,"/Data/","corrDistances.RData"))
load(paste0(resultsDirectory,"/Data/","speciesData.RData"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Test for correlations / collinearity

pdf(file = paste0(resultsDirectory,"/Figures/pairsPlot.pdf"), width=14, height=14 )
pairs(rasterLayers)
dev.off()

vifRasters <- vif(rasterLayers)
write.csv(vifRasters,file=paste0(resultsDirectory,"/Data/pairsVIF.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

crossValidation <- dataPartitioning(speciesData,rasterLayers,type=cvType,k=cvKFolds,distance=meanCorrDistance)
crossValidationBlocks <- crossValidation$crossValidation$blocks
crossValidationBlocksCentroids <- gCentroid(crossValidationBlocks,byid=TRUE)
crossValidationBlocksCentroids$ID <- crossValidationBlocks$folds
crossValidation <- crossValidation$crossValidation$folds

if( ! TRUE %in% sapply(1:length(crossValidation), function(x) { ( 1 %in% speciesData[crossValidation[[x]][[1]],1]) & (1 %in% speciesData[crossValidation[[x]][[2]],1]) } ) ) { next }

crossValidationPlot <- ggplot(data = worldMap) + geom_sf(color = NA, fill = "Gray") +
  geom_polygon(data = crossValidationBlocks , fill = NA, colour = "Black" , size = 0.5,  aes(long, lat, group = group)) +
  geom_text(data=as.data.frame(crossValidationBlocksCentroids), aes(label = ID, x = x, y = y), size = 3, fontface = "bold") +
  coord_sf(xlim = c(min(speciesData[,2]) - 7.5 , max(speciesData[,2]) + 7.5), ylim = c(min(speciesData[,3]) - 7.5, max(speciesData[,3]) + 7.5), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") + themeMap

pdf(file = paste0(resultsDirectory,"/Figures/crossValidationPlot.pdf"), width=14 )
print(crossValidationPlot)
dev.off()

# -----------------------

for( modelType in algorithms ) {

  # Perform full model
  
  source("Dependencies/modelMe.R")
  
  if( is.null(modelFull$model) ) { next }
  
  prediction <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
  accur <- accuracyPredicted( prediction ,speciesData ,cvIndex)
  accur <- data.frame( species=species, accur , modelFull$cvResultMatrix )

  accurThresholdBased <- accuracyPredicted( prediction ,speciesData , reclassificationIndex )
  accur$threshold <- accurThresholdBased$threshold
  accur$specificity <- accurThresholdBased$specificity
  accur$sensitivity <- accurThresholdBased$sensitivity
  accur$tss <- accurThresholdBased$tss

  rangePrediction <- c(cellStats(prediction,min),cellStats(prediction,max))

  save( modelFull , file=paste0( resultsDirectory, "/Models/",modelType,".RData"))
  save( accur , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"Performance.RData"))
  
  # Variable contribution
  
  source("Dependencies/getContribution.R")
  
  pdf(file = paste0(resultsDirectory,"/Figures/variableContribution",modelType,".pdf"), onefile=FALSE, width=10, height=7 )
  print(relativeContribution$plot)
  dev.off()
  
  relativeContributionDF <- relativeContribution$dataFrame
  save( relativeContributionDF , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"Contribution.RData"))

  # Partial Plots
  
  summaryModel <- relativeContributionDF # data.frame(Predictor=summary(modelFull$model)[,1],relImportance=summary(modelFull$model)[,2])
  summaryModel <- summaryModel[sort(summaryModel$relImportance,decreasing = TRUE,index.return=TRUE)$ix,]
  varsToPlot <- as.character(summaryModel[summaryModel$relImportance > 0,1])
  varsToPlot.n <- 0
  varsToPlotComposite <- character(0)
  varsToPlotTippingPoint <- data.frame()

  for( varsToPlot.i in varsToPlot ) {
    varsToPlot.n <- varsToPlot.n + 1
    variable.to.plot <- which(names(rasterLayers) == varsToPlot.i)
    contribution <- summaryModel[summaryModel$Predictor == varsToPlot.i,2]
    if(contribution <= 0 ) { contribution <- 0.01 }
    partialPlotVar <- modelPlot(model = modelFull$model,rasters = rasterLayers,distribution="binomial",variable.to.plot=variable.to.plot,speciesData,extrapolate = TRUE,plotExtremeThreshold= FALSE,plotRateChangeThreshold= FALSE,yLimits=rangePrediction)
    assign( paste0("partialPlotVar",varsToPlot.n) , partialPlotVar$partialPlot )
    varsToPlotTippingPoint <- rbind(varsToPlotTippingPoint,data.frame(variable=varsToPlot.i,partialPlotVar$tippingPoints))
    
    if(!is.na(partialPlotVar$tippingPoints[1])) { varsToPlotComposite <- paste0( varsToPlotComposite ,  paste0("partialPlotVar",varsToPlot.n ) , ifelse( (varsToPlot.n %% 2) != 0 , "" , " + ylab('')" ) , " , " )  }
    
  }
  
  pdf(file = paste0(resultsDirectory,"/Figures/partialPlots",modelType,".pdf"), onefile=FALSE, width=9 )
  print( eval(parse(text = paste0("grid.arrange(",substr(varsToPlotComposite,1,nchar(varsToPlotComposite)-3),", ncol=",ifelse(length(varsToPlot) > 1 , 2 , 1),")"))) )
  dev.off()

  # Export tipping points
  save(varsToPlotTippingPoint,file=paste0(resultsDirectory,"/SummaryModels/",modelType,"TippingPoints.RData"))
  
  # -------------------------------
  # Simplify / reduce model complexity
  
  if( reduceModelcomplexity ) {
    
    source("Dependencies/simplifyModel.R")
    
    modelFull <- modelReduction
    keptVariables <- modelReduction$keptVariables
    removedVariables <- modelReduction$removedVariables
    
    if( is.null(modelFull$model) ) { modelFull <- loadRData( paste0( resultsDirectory, "/Models/fullModel",modelType,".RData")) }
    
    prediction <- predictDistribution(rasterLayers,modelFull$model, reclassToOne=FALSE)
    accur <- accuracyPredicted( prediction ,speciesData , cvIndex )
    accur <- data.frame( species=species, accur  )
    
    accurThresholdBased <- accuracyPredicted( prediction ,speciesData , reclassificationIndex )
    accur$threshold <- accurThresholdBased$threshold
    accur$specificity <- accurThresholdBased$specificity
    accur$sensitivity <- accurThresholdBased$sensitivity
    accur$tss <- accurThresholdBased$tss
    
    save( modelFull , file=paste0( resultsDirectory, "/Models/reducedModel",modelType,".RData"))
    save( accur , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"ReducedPerformance.RData"))
    save( keptVariables,removedVariables , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"Reduction.RData"))
    
  }
  
}

# ---------------------------
# Ensemble variable contribution

relativeContribFiles <- list.files( paste0(resultsDirectory,"/SummaryModels/"), pattern="Contribution.RData", full.names = TRUE)
relativeImportance <- data.frame(Predictor=dataLayersName,matrix(NA,ncol=length(relativeContribFiles),nrow=length(dataLayersName)))

if( length(relativeContribFiles) > 0 ) { 
  
  for( file.i in 1:length(relativeContribFiles)) {
    relativeImportance.i <- loadRData(relativeContribFiles[file.i])
    relativeImportance[match(relativeImportance.i$Predictor,relativeImportance$Predictor),file.i+1] <- relativeImportance.i$relImportance
  }
  
  if(ncol(relativeImportance) > 2) {
    relativeImportance <- data.frame(Predictor=relativeImportance$Predictor,relImportance=apply(relativeImportance[,2:(length(relativeContribFiles)+1)],1,mean),relImportanceSD=apply(relativeImportance[,2:(length(relativeContribFiles)+1)],1,sd))
  }
  
  if( ncol(relativeImportance) == 2) {
    relativeImportance <- data.frame(Predictor=relativeImportance$Predictor,relImportance=relativeImportance[,2],relImportanceSD=0)
  }
  
  relativeImportance[which(relativeImportance$relImportance == 0), "relImportanceSD"] <- NA
  relativeImportance[is.na(relativeImportance)] <- 0
  
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
  
  relativeImportance <- data.frame( relativeImportance , significative= sapply(relativeImportance$Predictor , function(x) { ifelse(length(relativeImportanceStats[relativeImportanceStats$Predictor == x,2]) > 1 , relativeImportanceStats[relativeImportanceStats$Predictor == x,2] , NA) } ))
  row.names(relativeImportance) <- NULL
  
  if (length(relativeContribFiles) > 1 ) { relativeImportancePlot + geom_errorbar( aes(x=reorder(Predictor, relImportance), ymin=relImportance-relImportanceSD, ymax=relImportance+relImportanceSD), width=0.5, colour="Black", alpha=1, size=0.5) }
  
  pdf(file = paste0(resultsDirectory,"/Figures/variableContributionEnsemble.pdf"), onefile=FALSE, width=10, height=7 )
  print(relativeImportancePlot)
  dev.off()
  
  save( relativeImportance , file=paste0(resultsDirectory,"/SummaryModels/EnsembleContribution.RData"))

}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# End of Code