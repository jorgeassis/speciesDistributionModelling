
summaryModel <- relativeContributionDF
varsToPlot <- as.character(summaryModel[summaryModel$relImportance >= reduceModelcomplexityThreshold,1])
varsToPlot.n <- 0
varsToPlotComposite <- character(0)
varsToPlotTippingPoint <- data.frame()

for( varsToPlot.i in varsToPlot ) {
  
  varsToPlot.n <- varsToPlot.n + 1
  contribution <- summaryModel[summaryModel$Predictor == varsToPlot.i,2]
  partialPlotVar <- SDMmodelPlot(model = model,rasters = rasterLayers,variable.to.plot=varsToPlot.i,modelData,extrapolate = TRUE,plotExtremeThreshold= FALSE,plotRateChangeThreshold= FALSE,yLimits=NULL)
  assign( paste0("partialPlotVar",varsToPlot.n) , partialPlotVar$partialPlot )
  varsToPlotTippingPoint <- rbind(varsToPlotTippingPoint,data.frame(variable=varsToPlot.i,partialPlotVar$tippingPoints))
  if(!is.na(partialPlotVar$tippingPoints[1])) { varsToPlotComposite <- paste0( varsToPlotComposite ,  paste0("partialPlotVar",varsToPlot.n ) , ifelse( (varsToPlot.n %% 2) != 0 , "" , " + ylab('')" ) , " , " )  }
  
  partialPlotData <- partialPlotVar$partialPlotData
  save(partialPlotData,file=paste0(resultsDirectory,"/SummaryModels/",modelType,"partialPlotData",varsToPlot.i,".RData"))
  
}

partialPlots <- parse(text = paste0("grid.arrange(",substr(varsToPlotComposite,1,nchar(varsToPlotComposite)-3),", ncol=",ifelse(length(varsToPlot) > 1 , 2 , 1),")"))