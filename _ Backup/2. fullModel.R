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

cvFolds <- getBlocks(modelData,rasterLayers,corrDistance=minCorrDistance*2,k=cvKFolds)

crossValidationPlot <- ggplot() + 
                       geom_sf(data = worldMap, color = NA, fill = "Gray") +
                       geom_polygon(data=as(cvFolds$figureBlocks,"Spatial"), aes(x=long, y=lat, group=group),color="red", size=0.4, fill=NA ) +
                       geom_text(data.frame(Lon=getSpPPolygonsLabptSlots(as(cvFolds$figureBlocks,"Spatial"))[,1],Lat=getSpPPolygonsLabptSlots(as(cvFolds$figureBlocks,"Spatial"))[,2],id=cvFolds$figureBlocksK),mapping=aes(label = id, x = Lon, y = Lat), size=3) +
                       scale_color_identity() +
                       coord_sf(xlim = c(min(modelData@coords[,1]) - 7.5 , max(modelData@coords[,1]) + 7.5), ylim = c(min(modelData@coords[,2]) - 7.5, max(modelData@coords[,2]) + 7.5), expand = FALSE) +
                       xlab("Longitude") + ylab("Latitude") + themeMap + theme(legend.position = "none")

pdf(file = paste0(resultsDirectory,"/Figures/crossValidationPlot.pdf"), width=14 )
print(crossValidationPlot)
dev.off()

pdf(file = paste0(resultsDirectory,"/Figures/crossValidationPlotPart.pdf"), width=14 )
print(cvFolds$figurePartitioning)
dev.off()

# -----------------------

for( modelType in algorithms ) {

  # Perform full model
  
  source("Dependencies/modelMe.R", local = TRUE)
  save( model , file=paste0( resultsDirectory, "/Models/",modelType,".RData"))

  source("Dependencies/getContribution.R", local = TRUE)
  
  pdf(file = paste0(resultsDirectory,"/Figures/variableContribution",modelType,".pdf"), onefile=FALSE, width=10, height=7 )
  print(relativeContribution$plot)
  dev.off()
  
  relativeContributionDF <- relativeContribution$dataFrame
  save( relativeContributionDF , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"Contribution.RData"))

  source("Dependencies/getPartialPlots.R", local = TRUE)
  
  pdf(file = paste0(resultsDirectory,"/Figures/partialPlots",modelType,".pdf"), onefile=FALSE, width=9 )
  print( eval(partialPlots) )
  dev.off()

  # Export tipping points
  save(varsToPlotTippingPoint,file=paste0(resultsDirectory,"/SummaryModels/",modelType,"TippingPoints.RData"))
  
  # -------------------------------
  # Simplify / reduce model complexity
  
  if( reduceModelcomplexity ) {
    
    source("Dependencies/simplifyModel.R", local = TRUE)
    
    keptVariables <- modelReduced$keptVariables
    removedVariables <- modelReduced$removedVariables
    model <- modelReduced$model
    
    save( model , file=paste0( resultsDirectory, "/Models/reducedModel",modelType,".RData"))
    save( keptVariables,removedVariables , file=paste0(resultsDirectory,"/SummaryModels/",modelType,"ReductionPredictors.RData"))
    
  }

  # -------------------------------
  
}

# ---------------------------
# Ensemble variable contribution

relativeContribFiles <- list.files( paste0(resultsDirectory,"/SummaryModels/"), pattern="Contribution.RData", full.names = TRUE)
relativeContribFiles <- relativeContribFiles[!grepl("Ensemble",relativeContribFiles)]
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
  
  relativeImportancePlot <- ggplot(data=relativeImportance[sort(relativeImportance[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
    geom_bar( aes(x= reorder(Predictor, relImportance) , y=relImportance), stat="identity", fill="#ffd45b", alpha=0.85) +
    coord_flip() + 
    theme_bw() +
    theme(
      axis.text=element_text(size=10),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=13.5),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=13.5),
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#EFEFEF"), 
      panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
    ) + labs(x = "Predictor") + 
    labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
    annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)
  
  if (length(relativeContribFiles) > 1 ) { relativeImportancePlot <- relativeImportancePlot + geom_errorbar( aes(x= reorder(Predictor, relImportance), ymin=relImportance-ifelse(relImportanceSD == 0, NA, relImportanceSD), ymax=relImportance+ifelse(relImportanceSD == 0, NA, relImportanceSD)), width=.1, alpha=0.8, size=0.4, position=position_dodge(.9)) }
  
  pdf(file = paste0(resultsDirectory,"/Figures/variableContributionEnsemble.pdf"), onefile=FALSE, width=10, height=7 )
  print(relativeImportancePlot)
  dev.off()
  
  save( relativeImportance , file=paste0(resultsDirectory,"/SummaryModels/ensembleContribution.RData"))

}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# End of Code