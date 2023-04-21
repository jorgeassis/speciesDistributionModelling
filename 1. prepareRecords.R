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

# Rasterize occurrence records
occurrenceRecords <- data.frame(Lon=as.numeric(as.character(occurrenceRecords[,1])),Lat=as.numeric(as.character(occurrenceRecords[,2])))
occurrenceRecords <- occurrenceRecords[which(!duplicated(occurrenceRecords[,1])),]
occurrenceRecords <- xyFromCell(subset(rasterLayers,1),unique(cellFromXY(subset(rasterLayers,1),occurrenceRecords)))
colnames(occurrenceRecords) <- dataRecordsNames

# Determine study extent
# optimalRegionExtent <- getBackgroundExtent(occurrenceRecords,rasterLayers,regionBuffer, bufferStep=1 )
# optimalRegionExtentDist <- optimalRegionExtent$bestResult
# pdf(file = paste0(resultsDirectory,"/Figures/optimalRegionExtent.pdf"), width=12, height=8 )
# print(optimalRegionExtent$figure)
# dev.off()

# rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer=optimalRegionExtentDist, minDepth=ifelse( is.null(minDepth), "NULL" , ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer)) , maxDepth=ifelse( is.null(maxDepth), "NULL" , maxDepth+maxDepthBuffer),intertidal)
# names(rasterLayers) <- dataLayersName
# rasterLayers <- correctLayer(rasterLayers,"Salinity","he",28,28)
# rasterLayers <- stack(rasterLayers)

tryCatch( rasterLayers <- subset(rasterLayers,which(cellStats(rasterLayers,var) > 0)) , error = function(e) e <<- 0 )

# Determine autocorrelation
spatialAutocorr <- spatialAutocorrelation(rasterLayers)
minCorrDistance <- spatialAutocorr$minDistance
meanCorrDistance <- spatialAutocorr$meanDistance
medianCorrDistance <- spatialAutocorr$medianDistance

pdf(file = paste0(resultsDirectory,"/Figures/spatialAutocorrelation.pdf"), width=12, height=8 )
print(spatialAutocorr$figure)
dev.off()
save(spatialAutocorr,file=paste0(resultsDirectory,"/Data/","spatialAutocorr.RData"))

# Relocate records to closest cell (those falling on land / unlikely depth)
occurrenceRecords <- relocateNACoords(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth)

if( nrow(occurrenceRecords) < minOccurrenceRecords ) { next }

# Generate Pseudo-absences
pseudoAbsences <- generatePseudoAbsences(occurrenceRecords,rasterLayers,paType)
biasSurface <- pseudoAbsences$biasSurface
pseudoAbsences <- pseudoAbsences$records

# Spatial thinning
occurrenceRecords <- spatialThinning(occurrenceRecords,min(minCorrDistance,100),verbose=FALSE)
pseudoAbsences <- spatialThinning(pseudoAbsences,min(minCorrDistance,100),verbose=FALSE)
pseudoAbsences <- resampleRecords(pseudoAbsences,rasterLayers,max(paMinimum,nrow(occurrenceRecords)),EnvironmentStrat=TRUE,biasSurface)

if( nrow(occurrenceRecords) < minOccurrenceRecords ) { next }

# Drop no var Layers
rasterLayers <- subset(rasterLayers, which(apply(raster::extract(rasterLayers,rbind(occurrenceRecords,pseudoAbsences)),2, FUN= function(x) { length(unique(x)) } ) != 1))

# Generate modelData
modelData <- prepareModelData(occurrenceRecords, pseudoAbsences, rasterLayers)
save(modelData,file=paste0(resultsDirectory,"/Data/","modelData.RData"))

# Test for correlations / collinearity
pairsCorr <- corVar(modelData,method = "pearson")
write.csv(pairsCorr,file=paste0(resultsDirectory,"/Data/pairsCorrPearson.csv"), row.names = FALSE)

vifRasters <- vif(rasterLayers)
write.csv(vifRasters,file=paste0(resultsDirectory,"/Data/pairsVIF.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Maps and figures

sf_use_s2(FALSE)

occurrenceRecordsMapAbs <- ggplot() + geom_sf(data = worldMap, color = NA, fill = "Gray") +
  geom_point(data = modelData@coords[modelData@pa == 0,], aes(x = X, y =Y), size = 1, color="#707070") +
  geom_point(data = modelData@coords[modelData@pa == 1,], aes(x = X, y =Y), size = 1, color="#AC1919") + 
  xlab("Longitude") + ylab("Latitude") + themeMap +
  coord_sf(xlim = c(min(modelData@coords[,1]) - 5 , max(modelData@coords[,1]) + 5), ylim = c(min(modelData@coords[,2]) - 5, max(modelData@coords[,2]) + 5), expand = FALSE)

pdf(file = paste0(resultsDirectory,"/Figures/occurrenceRecordsAbsencesMap.pdf"), width=12 )
print(occurrenceRecordsMapAbs)
dev.off()

# --------------------------

occurrenceRecordsDepths <- data.frame(Depth=abs( raster::extract(raster(bathymetryDataLayer),occurrenceRecords) ))
occurrenceRecordsDepthsPlot <- ggplot(occurrenceRecordsDepths, aes(x=Depth)) +
  geom_histogram(color="black", fill="#CACACA", bins=ifelse(nrow(occurrenceRecordsDepths) > 50,50,nrow(occurrenceRecordsDepths)), size=0.25) + 
  geom_vline(aes(xintercept=quantile(Depth,probs=0.95,na.rm=T)  ),color="Black", linetype="dashed", size=0.3) + 
  xlab("Depth (m)") + ylab("Number of records")

occurrenceRecordsDepthsPlot <- occurrenceRecordsDepthsPlot + 
  annotate("text", y = max(ggplot_build(occurrenceRecordsDepthsPlot)$data[[1]]$count) , x = 10 + quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T) , label = paste0(round(quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T),digits=2)," m"), hjust = 0)

pdf(file = paste0(resultsDirectory,"/Figures/occurrenceRecordsDepthsPlot.pdf"), width=12, height=8 )
print(occurrenceRecordsDepthsPlot)
dev.off()

depthData <- data.frame(minKnownDepth=minDepth,
                        maxKnownDepth=maxDepth,
                        minObservedDepth=min(occurrenceRecordsDepths,na.rm=T),
                        maxObservedDepth=max(occurrenceRecordsDepths,na.rm=T),
                        q95ObservedDepth=quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T))

save(depthData,file=paste0(resultsDirectory,"/SummaryModels/","depthData.RData"))
gc(reset=TRUE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
