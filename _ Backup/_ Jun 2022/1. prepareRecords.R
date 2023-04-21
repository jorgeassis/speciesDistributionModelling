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

# Relocate records to closest cell (those falling on land / unlikely depth)
occurrenceRecords <- relocateNACoords(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth)

# Rasterize records
occurrenceRecords <- xyFromCell(subset(rasterLayers,1),unique(cellFromXY(subset(rasterLayers,1),occurrenceRecords)))
colnames(occurrenceRecords) <- c("Lon","Lat")

# Generate Pseudo-absence
absences <- pseudoAbsences(occurrenceRecords,rasterLayers,patchContinuity=FALSE,polygonPA,polygonPAType,paMindistance,paType,paBiasKernelSurface,paBiasKernelSurfaceProb,paRatio,paEnvironmentStratification,plotFigure=FALSE)

# Determine autocorrelation 
distanceUncorr <- spatialAutocorrelation(records=rbind(occurrenceRecords,absences),rasterLayers,autocorrelationClassDistance,autocorrelationSubsetRecords,autocorrelationMaxDistance,autocorrelationSignif)

pdf(file = paste0(resultsDirectory,"/spatialAutocorrelation.pdf"), width=12, height=8 )
print(distanceUncorr$figure)
dev.off()

meanCorrDistance <- distanceUncorr$distance
save(meanCorrDistance,file=paste0(resultsDirectory,"/RData/","corrDistances.RData"))

## ------------------------

occurrenceRecords <- spatialThinning(occurrenceRecords,rasterLayers,ifelse(meanCorrDistance <= 20 , meanCorrDistance , 20),verbose=FALSE)
absences <- spatialThinning(absences,rasterLayers,ifelse(meanCorrDistance <= 20 , meanCorrDistance , 20),verbose=FALSE)

# Generate speciesData
speciesData <- data.frame(PA=c(rep(1,nrow(occurrenceRecords)),rep(0,nrow(absences))) , rbind(occurrenceRecords,absences)  )
speciesData <- speciesData[sample(1:nrow(speciesData),nrow(speciesData), replace = FALSE),]

# Save objects
save(meanCorrDistance,file=paste0(resultsDirectory,"/RData/","corrDistances.RData"))
save(speciesData,file=paste0(resultsDirectory,"/RData/","speciesData.RData"))
write.csv(speciesData,paste0(resultsDirectory,"/RData/","speciesData.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Maps and figures

sf_use_s2(FALSE)

occurrenceRecordsMap <- ggplot(data = worldMap) + geom_sf(color = NA, fill = "Gray") +
  geom_point(data = occurrenceRecords, aes(x = Lon, y = Lat), size = 0.75, color="Black") + 
  xlab("Longitude") + ylab("Latitude") + themeMap +
  coord_sf(xlim = c(min(speciesData$Lon) - 5 , max(speciesData$Lon) + 5), ylim = c(min(speciesData$Lat) - 5, max(speciesData$Lat) + 5), expand = FALSE)

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsMap.pdf"), width=12 )
print(occurrenceRecordsMap)
dev.off()

occurrenceRecordsMapAbs <- occurrenceRecordsMap + 
  geom_point(data = absences, aes(x = Lon, y = Lat), size = 0.25, color="#656565")

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsAbsencesMap.pdf"), width=12 )
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

pdf(file = paste0(resultsDirectory,"/occurrenceRecordsDepthsPlot.pdf"), width=12, height=8 )
print(occurrenceRecordsDepthsPlot)
dev.off()

depthData <- data.frame(minKnownDepth=minDepth,
                        maxKnownDepth=maxDepth,
                        minObservedDepth=min(occurrenceRecordsDepths,na.rm=T),
                        maxObservedDepth=max(occurrenceRecordsDepths,na.rm=T),
                        q95ObservedDepth=quantile(occurrenceRecordsDepths$Depth,probs=0.95,na.rm=T))

save(depthData,file=paste0(resultsDirectory,"/summaryModel/","depthData.RData"))

gc()

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------