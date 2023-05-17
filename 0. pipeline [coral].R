## ------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------
##
##  Machine Learning Species Distribution Modelling [ Ver.301 ]
##
## ------------------------------------------------------------------

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)

mainfunctionsFile <- "Dependencies/mainFunctions.R"
mainConfigFile <- "0. config [coral].R"

# credentials::set_github_pat()
# exportRequirements()
# exportPackages()

## ------------------------------------------------------------------
## ------------------------------------------------------------------

source(mainfunctionsFile)
source(mainConfigFile)

## ---------

## For extent calculations

if( ! is.null(bathymetryDataLayerHR) ) {
  
  bathymetryDataLayerHR <- raster(bathymetryDataLayerHR)
  bathymetryDataLayerHR <- aggregate(bathymetryDataLayerHR,2,mean,na.rm=T)
  bathymetryDataLayerHR[bathymetryDataLayerHR <= maxDepth * (-1)] <- NA
  bathymetryDataLayerHR[bathymetryDataLayerHR >= minDepth * (-1)] <- NA
  bathymetryDataLayerHR[!is.na(bathymetryDataLayerHR)] <- 1
  
  bathymetryLayerFraction <- aggregate(bathymetryDataLayerHR,res(raster(bathymetryDataLayer))[1] / res(bathymetryDataLayerHR)[1] ,sum,na.rm=T)
  bathymetryLayerFraction <- bathymetryLayerFraction / cellStats(bathymetryLayerFraction, max, na.rm=T)
  bathymetryLayerFraction <- raster::resample(bathymetryLayerFraction,raster(bathymetryDataLayer), method="ngb")
  
}

## ---------

dataRecords <- read.csv(dataRecordsFile)
dataRecords <- dataRecords[complete.cases(dataRecords),]
colnames(dataRecords) <- c("speciesName","Lon","Lat")
speciesList <- sort(unique(dataRecords$speciesName))

## ------------------------------------------------------------------
## Model species

overwrite <- FALSE
if( ! overwrite) { 
  predictedSpecies <- character(0)
  for( species in list.files(mainResultsDirectory) ) {
    if( sum(sapply(algorithms, function(x) { length(list.files( paste0(mainResultsDirectory,"/",species,"/Models/"), pattern = x, recursive=T)) } )) >= length(algorithms) )  { predictedSpecies <- c(predictedSpecies,species) }
  }
  speciesList <- speciesList[ ! speciesList %in% predictedSpecies] 
}
length(speciesList)

cl <- parallel::makeCluster(nCores)
registerDoParallel(cl)

modelParallel <- foreach(species = speciesList, .export="monotonicity") %dopar% {
  
  source(mainfunctionsFile, local = TRUE)
  source(mainConfigFile, local = TRUE)
  
  ## ---------------------
  
  resultsDirectory <- paste0(mainResultsDirectory,"/",species,"/")
  
  ## ---------------------
  
  if( ! dir.exists(resultsDirectory) ) { dir.create(resultsDirectory, recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/Models/")) ) { dir.create(paste0(resultsDirectory,"/Models/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/SummaryModels/")) ) { dir.create(paste0(resultsDirectory,"/SummaryModels/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/Predictions/")) ) { dir.create(paste0(resultsDirectory,"/Predictions/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/Figures/")) ) { dir.create(paste0(resultsDirectory,"/Figures/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/Data/")) ) { dir.create(paste0(resultsDirectory,"/Data/"), recursive = TRUE) }
  
  ## ---------------------
  
  occurrenceRecords <- dataRecords[which(dataRecords$speciesName == species),]
  occurrenceRecords <- data.frame(Lon=occurrenceRecords[,dataRecordsNames[1]],Lat=occurrenceRecords[,dataRecordsNames[2]])
  # plotMap(occurrenceRecords,4,"black")
  
  ## ---------------------
  
  if( ! is.null(depthTraits) ) {
    
    depthTraitsData <- read.csv(depthTraits)
    depthTraitsData <- depthTraitsData[depthTraitsData$speciesName == species,]
    
    minDepth <- mean(depthTraitsData$DepthMin,na.rm=T)
    maxDepth <- mean(depthTraitsData$DepthMax,na.rm=T)
    
    if( is.na(minDepth) | is.na(maxDepth) ) {
      
      depthTraitsData <- read.csv(depthTraits)
      depthTraitsData <- depthTraitsData[which(grepl(strsplit(species, " ")[[1]][1],depthTraitsData$speciesName)),]
      minDepth <- mean(depthTraitsData$DepthMin,na.rm=T)
      maxDepth <- mean(depthTraitsData$DepthMax,na.rm=T)
      
    }
    
    if( is.na(minDepth) | is.na(maxDepth) ) { unlink(resultsDirectory, recursive=TRUE); occurrenceRecords <- data.frame() }

    if( ! is.na(minDepth) & ! is.na(maxDepth) ) { 
      if( abs(diff(c(minDepth,maxDepth))) < 25 ) { minDepth <- ifelse((minDepth - 25) < 0 , 0 , minDepth - 25); maxDepth <- maxDepth + maxDepth + 25 }
    }
    
  }
  
  ## ---------------------
  # Environmental layers
  
  if( nrow(occurrenceRecords) >= minOccurrenceRecords ) { 
    
    rasterLayers <- list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
    rasterLayers <- rasterLayers[!grepl("@",rasterLayers)]
    rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
    rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer, minDepth=ifelse( is.null(minDepth), "NULL" , ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer)) , maxDepth=ifelse( is.null(maxDepth), "NULL" , maxDepth+maxDepthBuffer),intertidal)
    names(rasterLayers) <- dataLayersName
    rasterLayers <- correctLayer(rasterLayers,"Salinity","he",28,28)
    
    worldMap <- ne_countries(scale = "medium", returnclass = "sf")
    
    source("1. prepareRecords.R", local = TRUE)
    
  }
  
  ## ---------------------
  
  if( nrow(occurrenceRecords) >= minOccurrenceRecords ) { 
    
    source("2. fullModel.R", local = TRUE)
    
  }
  
  ## ---------------------
  
  gc(reset=TRUE)
  return(NULL)
  
  ## ---------------------
  
}

stopCluster(cl); rm(cl)
closeAllConnections()
gc(reset=TRUE)

## ------------------------------------------------------------------
## Predict species

overwrite <- TRUE
speciesList <- sort(unique(dataRecords$speciesName))
speciesToPredict <- character(0)
for( species in list.files(mainResultsDirectory) ) {
  if( sum(sapply(algorithms, function(x) { length(list.files( paste0(mainResultsDirectory,"/",species,"/Models/"), pattern = paste0(x,".RData"), recursive=T)) } )) == length(algorithms) )  { speciesToPredict <- c(speciesToPredict,species) }
  if( ! overwrite) { 
    if( sum(sapply(c("Baseline","Reclass","Global"), function(x) { length(list.files( paste0(mainResultsDirectory,"/",species,"/Predictions/",scenariosToPredict[length(scenariosToPredict)],"/"), pattern = x, recursive=T)) > 1 } )) >= 2 )  { speciesToPredict <- speciesToPredict[which(speciesToPredict != species)] }
  }
}
length(speciesToPredict)

cl <- parallel::makeCluster(nCores)
registerDoParallel(cl)

predictParallel <- foreach(species = speciesToPredict ) %dopar% {
  
  source(mainfunctionsFile, local = TRUE)
  source(mainConfigFile, local = TRUE)
  
  ## ---------------------
  
  resultsDirectory <- paste0(mainResultsDirectory,"/",species,"/")
  modelData <- loadRData(paste0(resultsDirectory,"/Data/","modelData.RData"))
  
  if( ! is.null(depthTraits) ) {
    depthTraitsData <- read.csv(depthTraits)
    depthTraitsData <- depthTraitsData[depthTraitsData$speciesName == species,]
    minDepth <- mean(depthTraitsData$DepthMin)
    maxDepth <- mean(depthTraitsData$DepthMax)
    if( abs(diff(c(minDepth,maxDepth))) < 25 ) { minDepth <- ifelse((minDepth - 25) < 0 , 0 , minDepth - 25); maxDepth <- maxDepth + maxDepth + 25 }
  }
  
  ## ---------------------
  
  source("3. modelPredict.R", local = TRUE)
  
  ## ---------------------
  
  gc(reset=TRUE)
  return(NULL)
  
}

stopCluster(cl); rm(cl)
closeAllConnections()
gc(reset=TRUE)

## ------------------------------------------------------------------
## Estimate range shifts

speciesList <- sort(unique(dataRecords$speciesName))
speciesToPredict <- character(0)
for( species in list.files(mainResultsDirectory) ) {
  if( overwrite & sum(sapply(c("Baseline","Reclass","Global"), function(x) { length(list.files( paste0(mainResultsDirectory,"/",species,"/Predictions/",scenariosToPredict[length(scenariosToPredict)],"/"), pattern = x, recursive=T)) > 1 } )) >= 2 )  { speciesToPredict <- c(speciesToPredict,species) }
  if( ! overwrite & sum(sapply(c("Baseline","Reclass","Global"), function(x) { length(list.files( paste0(mainResultsDirectory,"/",species,"/Predictions/",scenariosToPredict[length(scenariosToPredict)],"/"), pattern = x, recursive=T)) > 1 } )) >= 2 )  { speciesToPredict <- c(speciesToPredict,species) }
}
length(speciesToPredict)

cl <- parallel::makeCluster(nCores)
registerDoParallel(cl)

estimateParallel <- foreach(species = speciesToPredict) %dopar% {
  
  source(mainfunctionsFile, local = TRUE)
  source(mainConfigFile, local = TRUE)
  
  ## ---------------------
  
  resultsDirectory <- paste0(mainResultsDirectory,"/",species,"/")
  modelData <- loadRData(paste0(resultsDirectory,"/Data/","modelData.RData"))
  occurrenceRecords <- modelData@coords[modelData@pa == 1,]
  
  if( ! is.null(depthTraits) ) {
    
    depthTraitsData <- read.csv(depthTraits)
    depthTraitsData <- depthTraitsData[depthTraitsData$speciesName == species,]
    minDepth <- mean(depthTraitsData$DepthMin)
    maxDepth <- mean(depthTraitsData$DepthMax)
    
    if( abs(diff(c(minDepth,maxDepth))) < 25 ) { minDepth <- ifelse((minDepth - 25) < 0 , 0 , minDepth - 25); maxDepth <- maxDepth + maxDepth + 25 }
    
  }
  
  ## ---------------------
  
  source("4. rangeShiftEstimates.R", local = TRUE)
  
  ## ---------------------
  
  return(NULL)
  
}

stopCluster(cl); rm(cl)
closeAllConnections()
gc(reset=TRUE)


## -----------------------------------------
## -----------------------------------------
# Summary predictions species

if(! dir.exists(paste0(stackResultsFolder))) { dir.create(paste0(stackResultsFolder), recursive = T) }
if(! dir.exists(paste0(stackResultsFolder,"/Summary"))) { dir.create(paste0(stackResultsFolder,"/Summary"), recursive = T) }
if(! dir.exists(paste0(stackResultsFolder,"/Maps"))) { dir.create(paste0(stackResultsFolder,"/Maps"), recursive = T) }
if(! dir.exists(paste0(stackResultsFolder,"/Maps/perSpecies"))) { dir.create(paste0(stackResultsFolder,"/Maps/perSpecies"), recursive = T) }

## ----------

speciesPredicted <- character(0)
for( species in list.files(mainResultsDirectory) ) {
  if( length(list.files( paste0(mainResultsDirectory,"/",species,"/SummaryModels/"), pattern = "rangeShiftEstimates" , recursive=T)) > 0 )  { speciesPredicted <- c(speciesPredicted,species) }
}
length(speciesPredicted)

## ----------

relativeContributionAll <- data.frame()
performanceAll <- data.frame()
tippingPointsAll <- data.frame()
rangeShiftsEstimatesAll <- data.frame()
depthRangeAll <- data.frame()

for( species in speciesPredicted ) {
  
  ## ---------------------
  
  cat("\014")
  cat("## --------------------------------- \n")
  cat(species," | ",which(speciesPredicted == species),"out of",length(speciesPredicted),"\n")
  cat("## --------------------------------- \n")
  cat("\n")
  
  ## ---------------------
  
  dataDirectory <- paste0(mainResultsDirectory,"/",species,"/")
  
  ## ---------------------
  # Does not have model or has been modelled
  
  if( ! TRUE %in% grepl("Contrib",list.files(dataDirectory, recursive = T)) ) { next }
  if( ! TRUE %in% grepl("Performance",list.files(dataDirectory, recursive = T)) ) { next }
  if( ! TRUE %in% grepl("TippingPoints",list.files(dataDirectory, recursive = T)) ) { next }
  if( ! TRUE %in% grepl("rangeShiftEstimates", list.files(dataDirectory, recursive = T)) ) { next }
  
  ## ---------------------
  
  source("6. summaryModels.R")
  
  ## ---------------------
  
  relativeContributionAll <- rbind.fill(relativeContributionAll,contributionDF)
  performanceAll <- rbind.fill(performanceAll,performanceDF)
  depthRangeAll <- rbind.fill(depthRangeAll,depthRange)
  tippingPointsAll <- rbind.fill(tippingPointsAll,tippingPointsDF)
  rangeShiftsEstimatesAll <- rbind.fill(rangeShiftsEstimatesAll,rangeShiftsEstimatesDF)
  
}

relativeContributionAll[is.na(relativeContributionAll)] <- 0
rangeShiftsEstimatesAll[ rangeShiftsEstimatesAll == -Inf] <- NA
rangeShiftsEstimatesAll[ rangeShiftsEstimatesAll == Inf] <- NA

## ----------

names(performanceAll)
speciesDiscardedNames <- speciesPredicted[which(performanceAll$auc.ensemble < 0.75)]
speciesPredicted <- speciesPredicted[which(performanceAll$auc.ensemble >= 0.75)]

write.csv(data.frame(species=speciesDiscardedNames), file=paste0(stackResultsFolder,"/speciesDiscarded.csv"), row.names = FALSE)
write.csv(data.frame(species=speciesPredicted), file=paste0(stackResultsFolder,"/speciesPredictedList.csv"), row.names = FALSE)

## ---------

taxaCrossChecker <- "Alismatales"
source("Dependencies/getWormsInfo.R")
write.csv(speciesListWorms, file=paste0(stackResultsFolder,"/speciesListWorms.csv"), row.names = FALSE)

## ----------

write.csv(relativeContributionAll[which(relativeContributionAll$name %in% speciesPredicted),], file=paste0(stackResultsFolder,"/relativeContribution.csv"), row.names = FALSE)
write.csv(performanceAll[which(performanceAll$name %in% speciesPredicted),], file=paste0(stackResultsFolder,"/performance.csv"), row.names = FALSE)
write.csv(tippingPointsAll[which(tippingPointsAll$name %in% speciesPredicted),], file=paste0(stackResultsFolder,"/tippingPoints.csv"), row.names = FALSE)
write.csv(depthRangeAll[which(rangeShiftsEstimatesAll$name %in% speciesPredicted),], file=paste0(stackResultsFolder,"/depthRanges.csv"), row.names = FALSE)
write.csv(rangeShiftsEstimatesAll[which(rangeShiftsEstimatesAll$name %in% speciesPredicted),], file=paste0(stackResultsFolder,"/rangeShiftsEstimates.csv"), row.names = FALSE)

## -----------------------------------------
## -----------------------------------------
# Summary of range extent and overlap between taxa (e.g., species)

speciesPredicted <- read.csv(paste0(stackResultsFolder,"/speciesPredictedList.csv"))[,1]
taxaLevel <- "Genus" # Species Genus

source("5. summaryRangeExtentOverlap.R")

## -----------------------------------------
## -----------------------------------------
# Summary presence-absence of species in polygon-Regions

resultsName <- "MarineRealm" # MarineEcoRegion MarineRealm
polygonPath <- "Dependencies/Data/Shapefiles/marine_ecoregions.shp"
polygonFeature <- "REALM" # ECOREGION PROVINCE REALM

resultsName <- "EEZ"
polygon <- load("Dependences/Data/Shapefiles/EEZGlobal.RData")
polygonFeature <- "EEZ"

source("6. summaryPresenceInsideRegion.R")

## -----------------------------------------
## -----------------------------------------
# Stack predictions species

source("7. stackingMultiplePredictions.R")

## -----------------------------------------
## -----------------------------------------
# Summary stack predictions

resultsName <- "Global" # MarineRealm MarineEcoRegion Global
polygon <- NULL
polygonFeature <- NULL # NULL ECOREGION PROVINCE REALM

resultsName <- "MarineRealm" # MarineRealm MarineEcoRegion
polygonPath <- "Dependencies/Data/Shapefiles/marine_ecoregions.shp"
polygonFeature <- "REALM" # NULL ECOREGION PROVINCE REALM

resultsName <- "EEZ" # MarineRealm MarineEcoRegion
polygon <- load("/Volumes/StingRay/Dropbox/Data/Spatial information/Shapefiles/EEZ/EEZGlobal.RData")
polygonFeature <- "EEZ" # EEZ

source("8. summaryStackingMultiplePredictions.R")

## ---------------
# Per latitude

source("8. summaryStackingMultiplePredictionsLat.R")

resultsDFTrunc <- resultsDF[resultsDF$binFrom >= -60 & resultsDF$binFrom >= -90 , ]
resultsDFTrunc[resultsDFTrunc == 0] <- NA

fig1 <- ggplot() +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessBaselineMean) , color="#575757", size=0.4) +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessssp119Mean) , color="#6EA5DC", size=0.4) +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessssp585Mean) , color="#DC6E6E", size=0.4) +
  theme_minimal() +
  xlab("Latitude (degree)") + 
  ylab("Species richness (average number of species)") + geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.7) + coord_flip()

pdf(paste0(stackResultsFolder,"/Summary","/speciesRichnessLatitude",typePrediction,".pdf"), width = 8, height=16 )
print(fig1)
dev.off()

fig1 <- ggplot() +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessBaselineArea) , color="#575757", size=0.4) +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessssp119Area) , color="#6EA5DC", size=0.4) +
  geom_line( data=resultsDFTrunc,aes(x=binFrom, y=speciesRichnessssp585Area) , color="#DC6E6E", size=0.4) +
  theme_minimal() +
  xlab("Latitude (degree)") + 
  ylab("Extent of suitable habitats (total area in km2)") + geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.7) + coord_flip()

pdf(paste0(stackResultsFolder,"/Summary","/areaLatitude",typePrediction,".pdf"), width = 8, height=16 )
print(fig1)
dev.off()

## ------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------
# End of Code