## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

if(! dir.exists(paste0(stackResultsFolder,"/Summary"))) { dir.create(paste0(stackResultsFolder,"/Summary"), recursive = T) }

## ----------------
## Get layers

layersToCalc <- list.files(stackResultsFolder, full.names = T, recursive = T, pattern = ".RData")
layersToCalcNames <- list.files(stackResultsFolder, full.names = F, recursive = T, pattern = ".RData")
layersToCalcNames <- gsub("\\.RData","",layersToCalcNames)
layersToCalcNames <- gsub("Maps/","",layersToCalcNames)
if(typePrediction == "Reachable") { 
  layersToCalc <- layersToCalc[grepl("Reachable",layersToCalc)] 
  layersToCalcNames <- layersToCalcNames[grepl("Reachable",layersToCalcNames)] 
}
if(typePrediction == "unConstrained") { 
  layersToCalc <- layersToCalc[!grepl("Reachable",layersToCalc)] 
  layersToCalcNames <- layersToCalcNames[!grepl("Reachable",layersToCalcNames)] 
}
layersToCalc <- layersToCalc[ sort(unique(unlist(sapply(scenariosToPredict, function(x) { which(grepl(x,layersToCalc)) } )))) ] 
layersToCalcNames <- layersToCalcNames[ sort(unique(unlist(sapply(scenariosToPredict, function(x) { which(grepl(x,layersToCalcNames)) } )))) ] 
layersToCalc <- layersToCalc[!grepl("speciesRichnessChanges",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("speciesRichnessChanges",layersToCalcNames)] 
layersToCalc <- layersToCalc[grepl("speciesRichness",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[grepl("speciesRichness",layersToCalcNames)] 

presentDayLayer <- layersToCalc[grepl("speciesRichnessBaseline",layersToCalc)] 
presentDayLayer <- presentDayLayer[!grepl("Uncertainty",presentDayLayer)] 

layersToCalc <- layersToCalc[!grepl("speciesRichnessBaseline",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("speciesRichnessBaseline",layersToCalcNames)] 
layersToCalc <- layersToCalc[!grepl("Uncertainty",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("Uncertainty",layersToCalcNames)] 

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

areaLayer <- raster::area(loadRData(presentDayLayer)) * bathymetryLayerFraction
presentDayLayer <- loadRData(presentDayLayer)
presentDayLayer[presentDayLayer == 0] <- NA
presentDayLayer[presentDayLayer > 1] <- 1
presentDayLayerArea <- presentDayLayer * areaLayer

if( resultsName == "Global" ) { 
  
  resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*2)+2,nrow=1))
  colnames(resultsDF) <- c("Region","Baseline",paste0(layersToCalcNames,"Area"),paste0(layersToCalcNames,"AreaChange"))
  resultsDF[,1] <- "Global"
  polygonFeatureNames <- "Global"
  polygon <- extent(areaLayer)
  polygon <- as(polygon, 'SpatialPolygons')  
  
}

if( resultsName != "Global" ) {
  
  if(class(polygonPath) == "character") { polygon <- shapefile(polygonPath) }
  
  if( ! is.null(polygonFeature) ) { 
    
    polygonFeatureNames <- polygon@data[which(names(polygon@data) == polygonFeature)][,polygonFeature] 
    polygon <- gUnaryUnion(polygon, id = polygonFeatureNames, checkValidity=2L)
    polygonFeatureNames <- sapply(1:length(polygon), function(x) slot(slot(polygon, "polygons")[[x]],"ID") )
    polygon$names <- polygonFeatureNames
    
  }
  
  resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*2)+2,nrow=length(unique(polygonFeatureNames))))
  colnames(resultsDF) <- c("Region","Baseline",paste0(layersToCalcNames,"Area"),paste0(layersToCalcNames,"AreaChange"))
  resultsDF[,1] <- unique(polygonFeatureNames)
  
}

## ------------

baselineAreas <- extract(presentDayLayerArea,polygon, small=TRUE)

for( i in 1:length(baselineAreas) ) { resultsDF[i,2] <- sum(baselineAreas[[i]],na.rm=T) }

## ------------

for( i in 1:length(layersToCalc) ) {
  
  rasterLayer <- loadRData(layersToCalc[i])
  rasterLayer[rasterLayer == 0] <- NA
  
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer
  
  for( p.i in 1:length(polygonFeatureNames) ){
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(i,"out of",length(layersToCalcNames),"\n")
    cat(p.i,"out of",length(polygonFeatureNames),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    polygon.pi <- polygon[p.i,]
    
    rasterLayer.pi <- crop(rasterLayer,polygon.pi)
    rasterLayer.pi <- raster::mask(rasterLayer.pi,polygon.pi)
    rasterLayerArea.pi <- crop(rasterLayerArea,rasterLayer.pi)
    rasterLayerArea.pi <- raster::mask(rasterLayerArea.pi,rasterLayer.pi)
    
    rasterLayerVals <- unlist(getValues(rasterLayer.pi))
    rasterLayerAreaVals <- unlist(getValues(rasterLayerArea.pi))
    
    if( ! is.null(rasterLayerAreaVals) ) { 
      
      resultsDF[p.i,i+2] <- sum( rasterLayerAreaVals ,na.rm=T )
      resultsDF[p.i,i+2+(length(layersToCalcNames))] <- (( sum( rasterLayerAreaVals , na.rm=T) - resultsDF[p.i,2] )  / resultsDF[p.i,2] ) * 100
      
    }
  }
  
}

resultsDF[resultsDF == -Inf] <- NA
resultsDF[resultsDF == Inf] <- NA

if( resultsName == "MarineEcoRegion" ) {
  
  polygon <- shapefile(polygonPath)
  
  resultsDF <- data.frame(realmRegion=polygon$REALM[match(resultsDF$Region,polygon$ECOREGION)], 
                          provinceRegion=polygon$PROVINCE[match(resultsDF$Region,polygon$ECOREGION)],
                          resultsDF)
  
}

names(resultsDF) <- gsub("Reachable","",names(resultsDF))
names(resultsDF) <- gsub("speciesRichness","",names(resultsDF))

write.csv(resultsDF, file=paste0(stackResultsFolder,"/Summary","/summaryExtent",resultsName,typePrediction,".csv"), row.names = FALSE)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

layersToCalc <- list.files(stackResultsFolder, full.names = T, recursive = T, pattern = ".RData")
layersToCalcNames <- list.files(stackResultsFolder, full.names = F, recursive = T, pattern = ".RData")
layersToCalcNames <- gsub("\\.RData","",layersToCalcNames)
layersToCalcNames <- gsub("Maps/","",layersToCalcNames)
if(typePrediction == "Reachable") { 
  layersToCalc <- layersToCalc[grepl("Reachable",layersToCalc)] 
  layersToCalcNames <- layersToCalcNames[grepl("Reachable",layersToCalcNames)] 
}
if(typePrediction == "unConstrained") { 
  layersToCalc <- layersToCalc[!grepl("Reachable",layersToCalc)] 
  layersToCalcNames <- layersToCalcNames[!grepl("Reachable",layersToCalcNames)] 
}
layersToCalc <- layersToCalc[ sort(unique(unlist(sapply(scenariosToPredict, function(x) { which(grepl(x,layersToCalc)) } )))) ] 
layersToCalcNames <- layersToCalcNames[ sort(unique(unlist(sapply(scenariosToPredict, function(x) { which(grepl(x,layersToCalcNames)) } )))) ] 
layersToCalc <- layersToCalc[!grepl("Uncertainty",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("Uncertainty",layersToCalcNames)] 

presentDayLayer <- layersToCalc[grepl("speciesRichnessBaseline",layersToCalc)] 
presentDayLayer <- presentDayLayer[!grepl("Uncertainty",presentDayLayer)] 

layersToCalc <- layersToCalc[grepl("rangeRefugia",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[grepl("rangeRefugia",layersToCalcNames)] 

layersToCalc <- layersToCalc[!grepl("StandPerBaseline",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("StandPerBaseline",layersToCalcNames)] 

areaLayer <- raster::area(loadRData(presentDayLayer)) * bathymetryLayerFraction
presentDayLayer <- loadRData(presentDayLayer)
presentDayLayer[presentDayLayer == 0] <- NA
presentDayLayer[presentDayLayer > 1] <- 1
presentDayLayerArea <- presentDayLayer * areaLayer

if( resultsName == "Global" ) { 
  
  resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*2)+2,nrow=1))
  colnames(resultsDF) <- c("Region","Baseline",paste0(layersToCalcNames,"Area"),paste0(layersToCalcNames,"ProportionBaseline"))
  resultsDF[,1] <- "Global"
  polygonFeatureNames <- "Global"
  polygon <- extent(areaLayer)
  polygon <- as(polygon, 'SpatialPolygons')  
  
}

if( resultsName != "Global" ) {
  
  if(class(polygonPath) == "character") { polygon <- shapefile(polygonPath) }
  
  if( ! is.null(polygonFeature) ) { 
    
    polygonFeatureNames <- polygon@data[which(names(polygon@data) == polygonFeature)][,polygonFeature] 
    polygon <- gUnaryUnion(polygon, id = polygonFeatureNames, checkValidity=2L)
    polygonFeatureNames <- sapply(1:length(polygon), function(x) slot(slot(polygon, "polygons")[[x]],"ID") )
    polygon$names <- polygonFeatureNames
    
  }
  
  resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*2)+2,nrow=length(unique(polygonFeatureNames))))
  colnames(resultsDF) <- c("Region","Baseline",paste0(layersToCalcNames,"Area"),paste0(layersToCalcNames,"ProportionBaseline"))
  resultsDF[,1] <- unique(polygonFeatureNames)
  
}

## ------------

baselineAreas <- extract(presentDayLayerArea,polygon, small=TRUE)

for( i in 1:length(baselineAreas) ) { resultsDF[i,2] <- sum(baselineAreas[[i]],na.rm=T) }

## ------------

for( i in 1:length(layersToCalc) ) {
  
  rasterLayer <- loadRData(layersToCalc[i])
  rasterLayer[rasterLayer == 0] <- NA
  
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer
  
  for( p.i in 1:length(polygonFeatureNames) ){
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(i,"out of",length(layersToCalcNames),"\n")
    cat(p.i,"out of",length(polygonFeatureNames),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    polygon.pi <- polygon[p.i,]
    
    rasterLayer.pi <- crop(rasterLayer,polygon.pi)
    rasterLayer.pi <- raster::mask(rasterLayer.pi,polygon.pi)
    rasterLayerArea.pi <- crop(rasterLayerArea,rasterLayer.pi)
    rasterLayerArea.pi <- raster::mask(rasterLayerArea.pi,rasterLayer.pi)
    
    rasterLayerVals <- unlist(getValues(rasterLayer.pi))
    rasterLayerAreaVals <- unlist(getValues(rasterLayerArea.pi))
    
    if( ! is.null(rasterLayerAreaVals) ) { 
      
      resultsDF[p.i,i+2] <- sum( rasterLayerAreaVals ,na.rm=T )
      resultsDF[p.i,i+2+(length(layersToCalcNames))] <- ( sum( rasterLayerAreaVals , na.rm=T) / resultsDF[p.i,2] ) * 100
      
    }
  }
  
}

resultsDF[resultsDF == -Inf] <- NA
resultsDF[resultsDF == Inf] <- NA

if( resultsName == "MarineEcoRegion" ) {
  
  polygon <- shapefile(polygonPath)
  
  resultsDF <- data.frame(realmRegion=polygon$REALM[match(resultsDF$Region,polygon$ECOREGION)], 
                          provinceRegion=polygon$PROVINCE[match(resultsDF$Region,polygon$ECOREGION)],
                          resultsDF)
  
}

names(resultsDF) <- gsub("Reachable","",names(resultsDF))
names(resultsDF) <- gsub("angeR","",names(resultsDF))

write.csv(resultsDF, file=paste0(stackResultsFolder,"/Summary","/summaryExtentRefugia",resultsName,typePrediction,".csv"), row.names = FALSE)
rm(polygon)
