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
layersToCalc <- layersToCalc[grepl("speciesRichness",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[grepl("speciesRichness",layersToCalcNames)] 
layersToCalc <- layersToCalc[!grepl("speciesRichnessChanges",layersToCalc)] 
layersToCalcNames <- layersToCalcNames[!grepl("speciesRichnessChanges",layersToCalcNames)] 

## Per latitude

regionOfInterest <- extent(loadRData(presentDayLayer))
bins <- seq(regionOfInterest[3],regionOfInterest[4], by=1)
bins <- data.frame(binFrom=bins[-length(bins)],binTo=bins[-1])
bins <- bins[ sort(bins$binFrom, decreasing = TRUE, index.return=T)$ix , ]

resultsDF <- as.data.frame(matrix(NA, ncol=(length(layersToCalcNames)*2),nrow=nrow(bins)))
colnames(resultsDF) <- c(paste0(layersToCalcNames,"Mean"),paste0(layersToCalcNames,"Area"))
resultsDF <- data.frame(bins,resultsDF)

## ----------------

for(i in 1:length(layersToCalc)) {
  
  cat("\014")
  cat("## --------------------------------- \n")
  cat(i,"out of",length(layersToCalcNames),"\n")
  cat("## --------------------------------- \n")
  cat("\n")
  
  rasterLayer <- loadRData(layersToCalc[i])
  rasterLayer[rasterLayer == 0] <- NA
  
  rasterLayerBinomial <- rasterLayer
  rasterLayerBinomial[rasterLayerBinomial > 1] <- 1
  rasterLayerArea <- rasterLayerBinomial * areaLayer
  
  for(b in 1:nrow(resultsDF)) {
    
    rasterLayer.bin <- crop(rasterLayer,extent(extent(regionOfInterest)[1],extent(regionOfInterest)[2],bins[b,1],bins[b,2]))
    rasterLayerArea.bin <- crop(rasterLayerArea,extent(extent(regionOfInterest)[1],extent(regionOfInterest)[2],bins[b,1],bins[b,2]))
    presentDayLayerArea.bin <- crop(presentDayLayerArea,extent(extent(regionOfInterest)[1],extent(regionOfInterest)[2],bins[b,1],bins[b,2]))
    
    resultsDF[b,i+2] <- cellStats(rasterLayer.bin,mean,na.rm=T)
    resultsDF[b,i+2+(length(layersToCalcNames))] <- sum( unlist(getValues(rasterLayerArea.bin)) ,na.rm=T )
    
  }
  
}

resultsDF[resultsDF == -Inf] <- NA
resultsDF[resultsDF == Inf] <- NA
resultsDF[is.na(resultsDF)] <- 0

names(resultsDF) <- gsub("Reachable","",names(resultsDF))
write.csv(resultsDF, file=paste0(stackResultsFolder,"/Summary","/summaryPerLatitudinalBin",resultsName,typePrediction,".csv"), row.names = FALSE)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------