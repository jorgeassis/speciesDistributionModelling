## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

if( ! dir.exists(paste0(stackResultsFolder,"/Summary/")) ) { dir.create(paste0(stackResultsFolder,"/Summary/"), recursive = TRUE) }

polygon <- rgdal::readOGR(dsn = DescTools::SplitPath(polygonPath)$dirname, layer = DescTools::SplitPath(polygonPath)$filename)

if( is.null(polygonFeature) ) { stop("Missing polygonFeature") }
if( ! polygonFeature %in% names(polygon@data) ) { stop(paste0(polygonFeature," not in polygon")) }

names(polygon)[which(names(polygon) == polygonFeature )] <- "regionName"
polygon <- as_Spatial(st_make_valid(st_as_sf(polygon)))
polygonFeatureNames <- unique(as.data.frame(polygon[,"regionName"])[,1])

for(scenario in scenariosToPredict ) {

  sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)
  
  sFiles <- sFiles[grepl("Predictions",sFiles)]
  sFiles <- sFiles[grepl("ensemble",sFiles)]
  sFiles <- sFiles[grepl("Reclass",sFiles)]
  
  if(typePrediction == "Reachable") { sFiles <- sFiles[grepl("Reachable",sFiles)] }
  if(typePrediction == "unConstrained") { sFiles <- sFiles[!grepl("Reachable",sFiles)] }
  
  sFiles <- sFiles[grepl(scenario,sFiles)]
  sFiles <- sFiles[grepl(paste0("Global.RData"),sFiles)]
  
  sFiles <- sFiles[!grepl("Gain",sFiles)]
  sFiles <- sFiles[!grepl("Loss",sFiles)]
  sFiles <- sFiles[!grepl("Refugia",sFiles)]
  sFiles <- sFiles[!grepl("RangeShifts",sFiles)]
  
  sFiles <- sFiles[unlist(sapply(speciesPredicted,function(x) { agrep(x, sFiles) } ))]
  sFiles <- unique(sFiles)
  
  if(length(sFiles) > length(speciesPredicted)) { stop("Error :: 332")}
  
  resMatrix <- matrix(0,nrow=length(polygonFeatureNames),ncol=length(sFiles))
  rownames(resMatrix) <- polygonFeatureNames
  colnames(resMatrix) <- speciesPredicted

  for( i in 1:length(sFiles)) {
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(scenario," || Species",i,"out of",length(sFiles),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    spRaster <- loadRData(sFiles[i])

    for( j in 1:length(polygonFeatureNames) ) {
      
      polygon.j <- polygon[which(polygon$regionName == polygonFeatureNames[j]),]
      spRaster.j <- crop(spRaster,polygon.j)
      spRaster.j <- mask(spRaster.j,polygon.j)
      
      if( cellStats(spRaster.j,max,na.rm=T) != 0 ) { resMatrix[j,i] <- 1 }
    
    }
    
  }
  
  ## --------
  
  resMatrix <- data.frame(region=polygonFeatureNames,resMatrix)
  rownames(resMatrix) <- NULL
  
  if( resultsName == "MarineEcoRegion" ) {
    polygonNames <- rgdal::readOGR(dsn = DescTools::SplitPath(polygonPath)$dirname, layer = DescTools::SplitPath(polygonPath)$filename)
    resMatrix <- data.frame(Realm=polygonNames$REALM,  Province=polygonNames$PROVINCE, resMatrix)
  }
  
  if( resultsName == "EEZOceans" ) {
    polygonNames <- rgdal::readOGR(dsn = DescTools::SplitPath(polygonPath)$dirname, layer = DescTools::SplitPath(polygonPath)$filename)
    resMatrix <- data.frame(oceanBasin=polygonNames$name,  EEZ=polygonNames$EEZ, resMatrix)
  }
  
  write.csv(resMatrix, file=paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,".csv"), row.names = FALSE)
  
}

## --------

summaryPresenceInsideBaseline <- read.csv(paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,"Baseline",typePrediction,".csv"), header=TRUE)
summaryPresenceInsideBaseline <- summaryPresenceInsideBaseline[,sapply(1:ncol(summaryPresenceInsideBaseline),function(x) is.numeric(summaryPresenceInsideBaseline[1,x]))]

for( scenario in scenariosToPredict[ scenariosToPredict != "Baseline"] ) {
  
  summaryPresenceInsideScenario <- read.csv(paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,".csv"), header=TRUE)
  summaryPresenceInsideScenario <- summaryPresenceInsideScenario[,sapply(1:ncol(summaryPresenceInsideScenario),function(x) is.numeric(summaryPresenceInsideScenario[1,x]))]
  
  resMatrix <- data.frame(baseline=apply(summaryPresenceInsideBaseline,1,sum),loss=NA,gain=NA, refugia=NA)
  
  for( region in 1:nrow(summaryPresenceInsideScenario))   {
    
    resultRegionDiff <- apply(data.frame(Baseline=unlist(summaryPresenceInsideBaseline[region,]),Scenario=unlist(summaryPresenceInsideScenario[region,])),1,diff)
    resultRegionSum <- apply(data.frame(Baseline=unlist(summaryPresenceInsideBaseline[region,]),Scenario=unlist(summaryPresenceInsideScenario[region,])),1,sum)
    
    resMatrix[region,"loss"] <- sum(resultRegionDiff == -1)
    resMatrix[region,"gain"] <- sum(resultRegionDiff == 1)
    resMatrix[region,"refugia"] <- sum(resultRegionSum == 2)
    
  }
  
  resMatrix <- data.frame(region=polygonFeatureNames,resMatrix)
  
  if( resultsName == "MarineEcoRegion" ) {
    polygonNames <- rgdal::readOGR(dsn = DescTools::SplitPath(polygonPath)$dirname, layer = DescTools::SplitPath(polygonPath)$filename)
    resMatrix <- data.frame(Realm=polygonNames$REALM,  Province=polygonNames$PROVINCE, resMatrix)
  }
  
  if( resultsName == "EEZOceans" ) {
    polygonNames <- rgdal::readOGR(dsn = DescTools::SplitPath(polygonPath)$dirname, layer = DescTools::SplitPath(polygonPath)$filename)
    resMatrix <- data.frame(oceanBasin=polygonNames$name,  EEZ=polygonNames$EEZ, resMatrix)
  }
  
  write.csv(resMatrix, file=paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,"Changes.csv"), row.names = FALSE)
  
}

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------