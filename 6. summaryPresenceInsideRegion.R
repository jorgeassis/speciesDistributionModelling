## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

if( ! dir.exists(paste0(stackResultsFolder,"/Summary/")) ) { dir.create(paste0(stackResultsFolder,"/Summary/"), recursive = TRUE) }

if(class(polygon) == "character") { polygon <- shapefile(polygon) }

if( ! is.null(polygonFeature) ) { 
  
  polygonFeatureNames <- polygon@data[which(names(polygon@data) == polygonFeature)][,polygonFeature] 
  polygon <- gUnaryUnion(polygon, id = polygonFeatureNames, checkValidity=2L)
  polygonFeatureNames <- sapply(1:length(polygon), function(x) slot(slot(polygon, "polygons")[[x]],"ID") )
  polygon$names <- polygonFeatureNames
  
}

for(scenario in scenariosToPredict ) {
  
  # rasterOptions(todisk = TRUE)
  # rasterOptions(tmpdir=tempFolder)
  
  sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)
  sFiles <- sFiles[unlist(sapply(speciesPredicted,function(x) { which(unlist(grepl(x,sFiles))) } ))]
  
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
  
  sFiles <- unique(sFiles)
  
  if(length(sFiles) > length(speciesPredicted)) { stop("Error :: 332")}
  
  sFilesNames <- sapply(sFiles, function(x) {  substr(x,nchar(mainResultsDirectory) + 1 ,  unlist(gregexpr("/Predictions",x)) - 1)   } )
  names(sFilesNames) <- NULL
  sFilesNames <- gsub("/","",sFilesNames)
  
  resMatrix <- matrix(0,nrow=nrow(polygon),ncol=length(sFiles))
  rownames(resMatrix) <- polygon$names
  colnames(resMatrix) <- sFilesNames
  
  for( i in 1:length(sFiles)) {
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(scenario," || Species",i,"out of",length(sFiles),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    spRaster <- loadRData(sFiles[i])
    
    for( j in 1:nrow(polygon)) {
      
      polygon.j <- polygon[j,]
      spRaster.j <- crop(spRaster,polygon.j)
      
      if( cellStats(spRaster.j,max,na.rm=T) == 0 ) { next }
      
      cells <- sfraster::cellFromPolygon(spRaster.j,polygon.j)
      
      if( 1 %in% spRaster.j[unlist(cells)] ) { resMatrix[j,i] <- 1 }
      
    }
    
  }
  
  ## --------
  
  resMatrix <- data.frame(region=polygon$names,resMatrix)
  write.csv(resMatrix, file=paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,".csv"), row.names = FALSE)
  
}

## --------

summaryPresenceInsideBaseline <- read.csv(paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,"Baseline",typePrediction,".csv"), header=TRUE)[,-1]

for( scenario in scenariosToPredict[ scenariosToPredict != "Baseline"] ) {
  
  summaryPresenceInsideScenario <- read.csv(paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,".csv"), header=TRUE)[,-1]
  resMatrix <- data.frame(region=polygon$names,baseline=apply(summaryPresenceInsideBaseline,1,sum),loss=NA,gain=NA, refugia=NA)
  
  for( region in 1:nrow(summaryPresenceInsideScenario))   {
    
    resultRegionDiff <- apply(data.frame(Baseline=unlist(summaryPresenceInsideBaseline[region,]),Scenario=unlist(summaryPresenceInsideScenario[region,])),1,diff)
    resultRegionSum <- apply(data.frame(Baseline=unlist(summaryPresenceInsideBaseline[region,]),Scenario=unlist(summaryPresenceInsideScenario[region,])),1,sum)
    
    resMatrix[region,"loss"] <- sum(resultRegionDiff == -1)
    resMatrix[region,"gain"] <- sum(resultRegionDiff == 1)
    resMatrix[region,"refugia"] <- sum(resultRegionSum == 2)
    
  }
  
  write.csv(resMatrix, file=paste0(stackResultsFolder,"/Summary","/summarySpeciesIn",resultsName,scenario,typePrediction,"Changes.csv"), row.names = FALSE)
  
}

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------