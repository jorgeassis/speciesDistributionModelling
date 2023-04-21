## -------------------------------------------------------------------------------
## --------------------------------------------------
## --------------------------------------------------

speciesListWorms.i <- speciesListWorms[speciesListWorms$scientificname %in% speciesPredicted,]
speciesListWorms.i.agg <- unique(speciesListWorms[speciesListWorms$scientificname %in% speciesPredicted,tolower(taxaLevel)])

for( scenario in scenariosToPredict ) {
  
    results.i <- as.data.frame(matrix(NA,ncol=length(speciesListWorms.i.agg),nrow=length(speciesListWorms.i.agg)))
    colnames(results.i) <- speciesListWorms.i.agg
    rownames(results.i) <- speciesListWorms.i.agg
      
    for( taxa.i in 1:length(speciesListWorms.i.agg) ) {
      
      for( taxa.j in 1:length(speciesListWorms.i.agg) ) {
        
        if(taxa.i == taxa.j) { next }
        
        if( ! exists("globalArea")) { globalArea <- raster::area(bathymetryLayerFraction) }
        
        sp.i <- speciesListWorms.i[which(speciesListWorms.i[,tolower(taxaLevel)] == speciesListWorms.i.agg[taxa.i]),"scientificname"]
    
        sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)
        sFiles <- sFiles[unlist(sapply(sp.i,function(x) { which(unlist(grepl(x,sFiles))) } ))]
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
        
        if( length(sFiles) > 1) {
          sp.i.range <- calc(stack(sapply(sFiles,function(x) { loadRData(x)})),sum,na.rm=T)
        }
        if( length(sFiles) == 1) {
          sp.i.range <- loadRData(sFiles)
        }
        
        sp.i.range[sp.i.range >= 1] <- 1
        sp.i.range.area <- cellStats( sp.i.range * globalArea * bathymetryLayerFraction, stat='sum', na.rm=T)
        
        sp.j <- speciesListWorms.i[which(speciesListWorms.i[,tolower(taxaLevel)] == speciesListWorms.i.agg[taxa.j]),"scientificname"]
        
        sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)
        sFiles <- sFiles[unlist(sapply(sp.j,function(x) { which(unlist(grepl(x,sFiles))) } ))]
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
        
        if( length(sFiles) > 1) {
          sp.j.range <- calc(stack(sapply(sFiles,function(x) { loadRData(x)})),sum,na.rm=T)
        }
        if( length(sFiles) == 1) {
          sp.j.range <- loadRData(sFiles)
        }
        
        sp.j.range[sp.j.range >= 1] <- 1
        sp.j.range.area <- cellStats( sp.j.range * globalArea * bathymetryLayerFraction, stat='sum', na.rm=T)
        
        overlapArea <- cellStats( (sp.i.range == 1 & sp.j.range == 1) * globalArea * bathymetryLayerFraction, stat='sum', na.rm=T)
    
        results.i[taxa.i,taxa.i] <- sp.i.range.area 
        results.i[taxa.j,taxa.j] <- sp.j.range.area 
        results.i[taxa.i,taxa.j] <- (1 - ((sp.i.range.area - overlapArea) / sp.i.range.area )) * 100 
        results.i[taxa.j,taxa.i] <- (1 - ((sp.j.range.area - overlapArea) / sp.j.range.area )) * 100 
        
      }
    }
    
    ## --------

    write.csv(results.i, file=paste0(stackResultsFolder,"/Summary","/summaryRangeExtentOverlap.",scenario,typePrediction,".csv"), row.names = FALSE)
    
}

## --------

results.baseline <- read.csv(paste0(stackResultsFolder,"/Summary","/summaryRangeExtentOverlap.","Baseline",typePrediction,".csv"))

for( scenario in scenariosToPredict[scenariosToPredict != "Baseline"] ) {
 
  results.i <- read.csv(paste0(stackResultsFolder,"/Summary","/summaryRangeExtentOverlap.",scenario,typePrediction,".csv"))
  results.i <- ((results.i - results.baseline) / results.baseline ) * 100
  
  write.csv(results.i, file=paste0(stackResultsFolder,"/Summary","/summaryRangeExtentOverlap.",scenario,typePrediction,"Diff.csv"), row.names = FALSE)
  
}

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------