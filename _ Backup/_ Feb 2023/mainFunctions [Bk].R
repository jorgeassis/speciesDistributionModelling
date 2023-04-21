# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

packages.to.use <- c(  
   "ecospat" 
  ,"plyr"
  ,"raster"
  ,"h3js"
  ,"sf"
  ,"maptools"
  ,"gridExtra"
  ,"xgboost"
  , "BiodiversityR"
  ,"devtools"
  ,"patchwork"
  ,"usdm"
  ,"fasterize"
  ,"sm"
  ,"sf"
  ,"SDMTools"
  , "blockCV"
  , "bcp"
  , "spThin"
  , "ecodist"
  , "credentials"
  , "rnaturalearth"
  , "ggplot2"
  , "mboost"
  , "blockCV"
  , "raster"
  , "gdata"
  , "dismo"
  , "gbm"
  , "sp"
  , "parallel"
  , "doParallel"
  , "biganalytics"
  , "nicheROVER"
  , "rgeos"
  , "rgdal"
  , "devtools"
  , "ENMeval"
  , "remotes"
  ,"leaflet"
  , "monmlp")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  
  if( ! package %in% rownames(installed.packages()) & package == "h3js" ) { devtools::install_github("saurfang/h3js") }
  if( ! package %in% rownames(installed.packages()) & package == "h3" ) { devtools::install_github("crazycapivara/h3-r") }
  if( ! package %in% rownames(installed.packages()) & package == "h3r" ) { devtools::install_github("scottmmjackson/h3r") }
  if( ! package %in% rownames(installed.packages()) & package == "h3r" ) { devtools::install_github("harryprince/h3r", ref="bug-fix/Makefile") }
  if( ! package %in% rownames(installed.packages()) & package == "h3jsr" ) { remotes::install_github("obrl-soil/h3jsr") }
  if( ! package %in% rownames(installed.packages()) & package == "rnaturalearthhires" ) { devtools::install_github("ropensci/rnaturalearthhires")  }
  if( ! package %in% rownames(installed.packages()) & package == "spatialEco" ) { remotes::install_github("jeffreyevans/spatialEco")  }
  
  if( ! package %in% rownames(installed.packages()) & package == "patchwork" ) { devtools::install_github("thomasp85/patchwork") }
  if( ! package %in% rownames(installed.packages()) & package == "SDMTools" ) { install_version("SDMTools", "1.1-221") }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source") }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

modelPerformance <- function(observed,predicted,index) {
  
    predicted.accuracy <- SDMTools::accuracy( observed , predicted , threshold = 100 )
    predicted.accuracy$deviance <- modEvA::Dsquared(obs = observed, pred=predicted, family = "binomial")
    predicted.accuracy$aicc <- AIC(glm(observed~predicted))
    
    # boyce
    boyce.p <- predicted[observed==1]
    boyce.p <- boyce.p[!is.na(boyce.p)]
    boyce.a <- predicted
    boyce.a <- boyce.a[!is.na(boyce.a)]
    
    if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01) }
    if( length(unique(boyce.p)) == 1) { boyce.p <- sapply(boyce.p, function(x) {   x + sample(c(0.001,-0.001), 1)}) }
    
    boyce.value <- ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE,rm.duplicate = TRUE)
    predicted.accuracy$boyce <- boyce.value$cor
    
    if(index == "boyce" ) { predicted.accuracy <- predicted.accuracy[which.min(abs(predicted.accuracy$threshold - min(  boyce.value$HS[ max( ifelse( length(which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))]) > 0 , which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))] , 1)) ]  ))),] }
    if(index == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
    if(index == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    if(index == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    
    predicted.accuracy$tss <- predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1
    names(predicted.accuracy)[2] <- "auc"
    
    ## -------

    return(predicted.accuracy)
  
}

## -----------------------

modelGridSearch <- function(modelAlgorithm,hyperParam,modelData,cvFolds,monotonicity,cvIndex){
  
  monotonicity.i <- monotonicity[match( names(modelData@data) , names(monotonicity))]
  
  combP <- expand.grid(hyperParam)
  combPResults <- data.frame()
  
  for( c in 1:nrow(combP) ) {
    
    performance <- data.frame()
    
    for( k in 1:ncol(cvFolds$train) ) {
      
      trainData <- data.frame(pa=modelData@pa[cvFolds$train[,k]],modelData@data[cvFolds$train[,k],])
      testData <- data.frame(pa=modelData@pa[cvFolds$test[,k]],modelData@data[cvFolds$test[,k],])
      
      if(modelAlgorithm == "BRT"){
        
        if( "n.trees" %in% names(combP)) { n.trees = combP[c,"n.trees"] } else { n.trees = 100 }
        if( "interaction.depth" %in% names(combP)) { interaction.depth = combP[c,"interaction.depth"] } else { interaction.depth = 2 }
        if( "shrinkage" %in% names(combP)) { shrinkage = combP[c,"shrinkage"] } else { shrinkage = 0.1 }
        
        model <- gbm(pa ~ ., distribution = "bernoulli", data=trainData, var.monotone=monotonicity.i, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, cv.folds = 0, verbose = FALSE)
        predTest <- predict(model,testData, type="response")
        obsTest <- testData$pa
        # predTrain <- predict(model,trainData, type="response")
        # obsTrain <- trainData$pa
        
      }

      performance <- rbind(performance,modelPerformance(obsTest,predTest,cvIndex))
      
    }
    
    combPResults <- rbind(combPResults,data.frame(combP[c,],t(data.frame(apply(performance,2,mean)))))
    
  }
  
  return(combPResults)
  
}

## -----------------------

hyperParam <- list(n.trees=100,interaction.depth = 2, shrinkage = 0.01)

modelTrainCV <- function(modelAlgorithm,hyperParam,modelData,cvFolds,monotonicity){
  
  monotonicity.i <- monotonicity[match( names(modelData@data) , names(monotonicity))]
  combP <- expand.grid(hyperParam)
  combPResults <- data.frame()

  performance <- data.frame()
  models <- list()
  
  for( k in 1:ncol(cvFolds$train) ) {
    
    trainData <- data.frame(pa=modelData@pa[cvFolds$train[,k]],modelData@data[cvFolds$train[,k],])
    testData <- data.frame(pa=modelData@pa[cvFolds$test[,k]],modelData@data[cvFolds$test[,k],])
    
    if(modelAlgorithm == "BRT"){
      
      if( "n.trees" %in% names(hyperParam)) { n.trees = hyperParam$n.trees } else { n.trees = 100 }
      if( "interaction.depth" %in% names(hyperParam)) { interaction.depth = hyperParam$interaction.depth } else { interaction.depth = 2 }
      if( "shrinkage" %in% names(hyperParam)) { shrinkage = hyperParam$shrinkage } else { shrinkage = 0.1 }
      
      model <- gbm(pa ~ ., distribution = "bernoulli", data=trainData, var.monotone=monotonicity.i, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, cv.folds = 0, verbose = FALSE)
      models <- c(models,list(model))
      predTest <- predict(model,testData, type="response")
      obsTest <- testData$pa
      # predTrain <- predict(model,trainData, type="response")
      # obsTrain <- trainData$pa
      
    }

    performance <- rbind(performance,modelPerformance(obsTest,predTest,cvIndex))
    
  }
  
  combPResults <- list(models=models,performance)
  
  return(combPResults)
  
}

## -----------------------








correctDateLineMap <- function(spatialObjectSF) {
  
  if(! TRUE %in% (class(spatialObjectSF) == "sf")) { stop("Object must be of class sf")}
  
  dggsOcean <- spatialObjectSF
  dggsOcean <- as(dggsOcean, "Spatial")
  listToRemove <- numeric()
  listToAdd <- list()
  
  for(i in 1:length(dggsOcean)) {
    
    if( extent(dggsOcean[i,])[1] < 0 & extent(dggsOcean[i,])[2] > 0 & max(extent(dggsOcean[i,])[1:2])-min(extent(dggsOcean[i,])[1:2]) > 180 ) { 
      
      listToRemove <- c(listToRemove,i)
      
      coords1 <- fortify(dggsOcean[i,])[,1:2]
      coords2 <- fortify(dggsOcean[i,])[,1:2]
      coords1[coords1[,1] < 0 ,1] <- 180
      coords2[coords2[,1] > 0 ,1] <- -180
      coords2 <- rbind(coords2,data.frame(long=-180,lat=max(coords2[,2])))
      coords2 <- rbind(coords2,data.frame(long=-180,lat=min(coords2[,2])))
      coords2 <- coords2[chull(coords2),] 
      
      coords1 <- spPolygons(as.matrix(coords1))
      coords2 <- spPolygons(as.matrix(coords2))
      
      listToAdd.i <- unionSpatialPolygons(union(coords1, coords2), c(1,1))
      listToAdd.i <- SpatialPolygonsDataFrame(listToAdd.i,as.data.frame(dggsOcean[i,]), match.ID = FALSE)
      
      if( length(listToAdd) >  0 ) { listToAdd <- bind(listToAdd, listToAdd.i) }
      if( length(listToAdd) == 0 ) { listToAdd <- listToAdd.i }
      
    }
    
    if( (extent(dggsOcean[i,])[1] > 0 & extent(dggsOcean[i,])[2] < 0) ) { stop(i) }
    
  }
  
  if( length(listToRemove) > 0) {
    
    dggsOcean <- dggsOcean[-listToRemove,]
    dggsOcean <- bind(dggsOcean,listToAdd)
    
  }
  
  dggsOcean <- st_as_sf(dggsOcean)
  
  return(dggsOcean)      
}

# ---------------------

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------------

plotMap <- function(coordinates,radius,color) {
  
  set.seed(42)
  
  m <- leaflet()
  m <- addTiles(m)
  # m <- addMarkers(m, lng=coordinates[,1], lat=coordinates[,2], popup=paste0( "Species record ") , icon = greenLeafIcon)
  
  m <- addCircleMarkers(m, lng=coordinates[,1], lat=coordinates[,2], 
                        popup=paste0( "Species record ") , 
                        radius = radius, color = color , 
                        stroke = FALSE, fillOpacity = 0.5 )
  
  return(m)
  
}

## -----------------------

processLayers <- function(rasterLayers,occurrenceRecords,regionBuffer,minDepth,maxDepth,intertidal) {
  
  if( ! is.null(minDepth) ) { if( minDepth == "NULL" ) { minDepth <- NULL } }
  if( ! is.null(maxDepth) ) { if( maxDepth == "NULL" ) { maxDepth <- NULL } }

  rastersBathymetryID <- NULL
  rastersIntertidalID <- NULL
  
  rasters <- rasterLayers
  
  if( ! is.null(regionBuffer) ) {
    
    if (length(regionBuffer) == 1) { regionBuffer <- rep(regionBuffer,4)}
    
    final.extent.rasters <- c( min(occurrenceRecords[,1] , na.rm=T) - regionBuffer[1],
                               max(occurrenceRecords[,1] , na.rm=T) + regionBuffer[2],
                               min(occurrenceRecords[,2] , na.rm=T) - regionBuffer[3],
                               max(occurrenceRecords[,2] , na.rm=T) + regionBuffer[4] )
    
    final.extent.rasters <- c(  ifelse(final.extent.rasters[1] < -180 , -180 , final.extent.rasters[1]),
                                ifelse(final.extent.rasters[2] > 180 , 180 , final.extent.rasters[2]),
                                ifelse(final.extent.rasters[3] < -90 , -90 , final.extent.rasters[3]),
                                ifelse(final.extent.rasters[4] > 90 , 90 , final.extent.rasters[4]) )
    
    rasters <- raster::crop(rasters, extent(final.extent.rasters))
    
  } else { final.extent.rasters <- c(-180,180,-90,90) }
  
  ## -----------------------
  
  regions <- subset(rasters,1)
  regions[!is.na(regions) ] <- 1
  regions <- aggregate(regions,5, fun=max, na.rm=TRUE)
  regions <- clump(regions)
  regionsUsed <- unique(raster::extract(regions,occurrenceRecords))
  regionsUsed <- regionsUsed[ !is.na(regionsUsed) ]
  regions[ ! regions %in% regionsUsed ] <- NA
  rastersLoc <- xyFromCell( subset(rasters,1) , Which( !is.na(subset(rasters,1)), cell=T))
  rastersLocNA <- which(is.na(raster::extract(regions,rastersLoc)))
  rasters[ Which( !is.na(subset(rasters,1)), cell=T)[rastersLocNA] ] <- NA
  
  ## -----------------------
  
  if( ! is.null(intertidal) ) {
    
    intertidalLayer <- raster(intertidal)
    intertidalLayer <- raster::crop(intertidalLayer, rasters)
    rastersIntertidalID <- Which( is.na(intertidalLayer), cells=TRUE)
    
  }
  
  ## -----------------------
  
  if( ! is.null(minDepth) | ! is.null(maxDepth) ) {
    
    if( is.null(minDepth) ) { minDepth <- 0 }
    if( is.null(maxDepth) ) { maxDepth <- 99999 }
    
    minDepth <- minDepth * (-1)
    maxDepth <- maxDepth * (-1)
    
    rclmat <- data.frame(from=c(-99999,maxDepth,minDepth),to=c(maxDepth,minDepth,0),value=c(NA,1,NA))
    rclmat <- rclmat[rclmat$from != rclmat$to,]
    
    bathymetry <- raster(bathymetryDataLayer)
    bathymetry <- crop(bathymetry, rasters)
    bathymetry <- reclassify(bathymetry, as.matrix(rclmat))
    rastersBathymetryID <- Which( is.na(bathymetry), cells=TRUE)
    
  }
  
  ## -----------------------
  
  rastersID <- unique(c(rastersIntertidalID,rastersBathymetryID))
  rasters[rastersID] <- NA
  
  if(class(rasters) == "RasterBrick") { rasters <- stack(rasters)}
  
  return(rasters)
  
}

## -----------------------

dropNoVariationLayers <- function(rasterLayers,records) {
  
  cave <- function(x) { length(unique(x)) / length(x) }
  varRasterLayers <- which( apply( raster::extract(rasterLayers,records) ,2,cave) > 0.025)
  rasterLayers <- subset(rasterLayers,varRasterLayers)
  return(rasterLayers)
  
}

## -----------------------

relocateNACoords <- function(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth) {
  
  set.seed(42)
  
  if(relocateType == "nonDistance") {
    
    vals <- raster::extract(subset(rasterLayers,1),occurrenceRecords)
    occurrenceRecords <- occurrenceRecords[which(!is.na(vals)),]
    
  }
  
  if(relocateType == "distance") {
    
    if( relocateSpeciesDepth ) { 
      
      if( is.null(minDepth) ) { minDepth <- 0 }
      if( is.null(maxDepth) ) { maxDepth <- 99998 }
      
      minDepth.i <- minDepth * (-1)
      maxDepth.i <- maxDepth * (-1)
      
      rclmat <- data.frame(from=c(-99999,maxDepth.i,minDepth.i),to=c(maxDepth.i,minDepth.i,0),value=c(NA,1,NA))
      rclmat <- rclmat[rclmat$from != rclmat$to,]
      
      shape <- subset(rasterLayers,1)
      
      bathymetry <- raster(bathymetryDataLayer)
      bathymetry <- crop(bathymetry, shape )
      bathymetry <- mask(bathymetry,shape )
      
      shape <- reclassify(bathymetry, as.matrix(rclmat))
      
    }
    
    if( ! relocateSpeciesDepth ) { shape <- subset(rasterLayers,1) } 
    
    to.relocate <- unique(which(is.na(raster::extract(shape,occurrenceRecords))))
    coordinates.to.relocate <- as.matrix(occurrenceRecords[to.relocate,],ncol=2)
    
    if( nrow(coordinates.to.relocate) > 0 ) { 
      
      old.presences <- occurrenceRecords[ (1:nrow(occurrenceRecords))[! 1:nrow(occurrenceRecords) %in% to.relocate] ,]
      
      correct.points <- xyFromCell(shape,  Which(!is.na(shape), cells=TRUE) ) # result
      
      cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
      cat( paste0("\n"))
      
      near.cells <- numeric(nrow(coordinates.to.relocate))
      
      for(p in 1:nrow(coordinates.to.relocate)) {
        
        near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,]),longlat=TRUE)
        
        if( near.cell.p[which.min(near.cell.p)] <= sqrt(sum(relocateSpeciesDistance^2,relocateSpeciesDistance^2)) ) {
          
          near.cell.p <- which.min(near.cell.p)
          
        } else {   near.cell.p <- NA }
        
        near.cells[p] <- near.cell.p
        
      }
      
      relocated <- which(!is.na(near.cells))
      
      if( length(relocated) > 0) {
        
        near.cells <- data.frame(Lon=correct.points[near.cells[relocated],1],Lat=correct.points[near.cells[relocated],2])
        colnames(near.cells) <- colnames(old.presences)
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
      cat( paste0("Done relocating"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
    if( nrow(coordinates.to.relocate) == 0) { 
      
      cat( paste0("None to Relocate"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
  }
  
  toRemove <- raster::extract(shape,occurrenceRecords)
  toRemove <- which(is.na(toRemove))
  
  if( length(toRemove) > 0 ) { occurrenceRecords <- occurrenceRecords[-toRemove,] }
  
  return( occurrenceRecords )
  
}

## -----------------------

spatialAutocorrelation <- function(occurrenceRecords,rasterLayers,autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif) {
  
  set.seed(42)
  
  if(nrow(occurrenceRecords) > 1000) { 
    
    cat( paste0("\n"))
    cat( paste0("\n"))
    cat("More than 1000 occurrence records.","\n")  
    cat("Using a maximum of 1000 random records.","\n")
    cat( paste0("\n"))
    
    set.seed(1)
    occurrenceRecords <- occurrenceRecords[sample(1:nrow(occurrenceRecords),1000,replace = FALSE),]
    
  }
  
  # --------
  
  presences.environment <- data.frame(raster::extract(rasterLayers,occurrenceRecords,stringsAsFactors=FALSE))
  
  occurrenceRecords <- occurrenceRecords[which(complete.cases(presences.environment)),]
  presences.environment <- presences.environment[which(complete.cases(presences.environment)),]
  
  # --------
  
  space <- spDists(as.matrix(occurrenceRecords),as.matrix(occurrenceRecords),longlat=TRUE)
  data <- ecodist::distance( presences.environment , method = "euclidean")
  data <- as.matrix(data)
  
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)
  
  resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,autocorrelationClassDistance),
                              classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationClassDistance),
                              R=NA,
                              pVal=NA)
  
  for( i in 1:nrow(resultsMatrix)) {
    
    d1 = resultsMatrix[i,1]
    d2 = resultsMatrix[i,2]
    
    data.d <- as.vector(data)
    space.d <- as.vector(space)
    
    remove <- which(space.d < d1 | space.d > d2)
    data.d <- data.d[-remove]
    space.d <- space.d[-remove]
    
    if( length(data.d) == 0 | length(space.d) == 0) { next}
    
    modelobject <- lm(space.d~data.d)
    
    f <- summary(modelobject)$fstatistic
    p <-0
    
    tryCatch( p <- pf(f[1],f[2],f[3],lower.tail=FALSE) , error=function(e) { Error <<- TRUE })
    
    #p <- summary.lm(modelobject)$coefficients[2,"Pr(>|t|)"] 
    resultsMatrix[i,3] <- summary(modelobject)$adj.r.squared
    resultsMatrix[i,4] <- p
    
  }
  
  vect.signif <- as.numeric(resultsMatrix[,4] > autocorrelationSignif)
  if(vect.signif[1] == 1) { vect.signif[1] <- 0}
  vect.signif.pch <- c(19,1)
  
  minDistance = min(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  meanDistance = mean(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  sdDistance = sd(resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ])
  
  if( is.na(minDistance)) { minDistance <- autocorrelationMaxDistance }
  
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Non-correlated distances | Min: ",minDistance,"; Mean: ",meanDistance," km"))
  cat( paste0("\n"))
  
  sac <- resultsMatrix[,c(1,3,4)]
  colnames(sac) <- c("Distance","R","signif")
  
  figure <- ggplot() +
    geom_line(data=sac, aes(x=Distance, y=R), linetype = "dashed", size=0.25)+
    geom_point(data=sac[sac$signif >= autocorrelationSignif,], aes(x=Distance, y=R), size=3) +
    geom_point(data=sac[sac$signif < autocorrelationSignif,], aes(x=Distance, y=R), size=3, shape=1) + 
    labs(x = "Distance (km)") + 
    labs(y = "Pearson' correlation")
  
  return( list(figure=figure, minDistance = minDistance , meanDistance = meanDistance , sdDistance = sdDistance  ))
  
}

## -----------------------		

correlatedPairs <- function(rasterLayers,speciesData,threhold,dataLayersMonotonocity.i) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasterLayers))-1)) {
    
    for( j in (i+1):length(names(rasterLayers))) {
      
      corr.val <- abs(cor( raster::extract(subset(rasterLayers,i),speciesData[,dataRecordsNames]) , raster::extract(subset(rasterLayers,j),speciesData[,dataRecordsNames]) ,use = "pairwise.complete.obs"))
      
      if(is.na(corr.val)) { corr.val <- 0 }
      
      if( !is.null(dataLayersMonotonocity.i)) {
        
        if( corr.val >= threhold & (sum(c(1,-1) %in% dataLayersMonotonocity.i[c(i,j)]) != 2)  ) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasterLayers,i)),
                                                                                                                                               Var.2=names(subset(rasterLayers,j)),
                                                                                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      }
      
      if( is.null(dataLayersMonotonocity.i)) {
        
        if( corr.val >= threhold  ) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasterLayers,i)),
                                                                                   Var.2=names(subset(rasterLayers,j)),
                                                                                   Cor=corr.val,stringsAsFactors = FALSE)) }
      }
      
    }
    
  }
  
  return(list.of.cor)
  
}

## -----------------------

spatialThinning <- function(occurrenceRecords,minDistance,verbose) {
  
  occurrenceRecords <- occurrenceRecords[which(!duplicated(occurrenceRecords)),]
  dataThinning <- data.frame(Name="Sp",occurrenceRecords)
  
  thinError <- FALSE
  
  tryCatch( 
    coordinates.t <- thin(dataThinning,
                          lat.col = dataRecordsNames[2],
                          long.col = dataRecordsNames[1],
                          spec.col = "Name",
                          thin.par = minDistance,
                          reps = 1,
                          write.files = FALSE,
                          locs.thinned.list.return = TRUE,
                          verbose = FALSE)[[1]]
    , error=function(e) { thinError <<- TRUE })
  
  if( thinError ) { coordinates.t <- occurrenceRecords }
  
  if(verbose) {
    
    cat( paste0("\n"))
    cat( paste0("\n"))    
    cat( paste0("Input Records: ",nrow(occurrenceRecords)))
    cat( paste0("\n"))
    cat( paste0("Final Records: ",nrow(coordinates.t)))
    
  }
  
  # Remove from main dataset of occurrences
  colnames( coordinates.t ) <- dataRecordsNames
  
  # Remove log file
  if(length(list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt"))){
    file.remove( list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt") ) 
  }
  
  return(coordinates.t)
  
}

## -----------------------

getBackgroundExtent <- function(occurrenceRecords,rasterLayers,maxBufferSize, bufferStep , type ) {
  
  extentOccurrenceRecords <- c(min(occurrenceRecords[,1]),max(occurrenceRecords[,1]),min(occurrenceRecords[,2]),max(occurrenceRecords[,2]))
  
  backgroundInformation <- subset(rasterLayers,1)
  backgroundInformation[!is.na(backgroundInformation)] <- 1
  
  cl <- parallel::makeCluster( nCores )
  registerDoParallel(cl)
  
  results <- foreach(buff = seq(3,maxBufferSize+2,by=bufferStep), .export=c("paMinimum","paRatio","cvIndex","accuracyEstimate","dataRecordsNames","accuracyEstimateCalc","generatePseudoAbsences"), .combine=rbind, .packages = c("rgeos","ecospat","raster","blockCV","SDMTools")) %dopar% {
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    sink.points.poly.list <- list()
    for( l in 1:length(sink.points.poly) ) {
      sink.points.poly.list <- c(sink.points.poly.list, list(gBuffer( sink.points.poly[l,], width=buff )))
    }
    
    sink.points.poly <- do.call("rbind",c(args = sink.points.poly.list, makeUniqueIDs = TRUE))
    row.names(sink.points.poly) <- as.character(1:length(sink.points.poly))
    sink.points.poly <- SpatialPolygonsDataFrame(sink.points.poly, data.frame(val=as.character(1:length(sink.points.poly))))
    sink.points.poly$val <- 1
    
    sink.points.poly <- gUnaryUnion(sink.points.poly, id = sink.points.poly$val)
    sink.points.poly <- rasterize(sink.points.poly, backgroundInformation)
    
    backgroundInformation.i <- raster::mask(backgroundInformation,sink.points.poly)
    backgroundInformation.i <- xyFromCell( backgroundInformation.i , Which(!is.na(backgroundInformation.i), cell=TRUE) )
    backgroundInformation.i <- backgroundInformation.i[sample(1:nrow(backgroundInformation.i),max(c(nrow(occurrenceRecords),nrow(backgroundInformation.i))),replace = F),]
    colnames(backgroundInformation.i) <- colnames(occurrenceRecords)
    
    absences.i <- generatePseudoAbsences(occurrenceRecords,rasterLayers,backgroundInformation.i,paRatio,type=type)
    
    dataModel <- data.frame(occ=1,raster::extract(rasterLayers,occurrenceRecords[,1:2]))
    dataModel <- rbind(dataModel,data.frame(occ=0,raster::extract(rasterLayers,absences.i)))
    
    observed <- dataModel$occ
    predicted <- predict( glm(dataModel, family=binomial(link = "logit")),dataModel,type="response")
    
    return(data.frame(buff=buff,performance= accuracyEstimate( observed ,predicted ,"tss")$tss))
    
  }
  
  stopCluster(cl); rm(cl)
  closeAllConnections()
  gc(reset=TRUE)
  
  ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2) }
  results[,2] <- as.numeric(ma(results[,2]))
  results <- results[complete.cases(results),]
  colnames(results) <- c("x","y")
  
  thresholdMax <- max(results$y) * (1-diff(range(results$y))) 
  threshold <- diff(range(results$y)) * 0.05  
  sorter <- sort(results$y, decreasing = T, index.return=T)$ix
  sorterTest <- 0
  
  for(t in (which.min(results[,2])):(nrow(results)-1)) {
    
    if( results[t+1,2] - results[t,2] <= threshold & results[t,2] >= thresholdMax ) { break }
    
  }
  
  bestResult <- results[t,"x"]
  bestResult <- ifelse(bestResult < maxBufferSize, bestResult , bestResult )
  
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Best buffer size distance: ",bestResult," degrees"))
  
  figure <- ggplot() +
    geom_vline(xintercept=bestResult, size=0.1) +
    geom_line(data = results, aes(x = x, y = y), size = 0.1, color="Black", linetype = "dashed") +
    geom_point(shape=19, data = results, aes(x = x, y = y), size = 1.75, color="Black") +
    ylab("Performance of models (TSS)") + xlab("Buffer size (degree)") + themePlot
  
  return( list(figure=figure, bestResult = bestResult) )
  
}

## -----------------------

generateBackgroundInformation <- function(rasterLayers,occurrenceRecords,n,spatialAutocorrelationVal) {
  
  options(warn=-1)
  
  shape <- as.data.frame(subset(rasterLayers,1), xy=T,na.rm=T)[,1:2]
  
  parallelChunk <- data.frame(from=seq(1,nrow(shape),by=5000)[-length(seq(1,nrow(shape),by=5000))],
                              to=seq(1,nrow(shape),by=5000)[-1]-1)
  parallelChunk[nrow(parallelChunk),2] <- nrow(shape)
  
  if(nrow(parallelChunk) == 0) { parallelChunk <- data.frame(from=1,to=nrow(shape))}
  
  shape.i <- data.frame()
  for(par in 1:nrow(parallelChunk)) {
    
    shape.i.temp <- spatialThinning(shape[parallelChunk[par,1]:parallelChunk[par,2],],spatialAutocorrelationVal,verbose=FALSE)
    shape.i <- rbind(shape.i,shape.i.temp)
    
  }
  
  n <- min(n,nrow(shape.i))
  
  pa.environment <- raster::extract( rasterLayers,shape.i )
  pa.environment <- scale(pa.environment)
  pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
  pa.environment <- kmeans(pa.environment,centers=n,iter.max = 10, nstart = 1)$cluster
  
  pa.clusters <- unique( pa.environment )
  finalKMClusters <- sapply(1:n,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
  shape.i <- shape.i[finalKMClusters,]
  
  return(shape.i)
  options(warn=0)
  
}

## -----------------------

generatePseudoAbsences <- function(occurrenceRecords,rasterLayers,paRatio,paType) {
  
  options(warn=-1)

  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) }
  if( paRatio > 1 ) { final.paRatio <- paRatio }
  
  if( nrow(occurrenceRecords) < 50 ) { final.paRatio <- 250 }

  if( paType != "sre" ) {
    
    if( paMinimum == 1) { paMinimum <- round(length(Which(!is.na(subset(rasterLayers,1)), cells=TRUE)) * 0.01) }
    
    final.paRatio <- min(nrow(backgroundInformation),max(final.paRatio,paMinimum))
    
    dataModel <- data.frame( Resp = c(rep(1,nrow(occurrenceRecords)) , rep(0,nrow(backgroundInformation) )) ,
                             Lon = c(occurrenceRecords[,"Lon"] , backgroundInformation[,"Lon"] ),
                             Lat = c(occurrenceRecords[,"Lat"] , backgroundInformation[,"Lat"]  ) )
    
    dataModel <- cbind(dataModel,raster::extract(rasterLayers,dataModel[,c("Lon","Lat")]))
    
  }
  
  if( paType == "random" ) {

    pseudoAbsences <- backgroundInformation
    final.paRatio <- min(nrow(pseudoAbsences) - 1,final.paRatio)
    
    pa.environment <- raster::extract( rasterLayers,pseudoAbsences )
    pa.environment <- scale(pa.environment)
    pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
    pa.environment <- kmeans(pa.environment,centers=final.paRatio,iter.max = 10, nstart = 1)$cluster
    
    pa.clusters <- unique( pa.environment )
    finalKMClusters <- sapply(1:final.paRatio,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
    pseudoAbsences <- pseudoAbsences[finalKMClusters,]
    
  }
  
  if( paType == "mess" ) {
    
    proj <- data.frame(long=backgroundInformation[,1],lat=backgroundInformation[,2],raster::extract(rasterLayers, backgroundInformation))
    cal <- data.frame(long=occurrenceRecords[,1],lat=occurrenceRecords[,2],raster::extract(rasterLayers, occurrenceRecords))
    toKeep <- intersect(which(apply(var(proj),1,sum) != 0) , which(apply(var(cal),1,sum) != 0))
    proj <- proj[,toKeep]
    cal <- cal[,toKeep]
    messSurface <- ecospat.mess( proj , cal , w="default")
    
    pseudoAbsences <- backgroundInformation[which(messSurface[,"MESSneg"] != 0),]
    final.paRatio <- min(nrow(pseudoAbsences) - 1,final.paRatio)
    
    pa.environment <- raster::extract( rasterLayers,pseudoAbsences )
    pa.environment <- pa.environment[,toKeep[-c(1,2)]-2]
    pa.environment <- scale(pa.environment)
    pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
    pa.environment <- kmeans(pa.environment,centers=final.paRatio,iter.max = 10, nstart = 1)$cluster
    
    pa.clusters <- unique( pa.environment )
    finalKMClusters <- sapply(1:final.paRatio,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
    pseudoAbsences <- pseudoAbsences[finalKMClusters,]
    
  }
  
  if( paType == "kernelDensity" ) {
    
    shape <- subset(rasterLayers,1)
    
    library(sm)
    bias <- cellFromXY(shape, as.matrix(occurrenceRecords))
    cells <- unique(sort(bias))
    kernelXY <- xyFromCell(shape, cells)
    samps <- as.numeric(table(bias))
    KDEsur <- sm.density(kernelXY, weights=samps, display="none", hmult=1, ylim=c(extent(shape)[3],extent(shape)[4]), xlim=c(extent(shape)[1],extent(shape)[2]), nbins=NA)
    KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
    KDErast <- SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
    KDErast <- raster(KDErast)
    KDErast <- raster::resample(KDErast, shape, method="bilinear")
    Bias.surface <- mask(KDErast,shape)
    Bias.surface[Bias.surface >= quantile(raster::extract(Bias.surface,occurrenceRecords), 0.1)] <- 0
    
    # plot(KDErastFigure)
    # plot(worldMap, add=T, col=NA)
    
    # histData <- raster::extract(Bias.surface,occurrenceRecords)
    # hist(histData, breaks=100)
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    sink.points.poly <- gBuffer( sink.points.poly, width=1 , byid=FALSE)
    Bias.surface[unlist(cellFromPolygon(Bias.surface,sink.points.poly))] <- 0
    
    pseudoAbsencesProbsBias <- raster::extract(Bias.surface,backgroundInformation)
    pseudoAbsencesVect <- sample( 1:nrow(backgroundInformation) , final.paRatio, replace = FALSE, prob = pseudoAbsencesProbsBias)
    pseudoAbsences <- backgroundInformation[pseudoAbsencesVect,]
    
  }
  
  if( paType == "sre") {
    
    library(biomod2)
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- dataRecordsNames
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    sink.points.poly$data <- rep(1,length(sink.points.poly))
    
    pseudoAbsences <- bm_PseudoAbsences(
      resp.var=terra::vect(sink.points.poly),
      expl.var= terra::rast(rasterLayers),
      nb.rep = 1,
      strategy = "sre",
      nb.absences = final.paRatio,
      sre.quant = 0 )
    
    pseudoAbsences <- pseudoAbsences$xy[which(is.na(pseudoAbsences$sp)),]
    
  }
  
  if( paType == "mahalanobis" ) {
    
    mahalanobis <- NULL
    
    library(adehabitatHS)
    rasterLayersSDF <- dataModel[complete.cases(dataModel),]
    toRemove <- c("None",names(sort(apply(rasterLayersSDF[,-c(1,2,3)],2,var, na.rm=T))))
    
    for(toRemove.i in 1:(length(toRemove)-1) ) {
      
      mahalanobis <- NULL
      rasterLayersSDF <- dataModel[complete.cases(dataModel),]
      
      if( toRemove.i == 1 ) { toKeep.col <- 1:ncol(rasterLayersSDF) }
      if( toRemove.i > 1 ) { toKeep.col <- which( ! names(rasterLayersSDF) %in% toRemove[2:toRemove.i] ) }
      
      rasterLayersSDF <- rasterLayersSDF[,toKeep.col]
      rasterLayersSDFClim <- rasterLayersSDF[,-c(1,2,3)]
      
      if( ncol(rasterLayersSDF) > 4 ){ rasterLayersSDFClim <- rasterLayersSDFClim[,which(apply(rasterLayersSDFClim,2,var) != 0)] }
      
      rasterLayersSDFClim <- scale(rasterLayersSDFClim)
      
      tryCatch( rasterLayersSDFClim <- dudi.pca(rasterLayersSDFClim, scannf=FALSE), error=function(e) { Error <<- TRUE })
      
      rasterLayersSDFClimContrib <- (rasterLayersSDFClim$eig - min(rasterLayersSDFClim$eig) ) / max((rasterLayersSDFClim$eig - min(rasterLayersSDFClim$eig) ) )
      rasterLayersSDFClimContrib <- rasterLayersSDFClimContrib / sum(rasterLayersSDFClimContrib)
      rasterLayersSDFClimContrib <- sapply(1:length(rasterLayersSDFClimContrib), function(x) { sum(rasterLayersSDFClimContrib[1:x]) })
      
      rasterLayersSDFClim <- rasterLayersSDFClim$tab[which(rasterLayersSDFClimContrib <= 0.95)]
      
      rasterLayersSDF <- SpatialPixelsDataFrame(rasterLayersSDF[,2:3], rasterLayersSDFClim)
      
      tryCatch( 
        mahalanobis <- raster(x=mahasuhab(rasterLayersSDF,SpatialPointsDataFrame(occurrenceRecords,data=data.frame(data=rep(1,nrow(occurrenceRecords)))), type="probability"))
        , error=function(e) { Error <<- TRUE })
      
      if( ! is.null(mahalanobis) ) { break }
      
    }
    
    if( ! is.null(mahalanobis)) {
      
      mahalanobisThresh <- raster::extract(mahalanobis,occurrenceRecords)
      mahalanobisThresh <- mahalanobisThresh[complete.cases(mahalanobisThresh)]
      mahalanobisThresh <- quantile(mahalanobisThresh, probs=0.025)
      mahalanobis[mahalanobis > mahalanobisThresh ] <- NA
      
      pseudoAbsencesProbs <- raster::extract(mahalanobis,backgroundInformation)
      pseudoAbsences <- backgroundInformation[which(!is.na(pseudoAbsencesProbs)),]
      
    }
    
    if( is.null(mahalanobis)) {
      
      pseudoAbsences <- pseudo.abs(coor=dataModel[,c("Lon","Lat")], status=dataModel$Resp, distance=1 , strategy = "circles", add.pres = FALSE, create.dataset=FALSE, species.name = "Species", plot = FALSE, acol = "grey80", pcol = "red")
      pseudoAbsences <- dataModel[pseudoAbsences,c("Lon","Lat")]
      
    }
    
    dataModel <- data.frame( Resp = c(rep(1,sum(dataModel$Resp == 1) ) , rep(0,nrow(pseudoAbsences) )) ,
                             Lon = c(dataModel[dataModel$Resp == 1,"Lon"] , pseudoAbsences[,1] ),
                             Lat = c(dataModel[dataModel$Resp == 1,"Lat"] , pseudoAbsences[,2] ) )
    
    final.paRatio <- min(nrow(pseudoAbsences) - 1,final.paRatio)
    pa.environment <- raster::extract( rasterLayers,pseudoAbsences )
    pa.environment <- scale(pa.environment)
    pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
    pa.environment <- kmeans(pa.environment,centers=final.paRatio,iter.max = 10, nstart = 1)$cluster
    
    pa.clusters <- unique( pa.environment )
    finalKMClusters <- sapply(1:final.paRatio,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
    pseudoAbsences <- pseudoAbsences[finalKMClusters,]
    
  }
  
  plot(pseudoAbsences, col="gray")
  points(occurrenceRecords, col="black")
  colnames(pseudoAbsences) <- colnames(occurrenceRecords)
  return(pseudoAbsences)
  
  options(warn=0)
  
}

## -----------------------

accuracyEstimate <- function(observed,predicted,cvIndex) {
  
  predicted.accuracy <- accuracyEstimateCalc( observed ,predicted ,cvIndex)

  return(predicted.accuracy)
  
}

## -----------------------

accuracyEstimateCalc <- function(observed,predicted,cvIndex) {
  
  options(warn=-1)
  
  for( rep in 1:3) {
    
    if( rep == 1) { sequenceT <- seq(0.01,1,length.out=100) }
    if( rep == 2) { sequenceT <- seq(predicted.accuracy$threshold - 0.01 ,predicted.accuracy$threshold + 0.01,length.out=100) }
    if( rep == 3) { sequenceT <- seq(predicted.accuracy$threshold - 0.001 ,predicted.accuracy$threshold + 0.001,length.out=100) }
    
    sequenceT[sequenceT <= 0] <- 0.000001
    sequenceT[sequenceT >= 1] <- 1 - 0.000001
    
    predicted.accuracy <- SDMTools::accuracy( observed , predicted , threshold = sequenceT )
    predicted.accuracy$deviance <- modEvA::Dsquared(obs = observed, pred=predicted, family = "binomial")
    predicted.accuracy$aicc <- AIC(glm(observed~predicted))
    
    # boyce
    boyce.p <- predicted[observed==1]
    boyce.p <- boyce.p[!is.na(boyce.p)]
    boyce.a <- predicted
    boyce.a <- boyce.a[!is.na(boyce.a)]
    
    if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01) }
    if( length(unique(boyce.p)) == 1) { boyce.p <- sapply(boyce.p, function(x) {   x + sample(c(0.001,-0.001), 1)}) }
    
    boyce.value <- ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE,rm.duplicate = TRUE)
    predicted.accuracy$boyce <- boyce.value$cor
    
    if(cvIndex == "boyce" ) { predicted.accuracy <- predicted.accuracy[which.min(abs(predicted.accuracy$threshold - min(  boyce.value$HS[ max( ifelse( length(which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))]) > 0 , which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))] , 1)) ]  ))),] }
    if(cvIndex == "area" ) { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity > 0.95)),] }
    if(cvIndex == "area10th" ) { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity > 0.9)),] }
    if(cvIndex == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
    if(cvIndex == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    if(cvIndex == "aicc") { predicted.accuracy <- predicted.accuracy[which.min(predicted.accuracy$AUC),] }
    if(cvIndex == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    
    if( is.na(predicted.accuracy$boyce ) & ( predicted.accuracy$sensitivity + predicted.accuracy$specificity ) == 2 ) { predicted.accuracy$boyce <- 1 }
    
  }
  
  ## -------
  
  if(nrow(predicted.accuracy) > 1) { stop("Error :: 922") }
  
  ## -------
  
  predicted.accuracy <- data.frame( boyce = predicted.accuracy$boyce,
                                    threshold = predicted.accuracy$threshold ,
                                    auc = predicted.accuracy$AUC ,
                                    specificity = predicted.accuracy$specificity ,
                                    sensitivity = predicted.accuracy$sensitivity ,
                                    tss = predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 ,
                                    area = predicted.accuracy$sensitivity ,
                                    area10th = predicted.accuracy$sensitivity ,
                                    aicc = predicted.accuracy$aicc,
                                    deviance = predicted.accuracy$deviance )
  
  options(warn=0)
  
  return(predicted.accuracy)
  
}

## -----------------------

accuracyPredicted <- function(predictedDistribution,speciesData,cvIndex) {
  
  observed <- speciesData$PA
  predicted <- raster::extract(predictedDistribution,speciesData[,dataRecordsNames])
  predicted.accuracy <- accuracyEstimate(observed,predicted,cvIndex)
  return(predicted.accuracy)
  
}

## -----------------------

predictDistribution <- function(rasters,model,reclassToOne) {
  
  if( class(model)[1] == "MaxEnt" ) {
    
    predicted.distribution <- predict( model , rasters )
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=TRUE)
    
  }   
  
  if( class(model)[1] == "gbm" ) {
    
    num.tress <- model$gbm.call$best.trees
    if(is.null(num.tress)) { num.tress <- length(model$trees) }
    predicted.distribution <- predict( rasters , model , n.trees=num.tress,type="response")
    
  }
  
  if( class(model)[1] == "xgb.Booster" ) {
    
    xgd_data <- as.data.frame(rasters, na.rm=T)
    xgd_data <- xgd_data[,model$feature_names]
    xgd_data = xgb.DMatrix(data = data.matrix(xgd_data), label = rep(0,nrow(xgd_data)))
    predicted.distribution <- subset(rasters,1)
    predicted.distribution[] <- NA
    predicted.distribution[cellFromXY( rasters, as.data.frame(rasters, xy=T, na.rm=T)[, c("x","y")])] <- predict(model,xgd_data, type = "response")
    
  }
  
  if( class(model)[1] == "mboost" ) {
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    options(warn=-1)
    predicted.distribution <- predict( rasters , model)
    predicted.distribution <- logit2prob(predicted.distribution)
    
    options(warn=0)
    
  }
  
  if( class(model)[1] == "list" ) {
    
    shape <- subset(rasters,1)
    shape[!is.na(shape)] <- 1
    cellsPredict <- Which(!is.na(shape),cells=TRUE)
    cellsPredictValues <- as.matrix(rasters[cellsPredict])
    
    correctLayers <- names(rasters)[which(dataLayersMonotonocity == -1)]
    for(correctLayers.i in correctLayers) {
      cellsPredictValues[,correctLayers.i] <- cellsPredictValues[,correctLayers.i] * (-1)
    }
    
    predicted.distribution <- monmlp.predict( cellsPredictValues , model)
    shape[cellsPredict] <- as.numeric(predicted.distribution)
    predicted.distribution <- shape
    
  }
  
  predicted.distribution[predicted.distribution < 0] <- 0
  predicted.distribution[predicted.distribution > 1] <- 1
  
  
  if(reclassToOne) {
    predicted.distribution <- predicted.distribution + ( min(getValues(predicted.distribution),na.rm=T) * (-1))
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=T)
  }
  
  return(predicted.distribution)
  
}

## -----------------------






SDMmodelPlot <- function(model,rasters,variable.to.plot,speciesData,extrapolate,plotExtremeThreshold,plotRateChangeThreshold,yLimits) {
  
  if( ! is.null(yLimits) ) { yLimits <- round(yLimits + c(-0.055,0.055), digits = 1) }
  if( is.null(yLimits) ) { yLimits <- c(0,1) }
  
  yLimits[yLimits < 0] <- 0
  yLimits[yLimits > 1] <- 1
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  data.to.plot.i <- raster::extract(subset(rasters,which(names(rasters) == variable.to.plot)),speciesData[,2:3])
  data.to.plot <- raster::extract(rasters,speciesData[,2:3])
  data.to.plot <- apply(data.to.plot,2,mean)
  data.to.plot <- do.call("rbind", replicate(100, data.to.plot, simplify = FALSE))
  data.to.plot[,variable.to.plot] <- seq(min(data.to.plot.i),max(data.to.plot.i),length.out=100)
  data.to.plot <- as.data.frame(data.to.plot)
  
  matrixEffect <- NULL
  
  for( m in 1:length(model@models)) {

    model.i <- model@models[[m]]@model
    if( "model" %in% slotNames(model.i) ) { model.i <- model.i@model }

    if ( class(model.i) == "xgb.Booster") {

      matrixEffect.m = xgb.DMatrix(data = data.matrix(data.to.plot), label = rep(0,nrow(data.to.plot)))
      matrixEffect.m <- data.frame(Effect= predict(model.i,matrixEffect.m, type = "response") ) 
      
    }
    
    if ( class(model.i) != "xgb.Booster") {
      
      matrixEffect.m <- data.frame(Effect= predict(model.i,data.to.plot)) 
      
    }
      
    matrixEffect.m <- logit2prob(matrixEffect.m)
      
    if(!is.null(matrixEffect)) { matrixEffect <- cbind(matrixEffect,data.frame(m=matrixEffect.m)) }
    if(is.null(matrixEffect)) { matrixEffect <- data.frame(m=matrixEffect.m) }
      
    }
    
  data.to.plot <- data.frame(Variable=data.to.plot[,variable.to.plot],value=apply(matrixEffect,1,mean,na.rm=T))
  data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
  
  if(data.to.plot[nrow(data.to.plot),2] > data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] < data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
  if(data.to.plot[nrow(data.to.plot),2] < data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] > data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
   
  ## -----------------------
  
  speciesDataUse <- raster::extract(subset(rasters,which(names(rasters) == variable.to.plot)),speciesData[speciesData[,1]==1,2:3])
  speciesDataUse <- data.frame(x=speciesDataUse,y=min(data.to.plot$value))
  speciesDataUse <- speciesDataUse[sort(speciesDataUse[,1], index.return =T)$ix,]
  
  ## -----------------------
  
  if( ! extrapolate ) { data.to.plot <- data.to.plot[data.to.plot$Variable >= min(speciesDataUse$x, na.rm=T) & data.to.plot$Variable <= max(speciesDataUse$x, na.rm=T) , ]  }
  
  ## -----------------------
  
  partialPlot <- ggplot() +
    geom_point( data=speciesDataUse,aes(x=x, y=y), size=2,shape=15, fill="gray", color="gray") +
    geom_line( data=data.to.plot,aes(x=Variable, y=value),color="Black", size=0.5) +
    themePlot +
    xlab(variable.to.plot) + ylab("Effect on response (probability)")
  
  ## -----------------------
  
  span.vals <- data.to.plot$value
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
  span.vals.var <- abs(span.vals.var)
  rateChangeThreshold <- which.max(span.vals.var)
  rateChangeThreshold <- data.to.plot[rateChangeThreshold,1]
  
  ## --
  
  if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min" }
  
  if( tail == "Min") {
    
    usedThreshold <- min(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.05,na.rm=T)
    extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value) + quantile(data.to.plot$value,0.05),"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value),"Variable"]) }
    
  }
  
  if( tail == "Max") {
    
    usedThreshold <- max(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.95,na.rm=T)
    extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) - quantile(data.to.plot$value,0.05) ,"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) ,"Variable"]) }
    
  }
  
  if( data.to.plot[1,2] == data.to.plot[nrow(data.to.plot),2]  ) { rateChangeThreshold <- NA }
  
  if( plotExtremeThreshold) {
    partialPlot <- partialPlot + geom_vline(xintercept = extremeThreshold,linetype = "dashed")
  }
  if( plotRateChangeThreshold ) {
    partialPlot <- partialPlot + geom_vline(xintercept = rateChangeThreshold,linetype = "dotted")
  }
  
  if(length(rateChangeThreshold) == 0) { rateChangeThreshold <- NA }
  if(length(extremeThreshold) == 0) { extremeThreshold <- NA }
  if(length(usedThreshold) == 0) { usedThreshold <- NA }
  if(length(usedThresholdQ95) == 0) { usedThresholdQ95 <- NA }
  
  tippingPoints <- data.frame( maxRateChangeThreshold=rateChangeThreshold,extremeThreshold=extremeThreshold,usedThreshold=usedThreshold,usedThresholdQ95=usedThresholdQ95 )
  
  return(list(partialPlot=partialPlot,tippingPoints=tippingPoints, partialPlotData=data.to.plot))
  
}

## -----------------------

modelPlot <- function(model,rasters,distribution,variable.to.plot,speciesData,extrapolate,plotExtremeThreshold,plotRateChangeThreshold,yLimits) {
  
  if( ! is.null(yLimits) ) { yLimits <- round(yLimits + c(-0.055,0.055), digits = 1) }
  if( is.null(yLimits) ) { yLimits <- c(0,1) }
  
  yLimits[yLimits < 0] <- 0
  yLimits[yLimits > 1] <- 1
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  if(class(model) ==  "gbm") {
    
    model.predictor <- model$var.names[variable.to.plot]
    
    data.to.plot <- plot(model,variable.to.plot,return.grid = TRUE)
    data.to.plot <- data.frame(Variable = data.to.plot[,1] , Effect=data.to.plot[,2] )
    data.to.plot[,2] <- logit2prob(data.to.plot[,2])
    
  }
  
  if(class(model) !=  "gbm") {
    
    model.predictor <- names(rasterLayers)[variable.to.plot]
    
    data.to.plot <- getValues(rasters)
    data.to.plot <- data.to.plot[complete.cases(data.to.plot),]
    data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot), ifelse( nrow(data.to.plot) > 10000 , 10000 , nrow(data.to.plot)  ) ,replace=FALSE),]
    
    means <- colMeans( data.to.plot[,-which( colnames(data.to.plot) == model.predictor )] )
    for( m in names(means) ) { data.to.plot[,m] <- means[m] }
    
    if ( class(model) == "xgb.Booster") {
      
      data.to.plot.xgd = xgb.DMatrix(data = data.to.plot, label = rep(0,nrow(data.to.plot)))
      data.to.plot <- data.frame(Variable = data.to.plot[,model.predictor] , Effect= predict(model,data.to.plot.xgd, type = "response") ) 
      
    } 
    
    if ( class(model) != "xgb.Booster") {
      
      data.to.plot <- as.data.frame(data.to.plot) 
      data.to.plot <- data.frame(Variable = data.to.plot[,model.predictor] , Effect= predict(model,data.to.plot)) 
      
    }
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
    data.to.plot[,2] <- logit2prob(data.to.plot[,2])
    
  }
  
  ## -----------------------
  
  speciesDataUse <- raster::extract(subset(rasters,which(names(rasters) == model.predictor)),speciesData[speciesData[,1]==1,2:3])
  speciesDataUse <- data.frame(x=speciesDataUse,y=min(data.to.plot$Effect))
  speciesDataUse <- speciesDataUse[sort(speciesDataUse[,1], index.return =T)$ix,]
  
  ## -----------------------
  
  if( ! extrapolate ) { data.to.plot <- data.to.plot[data.to.plot$Variable >= min(speciesDataUse$x, na.rm=T) & data.to.plot$Variable <= max(speciesDataUse$x, na.rm=T) , ]  }
  
  if( distribution == "binomial" ) { speciesDataUse$y <- min(data.to.plot$Effect) }
  
  ## -----------------------
  
  partialPlot <- ggplot() +
    geom_point( data=speciesDataUse,aes(x=x, y=y), size=2,shape=15, fill="gray", color="gray") +
    geom_line( data=data.to.plot,aes(x=Variable, y=Effect),color="Black", size=0.5) +
    themePlot +
    xlab(model.predictor) + ylab("Effect on response (probability)")
  
  # if(distribution == "binomial") { partialPlot <- partialPlot + ylim(min(yLimits), max(yLimits)) }
  
  ## -----------------------
  
  span.vals <- data.to.plot$Effect
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
  span.vals.var <- abs(span.vals.var)
  rateChangeThreshold <- which.max(span.vals.var)
  rateChangeThreshold <- data.to.plot[rateChangeThreshold,1]
  
  ## --
  
  if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min" }
  
  if( tail == "Min") {
    
    usedThreshold <- min(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.05,na.rm=T)
    extremeThreshold <- max(data.to.plot[data.to.plot$Effect <= min(data.to.plot$Effect),"Variable"])
    
  }
  
  if( tail == "Max") {
    
    usedThreshold <- max(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.95,na.rm=T)
    extremeThreshold <- min(data.to.plot[data.to.plot$Effect <= min(data.to.plot$Effect),"Variable"])
    
  }
  
  if( data.to.plot[1,2] == data.to.plot[nrow(data.to.plot),2]  ) { rateChangeThreshold <- NA }
  
  if( plotExtremeThreshold) {
    partialPlot <- partialPlot + geom_vline(xintercept = extremeThreshold,linetype = "dashed")
  }
  if( plotRateChangeThreshold ) {
    partialPlot <- partialPlot + geom_vline(xintercept = rateChangeThreshold,linetype = "dotted")
  }
  
  if(length(rateChangeThreshold) == 0) { rateChangeThreshold <- NA }
  if(length(extremeThreshold) == 0) { extremeThreshold <- NA }
  if(length(usedThreshold) == 0) { usedThreshold <- NA }
  if(length(usedThresholdQ95) == 0) { usedThresholdQ95 <- NA }
  
  tippingPoints <- data.frame( maxRateChangeThreshold=rateChangeThreshold,extremeThreshold=extremeThreshold,usedThreshold=usedThreshold,usedThresholdQ95=usedThresholdQ95 )
  
  return(list(partialPlot=partialPlot,tippingPoints=tippingPoints))
  
}

## -----------------------

correctLayer <- function(rasters,name,type,value.is,value.will) {		
  
  name.i <- which(names(rasters) == name)
  temp.r <- rasters[[name.i]]		
  rasters <- dropLayer(rasters,name.i)
  
  if(type =="he") { temp.r[temp.r >= value.is] <- value.will }		
  if(type =="le") { temp.r[temp.r <= value.is] <- value.will }		
  if(type =="h") { temp.r[temp.r > value.is] <- value.will }		
  if(type =="l") { temp.r[temp.r < value.is] <- value.will }		
  
  names(temp.r) <- name
  rasters <- stack( rasters , temp.r )
  
  return(rasters)		
  
}		

## -----------------------

reclassifyPredicted <- function(predicted.distribution,speciesData,method,reclassThreshold) {
  
  presences <- speciesData[speciesData$PA == 1,c("Lon","Lat")]  
  absences <- speciesData[speciesData$PA == 0,c("Lon","Lat")]  
  
  if(method == "directReclass") {
    
    
    predicted.distribution[predicted.distribution > reclassThreshold] <- 1
    predicted.distribution[predicted.distribution <= reclassThreshold] <- 0
    
  }
  
  if(method == "maxTSS") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold =100 )
    predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  
  if(method == "minAREA") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold = 1000 )
    predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity >= reclass.threshold,]
    predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  
  return(predicted.distribution)
  
}

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

library(SDMtune)


prepareModelData <- function(p,a,env) {
  
  m <- subset(env,1)
  
  p.i <- raster::extract(m,p)
  p <- p[which(!is.na(p.i)),]
  
  p <- xyFromCell(m,cellFromXY(m,p))
  p <- unique(p)
  
  a.i <- raster::extract(m,a)
  a <- a[which(!is.na(a.i)),]
  
  a <- xyFromCell(m,cellFromXY(m,a))
  a <- unique(a)
  
  return(prepareSWD(species = "Model species", p = p, a = a, env = env))
  
}

getBlocks <- function(modelData, rasters) {
  
  # get.block(modelData@coords[modelData@pa == 1,], modelData@coords[modelData@pa == 0,])
  # a <- randomFolds(modelData, 4)
  
  require(sp)
  require(sf)
  require(rgeos)
  
  occ <- modelData@coords
  occ <- SpatialPoints(coords = occ[, 1:2], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") )
  
  xy.sp <- st_as_sf(occ, crs = 4326)
  grd <- sf::st_make_grid(xy.sp, n=c(5,2))
  grd <- as_Spatial(grd)
  grd$id <- 1:length(grd)
  
  # plot(grd)
  # points(modelData@coords[which(modelData@pa == 1),], col="black")
  
  train <- matrix(FALSE,ncol=length(grd),nrow=length(occ))
  test <- matrix(FALSE,ncol=length(grd),nrow=length(occ))
  
  for( i in 1:length(grd)) {
    
    test.vect <- which(!is.na(over(occ,grd[i,])[,1]))
    test[ test.vect , i] <- TRUE
    train[ which(!test[ , i]), i] <- TRUE
    
  }
  
  train.vect.p <- which(apply(train[which(modelData@pa==1),],2,sum) > 0 )
  test.vect.p <- which(apply(test[which(modelData@pa==1),],2,sum) > 0 )
  train.vect.a <- which(apply(train[which(modelData@pa==0),],2,sum) > 0 )
  test.vect.a <- which(apply(test[which(modelData@pa==0),],2,sum) > 0 )
  
  folds <- intersect(intersect(intersect(train.vect.p,test.vect.p),train.vect.a),test.vect.a)
  folds <- list( grd=grd, train = train[,folds] , test = test[,folds] )
  
  if( length(intersect(train.vect,test.vect)) == 0) { 
    folds <- list( grd=grd, train = randomFolds(modelData, 5)$train , test = randomFolds(modelData, 5)$train )
  }
  
  
  return(folds)
  
}












.create_sdmtune_output <- function(models,
                                   metric,
                                   train_metric,
                                   val_metric) {
  
  tunable_hypers <- getTunableArgs(models[[1]])
  l <- length(tunable_hypers)
  labels <- c(tunable_hypers, .get_sdmtune_colnames(metric))
  
  res <- matrix(nrow = length(models), ncol = length(labels))
  colnames(res) <- labels
  fcs <- vector("character", length = length(models))
  distrs <- vector("character", length = length(models))
  
  for (i in seq_along(models)) {
    
    if (inherits(models[[i]], "SDMmodel")) {
      m <- models[[i]]
    } else {
      m <- models[[i]]@models[[1]]
    }
    
    for (j in 1:l) {
      if (tunable_hypers[j] == "distribution") {
        distrs[i] <- slot(m@model, tunable_hypers[j])
      } else if (tunable_hypers[j] == "fc") {
        fcs[i] <- m@model@fc
      } else {
        res[i, tunable_hypers[j]] <- slot(m@model, tunable_hypers[j])
      }
    }
    
    res[i, l + 1] <- train_metric[i, 2]
    
    if (metric != "aicc")
      res[i, l + 2] <- val_metric[i, 2]
  }
  
  if (metric != "aicc") {
    res[, l + 3] <- res[, l + 1] - res[, l + 2]
  } else {
    res[, l + 2] <- res[, l + 1] - min(res[, l + 1])
  }
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  if ("distribution" %in% tunable_hypers) {
    res$distribution <- distrs
  }
  
  if ("fc" %in% tunable_hypers) {
    res$fc <- fcs
  }
  
  SDMtune(results = res,
          models = models)
}

.get_train_args <- function(model) {
  
  args <- list(data = model@data)
  
  if (inherits(model, "SDMmodelCV")) {
    args$folds <- model@folds
    args$progress <- FALSE
    model <- model@models[[1]]@model
  } else {
    args$folds <- NULL
    model <- model@model
  }
  
  args$method <- class(model)
  
  if (args$method == "Maxent") {
    args$fc <- model@fc
    args$reg <- model@reg
    args$iter <- model@iter
  } else if (args$method == "Maxnet") {
    args$fc <- model@fc
    args$reg <- model@reg
  } else if (args$method == "ANN") {
    args$size <- model@size
    args$decay <- model@decay
    args$rang <- model@rang
    args$maxit <- model@maxit
  } else if (args$method == "RF") {
    args$mtry <- model@mtry
    args$ntree <- model@ntree
    args$nodesize <- model@nodesize
  } else if (args$method == "MBOOST") {
    args$distribution <- model@distribution
    args$df <- model@df
    args$shrinkage <- model@shrinkage
    args$mstop <- model@mstop
  } else if (args$method == "XGBOOST") {
    args$distribution <- model@distribution
    args$gamma <- model@gamma
    args$shrinkage <- model@shrinkage
    args$depth <- model@depth
    args$nrounds <- model@nrounds
  } else {
    args$distribution <- model@distribution
    args$n.trees <- model@n.trees
    args$interaction.depth <- model@interaction.depth
    args$shrinkage <- model@shrinkage
    args$bag.fraction <- model@bag.fraction
  }
  
  args
  
}

.get_hypers_grid <- function(model,hypers) {
  
  # Create data frame with all possible combinations of hyperparameters
  tunable_args <- .get_train_args(model)[getTunableArgs(model)]
  tunable_args[names(hypers)] <- hypers
  
  expand.grid(tunable_args, stringsAsFactors = FALSE)
}

.args_name <- function(x) {
  switch(x,
         "trainANN" = c("data", "size", "decay", "rang", "maxit"),
         "trainMBOOST" = c("data", "distribution", "df", "shrinkage", "mstop"),
         "trainXGBOOST" = c("data", "distribution", "gamma", "shrinkage", "depth", "nrounds"),
         "trainBRT" = c("data", "distribution", "n.trees", "interaction.depth", "shrinkage", "bag.fraction"),
         "trainMaxent" = c("data", "reg", "fc", "iter"),
         "trainMaxnet" = c("data", "reg", "fc"),
         "trainRF" = c("data", "mtry", "ntree", "nodesize")
  )
}


trainBRT <- function(data, distribution = "bernoulli", n.trees = 1000, interaction.depth = 1, shrinkage = 0.1, bag.fraction = 0.5 ) {
  
  result <- SDMmodel(data = data)
  
  df <- data@data
  df <- cbind(pa = data@pa, df)
  
  if( exists("monotonicity") ) {
    
    monotonicity.i <- as.numeric(monotonicity)[match(colnames(df)[-1], colnames(monotonicity) )]
    model <- gbm::gbm(pa ~ ., data = df, distribution = distribution,
                      n.trees = n.trees, interaction.depth = interaction.depth,
                      shrinkage = shrinkage, bag.fraction = bag.fraction, var.monotone=monotonicity.i)
    
  }
  if( ! exists("monotonicity")) {
    
    model <- gbm::gbm(pa ~ ., data = df, distribution = distribution,
                      n.trees = n.trees, interaction.depth = interaction.depth,
                      shrinkage = shrinkage, bag.fraction = bag.fraction)
    
  }
  
  model_object <- BRT(n.trees = n.trees, distribution = distribution,
                      interaction.depth = interaction.depth,
                      shrinkage = shrinkage, bag.fraction = bag.fraction,
                      model = model)
  
  result@model <- model_object
  
  return(result)
}

trainMBOOST <- function(data, distribution = "bernoulli", df = 5, shrinkage = 0.25, mstop = 50 ) {

  result <- SDMmodel(data = data)
  
  data.df <- data@data
  data.df <- cbind(pa = data@pa, data.df)

  cave <- function(x) { length(unique(x)) / length(x) }
  varLayers <- which( apply( data.df[-1] ,2,cave) > 0.025)
  data.df <- data.df[,c("pa",names(varLayers))]
  
  # ---------
  # Dummmy model
  
  model <- mboost(as.formula(paste("pa", paste(colnames(data.df)[-1] , collapse = " + "), sep = " ~ ")) , data= data.df)

  # ---------
  
  response <- "pa"
  data.df[,"pa"] <- as.factor(data.df[,"pa"])
  predictors <- colnames(data.df[-1])  
  
  if( exists("monotonicity") ) {
    
    monotonicity.i <- as.numeric(monotonicity)[match(predictors, colnames(monotonicity) )]
    monotonicity.i <- as.character(monotonicity.i)
    monotonicity.i[monotonicity.i == "-1"] <- "decreasing"
    monotonicity.i[monotonicity.i == "1"] <- "increasing"
    rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", monotonicity.i,"\", df = ",df,")", sep = "")),collapse = " + ")
  }
  
  if( ! exists("monotonicity")) {
    rhs <- paste(c(paste("bmono(", predictors, ", df = ",df,")", sep = "")),collapse = " + ")
  }

  fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
  ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
  
  tryCatch( suppressWarnings( model <- mboost(fm_mono, data = data.df, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" ))) , error=function(e) { a<<-1 })
  
  model_object <- MBOOST(df = df,
                         distribution = distribution,
                         shrinkage = shrinkage, 
                         mstop = mstop,
                         model = model)
  
  result@model <- model_object
  return(result)
  
}

trainXGBOOST <- function(data, distribution = "bernoulli", gamma = 1, shrinkage = 0.25, depth = 2, nrounds = 50 ) {
  
  result <- SDMmodel(data = data)
  
  data.df <- data@data
  data.df <- cbind(pa = data@pa, data.df)
  
  xgb_train = xgb.DMatrix(data = data.matrix(data.df[,-1]), label = data.df[,1])
  
  if( exists("monotonicity") ) {
    monotonicity.i <- as.numeric(monotonicity)[match(colnames(data.df[,-1]), colnames(monotonicity) )]
    model <- xgboost(data = xgb_train, monotone_constraints=monotonicity.i, max.depth = depth, gamma=gamma, nrounds = nrounds , verbose = 0, objective="binary:logistic")
  }
  
  if( ! exists("monotonicity")) {
    
    model <- xgboost(data = xgb_train, max.depth = depth, gamma=gamma, nrounds = nrounds , verbose = 0, objective="binary:logistic")
    
  }
  
  model_object <- XGBOOST(depth = depth,
                          distribution = distribution,
                          shrinkage = shrinkage, 
                          gamma = gamma, 
                          nrounds = nrounds,
                          model = model)
  
  result@model <- model_object
  return(result)
  
}

predict.xgboost <- function(model , data , type ) {
  
  if("model" %in% slotNames(model)) { model <- model@model }
  if("model" %in% slotNames(model)) { model <- model@model }

    if( class(data) == "RasterStack") { 
    xgd_data <- as.data.frame(data, na.rm=T)
    xgd_data <- xgd_data[,model$feature_names]
    xgd_data = xgb.DMatrix(data = data.matrix(xgd_data), label = rep(0,nrow(xgd_data)))
    predicted.distribution <- subset(data,1)
    predicted.distribution[] <- NA
    predicted.distribution[cellFromXY( data, as.data.frame(data, xy=T, na.rm=T)[, c("x","y")])] <- predict(model,xgd_data, type = type)
  } else {
    xgd_data <- as.data.frame(data, na.rm=T)
    xgd_data <- xgd_data[,model$feature_names]
    xgd_data = xgb.DMatrix(data = data.matrix(xgd_data), label = rep(0,nrow(xgd_data)))
    predicted.distribution <- predict(model,xgd_data, type = type)
  }
  return(predicted.distribution)
  
}

setOldClass("mboost")
MBOOST <- setClass("MBOOST",
                   slots = c(
                     distribution = "character",
                     df = "numeric",
                     shrinkage = "numeric",
                     mstop = "numeric",
                     model = "mboost")
)

#' @param object MBOOST object
#' @rdname MBOOST-class
setMethod("show",
          signature = "MBOOST",
          definition = function(object) {
            cat("Class            :", class(object), "\n")
            cat("distribution          :", object@distribution, "\n")
            cat("df          :", object@df, "\n")
            cat("mstop:", object@mstop, "\n")
            cat("shrinkage        :", object@shrinkage, "\n")
          })

setOldClass("xgb.Booster")
XGBOOST <- setClass("XGBOOST",
                    slots = c(
                      distribution = "character",
                      depth = "numeric",
                      shrinkage = "numeric",
                      gamma = "numeric",
                      nrounds = "numeric",
                      model = "xgb.Booster")
)

#' @param object XGBOOST object
#' @rdname XGBOOST-class
setMethod("show",
          signature = "XGBOOST",
          definition = function(object) {
            cat("Class            :", class(object), "\n")
            cat("distribution          :", object@distribution, "\n")
            cat("depth          :", object@depth, "\n")
            cat("nrounds          :", object@nrounds, "\n")
            cat("gamma:", object@gamma, "\n")
            cat("shrinkage        :", object@shrinkage, "\n")
          })

setGeneric("predict", function(object, ...)
  standardGeneric("predict")
)
setMethod(
  f = "predict",
  signature = "MBOOST",
  definition = function(object,
                        data,
                        type) {
    
    predict(object@model,
            newdata = data,
            type = "response")
  }
)


setGeneric("predict", function(object, ...)
  standardGeneric("predict")
)
setMethod(
  f = "predict",
  signature = "SDMmodel",
  definition = function(object,
                        data,
                        type = NULL,
                        clamp = TRUE,
                        filename = "",
                        overwrite = FALSE,
                        wopt = list(),
                        format = "",
                        extent = NULL,
                        progress = "",
                        ...) {

    model <- object
    if( "model" %in% slotNames(model) ) { model <- model@model }
    if( "model" %in% slotNames(model) ) { model <- model@model }
    
    # TODO: Remove with version 2.0.0
    if (!identical(progress, ""))
      cli::cli_warn(
        c("!" = "Argument {.field progress} is deprecated",
          "i" = "It will be removed in future releases")
      )
    
    # TODO: Remove with version 2.0.0
    if (format != "")
      cli::cli_warn(c(
        "!" = paste("The argument {.val format} is deprectated and will be",
                    "ignored. Use {.val wopt} instead and see {.val Details}",
                    "in {.fun terra::writeRaster}")
      ))
    
    if (inherits(data, "Raster") | inherits(data, "RasterStack") | inherits(data, "SpatRaster") ) {
      
      data <- stack(data)
      
      if(class(model) == "xgb.Booster") {
        pred <- predict.xgboost(model = model,data = data,type="response")
      }
      if(class(model) == "gbm") {
        pred <- predict(data,  model , n.trees=model$n.trees,type="response")
      }
      if(class(model) == "mboost") {
        pred <- predict(data,  model,type="response")
      }
      
    } else if (inherits(data, "SWD")) {
      data <- data@data
      
      if(class(model) == "xgb.Booster") {
        pred <- predict.xgboost(model = model,data = data,type="response")
      }
      if(class(model) == "gbm") {
        pred <- predict( model , data , n.trees=model$n.trees,type="response")
      }
      if(class(model) == "mboost") {
        pred <- predict(model, data,type="response")
      }
      
      pred <- as.vector(pred)
      
    } else if (inherits(data, "data.frame")) {
      data <- data
      
      if(class(model) == "xgb.Booster") {
        pred <- predict.xgboost(model = model,data = data,type="response")
      }
      if(class(model) == "gbm") {
        pred <- predict( model , data , n.trees=model$n.trees,type="response")
      }
      if(class(model) == "mboost") {
        pred <- predict(model, data,type="response")
      }
      
      pred <- as.vector(pred)
    } else {
      cli::cli_abort(c(
        "!" = paste("{.var data} must be an object of class",
                    "{.cls data.frame}, {.cls SWD} or {.cls Raster}"),
        "x" = "You have supplied a {.cls {class(data)}} instead."
      ))
    }
    
    return(pred)
  }
)


train <- function (method, data, folds = NULL, verbose = TRUE, ...) {
  
  l <- length(method)
  output <- vector("list", length = l)
  for (i in 1:l) {
    m <- match.arg(method[i], c("Maxent", "Maxnet", "ANN", "RF", "BRT","MBOOST","XGBOOST"))
    
    func <- paste0("train", m)
    ea <- list(...)
    if (is.null(folds)) {
      argus <- c(data = data, ea[names(ea) %in% .args_name(func)])
      output[[i]] <- do.call(func, args = argus)
    }
    else {
      folds <- .convert_folds(folds, data)
      k <- ncol(folds[[1]])
      if (verbose) {
        pb <- progress::progress_bar$new(format = "Cross Validation [:bar] :percent in :elapsedfull", 
                                         total = k, clear = FALSE, width = 60, show_after = 0)
        pb$tick(0)
      }
      models <- vector("list", length = k)
      for (j in 1:k) {
        train <- .subset_swd(data, folds$train[, j])
        argus <- c(data = train, ea[names(ea) %in% .args_name(func)])
        models[[j]] <- do.call(func, args = argus)
        if (verbose) 
          pb$tick(1)
      }
      output[[i]] <- SDMmodelCV(models = models, data = data, 
                                folds = folds)
    }
  }
  if (l == 1) {
    return(output[[1]])
  }
  else {
    names(output) <- method
    return(output)
  }
}

getTunableArgs <- function(model) {
  
  if (inherits(model, "SDMmodelCV")) {
    method <- class(model@models[[1]]@model)
  } else {
    method <- class(model@model)
  }
  
  if (method == "Maxent") {
    args <- c("fc", "reg", "iter")
  } else if (method == "Maxnet") {
    args <- c("fc", "reg")
  } else if (method == "ANN") {
    args <- c("size", "decay", "rang", "maxit")
  } else if (method == "RF") {
    args <- c("mtry", "ntree", "nodesize")
  } else if (method == "MBOOST") {
    args <- c("distribution", "df", "mstop", "shrinkage")
  } else if (method == "XGBOOST") {
    args <- c("distribution", "gamma", "shrinkage", "depth", "nrounds")
  } else {
    args <- c("distribution", "n.trees", "interaction.depth", "shrinkage", "bag.fraction")
  }
  
  return(args)
}

confMatrix <- function(model, test = NULL, th = NULL, type = NULL) {
  
  if (!inherits(model, "SDMmodel")) 
    stop("Function available only for SDMmodel objects.")
  if (is.null(test)) {
    data <- model@data
  }
  else {
    if (!inherits(test, "SWD")) 
      stop("\"test\" argument invalid, use an SWD object.")
    data <- test
  }
  n_p <- sum(data@pa == 1)
  n_a <- sum(data@pa == 0)
  
  pred <- predict(model, data, type = type)
  
  p_pred <- pred[1:n_p]
  a_pred <- pred[(n_p + 1):(n_p + n_a)]
  if (is.null(th)) {
    th <- sort(unique(pred))
    th <- c(0, th, 1)
  }
  tp <- fp <- vector(mode = "numeric", length = length(th))
  for (i in seq_along(th)) {
    tp[i] <- sum(p_pred >= th[i])
    fp[i] <- sum(a_pred >= th[i])
  }
  fn <- n_p - tp
  tn <- n_a - fp
  conf_matrix <- data.frame(th = th, tp = tp, fp = fp, fn = fn, 
                            tn = tn)
  return(conf_matrix)
}

environment(confMatrix) <- asNamespace('SDMtune')
assignInNamespace("confMatrix", confMatrix, ns = "SDMtune")

environment(trainBRT) <- asNamespace('SDMtune')
assignInNamespace("trainBRT", trainBRT, ns = "SDMtune")

environment(.create_sdmtune_output) <- asNamespace('SDMtune')
assignInNamespace(".create_sdmtune_output", .create_sdmtune_output, ns = "SDMtune")

environment(.get_hypers_grid) <- asNamespace('SDMtune')
assignInNamespace(".get_hypers_grid", .get_hypers_grid, ns = "SDMtune")

environment(.args_name ) <- asNamespace('SDMtune')
assignInNamespace(".args_name", .args_name, ns = "SDMtune")

environment(.get_train_args) <- asNamespace('SDMtune')
assignInNamespace(".get_train_args", .get_train_args, ns = "SDMtune")

environment(train) <- asNamespace('SDMtune')
assignInNamespace("train", train, ns = "SDMtune")

environment(getTunableArgs) <- asNamespace('SDMtune')
assignInNamespace("getTunableArgs", getTunableArgs, ns = "SDMtune")

environment(trainMBOOST) <- asNamespace('SDMtune')
environment(MBOOST) <- asNamespace('SDMtune')

environment(trainXGBOOST) <- asNamespace('SDMtune')
environment(XGBOOST) <- asNamespace('SDMtune')

setClassUnion("model", c("Maxent", "Maxnet", "ANN", "RF", "BRT","MBOOST","XGBOOST"))




