# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

# devtools::install_github('adamlilith/enmSdm')

packages.to.use <- c(  
  "credentials",
  "biomod2" ,
  "shiny",
  "rnaturalearth",
  "spatstat",
  "ggplot2",
  "mboost",
  "sm",
  "blockCV",
  "spThin",
  "pROC",
  "modEvA",
  "boot",
  "ecodist"
  , "raster"
  , "verification"
  , "gdata"
  , "leaflet"
  , "dismo"
  , "gbm"
  , "sp"
  ,"leaflet.extras"
  , "adehabitatHS"
  , "SDMTools"
  , "parallel"
  , "doParallel"
  , "biganalytics"
  , "nicheROVER"
  , "vegan"
  , "parallel"
  , "rgeos"
  , "rgdal"
  , "sdmpredictors"
  , "usdm"
 # , "doSNOW"
  , "ENMeval"
  , "maptools"
  , "sf"
  , "monmlp"
  , "gstat")

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source", dependencies=TRUE) }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

defineRegion <- function(records,lonName,latName) {
  
  records <- records[which(!is.na(records[,lonName])),c(lonName,latName)] 
  
  ui <- fluidPage(leafletOutput("mymap",height=500))
  
  server <- function(input, output) {
    
    output$mymap <- renderLeaflet(
      leaflet() %>%
        addProviderTiles("Esri.OceanBasemap",group = "Ocean Basemap") %>%
        
        addCircles(lng=records[,lonName], lat=records[,latName] , weight = 3,color="Black",radius = 6) %>%
        
        addDrawToolbar(
          targetGroup='draw')  
      
    )
    
    observeEvent(input$mymap_draw_new_feature,{
      feature <<- input$mymap_draw_new_feature
      
      print(feature)
      
    })
    
  }
  
  return(shinyApp(ui = ui, server = server))
  
}

## -----------------------

selectRecords <- function(records,lonName,latName) {
  
  records <- records[which(!is.na(records[,lonName])),] 
  
  nPoints <- length(feature$geometry$coordinates[[1]])
  sapply(1:nPoints,function(x) { unlist(feature$geometry$coordinates[[1]][[x]]) })
  poly <- spPolygons(t(sapply(1:nPoints,function(x) { unlist(feature$geometry$coordinates[[1]][[x]]) })))
  
  spobj1 <- SpatialPointsDataFrame(records[,c(lonName,latName)], data=records)
  crs(spobj1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  crs(poly) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  records <- records[as.numeric(which( ! is.na(over(spobj1,poly) ))),]
  
  return(records)
}

## -----------------------

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
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

readRecords <- function(dataRecordsFile,separator,excludeLowRes,excludeNA,excludeDuplic,extraColumns,maskRecordsShapefile) {
  
  if( grepl("xls",dataRecordsFile) ) {  presences <- read.xls( dataRecordsFile , sheet=1) }
  if( grepl("csv",dataRecordsFile) ) {  presences <- read.csv( dataRecordsFile , sep=";", header = T,row.names = NULL ) }
  
  extraColumns <- c("Lon","Lat",extraColumns)
  columnsUse <- sapply(extraColumns,function(x) { which( colnames(presences) %in% x) })
  presences <- presences[,columnsUse]
  
  cat("\n")
  cat("Initial: ",nrow(presences), "occurrence records")
  cat("\n")
  
  if(excludeNA) {
    records.to.keep <- which( ! is.na( presences$Lon ) & ! is.na( presences$Lat ) )
    presences <- presences[records.to.keep,]
  }
  
  if(excludeDuplic) {
    presences <- presences[  which( ! duplicated(presences[,c(which(colnames(presences) == "Lon"),which(colnames(presences) == "Lat"))] , fromLast = TRUE ) ) , ]
  }
  
  if(excludeLowRes) {
    presences <- presences[which( sapply( presences$Lon , function(x) { decimalplaces(x) } ) >= 2),]
  }
  
  presences <- presences[which(presences$Lon <= 180 & presences$Lon >= -180 & presences$Lat < 90 & presences$Lat > -90 ),]
  
  cat("Final: ",nrow(presences), "occurrence records")
  cat("\n")
  cat("\n")
  
  if( !is.null(maskRecordsShapefile)) {
    
    presence.records.mask.shp <- shapefile(maskRecordsShapefile)
    crs(presence.records.mask.shp) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    presences.pts <- presences[,1:2]
    coordinates(presences.pts) <- c("Lon", "Lat")
    crs(presences.pts) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    points.not.over <- which(is.na( sp::over(presences.pts,presence.records.mask.shp) ) )
    presences <- presences[-points.not.over,]
  }
  
  return(presences)
  
}

## -----------------------

list.sdm.predictors <- function() {
  
  raster.predictors <- list_layers()
  raster.predictors <- raster.predictors[raster.predictors$dataset_code == "Bio-ORACLE",2:3]
  return(raster.predictors)
}

## -----------------------

processLayers <- function(rasterLayers,occurrenceRecords,regionBuffer,minDepth,maxDepth,intertidal) {
  
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
    
    rasters <- crop(rasters, extent(final.extent.rasters))
    
  } else { final.extent.rasters <- c(-180,180,-90,90) }

  ## -----------------------
  
  if( intertidal ) {
    
    intertidal <- raster(coastLineDataLayer)
    intertidal <- crop(intertidal, rasters)
    rastersIntertidal <- mask(rasters,intertidal)
    
  }
  
  ## -----------------------
  
  regions <- clump(subset(rasters,1))
  regionsUsed <- unique(raster::extract(regions,occurrenceRecords))
  regionsUsed <- regionsUsed[ !is.na(regionsUsed) ]
  
  regions[ ! regions %in% regionsUsed ] <- NA
  rasters <- mask(rasters,regions)
  
  ## -----------------------
  
  if( ! is.null(minDepth) | ! is.null(maxDepth) ) {
    
    if( is.null(minDepth) ) { minDepth <- 0 }
    if( is.null(maxDepth) ) { maxDepth <- 99998 }
    
    minDepth <- minDepth * (-1)
    maxDepth <- maxDepth * (-1)
    
    rclmat <- data.frame(from=c(-99999,maxDepth,minDepth),to=c(maxDepth,minDepth,0),value=c(NA,1,NA))
    rclmat <- rclmat[rclmat$from != rclmat$to,]
    
    bathymetry <- raster(bathymetryDataLayer)
    bathymetry <- crop(bathymetry, rasters)
    
    bathymetry <- reclassify(bathymetry, as.matrix(rclmat))
    rastersBathymetry <- mask(rasters,bathymetry)
    
  }
  
  if( intertidal & exists("rastersBathymetry") ) {
    
    rasters.i <- stack( sapply(1:length(names(rasters)),function(x) {  calc(stack(subset(rastersBathymetry,x),subset(rastersIntertidal,x)) , mean,na.rm=TRUE ) }) )
    names(rasters.i) <- names(rastersBathymetry)
    rm(rasters)
    rasters <- rasters.i
    
  }
  
  if( intertidal & ! exists("rastersBathymetry") ) { rasters <- rastersIntertidal }
  
  if( ! intertidal & exists("rastersBathymetry") ) { rasters <- rastersBathymetry }
  
  return(rasters)
  
}

## -----------------------

dropNoVariationLayers <- function(rasterLayers) {
  
  cave <- function(x) { length(unique(x)) / length(x) }
  randomLocations <- Which(!is.na(subset(rasterLayers,1)),cells=TRUE)
  randomLocations <- xyFromCell( subset(rasterLayers,1) , sample(randomLocations,min(length(randomLocations),1000),replace=FALSE))
  varRasterLayers <- which( apply( raster::extract(rasterLayers,randomLocations) ,2,cave) > 0.025)
  rasterLayers <- subset(rasterLayers,varRasterLayers)
  return(rasterLayers)
  
}

## -----------------------

unCorrelatedPairs <- function(rasterLayers,speciesData,threhold,dataLayersMonotonocity) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasterLayers))-1)) {
    
    for( j in (i+1):length(names(rasterLayers))) {
      
      corr.val <- abs(cor( raster::extract(subset(rasterLayers,i),speciesData[,c("Lon","Lat")]) , raster::extract(subset(rasterLayers,j),speciesData[,c("Lon","Lat")]) ,use = "pairwise.complete.obs"))
      
      if( corr.val < threhold | (sum(c(1,-1) %in% dataLayersMonotonocity.i[c(i,j)]) == 2)  ) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasterLayers,i)),
                                                                                                                                          Var.2=names(subset(rasterLayers,j)),
                                                                                                                                          Cor=corr.val,stringsAsFactors = FALSE)) }
      
    }
    
  }
  
  return(list.of.cor)
  
  
}

## -----------------------

correlatedPairs <- function(rasterLayers,speciesData,threhold,dataLayersMonotonocity.i) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasterLayers))-1)) {
    
    for( j in (i+1):length(names(rasterLayers))) {
      
      corr.val <- abs(cor( raster::extract(subset(rasterLayers,i),speciesData[,c("Lon","Lat")]) , raster::extract(subset(rasterLayers,j),speciesData[,c("Lon","Lat")]) ,use = "pairwise.complete.obs"))
      
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
    coordinates.to.relocate <- occurrenceRecords[to.relocate,]
    
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
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
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

correctLayer <- function(rasters,name,type,value.is,value.will) {		
  
  name.i <- which(names(rasters) == name)
  temp.r <- rasters[[name.i]]		
  
  if(type =="he") { temp.r[temp.r >= value.is] <- value.will }		
  if(type =="le") { temp.r[temp.r <= value.is] <- value.will }		
  if(type =="h") { temp.r[temp.r > value.is] <- value.will }		
  if(type =="l") { temp.r[temp.r < value.is] <- value.will }		
  
  names(temp.r) <- name
  
  rasters <- stack(subset(rasters,(1:(name.i-1))) , temp.r , subset(rasters,((name.i+1):dim(rasters)[3])) )
  
  return(rasters)		
  
}		

## -----------------------		

correlatedPairsOld <- function(rasters,threhold) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasters))-1)) {
    
    for( j in (i+1):length(names(rasters))) {
      
      corr.val <- abs(cor(getValues(subset(rasters,i)),getValues(subset(rasters,j)),use = "pairwise.complete.obs"))
      
      if( corr.val >= threhold) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasters,i)),
                                                                               Var.2=names(subset(rasters,j)),
                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      
    }
    
  }
  
  return(list.of.cor)
  
}

## -----------------------		

spatialAutocorrelation <- function(records,rasterLayers,autocorrelationClassDistance,autocorrelationSubsetRecords,autocorrelationMaxDistance,autocorrelationSignif) {
  
  # --------
  
  if( length(names(rasterLayers)) > 1  ) { 
    
    corrPairs <- correlatedPairs(rasterLayers,speciesData=records,threhold=0.5,dataLayersMonotonocity=NULL)
    
    if( nrow(corrPairs) > 0) { rasters <- subset(rasterLayers,(1:length(names(rasterLayers)))[-which(names(rasterLayers) %in% corrPairs[,1])]) }
    if( nrow(corrPairs) == 0) { rasters <- rasterLayers }
    
  } else { rasters <- rasterLayers }
  
  # --------
  
  set.seed(1)
  presences.t <- records[sample(1:nrow(records), ifelse(autocorrelationSubsetRecords> nrow(records),nrow(records),autocorrelationSubsetRecords) ),c("Lon","Lat")]
  
  presences.environment <- data.frame(raster::extract(rasters,presences.t))
  to.remove <- unique(which(is.na(presences.environment),arr.ind = T)[,1])
  
  if(length(to.remove) > 0) {
    
    presences.t <- presences.t[-to.remove,]
    presences.environment <- raster::extract(rasterLayers,presences.t)
    
  }
  
  presences.environment <- presences.environment[,which(apply(presences.environment,2,var) != 0)]
  
  space <- spDists(as.matrix(presences.t),as.matrix(presences.t),longlat=TRUE)
  data <- ecodist::distance( presences.environment , method = "mahalanobis")
  
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)
  
  resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,by=autocorrelationClassDistance),
                              classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,by=autocorrelationClassDistance),
                              R=NA,
                              pVal=NA)
  
  for( i in 1:nrow(resultsMatrix)) {
    
    d1 = resultsMatrix[i,1]
    d2 = resultsMatrix[i,2]
    
    data.d <- c(as.matrix(data))
    space.d <- c(space)
    
    data.d[space.d < d1 | space.d > d2] <- NA
    space.d[space.d < d1 | space.d > d2] <- NA
    
    #data.d[space.d > d2] <- NA
    #space.d[space.d > d2] <- NA
    
    data.d <- data.d[!is.na(data.d)]
    space.d <- space.d[!is.na(space.d)]
    
    if(length(space.d) == 0) { next }
    if(length(unique(space.d)) == 1) { next }
    
    val <- cor.test(data.d, space.d, method=c("pearson"))
    val.corr <- val$estimate
    val.p.value <- val$p.value
    
    resultsMatrix[i,3] <- val.corr
    resultsMatrix[i,4] <- val.p.value
    
    # modelobject <- lm(space.d~data.d)
    # f <- summary(modelobject)$fstatistic
    # p <-0
    # tryCatch( p <- pf(f[1],f[2],f[3],lower.tail=FALSE) , error=function(e) { Error <<- TRUE })
    # resultsMatrix[i,3] <- summary(modelobject)$adj.r.squared
    # resultsMatrix[i,4] <- p
    
  }

  distance <- round( resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ][1] )
  if( is.na(distance)) { distance <- autocorrelationMaxDistance }
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("First non-correlated distance: ",distance," km"))
  
  figure <- ggplot() +
    geom_hline(yintercept=0, size=0.1) +
    geom_line(data = resultsMatrix, aes(x = classdistanceTo, y = R), size = 0.1, color="Black", linetype = "dashed") +
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal <= autocorrelationSignif,], aes(x = classdistanceTo, y = R), size = 1.5, color="#6E6E6E") + 
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal > autocorrelationSignif,], aes(x = classdistanceTo, y = R), size = 1.5, color="Black") +
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal > autocorrelationSignif,][1,], aes(x = classdistanceTo, y = R), size = 1.75, color="Red") +
    ylab("Correlation (R)") + xlab("Geographic distance (Km)") + mainTheme

  return( list(figure=figure, distance = distance) )
          
}

## -----------------------

spatialThinning <- function(occurrenceRecords,rasterLayers,minDistance,verbose) {
  
  if( missing(verbose) ) { verbose <- TRUE }
  
  shape <- subset(rasterLayers,1)
  
  occurrenceRecords.i <- rasterize(occurrenceRecords,shape)
  occurrenceRecords.i <- xyFromCell(occurrenceRecords.i, Which(!is.na(occurrenceRecords.i),cells=TRUE) )
  colnames( occurrenceRecords.i ) <- c( "Lon", "Lat" )
  
  dataThinning <- data.frame(Name="Sp",occurrenceRecords.i)
  
  coordinates.t <- thin(dataThinning,
                        lat.col = "Lat",
                        long.col = "Lon",
                        spec.col = "Name",
                        thin.par = minDistance,
                        reps = 1,
                        write.files = FALSE,
                        locs.thinned.list.return = TRUE,
                        verbose = FALSE)[[1]]
  
  if(verbose) {
    
    cat( paste0("\n"))
    cat( paste0("\n"))    
    cat( paste0("Input Records: ",nrow(occurrenceRecords)))
    cat( paste0("\n"))
    cat( paste0("Final Records: ",nrow(coordinates.t)))
    
  }
  
  # Remove from main dataset of occurrences
  colnames( coordinates.t ) <- c( "Lon", "Lat" )
  
  # Remove log file
  if(length(list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt"))){
    file.remove( list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt") ) 
  }
  
  return(coordinates.t)
  
}

## -----------------------

pseudoAbsences <- function(occurrenceRecords,rasterLayers,patchContinuity,polygonPA,polygonPAType,paMindistance,paType,paBiasKernelSurface,paBiasKernelSurfaceProb,paRatio,paEnvironmentStratification,plotFigure) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  
  if(patchContinuity) {
    
    shapeRegions <- raster::clump(shape, directions=8 )
    
    regionsUsed <- raster::extract(shapeRegions,occurrenceRecords)
    regionsUsed <- unique(regionsUsed[!is.na(regionsUsed)])
    shape <- shapeRegions %in% regionsUsed
    shape[shape == 0] <- NA
    shape[ !is.na(shape)] <- 1
    
    nonNACells <- Which(!is.na(shape), cells=TRUE) 
    sink.points <- xyFromCell(shape, nonNACells)
    
    sink.points <- as.data.frame(sink.points)
    coordinates( sink.points ) <- c( "x", "y" )
    proj4string( sink.points ) <- crs(regionsFinal)
    
    regionsFinal <- regionsFinal[unique(over(sink.points,regionsFinal)),]
    
  }
  
  ## ------------
  
  nonNACells <- Which(!is.na(shape), cells=TRUE) 
  sink.points <- xyFromCell(shape, nonNACells)
  
  ## ------------
  
  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) }
  if( paRatio > 1 ) { final.paRatio <- paRatio }
  
  final.paRatio <- max(final.paRatio,paMinimum)
  
  ## ------------
  
  if( paBiasKernelSurface ) {
    
    # Generates bias surface kernel
    
    library(sm)
    
    bias <- cellFromXY(shape, as.matrix(occurrenceRecords))
    cells <- unique(sort(bias))
    kernelXY <- xyFromCell(shape, cells)
    samps <- as.numeric(table(bias))
    KDEsur <- sm.density(kernelXY, weights=samps, display="none", hmult=2, ylim=c(extent(shape)[3]-2.5,extent(shape)[4]), xlim=c(extent(shape)[1]-2.5,extent(shape)[2]-2.5), nbins=NA)
    KDErast <- SpatialPoints(expand.grid(x=KDEsur$eval.points[,1], y=KDEsur$eval.points[,2]))
    KDErast <- SpatialPixelsDataFrame(KDErast, data.frame(kde = array(KDEsur$estimate, length(KDEsur$estimate))))
    KDErast <- raster(KDErast)
    
    KDErast <- raster::resample(KDErast, shape, method="bilinear")
    
    Bias.surface <- mask(KDErast,shape)
    projection(Bias.surface) <- CRS("+proj=longlat +datum=WGS84")
    Bias.surface <- Bias.surface - cellStats(Bias.surface,min)
    Bias.surface <- Bias.surface / cellStats(Bias.surface,max)
    
    Bias.surface[Bias.surface < paBiasKernelSurfaceProb ] <- NA
    Bias.surface[Bias.surface >= paBiasKernelSurfaceProb ] <- 1
    
    if(patchContinuity) {
      
      shapeRegions <- raster::clump(Bias.surface, directions=8 )
      regionsUsed <- raster::extract(shapeRegions,occurrenceRecords)
      regionsUsed <- unique(regionsUsed[!is.na(regionsUsed)])
      shape <- shapeRegions %in% regionsUsed
      shape[shape == 0] <- NA
      shape[ !is.na(shape)] <- 1
      Bias.surface <- mask(Bias.surface,shape)
      
    }
    
    sink.points.p <- raster::extract(Bias.surface,sink.points)
    sink.points.p <- sink.points[which( !is.na(sink.points.p) ),]
    
    if( nrow(sink.points.p) != 0 ) { sink.points <- sink.points.p }
    
  }
  
  ## ------------
  
  if( ! is.null(paMindistance) ) {
    
    # Removes those closer than paDist
    
    sink.points.poly <- as.data.frame(occurrenceRecords)
    coordinates( sink.points.poly ) <- c( "Lon", "Lat" )
    proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    sink.points.poly <- gBuffer( sink.points.poly, width=paMindistance / 111.699, byid=TRUE )
    # plot(sink.points.poly)
    
    sink.points.pts <- as.data.frame(sink.points)
    colnames( sink.points.pts ) <- c( "Lon", "Lat" )
    coordinates( sink.points.pts ) <- c( "Lon", "Lat" )
    proj4string( sink.points.pts ) <- CRS( "+proj=longlat +datum=WGS84" )
    
    to.remove.id <- sp::over(sink.points.pts,sink.points.poly)
    to.keep <- which(is.na(to.remove.id))
    sink.points <- sink.points[to.keep,]
    
    if( nrow(sink.points) == 0 ) { 
      
      shape <- subset(rasterLayers,1)
      shape[!is.na(shape)] <- 1
      nonNACells <- Which(!is.na(shape), cells=TRUE)
      sink.points <- xyFromCell(shape, nonNACells) 
      
      sink.points.poly <- as.data.frame(occurrenceRecords)
      coordinates( sink.points.poly ) <- c( "Lon", "Lat" )
      proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
      
      sink.points.poly <- gBuffer( sink.points.poly, width=paMindistance / 111.699, byid=TRUE )
      # plot(sink.points.poly)
      
      sink.points.pts <- as.data.frame(sink.points)
      colnames( sink.points.pts ) <- c( "Lon", "Lat" )
      coordinates( sink.points.pts ) <- c( "Lon", "Lat" )
      proj4string( sink.points.pts ) <- CRS( "+proj=longlat +datum=WGS84" )
      
      to.remove.id <- sp::over(sink.points.pts,sink.points.poly)
      to.keep <- which(is.na(to.remove.id))
      sink.points <- sink.points[to.keep,]
      
    }
    
  }

  ## ------------
  
  if( paEnvironmentStratification ) {
    
    # Generates environmental surface for stratification
    
    to.generate <- min( c( final.paRatio , nrow(sink.points) ) )
    
    if( final.paRatio < nrow(sink.points) - (nrow(sink.points) * 0.1) ) {
      
      pa.environment <- raster::extract( rasterLayers,sink.points[,1:2] )
      pa.environment <- scale(pa.environment)
      pa.environment <- pa.environment[,which(!is.na(sapply(1:ncol(pa.environment) , function(x) { var(pa.environment[,x]) })))]
      pa.environment <- kmeans(pa.environment,centers=to.generate,iter.max = 10, nstart = 1)$cluster
      
      pa.clusters <- unique( pa.environment )
      finalKMClusters <- sapply(1:to.generate,function(x){ ifelse( length(which(pa.environment == x)) == 1 , which(pa.environment == x) , sample(which(pa.environment == x),1,replace=F) ) } )
      sink.points <- sink.points[finalKMClusters,]
      
    }
    
    if( final.paRatio > nrow(sink.points) - (nrow(sink.points) * 0.1) ) {
      
      sink.points <- sink.points[sample(1:nrow(sink.points),to.generate,replace=F),]
      
    }
    
  }
  
  if( ! paEnvironmentStratification ) {
    
    to.generate <- min( c( final.paRatio , nrow(sink.points) ) )
    sink.points <- sink.points[sample(1:nrow(sink.points),to.generate,replace=F),]
    
  }
  
  ## ------------
  
  absences <- sink.points
  
  if(plotFigure) {
    
    plot(absences[,1:2],col="grey", pch=20 , cex=1 , main="Dataset : Presences and Pseudo-absences")
    points(occurrenceRecords,col="black", pch=20 , cex=1)
    
  }
  
  absences <- absences[,1:2]
  colnames(absences) <- c("Lon","Lat")
  return(absences)
  
  options(warn=0)
  
}

## -----------------------

pseudoAbsencesSRE <- function(occurrenceRecords,rasterLayers,paMindistance,paRatio,plotFigure) {
  
  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  shape[is.na(shape)] <- 1
  
  if( paRatio <= 1 ) { final.paRatio <- round(nrow(occurrenceRecords) * (1/paRatio)) }
  if( paRatio > 1 ) { final.paRatio <- paRatio }
  
  pa.number <- (1/paRatio) * nrow(occurrenceRecords)
  
  pa.object <- biomod2::BIOMOD_FormatingData(  resp.var = rep(1,nrow(occurrenceRecords)),
                                               expl.var = stack(rasters),
                                               resp.xy = occurrenceRecords[,c("lon","lat")],
                                               resp.name = paste0("Species"),
                                               PA.nb.rep = 1,
                                               PA.nb.absences = pa.number,
                                               PA.strategy = "sre",
                                               PA.dist.min = NULL,
                                               PA.sre.quant = 0.025 )  # plot(pa.object)
  
  sink.points <- pa.object@coord[ -(1:nrow(occurrenceRecords)) ,]
  return(sink.points)
}

## -----------------------

dataPartitioning <- function(speciesData,rasterLayers,type,cvRemoveEdges,maxCorrDistance,k) {
  
  maxCorrDistance.i <- min(maxCorrDistance,250)
  
  speciesData.i <- data.frame(Species="Sp",speciesData)
  coordinates(speciesData.i) <- ~Lon+Lat
  crs(speciesData.i) <- crs(rasterLayers)
  
  if( type=="randomBlocks" ) {
    
    shape <- subset(rasterLayers,1)
    shape <- as(extent(shape), "SpatialPolygons")
    crs(shape) <- crs(rasterLayers)
    shape <- st_as_sf(shape)

    grd_lrg <- st_make_grid(shape, cellsize = c(round(max(maxCorrDistance.i*1000,111325) / 111325), round(max(maxCorrDistance.i*1000,111325) / 111325)))
    grd_lrg <- as_Spatial(grd_lrg)
    
    grd_lrgOver <- sp::over( grd_lrg , speciesData.i)
    grd_lrgOver <- which(!is.na(grd_lrgOver[,1]))
    grd_lrg <- grd_lrg[grd_lrgOver,]
    grd_lrg <- SpatialPolygonsDataFrame(grd_lrg, data.frame(ID= sample(1:k,length(grd_lrg),replace=TRUE) ), match.ID = F)
    
    crossValidation <- list()
    
    for(i in 1:k){
      
      grd_lrg.Test <- grd_lrg[grd_lrg$ID == i,]
      grd_lrg.Train <- grd_lrg[grd_lrg$ID != i,]
      
      grd_lrg.Test <- which(!is.na(sp::over( speciesData.i , grd_lrg.Test)[,1]))
      grd_lrg.Train <- which(!is.na(sp::over( speciesData.i , grd_lrg.Train)[,1]))
      
      crossValidation <- c(crossValidation,list(list(grd_lrg.Train, grd_lrg.Test)))
      
    }
    
    crossValidationPlot <<- grd_lrg
    
  }
  
  if( type=="blocksLongitudinal" ) {
    
    stop("Fix here")
    if(cvRemoveEdges) { k <- k + 2}
    
    crossValidationPlot <<- grd_lrg$plots
    crossValidation <- crossValidation$folds
    
    if(cvRemoveEdges) { crossValidation <- crossValidation[c(-1,-length(crossValidation))]  }
    
  }
  
  if( type=="blocksLatitudinal" ) {
    
    stop("Fix here")
    if(cvRemoveEdges) { k <- k + 2}
    
    crossValidationPlot <<- crossValidation$plots
    crossValidation <- crossValidation$folds
    
    if(cvRemoveEdges) { crossValidation <- crossValidation[c(-1,-length(crossValidation))]  }
    
  }
  
  return( crossValidation )
  
}

## -----------------------

drawPolygon2 <- function(r1) {
  
  plot(r1)
  poly <- spatstat::clickpoly(add=TRUE)

    p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    plot(sps,add=TRUE,col="#727272")
    
    crs(sps) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    sps <- as(sps, "SpatialPolygons")

  return( sps )
  
}

## -----------------------

drawPolygon <- function(rasterLayers,occurrenceRecords,name) {
  
  x <- 0
  dataRecordsDirectory.i <- as.numeric(unlist(gregexpr("/",dataRecordsFile)))
  dataRecordsDirectory <- paste0(substr(dataRecordsFile,1,dataRecordsDirectory.i[length(dataRecordsDirectory.i)-1]),"Spatial/")
  
  if( ! dir.exists(dataRecordsDirectory) ) { dir.create(dataRecordsDirectory,recursive = TRUE) }
  
  if( paste0("Polygon.",name,".shp") %in% list.files(dataRecordsDirectory) ) {   
    
    cat("\n")
    cat("Polygon for region of interest alreay exists.")
    cat("\n")
    cat("\n")
    cat("1. Use available polygon")
    cat("\n")
    cat("2. Generate new polygon")
    cat("\n")
    
    while( x != 1 & x != 2) {
      
      x <- readline(":")  
      
    }
    
  }
  
  if( x == 2 | x == 0 ) {   
    
    r.1 <- subset(rasterLayers,1)
    r.1[!is.na(r.1)] <- 1
    
    plot(r.1,col="#000000")
    points(occurrenceRecords,col="red",pch=16)
    
    poly <- clickpoly(add=TRUE)
    
    p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    plot(sps,add=TRUE,col="#727272")
    
    crs(sps) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    sps <- as(sps, "SpatialPolygons")
    rgdal::writeOGR(as(sps, "SpatialPolygonsDataFrame"), dataRecordsDirectory ,paste0("Polygon.",name), driver="ESRI Shapefile",overwrite_layer=TRUE)
    
  }
  
  if( x == 1 ) {   
    
    r.1 <- subset(rasterLayers,1)
    r.1[!is.na(r.1)] <- 1
    
    plot(r.1,col="#000000")
    
    sps = shapefile(paste0(dataRecordsDirectory,"/","Polygon.",name,".shp"))
    plot(sps,add=TRUE,col="#727272")
    points(occurrenceRecords,col="red",pch=20)
    
    sps <- as(sps, "SpatialPolygons")
    
  }
  
  return( sps )
  
}

## -----------------------

nicheModel <- function(modelType,speciesData,rasterLayers,crossValidation,cvIndex,dataLayersMonotonocity,printResults) {
  
  ## -----------------------
  
  if( exists("model") ) { rm(model) }
  
  ## -----------------------
  
  predicted.accuracy <- NULL
  
  ## -----------------------
  
  sorensenIndex <- function(obs,pred) {
    
    tp <- sum( obs == 1 &  pred == 1 ) / sum( obs == 1 )
    fp <- sum( obs == 0 &  pred == 1 ) / sum( obs == 0 )
    fn <- sum( obs == 1 &  pred == 0 ) / sum( obs == 1 )
    
    sorensen.i <- (2 * tp) / ( fn + (2 * tp) + fp )
    return(sorensen.i)
    
  }
  
  ## -----------------------
  
  accuracy.area.c <- function(model,observed,predicted,predicted.map) {
    
    options(warn=-1)
    
    predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
    predicted.accuracy$auc <- SDMTools::auc(observed ,predicted)
    
    predicted.distribution.area <- NA
    model.deviance <- Dsquared(glm(predicted~observed,family=binomial()))
    aicc <- NA
    
    if( class(model)[1] == "MaxEnt") { n.param <- get.params(model) ; aicc <- calc.aicc(n.param, presences.test, predicted.map)$AICc }
    
    if( class(model)[1] == "gbm") { model.deviance <- 1 - model$cv.statistics$deviance.mean }
    
    if(cvIndex == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
    
    if(cvIndex == "deviance") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$deviance),] }
    
    if(cvIndex == "aicc") { predicted.accuracy <- predicted.accuracy[which.min(predicted.accuracy$aicc),] }
    
    if(cvIndex == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
    
    if(cvIndex == "area") { 
      
      predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity >= 0.95,] 
      predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),] 
      
      if( nrow(predicted.accuracy) == 0 ) { predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )[1,] }
      
      predicted.distribution <- crop(predicted.map,extent( min( occurrenceRecords[,"Lon"] ) , max( occurrenceRecords[,"Lon"] ) , min( occurrenceRecords[,"Lat"] ) , max( occurrenceRecords[,"Lat"] )))
      
      full.area <- predicted.distribution
      full.area[!is.na(full.area)] <- 1
      full.area <- raster::resample(r.area,full.area)
      full.area <- raster::mask(full.area,predicted.distribution)
      full.area <- sum(getValues(full.area),na.rm=T)
      
      predicted.distribution[ predicted.distribution > predicted.accuracy$threshold ] <- 1
      predicted.distribution[ predicted.distribution <= predicted.accuracy$threshold ] <- NA
      predicted.distribution.area <- raster::resample(r.area,predicted.distribution)
      predicted.distribution.area <- raster::mask(predicted.distribution.area,predicted.distribution)
      predicted.distribution.area <- 1 - (sum(getValues(predicted.distribution.area),na.rm=T) / full.area)
      
    }
    
    threshold <- predicted.accuracy$threshold 
    auc <- predicted.accuracy$AUC 
    specificity <- predicted.accuracy$specificity 
    sensitivity <- predicted.accuracy$sensitivity 
    tss <- predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 
    sorensen <- sorensenIndex( observed, ifelse(predicted >= threshold , 1 , 0) )
    
    if( length(threshold) == 0 ) { threshold <- NA }
    if( length(auc) == 0 ) { auc <- NA }
    if( length(specificity) == 0 ) { specificity <- NA }
    if( length(sensitivity) == 0 ) { sensitivity <- NA }
    if( length(tss) == 0 ) { tss <- NA }
    if( length(sorensen) == 0 ) { sorensen <- NA }
    
    predicted.accuracy <- data.frame( sorensen=sorensen,
                                      threshold = threshold ,
                                      auc = auc ,
                                      specificity = specificity ,
                                      sensitivity = sensitivity ,
                                      tss = tss ,
                                      area = predicted.distribution.area ,
                                      aicc = aicc,
                                      deviance = model.deviance )
    
    options(warn=0)
    
    return(predicted.accuracy)
    
  }
  
  ## ------------------------------------------------------------------------------
  ## Clip rasters for performance
  
  extent.t <- extent( c(min(speciesData[,"Lon"])-1,max(speciesData[,"Lon"])+1,min(speciesData[,"Lat"])-1,max(speciesData[,"Lat"])+1) )
  rasters <- crop(rasterLayers,extent.t)
  predictors <- names(rasters)
  
  ## -----------------------
  
  cv.k.span <- 1:cvKFolds
  
  ## -----------------------
  
  if( cvIndex == "area" ) { r.area <- raster::area(subset(rasters,1)) }
  
  ## -----------------------------------------------------------------------------------
  
  cross.validation.rounds <- list()
  
  ## ----------------
  
  if( modelType == "MaxEnt" ) {
    
    stop("Error 02: Check Section")
    
    comb = expand.grid(cv.k = cv.k.span, feature.class=sapply(1:length(maxent.feature.class.span),function(x) { paste0(t(maxent.feature.class.span[x]),collapse=",")} ),betamultiplier=maxent.betamultiplier.span )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: MaxEnt \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('Number of parameters:', length(maxent.feature.class.span), '+' , length(maxent.betamultiplier.span))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .packages = c("dismo", "rJava","SDMTools","ENMeval")) %dopar% {
      
      cv <- comb[c,1]
      feature <- as.character(unlist(strsplit(as.character(comb[c,2]),split=",")))
      beta <- comb[c,3]
      args <- c(paste0("betamultiplier=",beta),feature)
      
      # ----------------------
      
      if( TRUE %in% grepl("lat", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lat <= crossValidation[cv,1] | presences$Lat >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lat <= crossValidation[cv,1] | absences$Lat >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lat >= crossValidation[cv,1] & presences$Lat <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lat >= crossValidation[cv,1] & absences$Lat <= crossValidation[cv,2] , ]
        
      }
      
      if( TRUE %in% grepl("lon", colnames(crossValidation)) ) {
        
        presences.train <- presences[ presences$Lon <= crossValidation[cv,1] | presences$Lon >= crossValidation[cv,2] , ]
        absences.train <- absences[ absences$Lon <= crossValidation[cv,1] | absences$Lon >= crossValidation[cv,2] , ]
        presences.test <- presences[ presences$Lon >= crossValidation[cv,1] & presences$Lon <= crossValidation[cv,2] , ]
        absences.test <- absences[ absences$Lon >= crossValidation[cv,1] & absences$Lon <= crossValidation[cv,2] , ]
        
      }  
      
      presences.train <- presences.train[!is.na(presences.train[,1]),]
      absences.train <- absences.train[!is.na(absences.train[,1]),]
      
      presences.test <- presences.test[!is.na(presences.test[,1]),]
      absences.test <- absences.test[!is.na(absences.test[,1]),]
      
      if( nrow(presences.train) > 0 &nrow(absences.train) > 0 & nrow(presences.test) > 0 & nrow(absences.test) > 0 ) { 
        
        train.dataset <- data.frame( PA = c( rep(1,nrow(presences.train)) , rep(0,nrow(absences.train)) ) , raster::extract( rasters , rbind( presences.train, absences.train) ) )
        test.dataset <- data.frame( PA = c( rep(1,nrow(presences.test)) , rep(0,nrow(absences.test)) ) , raster::extract( rasters , rbind( presences.test, absences.test) ) )
        
        train.dataset <- train.dataset[complete.cases(train.dataset),]
        test.dataset <- test.dataset[complete.cases(test.dataset),]
        
        model <- dismo::maxent(x = rasters , p = presences.train, a=absences.train , args=c("-P","autorun=false",args) )
        
        observed <- test.dataset$PA
        predicted <- predict( model, test.dataset[,-1])
        
        if( cvIndex == "area" ) { 
          predicted.map <- predict( rasters , model )
        } else { predicted.map <- NULL}
        
        predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
        
        predicted.accuracy <- data.frame( cv.round=cv,
                                          feature=paste(feature,collapse="+"),
                                          beta=beta,
                                          predicted.accuracy)
        
        return(predicted.accuracy)
        
      }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- do.call(rbind,cv.accuracy)
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold),]
    
    # ------------------
    
    if(cvIndex == "aicc") { 
      best.model <- aggregate(list(aicc=cv.accuracy$aicc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.min(best.model[,"aicc"])
    }
    if(cvIndex == "tss") { 
      best.model <- aggregate(list(tss=cv.accuracy$tss), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"tss"])
    }
    if(cvIndex == "auc") { 
      best.model <- aggregate(list(auc=cv.accuracy$auc), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"auc"])
    }
    if(cvIndex == "area") { 
      best.model <- aggregate(list(area=cv.accuracy$area), by = list(feature=cv.accuracy$feature,beta=cv.accuracy$beta), mean)
      best.model.i <- which.max(best.model[,"area"])
    }
    
    # ------------------
    
    best.model <- best.model[best.model.i,]
    
    best.model.feature <- as.character(best.model$feature)
    best.model.beta <- as.numeric(best.model$beta)
    
    best.model.metrics <- cv.accuracy[cv.accuracy$feature == best.model.feature & cv.accuracy$beta == best.model.beta ,  ]
    best.cross.validation <- best.model.metrics[,cvIndex]
    
    model <- dismo::maxent(x = rasters , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
    
    if( simplifyModels ) {
      
      model.importance <- summary.model(model,print.data=FALSE)
      model.simplift.vars <- as.character(model.importance$Variable)[model.importance$Permutation <= 1]
      
      model <- dismo::maxent(x = dropLayer(rasters,which(names(rasters) %in% model.simplift.vars)) , p = presences, a=absences , args=c("-P","autorun=false",c(paste0("betamultiplier=",best.model.beta),unlist(strsplit(best.model.feature, "[+]")))) )
      model.predicted <- predict( rasters , model )
      
    }
    
  }
  
  ## ------------------------------------------------------------------------------
  
  if( modelType == "BRT" ) {
    
    dataLayersMonotonocity.i <- sapply(predictors,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )
    
    if( length(names(rasters)) == 1 ) {  
      
      dummy.raster <- rasters
      dummy.raster[1:length(rep(1:2,length(dummy.raster) / 2))] <- rep(1:2,length(dummy.raster) / 2)
      names(dummy.raster) <- "layer"
      rasters <- stack(rasters,dummy.raster)
      dataLayersMonotonocity.i <- data.frame(data.frame(t(dataLayersMonotonocity.i)),layer=+1)
      predictors <- c(predictors,"layer")
      
    }
    
    comb = expand.grid(cv.k = cv.k.span, learning.complex=brtLearning,tree.depth=brtTreeDepth , bag=brtBagFraction )
    
    if(printResults) {
      
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('\n')
      cat('Model: BRT \n')
      cat('Cross-validation rounds:', nrow(comb))
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      
    }
    
    cl <- parallel::makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .export=c("brtMaxTrees","Dsquared"), .packages = c("dismo","SDMTools","ENMeval")) %dopar% {
      
      cv <- comb[c,1]
      l.rate <- comb[c,2]
      tree.c <- comb[c,3]
      bag <- comb[c,4]
      
      train.dataset <- speciesData[crossValidation[[cv]][[1]],]
      test.dataset <- speciesData[crossValidation[[cv]][[2]],]
      
      if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2  ) { 
        
        if( sum(train.dataset$PA == 1) < 10 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  }
        
        train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasters , train.dataset[,c("Lon","Lat")] ) )
        train.dataset[train.dataset == "NaN"] <- NA
        test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasters , test.dataset[,c("Lon","Lat")] ) )
        test.dataset[test.dataset == "NaN"] <- NA
        
        
        tryCatch( model <- gbm.step( data=train.dataset, 
                                     gbm.x = which(colnames(train.dataset) %in% predictors),
                                     gbm.y = 1, 
                                     family = "bernoulli", 
                                     plot.main = FALSE,
                                     tree.complexity = tree.c, 
                                     learning.rate = l.rate, 
                                     bag.fraction = bag, 
                                     n.folds=10,
                                     step.size=50,
                                     max.trees=brtMaxTrees,
                                     silent=TRUE,
                                     var.monotone = dataLayersMonotonocity.i ,
                                     verbose=FALSE) , error=function(e) { model <<- NULL })
        
        if( ! is.null(model) ) {
          
          num.tress <- model$gbm.call$best.trees
          observed <- test.dataset$PA
          predicted <- predict( model , test.dataset[,-1] , n.trees=num.tress,type="response")
          
          if( cvIndex == "area" ) { 
            predicted.map <- predict( rasters , model , n.trees=num.tress,type="response")
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( is.null(model)) {
          
          predicted.accuracy <- data.frame( sorensen=NA, threshold=NA, auc=NA, specificity=NA, sensitivity=NA, tss=NA, area=NA, aicc=NA, deviance=NA )
          
        }
        
        predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,bag=bag,predicted.accuracy)
        return(predicted.accuracy)
        
      }
      
      if( length(unique(train.dataset$PA)) != 2 | length(unique(test.dataset$PA)) != 2 ) { return(NULL) }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    
    best.model <- cv.accuracy
    best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(bag=best.model$bag,tree.c=best.model$tree.c,l.rate=best.model$l.rate), mean)
    best.model.i <- best.model[sort(best.model$cvIndex, decreasing = TRUE, index.return = TRUE)$ix,]
    
    if( nrow(cv.accuracy) > 0 ) { 
      
      model <- NULL
      m.i <- 0
      
      while( is.null(model) ){
        
        m.i <- m.i + 1
        
        if( m.i > nrow(best.model.i) ) { stop("Error code :: INT001")}
        
        best.model.tc <- best.model.i[m.i,"tree.c"]
        best.model.lr <- best.model.i[m.i,"l.rate"]
        best.model.bag <- best.model.i[m.i,"bag"]
        if( m.i == nrow(best.model.i) ) { best.model.lr <- 0.00001 ; best.model.tc <- 2 }
        
        best.model.metrics <- cv.accuracy[cv.accuracy$tree.c == best.model.tc & cv.accuracy$l.rate == best.model.lr & cv.accuracy$bag == best.model.bag ,  ]
        best.model.metrics <- best.model.metrics[best.model.metrics$threshold != 0 ,]
        best.cross.validation <- best.model.metrics[,cvIndex]
        
        train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasters , speciesData[,c("Lon","Lat")] ) )
        train.dataset[train.dataset == "NaN"] <- NA
        if( sum(train.dataset$PA == 1) < 10 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  }
        
        tryCatch( model <- gbm.step( data=train.dataset, 
                                     gbm.x = which(colnames(train.dataset) %in% predictors),
                                     gbm.y = 1, 
                                     family = "bernoulli", 
                                     plot.main = FALSE,
                                     tree.complexity = best.model.tc, 
                                     learning.rate = best.model.lr, 
                                     bag.fraction = best.model.bag, 
                                     n.folds=10,
                                     step.size=50,
                                     max.trees=brtMaxTrees,
                                     silent=TRUE,
                                     var.monotone = dataLayersMonotonocity.i ,
                                     verbose=FALSE) , error=function(e) { Error <- TRUE })
        
      }
    }
    
    if( nrow(cv.accuracy) == 0 ) {
      
      model <- NULL 
      best.cross.validation <- rep(0,10)
      best.model.metrics <- data.frame(threshold = 1 ,
                                       cv <- rep(0,10) ,
                                       auc= 0 ,
                                       specificity = 0 ,
                                       sensitivity = 0 ,
                                       tss = 0 ,
                                       aicc = 0 ,
                                       area = 1 ,
                                       deviance = 0 )
      
    }
    
  }
  
  ## ------------------------------------------------------------------------------
  
  if( modelType == "MBOOST" ) {
    
    dataLayersMonotonocity.i <- sapply(predictors,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )
    
    comb = expand.grid(cv.k = cv.k.span, shrinkage=mboostShrinkage,df=mboostDF,mstop=mboostIterations )
    
    if(printResults) {
      
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
      cat('\n')
      cat('\n')
      cat('Model: MBOOST \n')
      cat('Cross-validation rounds:', nrow(comb))
      cat('\n')
      cat('\n')
      cat('----------------------------------------------------------------------')
    }
    
    cl <- parallel::makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine = rbind, .export=c("brtMaxTrees","Dsquared"), .packages = c("dismo","SDMTools","ENMeval","mboost")) %dopar% {
      
      cv <- comb[c,1]
      shrinkage <- comb[c,2]
      df <- comb[c,3]
      mstop <- comb[c,4]
      
      train.dataset <- speciesData[crossValidation[[cv]][[1]],]
      test.dataset <- speciesData[crossValidation[[cv]][[2]],]
      
      if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2  ) { 
        
        if( sum(train.dataset$PA == 1) < 10 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset, simplify = FALSE)))  }
        
        train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasters , train.dataset[,c("Lon","Lat")] ) )
        train.dataset[train.dataset == "NaN"] <- NA
        test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasters , test.dataset[,c("Lon","Lat")] ) )
        test.dataset[test.dataset == "NaN"] <- NA
        
        colnames(train.dataset) <- c("PA",names(rasters))
        colnames(test.dataset) <- c("PA",names(rasters))
        
        train.dataset[,1] <- as.factor(train.dataset[,1])
        
        predictors <- colnames(train.dataset)[-1]
        response <- colnames(train.dataset)[1]
        
        constr_mono <- as.character(dataLayersMonotonocity.i)
        constr_mono[constr_mono == "-1"] <- "decreasing"
        constr_mono[constr_mono == "1"] <- "increasing"
        
        rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",df,")", sep = "")),collapse = " + ")
        fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
        
        ctrl <- boost_control(mstop = mstop, trace = FALSE, nu=shrinkage,stopintern=TRUE)
        
        tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <<- TRUE })
        
        if( exists("model") ) { 
          
          observed <- test.dataset$PA
          predicted <- predict(model,test.dataset, type = "response")
          
          if( cvIndex == "area" ) { 
            predicted.map <- predict(rasters,model)
            predicted.map <- predicted.map + ( min(getValues(predicted.map),na.rm=T) * (-1))
            predicted.map <- predicted.map / max(getValues(predicted.map),na.rm=T)
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( ! exists("model") ) { 
          
          predicted.accuracy <- data.frame( sorensen=NA, threshold=NA, auc=NA, specificity=NA, sensitivity=NA, tss=NA, area=NA, aicc=NA, deviance=NA )
          
        }
        
        predicted.accuracy <- data.frame( cv.round=cv,
                                          shrinkage=shrinkage,
                                          df=df,
                                          n.iterations=mstop,
                                          predicted.accuracy)
        
        return(predicted.accuracy)
        
      }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy[cv.accuracy$threshold == 0 & !is.na(cv.accuracy$deviance) ,"threshold" ] <- 0.1
    
    best.model <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(shrinkage=best.model$shrinkage,df=best.model$df,n.iterations=best.model$n.iterations), mean)
    best.model.i <- best.model[sort(best.model$cvIndex, decreasing = TRUE, index.return = TRUE)$ix,]
    
    model <- NULL
    m.i <- 0
    
    while( is.null(model) ){
      
      m.i <- m.i + 1
      
      if(m.i > nrow(best.model.i) ) { break }

      best.model.shrinkage <- best.model.i[m.i,]$shrinkage
      best.model.df <- best.model.i[m.i,]$df
      best.model.iterations <- best.model.i[m.i,]$n.iterations
      
      best.model.metrics <- cv.accuracy[cv.accuracy$shrinkage == best.model.shrinkage & cv.accuracy$df == best.model.df & cv.accuracy$n.iterations == best.model.iterations ,  ]
      best.cross.validation <- best.model.metrics[,cvIndex]
      
      train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasters , speciesData[,c("Lon","Lat")] ) )
      train.dataset[train.dataset == "NaN"] <- NA
      
      if( sum(train.dataset$PA == 1) < 8 ) { train.dataset <- rbind(train.dataset,do.call("rbind", replicate( sum(train.dataset$PA == 1) * 5 , train.dataset[train.dataset$PA == 1 ,], simplify = FALSE)))  }
      
      
      train.dataset <- train.dataset[complete.cases(train.dataset),]
      train.dataset[,1] <- as.factor(train.dataset[,1])
      
      colnames(train.dataset) <- c("PA",names(rasters))
      
      predictors <- colnames(train.dataset)[-1]
      response <- colnames(train.dataset)[1]
      
      constr_mono <- as.character(dataLayersMonotonocity.i)
      constr_mono[constr_mono == "-1"] <- "decreasing"
      constr_mono[constr_mono == "1"] <- "increasing"
      
      rhs <- paste(c(paste("bmono(", predictors, ", constraint = \"", constr_mono,"\", df = ",best.model.df,")", sep = "")),collapse = " + ")
      fm_mono <- as.formula(paste(response, " ~ ", rhs, collapse = ""))
      
      ctrl <- boost_control(mstop = best.model.iterations, trace = FALSE, nu=best.model.shrinkage,stopintern=TRUE)
      
      tryCatch( model <- mboost(fm_mono, data = train.dataset, control = ctrl, family = Binomial(type = "adaboost" , link = "logit" )) , error=function(e) { Error <- TRUE })
      
    }
    
  }
  
  ## ------------------------------------------------------------------------------
  
  if( modelType == "MPNN" ) {
    
    dataLayersMonotonocity.i <- sapply(predictors,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )
    
    comb = expand.grid(cv.k = cv.k.span, hidden=mpnnHidden,itereractions=mpnnItereractions )
    
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    cat('\n')
    cat('\n')
    cat('Model: MPNN \n')
    cat('Cross-validation rounds:', nrow(comb))
    cat('\n')
    cat('\n')
    cat('----------------------------------------------------------------------')
    
    cl <- parallel::makeCluster(nCores)
    registerDoParallel(cl)
    
    cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("monmlp","dismo","SDMTools","ENMeval")) %dopar% {
      
      cv <- comb[c,1]
      hidden <- comb[c,2]
      itereractions <- comb[c,3]
      
      train.dataset <- speciesData[crossValidation[[cv]][[1]],]
      test.dataset <- speciesData[crossValidation[[cv]][[2]],]
      
      if( length(unique(train.dataset$PA)) == 2 & length(unique(test.dataset$PA)) == 2 ) { 
        
        train.dataset <- data.frame( PA = train.dataset$PA , raster::extract( rasters , train.dataset[,c("Lon","Lat")] ) )
        train.dataset[train.dataset == "NaN"] <- NA
        test.dataset <- data.frame( PA = test.dataset$PA , raster::extract( rasters , test.dataset[,c("Lon","Lat")] ) )
        test.dataset[test.dataset == "NaN"] <- NA
        
        correctLayers <- names(rasters)[which(dataLayersMonotonocity.i == -1)]
        
        for(correctLayers.i in correctLayers) {
          train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
          test.dataset[,correctLayers.i] <- test.dataset[,correctLayers.i] * (-1)
        }
        
        train.dataset <- train.dataset[complete.cases(train.dataset),]
        
        model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = hidden, bag = TRUE, iter.max = itereractions, iter.stopped = 10)
        
        if( ! is.null(model) ) {
          
          observed <- test.dataset$PA
          predicted <- as.numeric(monmlp.predict(as.matrix(test.dataset[,-1]), model))
          # odds <- exp(predicted)
          # predicted <- odds / (1 + odds)
          
          if( cvIndex == "area" ) { 
            predicted.map <- predictDistribution(rasters,model, reclassToOne=FALSE)
          } else { predicted.map <- NULL}
          
          predicted.accuracy <- accuracy.area.c(model,observed,predicted,predicted.map)
          
        }
        
        if( is.null(model)) {
          
          predicted.accuracy <- data.frame( sorensen=0,threshold=0,auc=0,specificity=0,sensitivity=0,tss=0,area=0,aicc=+Inf,deviance=0)
          
        }
        
        predicted.accuracy <- data.frame( cv.round=cv,hidden=hidden,itereractions=itereractions,predicted.accuracy)
        return(predicted.accuracy)
        
      } else { return(NULL) }
      
    }
    
    stopCluster(cl); rm(cl) ; gc(reset=TRUE)
    
    # ------------------
    
    cv.accuracy <- cv.accuracy[!is.na(cv.accuracy$threshold) & cv.accuracy$threshold > 0 ,]
    
    if( nrow(cv.accuracy) != 0 ) { 
      
      best.model <- cv.accuracy
      best.model <- aggregate(list(cvIndex=best.model[,cvIndex]), by = list(hidden=best.model$hidden,itereractions=best.model$itereractions), mean)
      best.model.i <- which.max(best.model[,"cvIndex"])
      best.model <- best.model[best.model.i,]
      
      best.model.hidden <- best.model$hidden
      best.model.itereractions <- best.model$itereractions
      
      best.model.metrics <- cv.accuracy[cv.accuracy$hidden == best.model.hidden & cv.accuracy$itereractions == best.model.itereractions  ,  ]
      best.cross.validation <- best.model.metrics[,cvIndex]
      
      train.dataset <- data.frame( PA = speciesData$PA , raster::extract( rasters , speciesData[,c("Lon","Lat")] ) )
      train.dataset[train.dataset == "NaN"] <- NA
      train.dataset <- train.dataset[complete.cases(train.dataset),]
      
      correctLayers <- names(rasters)[which(dataLayersMonotonocity.i == -1)]
      for(correctLayers.i in correctLayers) {
        train.dataset[,correctLayers.i] <- train.dataset[,correctLayers.i] * (-1)
      }
      
      model <- monmlp.fit(x = as.matrix(train.dataset[,-1]), y = as.matrix(train.dataset[,1]),  monotone = 1:(ncol(train.dataset)-1),  hidden1 = best.model.hidden, bag = TRUE, iter.max = best.model.itereractions, iter.stopped = 10)
      
    }
    
    if( nrow(cv.accuracy) == 0 ) {
      
      model <- NULL 
      best.cross.validation <- rep(0,10)
      best.model.metrics <- data.frame(threshold = 1 ,
                                       cv <- rep(0,10) ,
                                       auc= 0 ,
                                       specificity = 0 ,
                                       sensitivity = 0 ,
                                       tss = 0 ,
                                       aicc = 0 ,
                                       area = 1 ,
                                       deviance = 0 )
      
    }
    
  }
  
  ## -----------------------
  ## -----------------------
  
  if( printResults ) { 
    
    cat( paste0("\n"))
    cat( paste0("\n --------------------------------------------------------"))
    cat( paste0("\n"))
    cat( paste0("\n"))
    cat( paste0("Accuracy of model:"))
    cat( paste0("\n"))
    cat( paste0("AUC: ",round( mean(best.model.metrics[,"auc"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"auc"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("specificity: ",round( mean(best.model.metrics[,"specificity"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"specificity"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("sensitivity: ",round( mean(best.model.metrics[,"sensitivity"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"sensitivity"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("tss: ",round( mean(best.model.metrics[,"tss"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"tss"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("aicc: ",round( mean(best.model.metrics[,"aicc"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"aicc"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("deviance: ",round( mean(best.model.metrics[,"deviance"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"deviance"]) , digits = 4)) )
    cat( paste0("\n"))
    cat( paste0("area: ",round( mean(best.model.metrics[,"area"]) , digits = 4), " ± " , round( sd(best.model.metrics[,"area"]) , digits = 4)) )
    cat( paste0("\n"))
    
  }
  
  ## -----------------------------------------------
  
  return( list(model=model, 
               cv=best.cross.validation,
               auc=mean(best.model.metrics[,"auc"]),
               Specificity=mean(best.model.metrics[,"specificity"]),
               Sensitivity=mean(best.model.metrics[,"sensitivity"]),
               tss=mean(best.model.metrics[,"tss"]),
               area=mean(best.model.metrics[,"area"]),
               deviance=mean(best.model.metrics[,"deviance"]))
  )
  
}

## -----------------------

getContribution <- function(modelType,fullModel,rasterLayers,speciesData,varImpType) {
  
  if( varImpType == "deviance" ) {
    
    predictionModel <- predictDistribution(rasterLayers,fullModel$model, reclassToOne=FALSE)
    predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
    observed <- speciesData$PA
    observed <- observed[!is.na(predicted)]
    predicted <- predicted[!is.na(predicted)]
    fullModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
    
  }
  
  if( varImpType == "auc" ) {
    
    fullModelDE <- fullModel$auc
    
  }
  
  relativeImportance <- data.frame(Predictor=names(rasterLayers),relImportance=NA)
  
  for(v in 1:length(names(rasterLayers)) ) {
    
    rasterLayers.v <- subset(rasterLayers, (1:length(names(rasterLayers)))[-v] )
    dataLayersMonotonocity.v <- data.frame(t(sapply( names(rasterLayers.v) ,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )))
    
    model <- nicheModel(modelType=modelType,speciesData,rasterLayers.v,crossValidation,cvIndex,dataLayersMonotonocity.v,printResults=FALSE)
    
    if( varImpType == "deviance" ) {
      
      predictionModel <- predictDistribution(rasterLayers,model$model, reclassToOne=FALSE)
      predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
      observed <- speciesData$PA
      observed <- observed[!is.na(predicted)]
      predicted <- predicted[!is.na(predicted)]
      partialModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
      
    }
    
    if( varImpType == "auc" ) {
      
      partialModelDE <- model$auc
      
    }
    
    relativeImportance[v,2] <- fullModelDE - partialModelDE
    
  }
  
  relativeImportance[,2] <- (relativeImportance[,2] / sum( abs(relativeImportance[,2]) )) * 100
  
  relativeImportancePlot <- ggplot(relativeImportance[sort(relativeImportance[,2],decreasing = TRUE,index.return=TRUE)$ix,]) +
    geom_bar( aes(x= reorder(Predictor, relImportance) , y=relImportance), stat="identity", fill="black", alpha=0.5) +
    coord_flip() + theme(
      axis.text=element_text(size=10),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
      panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
      panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
    ) + labs(x = "Predictor") + 
    labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
    annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)
  
  return(list(dataFrame=relativeImportance,plot=relativeImportancePlot))
  
}

## -----------------------

simplifyModel <- function(modelType,fullModel,rasterLayers,speciesData,relativeImportance,allButOne,simplifyType) {
  
  if( simplifyType == "deviance" ) {
    
    predictionModel <- predictDistribution(rasterLayers,fullModel$model, reclassToOne=FALSE)
    predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
    observed <- speciesData$PA
    observed <- observed[!is.na(predicted)]
    predicted <- predicted[!is.na(predicted)]
    fullModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
    
  }
  
  if( simplifyType == "auc" ) {
    
    fullModelDE <- fullModel$auc
    
  }
  
  lossFunction <- data.frame(predictorLoss=rep(NA,length(names(rasterLayers))-1),loss=NA,loss0.01=NA,loss0.99=NA) 
  sortLossfunction <- as.numeric(sapply(as.character(relativeImportance[sort(relativeImportance[,2],decreasing = FALSE, index.return=TRUE)$ix,1]),function(x){ which(names(rasterLayers) == x )}))
  
  for(v in 1:(ifelse( length(names(rasterLayers)) >= 4 , length(names(rasterLayers)) - 1 , length(names(rasterLayers))-1 )) ) {
    
    rasterLayers.v <- subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] )
    dataLayersMonotonocity.v <- data.frame(t(sapply( names(rasterLayers.v) ,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )))
    
    model <- nicheModel(modelType=modelType,speciesData,rasterLayers.v,crossValidation,cvIndex,dataLayersMonotonocity.v,printResults=FALSE)
    
    if( simplifyType == "deviance" ) {
      
      lab <- "(Deviance Explained)"
      
      if( length(names(rasterLayers.v)) == 1 ) { 
        dummy.raster <- rasterLayers.v
        dummy.raster[1:length(rep(1:2,length(dummy.raster) / 2))] <- rep(1:2,length(dummy.raster) / 2)
        names(dummy.raster) <- "layer"
        rasterLayers.v <- stack(rasterLayers.v,dummy.raster) 
        }
      
      predictionModel <- predictDistribution(rasterLayers.v,model$model, reclassToOne=FALSE)
    
      predicted <- raster::extract(predictionModel,speciesData[,c("Lon","Lat")])
      observed <- speciesData$PA
      observed <- observed[!is.na(predicted)]
      predicted <- predicted[!is.na(predicted)]
      partialModelDE <- Dsquared(glm(predicted~observed,family=binomial()))
      partialModelDEConf <- 0
      partialModelDEConf$t <- c(0,0,0,0,0)
      
      # boot(dataset,function(data,indices) Dsquared(glm(predicted~observed,data=data[indices,],family=binomial())) ,R=9, parallel = c("snow"))
      
    }
    
    if( simplifyType == "auc" ) {
      
      lab <- "(Area Under the Curve)"
      
      partialModelDE <- mean(model$cv)
      partialModelDEConf <- list(t=c(model$cv,model$cv))
      
    }
    
    lossFunction[v,1] <- v
    lossFunction[v,2] <- fullModelDE - partialModelDE
    lossFunction[v,3] <- fullModelDE - partialModelDE - (sd(partialModelDEConf$t) / sqrt(length(partialModelDEConf$t) ))
    lossFunction[v,4] <- fullModelDE - partialModelDE + (sd(partialModelDEConf$t) / sqrt(length(partialModelDEConf$t) ))
    
  }
  
  keepPredictors <- which(round(lossFunction[,4],digits = 3) >= 0.01)

  if( length(keepPredictors) == 0) { keepPredictors <- 1 }
  if( ! 1 %in% keepPredictors) { removePredictors <- 1:(min(keepPredictors)-1) }
  if(   1 %in% keepPredictors) { removePredictors <- numeric(0) }
  
  colorsPredictors <- rep("Black",nrow(lossFunction))
  colorsPredictors[removePredictors] <- "Red"
  
  lossFunctionPlot <- ggplot(lossFunction, aes(x=predictorLoss, y=loss)) + 
    geom_line(linetype = "dashed",color="Gray") + #  ylim(-0.5,0.5) +
    geom_point(color=colorsPredictors, size=2.5) +
    geom_hline(aes(yintercept=0  ),color="Black", size=0.2) + 
    xlab("Variables removed (number)") + ylab(paste0("Loss function ",lab)) +
    geom_ribbon(aes(ymin = lossFunction[,3],ymax = lossFunction[,4]), alpha = 0.3)
  
  if( length(removePredictors) > 0 ) {
    if(allButOne) { if(length(removePredictors) > 1) { removePredictors <- removePredictors[-length(removePredictors)]} }
    v <- ifelse( length(keepPredictors) > 0 , min(keepPredictors) - 1 , max(removePredictors) ) 
    rasterLayers.v <- subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] )
    dataLayersMonotonocity.v <- data.frame(t( sapply( names(rasterLayers.v) ,function(x) { as.numeric(dataLayersMonotonocity[colnames(dataLayersMonotonocity) == x]) } )))
    modelRecuced <- nicheModel(modelType=modelType,speciesData,rasterLayers.v,crossValidation,cvIndex,dataLayersMonotonocity.v,printResults=FALSE)
    
    keptVariables <- names(subset(rasterLayers, (1:length(names(rasterLayers)))[-sortLossfunction[1:v]] ))
    removedVariables <- names(subset(rasterLayers, (1:length(names(rasterLayers)))[sortLossfunction[1:v]] ))
    
  }
  
  if( length(removePredictors) == 0 ) { 
    
    modelRecuced <- fullModel 
    
    keptVariables <- names(rasterLayers)
    removedVariables <- character(0)
    
    
  }
  
  
  return(list(model=modelRecuced$model,plot=lossFunctionPlot,keptVariables=keptVariables,removedVariables=removedVariables))
  
}

## -----------------------

summaryModel <- function(model,print.data) {
  
  if( class(model) == "MaxEnt" ) { 
    
    data.to.plot <- ENMeval::var.importance(model)
    colnames(data.to.plot) <- c("Variable","Percentage","Permutation")
    rownames(data.to.plot) <- NULL
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Percentage,index.return=T)$ix,]
    
    barplot(  t(data.to.plot[,2:3]), 
              main="Variable contribution", 
              horiz=TRUE,
              legend = TRUE,
              beside = TRUE, 
              names.arg=as.character(data.to.plot$Variable))
    
    if(print.data) { print(data.to.plot) }
    
    
  }
  
  if( class(model) == "gbm" ) { 
    
    data.to.plot <- summary(model)
    colnames(data.to.plot) <- c("Variable","Percentage")
    rownames(data.to.plot) <- NULL
    
    if(print.data) { print(data.to.plot) }
    
  }
  
  if( class(model) == "mboost" ) { 
    
    data.to.plot <- data.frame(mboost::varimp(model))
    data.to.plot <- data.frame(Variable=data.to.plot$variable,Percentage=(data.to.plot$reduction / sum(data.to.plot$reduction) ) * 100 )
    
    if(print.data) { print(data.to.plot) }
    
  }
  
  return(data.to.plot)
  
}

## -----------------------

modelPlot <- function(model,rasters,distribution,variable.to.plot,occurrenceRecords,print.limiting,auto.limiting,distribution.threshold,distribution.predicted,val.limiting,export.data.to.plot,plot.x.lim,plot.y.lim) {
  
  # model = modelBRTFull$model
  # rasters = rasterLayers
  # distribution="binomial"
  # variable.to.plot=var.n
  # print.limiting=FALSE
  # auto.limiting=FALSE
  # distribution.threshold=thresholdBRT
  # distribution.predicted=predictionBRT
  # val.limiting=NULL
  # export.data.to.plot=FALSE
  # plot.x.lim=NULL
  # plot.y.lim=NULL
  
  if( missing(plot.x.lim) ) { plot.x.lim <- NULL }
  if( missing(plot.y.lim) ) { plot.y.lim <- NULL }
  if( missing(distribution.threshold) ) { distribution.threshold <- NULL }
  if( missing(distribution.predicted) ) { distribution.predicted <- NULL }
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  if( class(model) == "MaxEnt" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- data.frame(Variable=data.frame(response(model,variable.to.plot,expand=0))[,1] , Effect= data.frame(response(model,variable.to.plot,expand=0))[,2] )
      
      data.to.plot <- data.to.plot[ data.to.plot[,1] >= min(getValues(subset(rasters,variable.to.plot)),na.rm=T) &
                                      data.to.plot[,1] <= max(getValues(subset(rasters,variable.to.plot)),na.rm=T) , ]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot]
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "gbm" ) { 
    
    model.predictor <- model$var.names[variable.to.plot]
    
    if(variable.to.plot > length(model$gbm.call$predictor.names) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(model$gbm.call$predictor.names) ) { 
      
      data.to.plot <- data.frame(Variable=plot(model,variable.to.plot, return.grid=TRUE)[,1] , Effect=scale(plot(model,variable.to.plot, return.grid=TRUE)[,2],scale=FALSE))
      
      if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
      
      variable.name <- model$var.names[variable.to.plot]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "mboost" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- getValues(rasters)
      data.to.plot <- data.to.plot[!is.na(data.to.plot[,1]),]
      data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot), ifelse( nrow(data.to.plot) > 100 , 100 , nrow(data.to.plot)  ) ,replace=FALSE),]
      
      data.to.plot.x <- data.to.plot[,variable.to.plot]
      
      means <- colMeans(matrix(data.to.plot[,-variable.to.plot],ncol=length(names(rasters) )-1))
      
      for(m in 1:length(means) ) {
        
        data.to.plot[,(1:length(names(rasters)))[-variable.to.plot][m]] <- means[m]
        
      }
      
      data.to.plot <- data.frame(Variable=data.to.plot.x , Effect= predict(model,as.data.frame(data.to.plot)) )
      data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
      data.to.plot[,2] <- logit2prob(data.to.plot[,2])
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot] 
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "list" ) { 
    
    model.predictor <- names(rasters)[variable.to.plot]
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    data.to.plot <- as.data.frame(rasters)
    data.to.plot <- data.to.plot[complete.cases(data.to.plot),]
    data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot),2500),]
    
    correctLayers <- names(rasters)[which(dataLayersMonotonocity.i == -1)]
    for(correctLayers.i in correctLayers) {
      data.to.plot[,correctLayers.i] <- data.to.plot[,correctLayers.i] * (-1)
    }
    
    for( var in (1:(length(names(rasters))))[-variable.to.plot] ) {
      data.to.plot[,var] <- mean(data.to.plot[,var])
    }
    
    predicted <- monmlp.predict( as.matrix(data.to.plot), model)
    data.to.plot <- data.frame(Variable = as.numeric(data.to.plot[,variable.to.plot]) , Effect=as.numeric(predicted) )
    
    if(names(rasters)[variable.to.plot] == correctLayers) {
      data.to.plot[,"Variable"] <- data.to.plot[,"Variable"] * (-1)
    }
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
    
    if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
    
    variable.name <- model.predictor
    
    if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
      plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
    }
    
    if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
      plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
    }
    
  }
  
  ## -----------------------
  
  if( length(unique(plot.y.lim)) < 2 ) { plot.y.lim <- c(plot.y.lim[1]-0.1,plot.y.lim[1]+0.1) }
  
  ## -----------------------
  
  lab <- ifelse( exists("contribution"),paste0(variable.name," (",round(contribution,digits=2),"%)"),variable.name)
  
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  plot(data.to.plot,ylim=plot.y.lim,xlim=plot.x.lim,lty=1,col="#8E8E8E",type="l",ylab="",xlab=lab,axes=FALSE)
  axis(2,las=2,col="White",col.ticks="Black")
  axis(1,las=0,col="White",col.ticks="Black")
  box()
  title(ylab="Effect on Response",mgp=c(4,1,0)) 
  
  span.vals <- data.to.plot[,2]
  
  if( !is.null(occurrenceRecords)) {
    
    pointsAdd <- data.frame(x=raster::extract(subset(rasters,model.predictor),occurrenceRecords),y=min(plot.y.lim))
    points(pointsAdd,pch=15,col="#5C5C5C")
    
  }
  
  
  if( length(unique(span.vals)) >= 2 ) { 
    
    span.vals.var <- numeric(length(span.vals)-1)
    for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
    span.vals.var <- abs(span.vals.var)
    
    if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min"}
    
    if( tail == "Min") {  
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[1],1]
      
    }
    
    if( tail == "Max" ) {   
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[length(which(span.vals.var >= 0.05))],1] 
      
    }
  }
  
  if( length(unique(span.vals)) < 2 ) { print.limiting <- FALSE ; tipping.point <- 0 }
  
  if( ! auto.limiting ) { tipping.point <- val.limiting }
  
  if( ! is.null(distribution.threshold) ) {
    
    distribution.predicted[distribution.predicted >= distribution.threshold] <- 1
    distribution.predicted[distribution.predicted < distribution.threshold] <- NA
    
    distribution.predicted <- raster::clump(distribution.predicted)
    
    regionsUsed <- raster::extract(distribution.predicted,occurrenceRecords)
    regionsUsed <- regionsUsed[!is.na(regionsUsed)]
    
    distribution.predicted <- distribution.predicted %in% regionsUsed
    distribution.predicted[distribution.predicted == 0] <- NA
    
    distribution.predicted <- distribution.predicted * subset(rasters,variable.to.plot)
    
    tipping.point <- getValues(distribution.predicted)
    tipping.point <- tipping.point[!is.na(tipping.point)]
    # tipping.point <- quantile(tipping.point,probs=c(0.01,0.99))
    tipping.point <- range(tipping.point)
    
    if(dataLayersMonotonocity[variable.to.plot] > 0) { tipping.point <- min(tipping.point)}
    if(dataLayersMonotonocity[variable.to.plot] < 0) { tipping.point <- max(tipping.point)}
    
  }
  
  
  cat( paste0("\n"))
  cat( paste0("Predictor: ",model.predictor ))
  cat( paste0("\n"))
  
  if( print.limiting ) {
    
    cat( paste0("Limiting point: ",tipping.point))
    cat( paste0("\n"))
    abline(v = tipping.point , lty=2 , col="#AC4F4F")
    
  }
  
  options(warn=0)
  
  tippingPoint <<- tipping.point
  if(export.data.to.plot) { return(data.to.plot)}
  
}

## -----------------------

modelPlotOld <- function(model,rasters,distribution,variable.to.plot,print.limiting,auto.limiting,distribution.threshold,distribution.predicted,val.limiting,export.data.to.plot,plot.x.lim,plot.y.lim) {
  
  if( missing(plot.x.lim) ) { plot.x.lim <- NULL }
  if( missing(plot.y.lim) ) { plot.y.lim <- NULL }
  if( missing(distribution.threshold) ) { distribution.threshold <- NULL }
  if( missing(distribution.predicted) ) { distribution.predicted <- NULL }
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  if( class(model) == "MaxEnt" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- data.frame(Variable=data.frame(response(model,variable.to.plot,expand=0))[,1] , Effect= data.frame(response(model,variable.to.plot,expand=0))[,2] )
      
      data.to.plot <- data.to.plot[ data.to.plot[,1] >= min(getValues(subset(rasters,variable.to.plot)),na.rm=T) &
                                      data.to.plot[,1] <= max(getValues(subset(rasters,variable.to.plot)),na.rm=T) , ]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot]
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "gbm" ) { 
    
    model.predictor <- model$var.names[variable.to.plot]
    
    if(variable.to.plot > length(model$gbm.call$predictor.names) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(model$gbm.call$predictor.names) ) { 
      
      data.to.plot <- data.frame(Variable=plot(model,variable.to.plot, return.grid=TRUE)[,1] , Effect=scale(plot(model,variable.to.plot, return.grid=TRUE)[,2],scale=FALSE))
      
      if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
      
      variable.name <- model$var.names[variable.to.plot]
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "mboost" ) { 
    
    names.on.model <- names(rasters)
    model.predictor <- names.on.model[variable.to.plot]
    
    if(variable.to.plot > length(names.on.model) ) { cat("Predictor does no exist in the model") }
    
    if(variable.to.plot <= length(names.on.model) ) { 
      
      data.to.plot <- getValues(rasters)
      data.to.plot <- data.to.plot[!is.na(data.to.plot[,1]),]
      data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot), ifelse( nrow(data.to.plot) > 100 , 100 , nrow(data.to.plot)  ) ,replace=FALSE),]
      
      data.to.plot.x <- data.to.plot[,variable.to.plot]
      
      means <- colMeans(matrix(data.to.plot[,-variable.to.plot],ncol=length(names(rasters) )-1))
      
      for(m in 1:length(means) ) {
        
        data.to.plot[,(1:length(names(rasters)))[-variable.to.plot][m]] <- means[m]
        
      }
      
      data.to.plot <- data.frame(Variable=data.to.plot.x , Effect= predict(model,as.data.frame(data.to.plot)) )
      data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
      data.to.plot[,2] <- logit2prob(data.to.plot[,2])
      
      if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
        plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
      }
      
      if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
        plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
      }
      
      variable.name <- names.on.model[variable.to.plot] 
      xlim <- plot.x.lim
      ylim <- plot.y.lim
      
    }
    
  }
  
  ## -----------------------
  
  if( class(model) == "list" ) { 
    
    model.predictor <- names(rasters)[variable.to.plot]
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    data.to.plot <- as.data.frame(rasters)
    data.to.plot <- data.to.plot[complete.cases(data.to.plot),]
    data.to.plot <- data.to.plot[sample(1:nrow(data.to.plot),2500),]
    
    correctLayers <- names(rasters)[which(dataLayersMonotonocity.i == -1)]
    for(correctLayers.i in correctLayers) {
      data.to.plot[,correctLayers.i] <- data.to.plot[,correctLayers.i] * (-1)
    }
    
    for( var in (1:(length(names(rasters))))[-variable.to.plot] ) {
      data.to.plot[,var] <- mean(data.to.plot[,var])
    }
    
    predicted <- monmlp.predict( as.matrix(data.to.plot), model)
    data.to.plot <- data.frame(Variable = as.numeric(data.to.plot[,variable.to.plot]) , Effect=as.numeric(predicted) )
    
    if(names(rasters)[variable.to.plot] == correctLayers) {
      data.to.plot[,"Variable"] <- data.to.plot[,"Variable"] * (-1)
    }
    
    data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
    
    if( distribution != "gaussian" ) {  data.to.plot[,2] <- logit2prob(data.to.plot[,2]) }
    
    variable.name <- model.predictor
    
    if( missing(plot.x.lim) | is.null(plot.x.lim) ) { 
      plot.x.lim <- c(min(data.to.plot[,1]),max(data.to.plot[,1]))
    }
    
    if( missing(plot.y.lim) | is.null(plot.y.lim) ) { 
      plot.y.lim <- c(min(data.to.plot[,2]),max(data.to.plot[,2]))
    }
    
  }
  
  ## -----------------------
  
  if( length(unique(plot.y.lim)) < 2 ) { plot.y.lim <- c(plot.y.lim[1]-0.1,plot.y.lim[1]+0.1) }
  
  ## -----------------------
  
  lab <- ifelse( exists("contribution"),paste0(variable.name," (",round(contribution,digits=2),"%)"),variable.name)
  
  par(mar = c(4.5, 5.5, 4.5, 4.5))
  plot(data.to.plot,ylim=plot.y.lim,xlim=plot.x.lim,lty=1,col="#8E8E8E",type="l",ylab="",xlab=lab,axes=FALSE)
  axis(2,las=2,col="White",col.ticks="Black")
  axis(1,las=0,col="White",col.ticks="Black")
  box()
  title(ylab="Effect on Response",mgp=c(4,1,0)) 
  
  span.vals <- data.to.plot[,2]
  
  if( length(unique(span.vals)) >= 2 ) { 
    
    span.vals.var <- numeric(length(span.vals)-1)
    for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) } 
    span.vals.var <- abs(span.vals.var)
    
    if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min"}
    
    if( tail == "Min") {  
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[1],1]
      
    }
    
    if( tail == "Max" ) {   
      
      tipping.point <- data.to.plot[which(span.vals.var >= 0.05)[length(which(span.vals.var >= 0.05))],1] 
      
    }
  }
  
  if( length(unique(span.vals)) < 2 ) { print.limiting <- FALSE ; tipping.point <- 0 }
  
  if( ! auto.limiting ) { tipping.point <- val.limiting }
  
  if( ! is.null(distribution.threshold) ) {
    
    distribution.predicted[distribution.predicted >= distribution.threshold] <- 1
    distribution.predicted[distribution.predicted < distribution.threshold] <- NA
    distribution.predicted <- distribution.predicted * subset(rasters,variable.to.plot)
    
    tipping.point <- getValues(distribution.predicted)
    tipping.point <- tipping.point[!is.na(tipping.point)]
    # tipping.point <- quantile(tipping.point,probs=c(0.01,0.99))
    tipping.point <- range(tipping.point)
    
    if(dataLayersMonotonocity[variable.to.plot] > 0) { tipping.point <- min(tipping.point)}
    if(dataLayersMonotonocity[variable.to.plot] < 0) { tipping.point <- max(tipping.point)}
    
  }
  
  
  cat( paste0("\n"))
  cat( paste0("Predictor: ",model.predictor ))
  cat( paste0("\n"))
  
  if( print.limiting ) {
    
    cat( paste0("Limiting point: ",tipping.point))
    cat( paste0("\n"))
    abline(v = tipping.point , lty=2 , col="#AC4F4F")
    
  }
  
  options(warn=0)
  
  tippingPoint <<- tipping.point
  if(export.data.to.plot) { return(data.to.plot)}
  
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
  
  if(reclassToOne) {
    predicted.distribution <- predicted.distribution + ( min(getValues(predicted.distribution),na.rm=T) * (-1))
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=T)
  }
  
  return(predicted.distribution)
  
}

## -----------------------

accuracyPredicted <- function(predicted.distribution,speciesData,type) {
  
  observed <- speciesData$PA
  predicted <- raster::extract(predicted.distribution,speciesData[,c("Lon","Lat")])
  
  predicted.accuracy <- accuracy( observed , predicted , threshold = 100 )
  predicted.accuracy$AUC <- SDMTools::auc(observed ,predicted)
  
  if(type == "tss") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
  if(type == "area") { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity >= 0.95)),] }
  if(type == "Kappa") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$Kappa),] }
  
  predicted.accuracy <- data.frame( threshold =  predicted.accuracy$threshold ,
                                    auc = predicted.accuracy$AUC ,
                                    specificity = predicted.accuracy$specificity ,
                                    sensitivity = predicted.accuracy$sensitivity ,
                                    tss = predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 )
  
  return(predicted.accuracy)
  
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
    predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity > reclass.threshold,]
    predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  return(predicted.distribution)
}

## -----------------------

generate.points.from.prediction <- function(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,threshold) {
  
  set.seed(42)
  
  # Best reclassification
  
  predictive.model.sp.1.v <- unique(getValues(predictive.model.sp.1))
  predictive.model.sp.2.v <- unique(getValues(predictive.model.sp.2))
  
  if( sum(!is.na(predictive.model.sp.1.v)) > 2 ) { stop("Wrong reclass file type") }
  if( sum(!is.na(predictive.model.sp.2.v)) > 2 ) { stop("Wrong reclass file type") }
  
  if( class(predictive.model.sp.1)[1] != "RasterLayer") { raster.predict.1 <- raster(predictive.model.sp.1) }
  if( class(predictive.model.sp.1)[1] == "RasterLayer") { raster.predict.1 <- predictive.model.sp.1 }
  
  predict.1 <- as.data.frame(raster.predict.1,xy=TRUE)
  presences.1 <- predict.1[which(predict.1[,3] == 1),1:2]
  
  if( class(predictive.model.sp.2)[1] != "RasterLayer") { raster.predict.2 <- raster(predictive.model.sp.2) }
  if( class(predictive.model.sp.2)[1] == "RasterLayer") { raster.predict.2 <- predictive.model.sp.2 }
  
  predict.2 <- as.data.frame(raster.predict.2,xy=TRUE)
  presences.2 <- predict.2[which(predict.2[,3] == 1),1:2]
  colnames(presences.2) <- colnames(presences.1)
  
  max.points <- min(nrow(presences.1),nrow(presences.2))
  
  set.seed(42)
  presences.1 <- cbind(presences.1[ sample(1:nrow(presences.1),max.points,replace=F) ,],rep(1, max.points )) ; colnames(presences.1) <- c("Lon","Lat","Sp")
  
  set.seed(42)
  presences.2 <- cbind(presences.2[ sample(1:nrow(presences.2),max.points,replace=F) ,],rep(2, max.points )) ; colnames(presences.2) <- c("Lon","Lat","Sp")
  
  m <- leaflet()
  m <- addTiles(m)
  m <- addCircleMarkers(m, lng=c(presences.1[,1],presences.2[,1]), lat=c(presences.1[,2],presences.2[,2]), popup=paste0( "Species record ") , radius = 3, color = c( rep("Red",nrow(presences.1)) , rep("Blue",nrow(presences.2))) )
  print(m)
  
  return(rbind(presences.1,presences.2))
  
}

## -----------------------

pca.ordination <- function(rasters.1,rasters.2,rnd.points.sp) {
  
  rnd.points.env.sp.1 <- as.data.frame(rasters.1,na.rm=TRUE)
  rnd.points.env.sp.1 <- rnd.points.env.sp.1[sample(1:nrow(rnd.points.env.sp.1),1000,replace=FALSE),]
  
  rnd.points.env.sp.2 <- as.data.frame(rasters.2,na.rm=TRUE)
  rnd.points.env.sp.2 <- rnd.points.env.sp.2[sample(1:nrow(rnd.points.env.sp.2),1000,replace=FALSE),]
  
  rnd.points.sp.1.env <- raster::extract(rasters.1,rnd.points.sp[rnd.points.sp$Sp==1,1:2])
  rnd.points.sp.1.env <- rnd.points.sp.1.env[complete.cases(rnd.points.sp.1.env),]
  rnd.points.sp.2.env <- raster::extract(rasters.2,rnd.points.sp[rnd.points.sp$Sp==2,1:2])
  rnd.points.sp.2.env <- rnd.points.sp.2.env[complete.cases(rnd.points.sp.2.env),]
  
  env.used <- rbind(rnd.points.env.sp.1,rnd.points.env.sp.2,rnd.points.sp.1.env,rnd.points.sp.2.env)
  # env.used <- rbind(rnd.points.sp.1.env,rnd.points.sp.2.env)
  pca.env.used <- dudi.pca(env.used,scannf=F,nf=2)
  ecospat.plot.contrib(contrib=pca.env.used$co, eigen=pca.env.used$eig)
  
  sp.vect <- c( rep("env.sp1",nrow(rnd.points.env.sp.1)) , 
                rep("env.sp2",nrow(rnd.points.env.sp.2)) , 
                rep("sp1",nrow(rnd.points.sp.1.env)) , 
                rep("sp2",nrow(rnd.points.sp.2.env)) )
  
  
  # PCA scores for the whole study area
  scores.globclim <- pca.env.used$li
  # PCA scores for the species native distribution
  scores.sp.1 <- suprow(pca.env.used,rnd.points.sp.1.env)$li
  scores.sp.2 <- suprow(pca.env.used,rnd.points.sp.2.env)$li
  
  scores.clim.sp.1 <- suprow(pca.env.used,rnd.points.env.sp.1)$li
  scores.clim.sp.2 <- suprow(pca.env.used,rnd.points.env.sp.2)$li
  
  grid.clim.sp.1 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.1,
                                          sp=scores.sp.1, R=100,
                                          th.sp=0)
  grid.clim.sp.2 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.2,
                                          sp=scores.sp.2, R=100,
                                          th.sp=0)
  
  newList <- list("grid.clim.sp.1" = grid.clim.sp.1, "grid.clim.sp.2" = grid.clim.sp.2)
  
  return(newList)
  
}

## -----------------------

prepare.data.niche.overlap <- function(rasters.1,rasters.2,presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,type) {
  
  if(type == 1) {
    
    rnd.points.sp <- data.frame( rbind(presences.sp.1,presences.sp.2) , Sp = c(rep(1,nrow(presences.sp.1)) , rep(2,nrow(presences.sp.2))))
    
  }
  
  if(type == 2) {
    
    rnd.points.sp <- generate.points.from.prediction(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,0.95)
    
  }
  
  environment.used.by.species.1 <- raster::extract(rasters.1,rnd.points.sp[rnd.points.sp$Sp==1,1:2])
  environment.used.by.species.1 <- environment.used.by.species.1[complete.cases(environment.used.by.species.1),]
  environment.used.by.species.2 <- raster::extract(rasters.2,rnd.points.sp[rnd.points.sp$Sp==2,1:2])
  environment.used.by.species.2 <- environment.used.by.species.2[complete.cases(environment.used.by.species.2),]
  
  data <- data.frame(  species = c(rep(1,nrow(environment.used.by.species.1)),rep(2,nrow(environment.used.by.species.2))) ,
                       rbind(environment.used.by.species.1,environment.used.by.species.2) )
  
  print(aggregate(data[2:ncol(data)], data[1], mean))
  return(data)
  
}

## -----------------------

niche.overlap.plot <- function(data,nsamples) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  clrs <- c("black", "red")  # colors for each species
  nsamples <- 10
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  species.par.data <- tapply(1:nrow(data), data$species, function(ii) X = data[ii, 2:ncol(data)])
  niche.plot(niche.par = species.par, ndens=1000, niche.data = species.par.data, pfrac = 0, col = clrs, xlab = expression("Niche overlap"))
  
}

## -----------------------

mean.niche.overlap.plot <- function(data,nsamples,colors.sp) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  # Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
  # accuracy.  the variable over.stat can be supplied directly to the
  # overlap.plot function
  
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 10000, alpha = c(0.95, 0.99))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("The mean overlap metrics calculated across iteratations for both niche region sizes"))
  cat( paste0("\n"))
  
  over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
  print(round(over.mean[, ,  1], 2))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Credible interval"))
  cat( paste0("\n"))
  
  over.cred <- apply(over.stat * 100, c(1:2, 4), quantile, prob = c(0.025, 0.975), na.rm = TRUE)
  print(round(over.cred[, , , 1]),2)  # display alpha = .95 niche region
  
  # Overlap plot.Before you run this, make sure that you have chosen your
  # alpha ecoregionsLevel.
  
  clrs <- colors.sp  # colors for each species
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 1000, alpha = 0.95)
  overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE, xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
  
}

## -----------------------

prepare.data.niche.overlap.2 <- function(rasters.1,rasters.2,presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,type=1) {
  
  
  if(type == 1) {
    
    rnd.points.sp <- data.frame( rbind(presences.sp.1,presences.sp.2) , Sp = c(rep(1,nrow(presences.sp.1)) , rep(2,nrow(presences.sp.2))))
    
  }
  
  if(type == 2) {
    
    rnd.points.sp <- generate.points.from.prediction(presences.sp.1,presences.sp.2,predictive.model.sp.1,predictive.model.sp.2,0.95)
    
  }
  
  grid.clim <- pca.ordination(rasters.1,rasters.2,rnd.points.sp)
  
  return(grid.clim)
  
}

## -----------------------

dataTranformation <- function(vector,new.range) {
  
  exp.lm <- data.frame( input = sort(unique(vector)),
                        outup = seq(0,1,length.out = length(unique(vector))))
  
  fit <- lm(outup ~ input, data=exp.lm)
  
  new.data <- data.frame(input=vector)
  new.data <- unlist(predict.lm(fit,new.data ))
  new.data <- as.vector(new.data)
  new.data[which(new.data == min(new.data))] <- 0
  new.data[which(new.data == max(new.data))] <- 1
  
  return(new.data)
  
}

# ---------------------------------------------------------------------------------------------------------------

getTaxonomyNewSpecies <- function(taxa,rank) {
  
  library(worms)
  library(worrms)
  
  if( rank != "species" ) {
    
    error <- FALSE
    
    tryCatch({ list.of.species <- wormsbynames(taxon_names=taxa)$AphiaID }, error = function(e) { error <<- TRUE })
    
    if( ! error ) {     
      
      final.list <- data.frame()
      higher.taxa <- TRUE
      
      while(higher.taxa){
        
        list.of.childen <- data.frame()
        
        for(sp in 1:length(list.of.species)) {
          
          for(offset in seq(1,1000000,50)) {
            
            if( exists("list.of.childen.t") ) {  rm(list.of.childen.t)  }
            
            tryCatch( list.of.childen.t <- wm_children(id=list.of.species[sp],offset=offset,marine_only=FALSE) , error=function(e) { error <- TRUE })
            
            if( ! exists("list.of.childen.t") ) {  break  }
            
            if( exists("list.of.childen.t") ) { 
              
              list.of.childen <- rbind(list.of.childen,list.of.childen.t)
              rm(list.of.childen.t) 
              
            }
            
          }
          
        }
        
        list.of.species <- list.of.childen$AphiaID
        list.of.rank <- list.of.childen$rank
        
        to.include <- which(list.of.childen$rank == "Species")
        
        if( "Species" %in% list.of.rank ) { 
          
          final.list <- rbind(final.list,list.of.childen[to.include,] )
          
        }
        
        if( "Species" %in% list.of.rank & length(unique(list.of.rank)) == 1  ) { higher.taxa <- FALSE }
        if( "Species" %in% list.of.rank & length(unique(list.of.rank)) > 1  ) { higher.taxa <- TRUE ; list.of.species <- list.of.species[-to.include] }
        
        # wormsbyid(list.of.species)$scientificname
        
      }
      
      taxa <- final.list$AphiaID
      
    }
    
    if( error ) { taxa <- NULL }
    
  }
  
  # --------------------------
  
  if( rank == "species" ) {
    
    error <- FALSE
    
    tryCatch({ taxa <- wormsbynames(taxon_names=taxa)$AphiaID }, error = function(e) { error <<- TRUE })
    
    if( error ) { taxa <- NULL }
    
  }
  
  # ------------------ 
  
  if( !is.null(taxa)) {
    
    taxa <- unique(taxa)
    
    new.species <- data.frame()
    
    for( t in 1:length(taxa) ) {
      
      taxa.to.search <- taxa[t]
      new.species <- rbind(new.species,wormsbyid(taxa.to.search))
      
    }
    
    return(new.species)
    
  }
  
  # ------------------ 
  
  if( is.null(taxa)) {
    
    return(NULL)
    
  }
}

# ---------------------------------------------------------------------------------------------------------------

getExternalDataObis <- function(taxa) {
  
  if( missing(taxa)) { errormessage("no taxa (Worms name) introduced.") }
  
  library(robis)
  
  for( i in 1:length(taxa) ) { 
    
    print(i)
    
    taxa.i <- taxa[i]
    error <- FALSE
    
    tryCatch( my_occs_obis <- occurrence(scientificname = taxa.i ) , error=function(e) { error <<- TRUE })
    
    if( ! error ) { if( nrow(my_occs_obis) == 0 ) { my_occs_obis <- data.frame() } }
    
    if( error ) { my_occs_obis <- data.frame() }
    
    if( nrow(my_occs_obis) > 0) {
      
      my_occs_obis <- subset(my_occs_obis, my_occs_obis$decimalLongitude !=0 & my_occs_obis$decimalLatitude !=0)
      
    }
    
    if( nrow(my_occs_obis) > 0) {
      
      my_occs_obisInfo <- my_occs_obis$dataset_id
      my_occs_obisInfo <- unique(my_occs_obis$dataset_id)
      
      for(z in 1:length(my_occs_obisInfo) ) {
        
        error <- TRUE
        errortrials <- 0
        
        while(error & errortrials < 10) {
          error <- FALSE
          errortrials <- errortrials + 1
          tryCatch(  z.Res <- RJSONIO::fromJSON(paste0("https://api.obis.org/v3/dataset/",my_occs_obisInfo[z])) , error=function(e) { error <- TRUE })
        }
        
        if(!error) {  
          
          institutionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"institutionCode"]
          collectionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"collectionCode"]
          
          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"accessRights"] <- z.Res$results[[1]]$intellectualrights
          
          z.Res <- paste0( z.Res$results[[1]]$citation,
                           " ",
                           ifelse(!is.na(institutionCode) | !is.null(institutionCode) , institutionCode , ""),
                           " ",
                           ifelse(!is.na(collectionCode) | !is.null(collectionCode) , collectionCode , ""),
                           " (Available: Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO. www.iobis.org. Accessed: ", Sys.time())
          
          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"bibliographicCitation"] <- z.Res
          
        }
        
      }
      
    }
    
  }
  
  return(my_occs_obis)
  
}

# ---------------------------------------------------------------------------------------------------------------

getExternalDataGbif <- function(taxa) {
  
  if( missing(taxa)) { errormessage("no taxa (Worms name) introduced.") }
  
  library(rgbif)
  
  nRecords <- gbif(strsplit(as.character(taxa.i), " ")[[1]][1], strsplit(as.character(taxa.i), " ")[[1]][2], geo=T, removeZeros=T , download=FALSE )
  
  if( nRecords > 300 ) {
    
    seqListing <- seq(0,nRecords,by =300)
    if(max(seqListing) < nRecords) { seqListing <- c(seqListing,nRecords) }
    parallelChunks <- data.frame(from = seqListing[-length(seqListing)], to = c(seqListing[-c(1,length(seqListing))] -1 , nRecords ) )
    
    my_occs_gbif <- data.frame()
    
    for(ch in 1:nrow(parallelChunks)) { 
      
      tmpfile <- paste(tempfile(), ".json", sep = "")
      download.file(paste0("http://api.gbif.org/v1/occurrence/search?scientificname=",strsplit(as.character(taxa.i), " ")[[1]][1],"+",strsplit(as.character(taxa.i), " ")[[1]][2],"&offset=",parallelChunks[ch,1],"&limit=300"), tmpfile, quiet = TRUE, method="curl") 
      json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
      json <- chartr("\a\v", "  ", json)
      x <- jsonlite::fromJSON(json)
      r <- x$results
      r <- r[, !sapply(r, class) %in% c("data.frame", "list")]
      rownames(r) <- NULL
      my_occs_gbif <- smartbind(my_occs_gbif,r)
      
    }
    
    if (length(my_occs_gbif) == 0) {
      my_occs_gbif <- data.frame()
    }
    if (length(my_occs_gbif) == 1) {
      my_occs_gbif <- my_occs_gbif[[1]]
    }
    
  }
  
  if( nRecords <= 300 ) {
    
    tmpfile <- paste(tempfile(), ".json", sep = "")
    test <- try(download.file(paste0("http://api.gbif.org/v1/occurrence/search?scientificname=",strsplit(as.character(taxa.i), " ")[[1]][1],"+",strsplit(as.character(taxa.i), " ")[[1]][2],"&offset=",0,"&limit=300"), tmpfile, quiet = TRUE))
    
    json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
    json <- chartr("\a\v", "  ", json)
    x <- jsonlite::fromJSON(json)
    r <- x$results
    
    if( length(r) > 0) {
      
      r <- r[, !sapply(r, class) %in% c("data.frame", "list")]
      rownames(r) <- NULL
      my_occs_gbif <- r 
      
    }
    
  }
  
  if( exists("my_occs_gbif") ) { if( is.null(my_occs_gbif) ) { my_occs_gbif <- data.frame() } }
  
  if( ! exists("my_occs_gbif") ) { my_occs_gbif <- data.frame() }
  
  if( ! is.null(my_occs_gbif$decimalLatitude) & nrow(my_occs_gbif) > 0 ) {
    
    my_occs_gbif <- subset(my_occs_gbif, decimalLatitude !=0 & decimalLongitude !=0)
    
  }
  
  if( nrow(my_occs_gbif) > 0 ) {
    
    my_occs_gbif_all <- unique(my_occs_gbif$datasetKey)
    
    for(z in 1:length(my_occs_gbif_all) ) {
      
      z.Res <- gbif_citation(x=my_occs_gbif_all[z])
      
      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z] ,"accessRights"] <- ifelse(!is.null(z.Res$rights),z.Res$rights,"")
      
      z.Res <- z.Res$citation$citation
      
      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z],"bibliographicCitation"] <- z.Res
      
    }
    
  }
  
  return(my_occs_gbif)
  
}

# ---------------------------------------------------------------------------------------------------------------

mySmartBind <- function(dfr1,dfr2) {
  
  common_cols <- intersect(colnames(dfr1), colnames(dfr2))
  common_cols <- rbind(
    subset(dfr1, select = common_cols), 
    subset(dfr2, select = common_cols)
  )
  return(common_cols)
}

