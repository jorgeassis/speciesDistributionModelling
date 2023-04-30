## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------
##
##  Machine Learning Species Distribution Modelling [ Ver.301 ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

closeAllConnections()
gc(reset=TRUE)

if( ! "h3js" %in% rownames(installed.packages()) ) { devtools::install_github("saurfang/h3js") }
if( ! "h3jsr" %in% rownames(installed.packages()) ) { remotes::install_github("obrl-soil/h3jsr") }

library( rnaturalearth )
library( h3js )
library( h3jsr )
library( sf )
library(maptools)

# source("Dependencies/mainFunctions.R")

themeMapEqual <- theme( text = element_text(family = "Helvetica", color = "#22211d"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_line(color = "black", size = 0.1), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#FFFFFF", color = NA), panel.background = element_rect(fill = "#FFFFFF", color = NA), panel.border = element_blank(), legend.background = element_rect(fill = "#FFFFFF", color = NA), legend.position="bottom", legend.box = "horizontal", legend.margin=margin(0,0,0,0), legend.box.margin=margin(-10,-10,-10,-10), legend.key.height= unit(0.25, 'cm'), legend.key.width= unit(0.75, 'cm') )

projection <- CRS("+proj=robin +over")

bb <- sf::st_union(sf::st_make_grid( st_bbox(c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)), n = 100))
bb <- st_transform(bb, projection)

worldMap <- ne_countries(scale = 10, returnclass = "sp")
worldMap <- spTransform(worldMap, CRSobj = projection)
worldMap <- gBuffer(worldMap, byid=TRUE, width=0.001)
worldMap <- crop(worldMap, as(bb, "Spatial"))

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.box.background = element_rect(fill='#FFFFFF'),
    panel.border = element_blank()
  )

## ------------------------

mainGlobalMap <- ggplot() + geom_sf(data = st_as_sf(ne_countries(scale = 10, returnclass = "sp")), color = "#CDCDCD", fill = "#CDCDCD", size=0.1) + themeMap

mainGlobalMapEqual <- ggplot() + 
  geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill="#CDCDCD", colour = "#CDCDCD" , size=0.1 ) +
  # coord_sf(xlim = c(-16808016, 16927051), ylim = c(-6525155,6242353), expand = FALSE) + 
  themeMapEqual

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Species by Species maps

directoryMaps <- paste0(stackResultsFolder,"/Maps/perSpecies/")
scenariosToPredict <- c("Baseline","ssp119","ssp370","ssp585")
modelType <- ""

if( ! dir.exists(directoryMaps) ) { dir.create(directoryMaps, recursive = TRUE) }

for( scenario in scenariosToPredict ) {

    allMaps <- list.files( mainResultsDirectory, pattern = "Global.RData", full.names = T, recursive=T )
    allMapsNames <- list.files( mainResultsDirectory, pattern = "Global.RData", full.names = F , recursive=T )
    allMaps <- allMaps[grepl("ensembleReclassReachable",allMaps)] # Reachable unConstrained
    allMapsNames <- allMapsNames[grepl("ensembleReclassReachable",allMapsNames)] # Reachable unConstrained
    allMaps <- allMaps[grepl(scenario,allMaps)] 
    allMapsNames <- allMapsNames[grepl(scenario,allMapsNames)] # Baseline ssp119 ssp585
    allMapsNames <- gsub(paste0("/Predictions/",scenario,"/ensembleReclassReachableGlobal.RData"),"",allMapsNames); allMapsNames
    
    # allMaps <- allMaps[which( allMapsNames == species)]
    # allMapsNames <- allMapsNames[which( allMapsNames == species)]
    
    cl <- parallel::makeCluster(nCores)
    registerDoParallel(cl)
    
    mappingParallel <- foreach(i = 1:length(allMaps) ) %dopar% {

      source(mainfunctionsFile, local = TRUE)
      source(mainConfigFile, local = TRUE)
      
      ## ---------------------
      
      m <- allMaps[i]
      name <- allMapsNames[i]
      
      if( scenario == "Baseline" ) {
        occurrecence <- loadRData(paste0(mainResultsDirectory,"/",name,"/Data/modelData.RData"))
        occurrecence <- occurrecence@coords[occurrecence@pa==1,]
      }
      
      rasterMap <- loadRData( m )
      
      resolutionH3 <- 3
      rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
      rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]

      if(nrow(rasterMapDF) == 0 ) { 
        rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
        rasterMapDF <- rasterMapDF[1:2,]
      }
      
      rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
      rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))

      rasterMapDF.polygons <- h3jsr::cell_to_polygon(input = rasterMapDF$hex, simple = FALSE)
      rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
      rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
      ref <- st_crs(rasterMapDF.polygons)
      rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
      st_crs(rasterMapDF.polygons) <- ref

      minLegend <- min(rasterMapDF.polygons$value)
      maxLegend <- max(rasterMapDF.polygons$value)
      
      nColors <- 7
      colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
      colorBreaks[1] <- min(rasterMapDF.polygons$value)
      colorBreaks[nColors] <- max(rasterMapDF.polygons$value)
      
      myColors <- c("#E4FAFF","#E8EF15","#ec7a06","#e31515", "#450751") # Blue Yellow Red Purple
      myColors <- colorRampPalette(myColors)(nColors)
      
      #----------------------
      
      plot1 <- mainGlobalMap +
        geom_sf(data = rasterMapDF.polygons, aes(fill="black"), colour ="gray", size = 0.05) +
        theme(legend.position="none",
              legend.margin=margin(0,0,0,0),
              legend.key.height= unit(0.25, 'cm'),
              legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
      
      if( scenario == "Baseline" ) {
        plot1 <- plot1 + geom_point(data = occurrecence, aes(x = X, y = Y), size = 0.5, color="Black")
      }
      
      # plot1
      
      pdf(file=paste0(directoryMaps,modelType,"/",name," ",scenario,".pdf"),width=12,useDingbats=FALSE)
      print(plot1)
      dev.off()
      
      return(NULL)
      
    }
    
    stopCluster(cl); rm(cl)
    closeAllConnections()
    gc(reset=TRUE)

}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Continuous Map 

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive = TRUE)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
files <- files[!grepl("ExchangeRatio1",files)] 
files <- files[!grepl("extinction",files)] 
files <- files[!grepl("Changes",files)] 

autoLegend <- FALSE

cl <- parallel::makeCluster(10)
registerDoParallel(cl)

mappingParallel <- foreach(m = 1:length(files) ) %dopar% {
  
  source(mainfunctionsFile, local = TRUE)
  source(mainConfigFile, local = TRUE)
  
  ## ---------------------
  
  rasterMap <- loadRData( files[m] )
  mapName <- gsub(paste0(stackResultsFolder,"/Maps"),"",files[m])
  mapName <- gsub("/","",mapName)
  mapName <- gsub(".RData","",mapName)

  if( grepl("Refugia", mapName ) | grepl("Loss", mapName ) ) {
    rasterMapPresent <- files[grepl("RichnessBaseline",files)]
    rasterMapPresent <- rasterMapPresent[!grepl("Uncertainty",rasterMapPresent)]
    rasterMapPresent <- loadRData(rasterMapPresent)
    rasterMapPresent[ rasterMapPresent == 0 ] <- NA
    rasterMap <- mask(rasterMap,rasterMapPresent)
  }

  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])

  if( ! grepl("Refugia", mapName ) & ! grepl("Loss", mapName ) & ! grepl("ExchangeRatio", mapName ) ) {
    rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]
  }
  
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  rasterMapDF.polygons <- h3jsr::cell_to_polygon(input = rasterMapDF$hex, simple = FALSE)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
  
  minLegend <- min(rasterMapDF.polygons$value)
  maxLegend <- max(rasterMapDF.polygons$value)
  
  autoLegendDF <- data.frame(mapName=mapName,minLegend=minLegend,maxLegend=maxLegend)
    
  if( ! autoLegend ) {
    
    if( grepl("Uncertainty", mapName ) ) {
      minLegend <- 0
      maxLegend <- 1
    }
    if( grepl("Gain", mapName ) ) {
      minLegend <- 0
      maxLegend <- 14
    }
    if( grepl("Loss", mapName ) ) {
      minLegend <- 0
      maxLegend <- 13
    }
    if( grepl("Refugia", mapName ) ) {
      minLegend <- 0
      maxLegend <- 19
    }
    if( grepl("speciesRichness", mapName ) ) {
      minLegend <- 1
      maxLegend <- 20
    }
    if( grepl("speciesExchangeRatio", mapName ) ) {
      minLegend <- 0
      maxLegend <- 1
    }
    if( grepl("PerBaselineRichness", mapName ) ) {
      minLegend <- 0
      maxLegend <- 1
    }
  }
  
  nColors <- 7
  colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
  colorBreaks[1] <- min(rasterMapDF.polygons$value)
  colorBreaks[nColors] <- max(rasterMapDF.polygons$value)
  
  myColors <- c("#E4FAFF","#E8EF15","#ec7a06","#e31515", "#450751") # Blue Yellow Red Purple
  myColors <- colorRampPalette(myColors)(nColors)
  
  #----------------------

  hexagons <- st_transform(rasterMapDF.polygons, projection)
  hexagons <- as_Spatial(hexagons)
  hexagons <- gBuffer(hexagons, byid=TRUE, width=0.001)
  hexagons <- crop(hexagons, as(bb, "Spatial"))
  hexagons <- st_as_sf(hexagons)

  plot1 <- mainGlobalMapEqual +
  geom_sf(data = hexagons, aes(fill=value), colour ="black", size = 0.05) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
    geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())

  pdf(file=gsub(".RData",".pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
  return(autoLegendDF)
  
}

stopCluster(cl); rm(cl)
closeAllConnections()
gc(reset=TRUE)

autoLegendValues <- do.call(rbind,mappingParallel)
autoLegendValues

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Zero centered map

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive = T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
files <- files[!grepl("extinction",files)] 
files <- files[!grepl("Changes",files)] 
files <- files[!grepl("Uncertainty",files)] 
files <- files[grepl("speciesRichness",files)] 

maskIntertidal <- raster("Dependencies/Data/Rasters/coastLineRes005.tif")

maskDistribution <- calc(stack(sapply(files, function(x) { loadRData( x ) } )), sum)
maskDistribution[maskDistribution > 0] <- 1
maskDistribution[maskDistribution == 0] <- NA

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive = T)
filesNames <- list.files( stackResultsFolder, pattern = "RData", full.names = F , recursive = T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
filesNames <- filesNames[grepl("Reachable",filesNames)] # Reachable unConstrained
files <- files[grepl("Changes",files)] 
filesNames <- filesNames[grepl("Changes",filesNames)]; filesNames

autoLegend <- FALSE

for( m in 1:length(files) ) {
  
  rasterMap <- loadRData( files[m] )
  mapName <- gsub(".RData","",filesNames[m])
  #rasterMap <- raster::mask(rasterMap,maskIntertidal)
  
  rasterMap[intersect(Which(is.na(maskDistribution), cells=TRUE) ,Which(!is.na(rasterMap), cells=TRUE))] <- NA
  
  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])

  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  rasterMapDF.polygons <- h3jsr::cell_to_polygon(input =rasterMapDF$hex, simple = FALSE)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
  
  minLegend <- min(rasterMapDF.polygons$value)
  maxLegend <- max(rasterMapDF.polygons$value)
  
  if ( ! autoLegend ) {
    
    minLegend <- -15
    maxLegend <- 12
    
  }

  plot2 <- mainGlobalMap +

    geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.05) +
    scale_colour_gradient2(low = "#450751",mid = "white",high = "#135107",midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar", aesthetics = "fill", limits=c(minLegend,maxLegend) ) + 
    
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())

  pdf(file=gsub(".RData",".pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot2)
  dev.off()
  
  # ---------
  
  hexagons <- st_transform(rasterMapDF.polygons, projection)
  hexagons <- as_Spatial(hexagons)
  hexagons <- gBuffer(hexagons, byid=TRUE, width=0.001)
  hexagons <- crop(hexagons, as(bb, "Spatial"))
  hexagons <- st_as_sf(hexagons)
  
  plot1 <- mainGlobalMapEqual +
    
    geom_sf(data = hexagons, aes(fill=value), colour ="black", size = 0.05) +
    scale_colour_gradient2(low = "#450751",mid = "white",high = "#135107",midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar", aesthetics = "fill", limits=c(minLegend,maxLegend) ) + 
    
    geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
  
  pdf(file=gsub(".RData","Equal.pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
  # Class based colors
  hexagons$categories <- as.factor(cut(hexagons$value, breaks = c(-15,-10,-5,-0.1,0.1,4,8,12), labels = c(1,2,3,4,5,6,7)))
  
  library(RColorBrewer)
  myColors <- rev(c("#164B0B","#3D982A","#7FF567","#FFFFFF","#D679E9","#9833AC","#4E0E5B"))
  names(myColors) <- levels(hexagons$categories)

  plot1 <- mainGlobalMapEqual +
    geom_sf(data = hexagons, aes(fill=categories), colour ="black", size = 0.05) +
    scale_color_manual(values = myColors, aesthetics="fill", na.value = "grey50") +
    geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
  
  pdf(file=gsub(".RData","EqualClasses.pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
}

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Extrapolation map

load(paste0(stackResultsFolder,"/RData/","speciesData.RData"))

dataLayersFileType <- "tif"
dataLayers <- c("CoastalExposureW_Pred_Max.tif","AirTemperature Surface Pred LtMax.tif","AirTemperature Surface Pred LtMin.tif","Precipitation Surface Pred Mean.tif","distanceToDelta.tif","Slope.tif")
dataLayersName <- c("CoastalExposure","TempMax","TempMin","PrecipMean","distanceToDelta","Slope") # 

climateLayersDirectory <- "../Data/Climate/Present/"
rasterLayers <- list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
rasterLayers <- processLayers(rasterLayers,occurrenceRecords=speciesData[speciesData$PA==1,2:3],regionBuffer=NULL, minDepth="NULL" , maxDepth="NULL",intertidal="Dependencies/Data/Rasters/BO2CoastLine.tif")
names(rasterLayers) <- dataLayersName

calibrationClimaticRegion <- raster::extract(rasterLayers,speciesData[,2:3])

climateLayersDirectory <- "../Data/Climate/RCP26/"
rasterLayers <- list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
rasterLayers <- processLayers(rasterLayers,occurrenceRecords=speciesData[speciesData$PA==1,2:3],regionBuffer=NULL, minDepth="NULL" , maxDepth="NULL",intertidal="Dependencies/Data/Rasters/BO2CoastLine.tif")
names(rasterLayers) <- dataLayersName
names(rasterLayers)

rasterMap <- subset(rasterLayers,2)
rasterMap[rasterMap <= max(calibrationClimaticRegion[,names(rasterMap)]) ] <- 0
rasterMap[rasterMap > 0] <- 1

rasterMap <- crop(rasterMap,extent(c(min(speciesData[,2]),max(speciesData[,2]),min(speciesData[,3]),max(speciesData[,3]))))

resolutionH3 <- 3
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { max(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)

minLegend <- min(rasterMapDF.polygons$value); minLegend
maxLegend <- max(rasterMapDF.polygons$value); maxLegend

myColors <- c("#C3E6EF","#e31515") # Blue Yellow Red Purple #450751

#----------------------

plot3 <- mainGlobalMap +
  geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.1) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        legend.key.height= unit(0.25, 'cm'),
        legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
plot3

pdf(file=paste0(stackResultsFolder,"/extrapolationMaxTempRCP26.pdf"),width=12,useDingbats=FALSE)
plot3
dev.off()

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Zero centered Climate map

list.files("../Data/Climate/Baseline/")
# Precipitation Surface Pred LtMin.tif

rasterMapP <- raster("../Data/Climate/Baseline/OceanTemperature BenthicMin Ltmax 2010-2020.tif")
rasterMapF <- raster("../Data/Climate/ssp585/OceanTemperature BenthicMin Ltmax 2090-2100.tif")
rasterMap <- rasterMapF - rasterMapP
rasterMap

maskOcc <- loadRData( files[17] )
maskOcc[maskOcc != 1] <- NA
rasterMap <- mask(rasterMap,maskOcc)
rasterMap <- crop(rasterMap,maskOcc)
rasterMap

resolutionH3 <- 3
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)

minLegend <- min(rasterMapDF.polygons$value); minLegend
maxLegend <- max(rasterMapDF.polygons$value); maxLegend

nColors <- 7

colfuncLower <- colorRampPalette(c("#9A0F0F","#FCB46D", "#FFFFFF")) # FFFEC7
colfuncUpper <- colorRampPalette(c("#FFFFFF","#82DC9D","#11702E"))
colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
myColors <- c(colfuncLower(sum(colorBreaks <= 0))[- which(colfuncLower(sum(colorBreaks <= 0)) == "#FFFFFF")],"#FFFFFF",colfuncUpper(sum(colorBreaks >= 0))[-1])

if( length(myColors) > 14 ) {
  myColors <- myColors[  which(sapply( colorBreaks , function(x) f(abs(x)) | x == 0)) ]
  colorBreaks <- colorBreaks[  which(sapply( colorBreaks, function(x) f(abs(x)) | x == 0)) ]
}

plot4 <- mainGlobalMap +
  geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.05) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        legend.key.height= unit(0.25, 'cm'),
        legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
plot4

pdf(file=paste0(stackResultsFolder,"/climateAnomalyMaxTemperaturessp585.pdf"),width=12,useDingbats=FALSE)
plot4
dev.off()

## ------------------------
