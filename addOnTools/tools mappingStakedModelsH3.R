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
rm(list=(ls()))
gc(reset=TRUE)

source("Dependencies/mainFunctions.R")
source("0. config.R")

mapType <- c("Reachable") # "Reachable | unConstrained"
scenariosToMap <- c("Baseline","ssp119","ssp585")

stackResultsFolder <-  "/Volumes/StingRay/Dropbox/Manuscripts/Major restructuring of marine forests diversity under projected climate change/Results/"
stackResultsFolderSpecies <- "/Volumes/StingRay/Dropbox/Data/Distribution Models/Global distribution of marine forests/Results/"

# ---------------------

worldMap <- ne_countries(scale = "medium", returnclass = "sf")

## ------------------------

mainGlobalMap <- ggplot() + geom_sf(data = worldMap, color = "#B0B0B0", fill = "#DFDFDF", size=0.1) + themeMap

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Species by Species maps

if(! dir.exists(paste0(stackResultsFolder,"/Maps/perSpecies"))) { dir.create(paste0(stackResultsFolder,"/Maps/perSpecies"), recursive = T) }

allMapsAll<- list.files( stackResultsFolderSpecies, pattern = "Global.RData", full.names = T, recursive=T )
allMapsNamesAll <- list.files( stackResultsFolderSpecies, pattern = "Global.RData", full.names = F , recursive=T )
allMapsAll <- allMapsAll[grepl(paste0("ensembleReclass",mapType),allMapsAll)] 
allMapsNamesAll <- allMapsNamesAll[grepl(paste0("ensembleReclass",mapType),allMapsNamesAll)]
allMapsAll <- unlist(sapply(scenariosToMap,function(x) { allMapsAll[grepl(x,allMapsAll)] }))
allMapsNamesAll <- unlist(sapply(scenariosToMap,function(x) { allMapsNamesAll[grepl(x,allMapsNamesAll)] }))
allMapsNamesAll <- sapply(1:length(allMapsNamesAll), function(x) {  substr(allMapsNamesAll[x],1,gregexec("/Predictions/",allMapsNamesAll)[[x]][1]-1) })

# which( allMapsNames == "Alaria crispa")
# which( allMapsNames == "Zostera marina")

cl <- makeCluster(nCores)
registerDoParallel(cl)

mapping <- foreach(i = 1:length(allMapsNamesAll), .export=c("unionSpatialPolygons"), .packages = c("sf", "sp", "ggplot2","raster")) %dopar% {
  
  m <- allMapsAll[i]
  name <- allMapsNamesAll[i]
  scenario <- scenariosToMap[which(sapply(scenariosToMap,function(x) { grepl(x,m)}))]
    
  if(length(scenario) > 1) { stop("Error :: 992")}
  
  if( scenario == "Baseline" ) {
    occurrecence <- loadRData(paste0(stackResultsFolderSpecies,"/",name,"/Data/speciesData.RData"))
    occurrecence <- occurrecence[occurrecence$PA == 1,2:3]
  }
  
  rasterMap <- loadRData( m )
  
  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]
  
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  
  if(nrow(rasterMapDF) != 0 ) {  
      
      rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
      rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
      rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
      rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
      
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
        geom_sf(data = rasterMapDF.polygons, fill="#1284B5", colour ="gray", size = 0.05) +
        theme(legend.position="none",
              legend.margin=margin(0,0,0,0),
              legend.key.height= unit(0.25, 'cm'),
              legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
    
      # plot1
      
      pdf(file=paste0(stackResultsFolder,"/Maps/perSpecies/",name," ",scenario,".pdf"),width=12,useDingbats=FALSE)
      print(plot1)
      dev.off()
    
      if( scenario == "Baseline" ) {
        
        plot1 <- plot1 + geom_point(data = occurrecence, aes(x = Lon, y = Lat), size = 0.5, color="Black")
        pdf(file=paste0(stackResultsFolder,"/Maps/perSpecies/",name," ",scenario,"Occurrences.pdf"),width=12,useDingbats=FALSE)
        print(plot1)
        dev.off()
        
      }
  
  }
  
  return(NULL)
  
}

stopCluster(cl); rm(cl) ; closeAllConnections(); gc(reset=TRUE)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Continuous Map 

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T, recursive=T )
filesNames <- list.files( stackResultsFolder, pattern = "RData", full.names = F, recursive=T )
files <- files[grepl("Reachable",files)] # Reachable unConstrained
filesNames <- filesNames[grepl("Reachable",filesNames)] # Reachable unConstrained
files <- files[!grepl("extinction",files)] 
files <- files[!grepl("Changes",files)] 
filesNames <- filesNames[!grepl("extinction",filesNames)]
filesNames <- filesNames[!grepl("Changes",filesNames)]; filesNames

for( m in 1:length(files) ) {
  
  print(files[m])
  rasterMap <- loadRData( files[m] )
  mapName <- gsub(".RData","",filesNames[m])

  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]
  
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
  
  minLegend <- min(rasterMapDF.polygons$value); minLegend
  maxLegend <- max(rasterMapDF.polygons$value); maxLegend
  
  nColors <- 7
  colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
  colorBreaks[1] <- min(rasterMapDF.polygons$value)
  colorBreaks[nColors] <- max(rasterMapDF.polygons$value)
  
  myColors <- c("#E4FAFF","#E8EF15","#ec7a06","#e31515", "#450751") # Blue Yellow Red Purple
  myColors <- colorRampPalette(myColors)(nColors)
  
  #----------------------
  
  plot1 <- mainGlobalMap +
    geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="gray", size = 0.05) + # round signif
    scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
  
  # plot1
  
  pdf(file=gsub(".RData",".pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
}

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Zero centered map

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive=T)
filesNames <- list.files( stackResultsFolder, pattern = "RData", full.names = F , recursive=T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
filesNames <- filesNames[grepl("Reachable",filesNames)] # Reachable unConstrained
files <- files[!grepl("extinction",files)] 
files <- files[!grepl("Changes",files)] 
files <- files[grepl("speciesRichness",files)] 

maskIntertidal <- raster("Dependencies/Data/Rasters/coastLineRes005.tif")
maskDistribution <- calc(stack(sapply(files, function(x) { loadRData( x ) } )), sum)
maskDistribution[maskDistribution > 0] <- 1
maskDistribution[maskDistribution == 0] <- NA

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive=T)
filesNames <- list.files( stackResultsFolder, pattern = "RData", full.names = F , recursive=T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
filesNames <- filesNames[grepl("Reachable",filesNames)] # Reachable unConstrained
files <- files[grepl("Changes",files)] 
filesNames <- filesNames[grepl("Changes",filesNames)]; filesNames

for( m in 1:length(files) ) {

  rasterMap <- loadRData( files[m] )
  mapName <- gsub(".RData","",filesNames[m])
  rasterMap <- raster::mask(rasterMap,maskIntertidal)
  
  rasterMap[intersect(Which(is.na(maskDistribution), cells=TRUE) ,Which(!is.na(rasterMap), cells=TRUE))] <- NA
  
  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  #rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]
  
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
  rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)

  minLegend <- min(rasterMapDF.polygons$value)
  maxLegend <- max(rasterMapDF.polygons$value)
  
  #nColors <- 7
  
  # colfuncLower <- colorRampPalette(c("#9A0F0F","#FCB46D", "#FFFFFF")) # FFFEC7
  # colfuncUpper <- colorRampPalette(c("#FFFFFF","#82DC9D","#11702E"))
  # colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
  # myColors <- c(colfuncLower(sum(colorBreaks <= 0))[- which(colfuncLower(sum(colorBreaks <= 0)) == "#FFFFFF")],"#FFFFFF",colfuncUpper(sum(colorBreaks >= 0))[-1])
  
  # if( length(myColors) > 14 ) {
  #  myColors <- myColors[  which(sapply( colorBreaks , function(x) f(abs(x)) | x == 0)) ]
  #  colorBreaks <- colorBreaks[  which(sapply( colorBreaks, function(x) f(abs(x)) | x == 0)) ]
  # }
  
  # rasterMapDF.polygons[rasterMapDF.polygons$value > 0 & rasterMapDF.polygons$value <= 1,"value"] <- 1
  # rasterMapDF.polygons[rasterMapDF.polygons$value >= -1 & rasterMapDF.polygons$value < 0,"value"] <- -1
  # rasterMapDF.polygons <- rasterMapDF.polygons[rasterMapDF.polygons$value != 0,]
  
  plot2 <- mainGlobalMap +
    
    geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.05) +
    scale_colour_gradient2(low = "#BC0000",mid = "white",high = "#0B6B28",midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar", aesthetics = "fill" ) + 
    
    # geom_sf(data = rasterMapDF.polygons[rasterMapDF.polygons$value < 0,], aes(fill=value), colour ="black", size = 0.05) +
    # scale_colour_gradient(low = "#A40404", high = "#FFFFFF", aesthetics = "fill") +
    # ggnewscale::new_scale_fill() +
    # geom_sf(data = rasterMapDF.polygons[rasterMapDF.polygons$value > 0,], aes(fill=value), colour ="black", size = 0.05) +
    # scale_colour_gradient(low = "#FFFFFF", high = "#097C18", aesthetics = "fill") +
    
    # scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
  # plot2
  
  pdf(file=gsub(".RData",".pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot2)
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
