## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##  
##  Centre of Marine Sciences [ ccmar.ualg.pt ]
##  Faro, Portugal
##
## ---------------------------------------------------------------------
##
##  Machine Learning Species Distribution Modelling [ Ver.301 ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

setwd("/Volumes/Jellyfish/Dropbox/theMarineDataScientist/gitRepositories/speciesDistributionModelling")

source("Dependencies/mainFunctions.R")
source("0. config.R")

# ---------------------

mainResultsDirectory <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/"
dataLayersDirectory <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Data/Climate/"
dataLayersFileType <- "tif"

# -------------------

dataRecords <- "/Volumes/Jellyfish2/Backup [Ongoing studies]/Future diversity of cold water corals in a changing climate/Data/Occurrence/biodiversityDataBaseMergedFlaggedPrunned.csv"

dataLayers <- c("DissolvedMolecularOxygen BenthicMean Mean 2010-2020.tif",
                "OceanTemperature BenthicMean Ltmax 2010-2020.tif",
                "OceanTemperature BenthicMean Ltmin 2010-2020.tif",
                "pH BenthicMean Mean 2010-2020.tif",
                "Salinity BenthicMean Mean 2010-2020.tif",
                "SeaWaterSpeed BenthicMean Mean 2010-2020.tif",
                "Slope BenthicMean.tif",
                "Terrain Ruggedness Index BenthicMean.tif",
                "TotalPrimaryProductionPhyto Surface Mean 2010-2020.tif")

dataLayersName <- c("Oxygen","TemperatureMax","TemperatureMin","pH","Salinity","WaterVelocity","Slope","Terrain","PrimaryProduction")
dataLayersMonotonocity <- c(+1,-1,+1,-1,+1,+1,+1,+1,+1)
dataLayersMonotonocity <- setNames(data.frame(t(dataLayersMonotonocity)), dataLayersName)
dataLayersMonotonocity

# -------------------

bathymetryDataLayer <- "/Volumes/Jellyfish2/Backup [Ongoing studies]/Future diversity of cold water corals in a changing climate/Data/Bathymetry Mean BenthicMean.tif"
coastLineDataLayer <-  "Dependencies/Data/Rasters/BO2CoastLine.tif"
intertidal <- FALSE

depthTraits <- "/Volumes/Jellyfish2/Backup [Ongoing studies]/Future diversity of cold water corals in a changing climate/Data/Traits/coldWaterCorals.csv" # NULL

minDepth <- 0
maxDepth <- 0

minDepthBuffer <- 100
maxDepthBuffer <- 100

## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

worldMap <- ne_countries(scale = "medium", returnclass = "sf")

dataRecords <- read.csv(dataRecords)
speciesToModel <- unique(dataRecords$acceptedName)

for(species in speciesToModel) {
  
  # ---------------------
  
  cat("\014")
  cat("# --------------------------------- \n")
  cat(species," | ",which(speciesToModel == species),"out of",length(speciesToModel),"\n")
  cat("# --------------------------------- \n")
  cat("\n")
  
  # ---------------------
  
  resultsDirectory <- paste0(mainResultsDirectory,"/",species,"/")
  
  if( dir.exists(resultsDirectory) ) { next }
  
  # ---------------------
  
  if( ! dir.exists(resultsDirectory) ) { dir.create(resultsDirectory, recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/RData/")) ) { dir.create(paste0(resultsDirectory,"/RData/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/summaryModel/")) ) { dir.create(paste0(resultsDirectory,"/summaryModel/"), recursive = TRUE) }
  if( ! dir.exists(paste0(resultsDirectory,"/Rasters/")) ) { dir.create(paste0(resultsDirectory,"/Rasters/"), recursive = TRUE) }
  
  # ---------------------
  
  occurrenceRecords <- dataRecords[which(dataRecords$acceptedName == species),]
  occurrenceRecords <- data.frame(Lon=occurrenceRecords$decimalLongitude,Lat=occurrenceRecords$decimalLatitude)
  
  # plotMap(occurrenceRecords,4,"black")
  
  # ---------------------
  
  if( !is.null(depthTraits) ) {
    
    depthTraitsData <- read.csv(depthTraits)
    depthTraitsData <- depthTraitsData[depthTraitsData$Species == species,]
    minDepth <- depthTraitsData$minDepth
    maxDepth <- depthTraitsData$maxDepth

  }
  
  # ---------------------
  # Environmental layers
  
  rasterLayers <- list.files(dataLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
  rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
  rasterLayers <- processLayers(rasterLayers,occurrenceRecords,regionBuffer, minDepth=ifelse(minDepth-minDepthBuffer < 0 , 0 , minDepth-minDepthBuffer) , maxDepth=maxDepth+maxDepthBuffer,intertidal)
  
  names(rasterLayers) <- dataLayersName
  
  rasterLayers <- dropNoVariationLayers(rasterLayers)
  rasterLayers <- correctLayer(rasterLayers,"Salinity","he",28,28)
  
  # ---------------------
  
  source("1. prepareRecords.R")
 
  # ---------------------
  
  source("2. fullModel.R")
  
  # ---------------------
  
}

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# End of Code