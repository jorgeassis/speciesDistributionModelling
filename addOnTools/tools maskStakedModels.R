## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")

nCores <- 40

tempFolder <- "~/Projects/SDM/tempFolder/"

modelsDirectory <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Results [Full models]/_ Per Species/"
resultsFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future diversity of marine forests in a changing climate/Results [Kelp]"

## ----------------

resultsFolder <- paste0(resultsFolder,"/Rasters/")
list.files(resultsFolder)

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
# Mask with bathymetry

filesToMask <- list.files(resultsFolder, full.names = T, pattern = "tif")
filesToMask <- filesToMask[!grepl("MaskedDepth",filesToMask)]
filesToMask

bathymetry <- raster("Dependencies/Data/Rasters/BO2BathymetryDepthMean.tif")
bathymetry <- crop( bathymetry , raster(filesToMask[1]))

thresholdMax <- 0
thresholdMin <- -30

bathymetry[bathymetry > thresholdMax] <- NA
bathymetry[bathymetry < thresholdMin] <- NA

# sum(getValues(( !  is.na(bathymetry) ) * raster::area(bathymetry)), na.rm=T)

for(f in filesToMask) {
  
  file <- raster(f)
  file <- mask(file,bathymetry)
  writeRaster(file,filename=gsub("\\.tif","MaskedDepth.tif",f),format="GTiff",overwrite=T)
  
}

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
# Mask ensembles with polygon [remove unwanted regions]

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

filesToMask <- list.files(resultsFolder, full.names = T, pattern = "tif")
filesToMask

mask.r <- calc(!is.na(stack(filesToMask)),sum)
mask <- drawPolygon2(mask.r)

# ---------------------------------------
# ---------------------------------------
