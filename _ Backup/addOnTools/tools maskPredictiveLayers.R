## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")

nCores <- 12

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

species <- "Zostera marina"

r1Files <- list.files(paste0(mainResultsDirectory,"/"),pattern = "ensembleReclassReachableGlobal", full.names = T, recursive = TRUE)
r1Files <- r1Files[grepl(species,r1Files)]
r1Files

r1 <- loadRData(r1Files[1])
plot(r1,col=c("white","black"))
regionMask <- drawPolygon(r1)
cells <- unlist(cellFromPolygon(r1,regionMask))

# 0. Decide query value
cells <- cells[ which(r1[cells] == 1) ]

# 1. Substitute over the same raster
r1[cells] <- 0
plot(r1)
save(r1,filename="")

# 2. Substitute over other rasters
r1Files
file.r <- 4
r2 <- loadRData(r1Files[file.r])
r2[cells] <- 1
save(r2,file=r1Files[file.r])

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
# Mask with bathymetry

filesToMask <- list.files(resultsFolder, full.names = T, pattern = "tif")

bathymetry <- raster("/Volumes/Jellyfish/Dropbox/Data/Spatial information/Rasters/BathymetryDepthMean.tif")
bathymetry <- crop( bathymetry , raster(filesToMask[1]))

thresholdMax <- 0
thresholdMin <- -30

bathymetry[bathymetry > thresholdMax] <- NA
bathymetry[bathymetry < thresholdMin] <- NA

for(f in filesToMask) {
  
  file <- raster(f)
  file <- mask(file,bathymetry)
  writeRaster(file,filename=f,format="GTiff",overwrite=T)
  
}

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

occurrenceRecords <- read.table(gsub(".csv"," Processed.csv",dataRecordsFile),sep=";",header=T)

resultsFolder <- resultsDirectory # "/media/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Results/SDM North"
r1Files <- list.files(resultsFolder,pattern = "EnsembleReclass", full.names = T)
r1Files <- r1Files[grepl(".Masked.",r1Files)]
r1Files

# file.i <- 4

r1Files[file.i]

r1 <- raster(r1Files[file.i])
plot(r1,col=c("white","black"))
removeMask <- drawPolygon(r1,occurrenceRecords,"rcp852") # maskPacific rcp85 rcp852

## --------------

cells <- unlist(cellFromPolygon(r1,removeMask))
cells <- cells[ which(r1[cells] == 1) ]
r1[cells] <- 0
plot(r1)

writeRaster(r1,filename=gsub(".tif",".Masked.tif",r1Files[file.i]),format="GTiff",overwrite=T)

## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
## Multiple

maskPrediction <- drawPolygon(predictLayers,occurrenceRecords,"predictMaskSurface")
maskPrediction <- shapefile("/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Spatial/maskPresent.shp")

# ------------------

files <- list.files("/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Results/SDM/South", full.names = T,pattern = "Present") 
filesProb <- files[!grepl(pattern = "Reclass",files)]
filesReclass <- files[grepl(pattern = "Reclass",files)]
filesProb <- filesProb[!grepl(pattern = "Surface",filesProb)]
filesReclass <- filesReclass[!grepl(pattern = "Surface",filesReclass)]
filesProb <- filesProb[!grepl(pattern = "SD[.]",filesProb)]
filesProb <- filesProb[!grepl(pattern = "SDMasked",filesProb)]
filesReclass <- filesReclass[!grepl(pattern = "SD[.]",filesReclass)]

cells <- unlist(cellFromPolygon(raster(filesReclass[1]),maskPrediction))
cells <- cells[ which(raster(filesReclass[1])[cells] == 1) ]

for(f in 1:length(filesProb)) {
  
  r <- raster(filesProb[f])
  r[cells] <- sample( seq(0.01,0.023,length.out = length(cells)) , length(cells))
  writeRaster(r,filename=filesProb[f],format="GTiff",overwrite=T)
  
  r <- raster(filesReclass[f])
  r[cells] <- 0
  writeRaster(r,filename=filesReclass[f],format="GTiff",overwrite=T)
  
}

# ------------------

files <- list.files("/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Results/SDM/North", full.names = T,pattern = "Present") 
files <- files[!grepl(pattern = "Masked",files)]
files <- files[!grepl(pattern = "Surface",files)]

maskPrediction <- shapefile("/Volumes/Jellyfish/Dropbox/Manuscripts/The Phylogeography of Macrocystis pyrifera/Data/Spatial/maskDistribution.shp")

for(f in 1:length(files)) {
  
  r <- raster(files[f])
  r <- mask(r,maskPrediction)
  writeRaster(r,filename=gsub(".tif","SpeciesDistribution.tif",files[f]),format="GTiff",overwrite=T)
  
}

# ------------------
