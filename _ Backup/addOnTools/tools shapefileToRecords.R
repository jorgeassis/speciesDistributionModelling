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

# setwd("/media/Jellyfish/Backup [Ongoing studies]/Future diversity of cold water corals in a changing climate/speciesDistributionModelling")

source("Dependencies/mainFunctions.R")

## -------------

mainResultsDirectory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/_ Under Revision/Mangrove productivity losses under contrasting scenarios of future climate change/"

climateLayersDirectory <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers BO3.0/Present/Benthic/LongTerm/"
dataLayersFileType <- "tif"

shape <- list.files(climateLayersDirectory,pattern=dataLayersFileType, full.names = TRUE, recursive=TRUE)[1]
shape <- raster(shape)
shape <- disaggregate(shape, fact=2)

occurrenceRecords <- shapefile(paste0(mainResultsDirectory,"/Data/Records/SimplifiedSingle2.shp"))
occurrenceRecordsSf <- st_as_sf(occurrenceRecords)
occurrenceRecordsRaster <- fasterize( occurrenceRecordsSf, shape, field = NULL, background = NA_real_)
occurrenceRecordsRasterPts <- Which(occurrenceRecordsRaster==1,cells=TRUE)
occurrenceRecordsRasterPts <- xyFromCell(occurrenceRecordsRaster,occurrenceRecordsRasterPts)

occurrenceRecords <- data.frame(speciesName="Mangrove",occurrenceRecordsRasterPts)
colnames(occurrenceRecords) <- c("speciesName","Lon","Lat")
write.csv(occurrenceRecords,file=paste0( mainResultsDirectory, "/Data/Records/occurrenceFinalVersion.csv"), row.names = FALSE)

## -------------
## -------------