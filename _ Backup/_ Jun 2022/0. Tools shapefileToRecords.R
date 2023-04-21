shape <- listAllFiles(dataLayersDirectory,dataLayersFileType)[1]
shape <- raster(shape)

library(fasterize)

occurrenceRecords <- shapefile("../Data/Records/SimplifiedSingle2")
occurrenceRecordsSf <- st_as_sf(occurrenceRecords)
occurrenceRecordsRaster <- fasterize( occurrenceRecordsSf, shape, field = NULL, background = NA_real_)
occurrenceRecordsRasterPts <- Which(occurrenceRecordsRaster==1,cells=TRUE)
occurrenceRecordsRasterPts <- xyFromCell(occurrenceRecordsRaster,occurrenceRecordsRasterPts)

occurrenceRecords <- occurrenceRecordsRasterPts
occurrenceRecords <- data.frame(occurrenceRecords)
colnames(occurrenceRecords) <- c("Lon","Lat")
save(occurrenceRecords,file=paste0( resultsDirectory, "/RData/occurrenceRecordsInit.RData"))
