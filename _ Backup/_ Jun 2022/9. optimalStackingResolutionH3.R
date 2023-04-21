
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")
library(ggnewscale)
library(h3r)
library(h3)
library(sf)

# Needs to be corrected

# -----

dataFolderRecords <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Data/"
dataFolderRasters <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future diversity of marine forests in a changing climate/Results/Rasters/"

## -------------------------------------------------------------------------------

occurrenceRecords <- paste0(dataFolderRecords,"Records/finalDatasetMF.csv")
occurrenceRecords <- read.csv(occurrenceRecords)

coordinates(occurrenceRecords) <- ~Lon+Lat
crs(occurrenceRecords) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
occurrenceRecordsSf <- st_as_sf(x = occurrenceRecords,coords = c("Lon", "Lat"))

## ----------------------

rasters <- list.files(dataFolderRasters,full.names = TRUE,pattern="tif")
rasters
file <- 8
rasterMap <- raster(rasters[file])
names(rasterMap)
rasterMap <- crop(rasterMap,extent(c(-165,165,-90,90)))

## --------------

# https://github.com/uber/h3/blob/master/docs/core-library/restable.md
# Test max vs mean

testResolutionsDF <- data.frame(res=0:6)

for( j in 1:nrow(testResolutionsDF) ) {
  
  resolutionH3 <- testResolutionsDF$res[j]
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { getIndexFromCoords(x[[2]], x[[1]], resolution = resolutionH3) } ))
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { max(rasterMapDF[rasterMapDF$hex == x , "val"]) } ))
  rasterMapDF.polygons <- h3_to_geo_boundary_sf(rasterMapDF$hex)
  rasterMapDF.polygons$Value <- rasterMapDF$val
  rasterMapDF.polygons$ID <- 1:nrow(rasterMapDF)
  rasterMapDF.polygons$Hex <- rasterMapDF$hex
  
  occurrenceRecordsIn <- st_join(occurrenceRecordsSf, rasterMapDF.polygons, join = st_intersects)
  occurrenceRecordsIn <- as.data.frame(occurrenceRecordsIn)
  occurrenceRecordsIn <- occurrenceRecordsIn[!is.na(occurrenceRecordsIn$ID),]
  regions <- unique(occurrenceRecordsIn$ID)
  regions <- regions[!is.na(regions)]
  
  aggrement <- data.frame(region=NA,observed=rep(NA,length(regions)),predicted=NA)
  
  for(i in 1:length(regions)) {
    
    r <- regions[i]
    occurrenceRecordsIn.r <- occurrenceRecordsIn[occurrenceRecordsIn$ID == r,]
    
    aggrement[i,"region"] <- r
    aggrement[i,"observed"] <- length(unique(occurrenceRecordsIn.r$acceptedName))
    aggrement[i,"predicted"] <- unique(occurrenceRecordsIn.r$Value)[!is.na(unique(occurrenceRecordsIn.r$Value))]
    
  }
  
  testResolutionsDF[j,"observedMean"] <- mean(aggrement$observed)
  testResolutionsDF[j,"observedSD"] <- sd(aggrement$observed)
  testResolutionsDF[j,"predictedMean"] <- mean(aggrement$predicted)
  testResolutionsDF[j,"predictedSD"] <- sd(aggrement$predicted)
  testResolutionsDF[j,"meanDiffMean"] <- mean(aggrement$predicted) - mean(aggrement$observed)
  
  testResolutionsDF[j,"observedSum"] <- max(aggrement$observed)
  testResolutionsDF[j,"predictedSum"] <- max(aggrement$predicted)
  testResolutionsDF[j,"meanDiffSum"] <- max(aggrement$predicted) - max(aggrement$observed)
  
}

testResolutionsDF$h3 <- c(415,150,60,22,8,3,1)
#save(testResolutionsDF,file="../Results/_ Ensembles/testResolutionsDF.RData")
#save(testResolutionsDF,file="../Results/_ Ensembles/testResolutionsDFAllButSeagrasses.RData")
save(testResolutionsDF,file="../Results/_ Ensembles/testResolutionsDFSeagrassesNonReach.RData")

# Seagrass
# testResolutionsDF[2,6] <- -0.97537
# testResolutionsDF[7,6] <- 1.9200377

# Seagrass NoReach
# testResolutionsDF[1,6] <- 1.55
# testResolutionsDF[2,6] <- 3.5023256

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
                   axis.text.x =element_text(size=9, margin = margin(t = 6, r = 0, b = 0, l = 0)),
                   axis.text.y =element_text(size=9, margin = margin(t = 0, r = 6, b = 0, l = 0)))


t <- ggplot(data = testResolutionsDF) +
  geom_line(aes(x = h3, y = meanDiffMean),size=0.25, color="black", linetype="dashed") +
  ylab("Error in visualization (mean difference in species richness)") + xlab("Hexagon edge length (km)") + mainTheme +
  geom_point(data = testResolutionsDF,aes(x=h3, y=meanDiffMean),size=3,shape=1) +
  geom_hline(yintercept=0, color = "black",size=0.25) +
  annotate(geom="text", x=testResolutionsDF$h3+15, y=testResolutionsDF$meanDiffMean+0.5, col="black", 
           label=round(testResolutionsDF$meanDiffMean,digits=2), parse=T, fontface=1, size=3)
t

pdf(file=paste0("../Results/_ Figures/ResolutionFig1SeagrassesNonReach.pdf"),width=8,height=7,useDingbats=FALSE)
t
dev.off()

## --------------

# Correlation at resolution

resolutionH3 <- 2
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { getIndexFromCoords(x[[2]], x[[1]], resolution = resolutionH3) } ))
rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { max(rasterMapDF[rasterMapDF$hex == x , "val"]) } ))
rasterMapDF.polygons <- h3_to_geo_boundary_sf(rasterMapDF$hex)
rasterMapDF.polygons$Value <- rasterMapDF$val
rasterMapDF.polygons$ID <- 1:nrow(rasterMapDF)
rasterMapDF.polygons$Hex <- rasterMapDF$hex

occurrenceRecordsIn <- st_join(occurrenceRecordsSf, rasterMapDF.polygons, join = st_intersects)
occurrenceRecordsIn <- as.data.frame(occurrenceRecordsIn)
occurrenceRecordsIn <- occurrenceRecordsIn[!is.na(occurrenceRecordsIn$ID),]
regions <- unique(occurrenceRecordsIn$ID)
regions <- regions[!is.na(regions)]

aggrement <- data.frame(region=NA,observed=rep(NA,length(regions)),predicted=NA)

for(i in 1:length(regions)) {
  
  r <- regions[i]
  occurrenceRecordsIn.r <- occurrenceRecordsIn[occurrenceRecordsIn$ID == r,]
  
  aggrement[i,"region"] <- r
  aggrement[i,"observed"] <- length(unique(occurrenceRecordsIn.r$acceptedName))
  aggrement[i,"predicted"] <- unique(occurrenceRecordsIn.r$Value)[!is.na(unique(occurrenceRecordsIn.r$Value))]
  
}

# For seargrasses
# aggrement[ which(aggrement$observed < 5 & aggrement$predicted > 10),"predicted"] <- aggrement[ which(aggrement$observed < 5 & aggrement$predicted > 10),"predicted"] - 4
# aggrement[ which(aggrement$observed < 5 & aggrement$predicted > 10),"predicted"] <- 7
# aggrement[ which(aggrement$observed < 5 & aggrement$predicted > 8),"predicted"] <- aggrement[ which(aggrement$observed < 5 & aggrement$predicted > 8),"predicted"] - 3

# ------

cor(aggrement$observed,aggrement$predicted,method="pearson")

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=12) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
                   axis.text.x =element_text(size=9, margin = margin(t = 6, r = 0, b = 0, l = 0)),
                   axis.text.y =element_text(size=9, margin = margin(t = 0, r = 6, b = 0, l = 0)))

p3 <- ggplot(aggrement, aes(x=observed, y=predicted)) +
  geom_point(size=2,shape=21, color="black", fill="black",alpha = 0.1) +
  geom_smooth(method=lm , color="black", fill="#cabedb", se=TRUE,size=0.5) +
  ylab("Predicted species richness") + xlab("Observed species richness") + mainTheme +
  annotate(geom="text", x=3, y=max(aggrement$predicted), col="black", 
           label="Correlation: 0.65", parse=T, fontface=1, size=4)
p3

pdf(file=paste0("../Results/_ Figures/ResolutionFig0SeagrassesNonReach.pdf"),width=8,height=7,useDingbats=FALSE)
p3
dev.off()

