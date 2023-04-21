## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)
source("Dependencies/mainFunctions.R")

nCores <- 10

tempFolder <- "/Volumes/Jellyfish/Temp"
modelsDirectory <- "/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Results [Full models]/_ Per Species/"
resultsFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Future diversity of marine forests in a changing climate/Results [Kelp]"

## ----------------

resultsFolder <- paste0(resultsFolder,"/Rasters/")
if(! dir.exists(paste0(resultsFolder))) { dir.create(paste0(resultsFolder), recursive = T) }

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

# Subset species

occRecords <- read.csv("/Volumes/Jellyfish/Dropbox/Data/Distribution Models/Global distribution of marine forests/Data/Records/finalDatasetMF.csv")
unique(occRecords$order)

# !
unique(occRecords[which(is.na(occRecords$order)),"name"])

spVect.1 <- unique(occRecords$name)[grepl("Desmarestia",unique(occRecords$name)) ]
spVect.2 <- unique(occRecords[occRecords$order %in% c("Laminariales","Tilopteridales"),"name"])
spVect <- unique(c(spVect.1,spVect.2))
length(spVect)

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

finalExtent <- c(-180,180,-90,90)

## ---------------------------------
## Diversity

for(reach in c("Reachable")) { # c("","Reachable")
  
  for(period in c("Present","RCP26","RCP85")) { # ,"RCP26","RCP45","RCP60","RCP85"
    
    print(paste(period,"::",reach))
    rasterOptions(todisk = TRUE)
    rasterOptions(tmpdir=tempFolder)
    
    sFiles <- list.files(modelsDirectory,pattern=paste0("Ensemble","Reclass",reach,"_",period,".tif"),full.names = TRUE,recursive=TRUE)
    # sFiles <- sFiles[grepl("/Maps/Global/",sFiles)]
    # sFiles <- sFiles[!grepl("/_ Backup/",sFiles)]
    sFiles <- sFiles[grepl("Global/EnsembleReclass",sFiles)]
    
    sFiles <- sFiles[unlist(sapply(spVect, function(x) { which(grepl(x,sFiles)) } ))]
    
    diversity <- raster(sFiles[1])
    
    if( ! all.equal(finalExtent,c(-180,180,-90,90)) ) {
      diversity <- crop(diversity,extent(finalExtent))
    }
    
    for( i in 2:length(sFiles)) {
      
      diversity.i <- raster(sFiles[i])
      
      if( ! all.equal(finalExtent,c(-180,180,-90,90)) ) {
        diversity.i <- crop(diversity.i,extent(finalExtent))
      }
      
      diversity <- diversity + diversity.i
      
    }
    
    ## --------
    
    if( period != "Present" & reach == "" ) {
      
      diversityP <- diversity
      diversityP[diversityP == 0] <- NA
      writeRaster(diversityP,filename=paste0(resultsFolder,"/speciesRichness",period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
      
      diversity <- diversity - diversityPresent
      diversity[diversity == 0] <- NA
      writeRaster(diversity,filename=paste0(resultsFolder,"/speciesRichnessChange",period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
      
    }
    
    ## --------
    
    if( period == "Present" ) { 
      
      diversityPresent <- diversity  
      diversity[diversity == 0] <- NA
      writeRaster(diversity,filename=paste0(resultsFolder,"/speciesRichness",period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
      
    }
    
  }
  
  ## --------
  
  gc()
  file.remove(list.files(tempFolder, full.names = T))
  
}

## ---------------------------------
## Gain and Loss and Turnover

reach <- "Reachable" # "Reachable" ""

for(period in c("RCP26","RCP85")) { # ,"RCP26","RCP45","RCP60","RCP85"
  
  for(trait in c("Gain","Loss","Refugia")) {
    
    print(paste(period,"::",reach))
    rasterOptions(todisk = TRUE)
    rasterOptions(tmpdir=tempFolder)
    
    sFiles <- list.files(modelsDirectory,pattern=paste0("EnsembleReclass",reach,"_",period,".tif"),full.names = TRUE,recursive=TRUE)
    sFiles <- sFiles[!grepl("GainLoss",sFiles)]
    sFiles <- sFiles[grepl(paste0("/",trait),sFiles)]

    sFiles <- sFiles[unlist(sapply(spVect, function(x) { which(grepl(x,sFiles)) } ))]
    
    metric <- raster(sFiles[1])
    
    if( ! all.equal(finalExtent,c(-180,180,-90,90)) ) {
      metric <- crop(metric,extent(finalExtent))
    }
    
    for( i in 2:length(sFiles)) {
      
      metric.i <- raster(sFiles[i])
      
      if( ! all.equal(finalExtent,c(-180,180,-90,90)) ) {
        metric.i <- crop(metric.i,extent(finalExtent))
      }
      
      metric <- metric + metric.i
      
    }
    
    if(trait  == "Gain") { metricGain <- metric }
    if(trait  == "Loss") { metricLoss <- metric }
    if(trait  == "Refugia") { metricRefugia <- metric }
    
    metric[metric == 0] <- NA
    writeRaster(metric,filename=paste0(resultsFolder,"/range",trait,period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
    
  }
  
  ## --------
  ## Turnover rasters
  
  # 1 - (refugia / RichnessFuture)
  # 0 means all species persist and 1 all species are exchanged
  
  speciesRichnessPresent <- raster(paste0(resultsFolder,"/speciesRichness","Present",ifelse(reach=="","","Reachable"),".tif"))
  speciesRichnessPresent[is.na(speciesRichnessPresent)] <- 0
  speciesRichnessPeriod <- raster(paste0(resultsFolder,"/speciesRichness",period,ifelse(reach=="","","Reachable"),".tif"))
  speciesRichnessPeriod[is.na(speciesRichnessPeriod)] <- 0
  
  metric <- 1 - (metricRefugia / speciesRichnessPeriod)
  writeRaster(metric,filename=paste0(resultsFolder,"/speciesExchangeRatio",period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
  
  gc()
  file.remove(list.files(tempFolder, full.names = T))  
  
}

## ---------------------------------
## Regions of extinction

reach <- "Reachable"

for(period in c("RCP26","RCP85")) { # ,"RCP26","RCP45","RCP60","RCP85"
  
  print(paste(period,"::",reach))
  rasterOptions(todisk = TRUE)
  rasterOptions(tmpdir=tempFolder)
  
  speciesRichnessPresent <- raster(paste0(resultsFolder,"/speciesRichness","Present",ifelse(reach=="","","Reachable"),".tif"))
  speciesRichnessPresent[is.na(speciesRichnessPresent)] <- 0
  speciesRichnessPeriod <- raster(paste0(resultsFolder,"/speciesRichness",period,ifelse(reach=="","","Reachable"),".tif"))
  speciesRichnessPeriod[is.na(speciesRichnessPeriod)] <- 0
  
  ## --------
  
  f1L <- function(x) { ifelse( x[[1]] > 0 & x[[2]] == 0 , 1 , 0)  }
  extinctionRegions <- stack(speciesRichnessPresent,speciesRichnessPeriod)
  extinctionRegions <- calc(extinctionRegions, fun=f1L)
  
  ## --------
  
  writeRaster(extinctionRegions,filename=paste0(resultsFolder,"/extinctionAreas",period,ifelse(reach=="","","Reachable"),".tif"),format="GTiff",overwrite=T)
  
  gc()
  file.remove(list.files(tempFolder, full.names = T))
  
}

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------