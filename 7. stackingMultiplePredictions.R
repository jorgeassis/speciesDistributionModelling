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

## ---------------------------------
## Diversity

for(scenario in scenariosToPredict ) {
  
  # rasterOptions(todisk = TRUE)
  # rasterOptions(tmpdir=tempFolder)
  
  sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)

  sFiles <- sFiles[grepl("Predictions",sFiles)]
  sFiles <- sFiles[grepl("ensemble",sFiles)]
  sFiles <- sFiles[grepl("Reclass",sFiles)]
  
  if(typePrediction == "Reachable") { sFiles <- sFiles[grepl("Reachable",sFiles)] }
  if(typePrediction == "unConstrained") { sFiles <- sFiles[!grepl("Reachable",sFiles)] }
  
  sFiles <- sFiles[ grepl(scenario,sFiles) ]
  sFiles <- sFiles[ grepl(paste0("Global.RData"),sFiles) ]
  
  sFiles <- sFiles[!grepl("Gain",sFiles)]
  sFiles <- sFiles[!grepl("Loss",sFiles)]
  sFiles <- sFiles[!grepl("Refugia",sFiles)]
  sFiles <- sFiles[!grepl("RangeShifts",sFiles)]
  sFiles <- sFiles[!grepl("ensembleSD",sFiles)]

  sFiles <- sFiles[unlist(sapply(speciesPredicted,function(x) { agrep(x, sFiles) } ))]
  sFiles <- unique(sFiles)
  
  if(length(sFiles) > length(speciesPredicted)) { stop("Error :: 333")}
  
  diversity <- loadRData(sFiles[1])
  
  for( i in 2:length(sFiles)) {
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(scenario," || ",i,"out of",length(sFiles),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    diversity <- diversity + loadRData(sFiles[i])
    
  }
  
  ## --------
  
  if( scenario != "Baseline" ) {
    
    save( diversity , file=paste0(stackResultsFolder,"/Maps/","/speciesRichness",scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)
    
    diversity <- diversity - diversityBaseline
    save( diversity , file=paste0(stackResultsFolder,"/Maps/","/speciesRichnessChange",scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)
    
  }
  
  ## --------
  
  if( scenario == "Baseline" ) { 
    
    diversityBaseline <- diversity  
    save( diversity , file=paste0(stackResultsFolder,"/Maps/","/speciesRichness",scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)
    
  }
  
}

## ---------------------------------
## Gain, Loss and Turnover

for(scenario in scenariosToPredict[scenariosToPredict != "Baseline"] ) {
  
  for(trait in c("Gain","Loss","Refugia")) {

    sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)

    sFiles <- sFiles[grepl("Predictions",sFiles)]
    sFiles <- sFiles[grepl("ensemble",sFiles)]
    sFiles <- sFiles[grepl("Reclass",sFiles)]
    sFiles <- sFiles[!grepl("ensembleSD",sFiles)]
    
    if(typePrediction == "Reachable") { sFiles <- sFiles[grepl("Reachable",sFiles)] }
    if(typePrediction == "unConstrained") { sFiles <- sFiles[!grepl("Reachable",sFiles)] }
    
    sFiles <- sFiles[grepl(paste0(trait,"Reclass"),sFiles)]
    sFiles <- sFiles[grepl(scenario,sFiles)]
    sFiles <- sFiles[grepl("Global.RData",sFiles)]

    sFiles <- sFiles[unlist(sapply(speciesPredicted,function(x) { agrep(x, sFiles) } ))]
    sFiles <- unique(sFiles)
    
    if(length(sFiles) > length(speciesPredicted)) { stop("Error :: 334")}
    
    metric <- loadRData(sFiles[1])
    
    for( i in 2:length(sFiles)) {
      
      cat("\014")
      cat("## --------------------------------- \n")
      cat(trait," || ",scenario," || ",typePrediction," || ",i,"out of",length(sFiles),"\n")
      cat("## --------------------------------- \n")
      cat("\n")
      
      metric <- metric + loadRData(sFiles[i])
      
    }
    
    if( trait  == "Gain" ) { metricGain <- metric }
    if( trait  == "Loss" ) { metricLoss <- metric }
    if( trait  == "Refugia" ) { metricRefugia <- metric }
    
    save( metric , file=paste0(stackResultsFolder,"/Maps/","/range",trait,scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)
    
    ## --------
    # Standardize per richness
    # e.g. Proportion of extiniton
    # (Thuiller et al. 2005) # Thuiller W, Lavorel S, Araújo MB, Sykes MT, Prentice IC (2005) Climate change threats to plant diversityin Europe. Proc Natl Acad Sci USA 102:8245–8250. https://doi.org/10.1073/pnas.0409902102 
    
    richnessBaseline <- loadRData(paste0(stackResultsFolder,"/Maps/","/speciesRichness","Baseline",typePrediction,".RData"))
    metric.std <- metric / richnessBaseline
    metric.std[metric.std == Inf] <- NA
    metric.std[metric.std == -Inf] <- NA
    
    if( trait  != "Gain" ) { metric.std[metric.std > 1] <- 1 }
    
    save( metric.std , file=paste0(stackResultsFolder,"/Maps/","/range",trait,scenario,typePrediction,"StandPerBaselineRichness.RData"), compress=TRUE, compression_level=6)
    
  }
  
  ## --------
  ## Turnover rasters
  
  # 1 - (refugia / RichnessFuture)
  # 0 means all species persist and 1 all species are exchanged
  # metric <- 1 - (metricRefugia / speciesRichnessscenario)
  
  # ((G + L)/(S + G))
  # (Thuiller et al. 2005) # Thuiller W, Lavorel S, Araújo MB, Sykes MT, Prentice IC (2005) Climate change threats to plant diversityin Europe. Proc Natl Acad Sci USA 102:8245–8250. https://doi.org/10.1073/pnas.0409902102 
  
  speciesRichnessscenario <- loadRData( file=paste0(stackResultsFolder,"/Maps/","/speciesRichness",scenario,typePrediction,".RData"))
  speciesRichnesBaseline <- loadRData( file=paste0(stackResultsFolder,"/Maps/","/speciesRichness","Baseline",typePrediction,".RData"))
  
  metric <- 1 - ( metricRefugia / speciesRichnessscenario) 
  metric[metric > 1] <- 1
  metric[metric < 0] <- 0
  
  save( metric , file=paste0(stackResultsFolder,"/Maps/","/speciesExchangeRatio1",scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)
  
  metric <- (metricGain + metricLoss) / (speciesRichnesBaseline + metricGain)
  metric[metric > 1] <- 1
  metric[metric < 0] <- 0
  
  save( metric , file=paste0(stackResultsFolder,"/Maps/","/speciesExchangeRatio2",scenario,typePrediction,".RData"), compress=TRUE, compression_level=6)

}

## ---------------------------------
## Uncertainty

for(scenario in scenariosToPredict ) {

  sFiles <- list.files(mainResultsDirectory,pattern="RData",full.names = TRUE,recursive=TRUE)
  
  sFiles <- sFiles[ grepl("Predictions",sFiles) ]
  sFiles <- sFiles[ grepl("ensembleSDGlobal",sFiles) ]
  sFiles <- sFiles[ grepl(scenario,sFiles) ]
  
  sFiles <- sFiles[unlist(sapply(speciesPredicted,function(x) { agrep(x, sFiles) } ))]
  sFiles <- unique(sFiles)
  
  uncertainty <- loadRData(sFiles[1])
  
  for( i in 2:length(sFiles)) {
    
    cat("\014")
    cat("## --------------------------------- \n")
    cat(scenario," || ",i,"out of",length(sFiles),"\n")
    cat("## --------------------------------- \n")
    cat("\n")
    
    uncertainty <- uncertainty + loadRData(sFiles[i])
    
  }
  
  uncertainty <- uncertainty / length(sFiles)
  uncertainty[uncertainty == 0] <- NA
  
  ## --------
  
  save( uncertainty , file=paste0(stackResultsFolder,"/Maps/","/speciesRichness",scenario,typePrediction,"Uncertainty.RData"), compress=TRUE, compression_level=6)
  
}

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------