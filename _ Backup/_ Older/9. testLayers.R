## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##                                      #####
##                                ####  #####
##                                ####       
##          ####                         
##         ##################             
##           ##################           
##       #######################
##   ##################################   
##  #######################################
##  ######################################
##  ###################################### 
##  ####################################
##  ##################################     
##  ####################                   
##  ###################                    
##  ##################                     
##  #################                      
##  ###############                                     
##      
##  theMarineDataScientist
##
##  github.com/jorgeassis
##  medium.com/themarinedatascientist
##  medium.com/@jorgemfa
##
## -------------------------------------------------------------------------------
##
##  SDM 3.0
##  R Pipelines for Marine Species Distribution Modelling
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

library(raster)
library(rasterVis)

dumpFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers/Test Layers/Images/"
layersFolder <- "/Volumes/Jellyfish/Dropbox/Manuscripts/Bio-ORACLE Across Climate Changes/Bioclimatic Layers"

# ------------------------

files <- list.files( layersFolder, pattern="tif" , full.names = T, recursive = T)
files.names <- list.files( layersFolder, pattern="tif" , full.names = F, recursive = T)

repeat { dev.off() }

for(f in 1:length(files)) {

  if( grepl("PerYear",files[f]) )  { next() }
  if( grepl("Anomaly",files[f]) )  { next() }
  if( grepl("Backup",files[f]) )  { next() }
  if( grepl("Spatial Data",files[f]) )  { next() }
  
  r <- raster(files[f])
  r.val <- range(r[],na.rm=T)

  name <- substr(files[f], as.numeric(gregexpr("Bioclimatic Layers",files[f])[[1]]), nchar(files[f]))
  name.i <- gregexpr("[/]",name)[[1]][1]
  name <- substr(name, name.i+1, nchar(name))
  name <- gsub("[/]"," - ",name)
  name <- gsub("[.]"," ",name)
  name <- gsub(" tif","",name)
  
  jpeg(width=4320, height=2160,units = "px",quality = 75,file=paste0(dumpFolder,"/",name,".jpg"))

  par(mar = c(0, 2, 4 ,12))
  
  pal <- colorRampPalette(c("#9ad1ee","#EADB54","#b10303"))
  breakpoints <- seq(r.val[1], r.val[2], length.out=25)
  breakpoints.l <- seq(r.val[1], r.val[2], length.out=5)
  plot(r, breaks=breakpoints,col=pal(25), axes=FALSE,box=FALSE,legend=FALSE)
  title(name, line = 1,cex.main=2)
  
  arg <- list(at=breakpoints.l, labels=round(breakpoints.l, 3),cex.axis=2)
  
  plot(r, legend.only=TRUE, breaks=breakpoints,col=pal(25), axis.args=arg,
       legend.width=1,
       smallplot=c(0.96,0.965, .075,.962),
       legend.args=list("",side=4, line=-40, cex=0.8))
  
  dev.off()
  
}

