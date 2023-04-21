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

## library(credentials)
## set_github_pat()

# ---------------------
# ---------------------

nCores <- 6

# ---------------------

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.box.background = element_rect(fill='#FFFFFF'),
    panel.border = element_blank()
  )

themePlot <- theme(
  text = element_text(size=12) ,
  panel.background = element_rect(fill = "#F0F0F0", colour = "#F0F0F0",size = 0, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',colour = "#EFEFEF"), 
  axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
  axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
  axis.text.x = element_text(size=11, margin = margin(t = 8, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(size=11, margin = margin(t = 0, r = 8, b = 0, l = 0)))

# -------------------------------------------------
# Raster layers and autocorrelation

regionBuffer <- 25
regionBufferPoly <- NULL
regionBufferPolyLevel <- NULL

autocorrelationClassDistance <- 10
autocorrelationMaxDistance <- 500

relocateType <- "distance" # distance / nonDistance
relocateSpeciesDistance <- 100
relocateSpeciesDepth <- TRUE # only if relocateType == "distance"

minOccurrenceRecords <- 3

# -------------------------------------------------
# Cross-validation
# References: blockCV. Methods Ecol Evol. 2019; 10:225–232. https://doi.org/10.1111/2041-210X.13107

cvKFolds <- 10
cvType <- "blocks" # blocks randomBlocks" blocksLatitudinal blocksLongitudinal
cvIndex <- "tss" # tss auc aicc
reclassificationIndex <- "tss" # tss auc minimumTrain minimumTrain95

# -------------------------------------------------
# Pseudo-absences

paRatio <- 1 # if ratio.p.pa < 1 is prevalence ( 1/ratio.p.pa ); If > 1 is absolute
paMinimum <- 1000 #
paType <- "mess" # use mess sre kernelDensity random mahalanobis

# ------------------------------------------------------------------------------------
# Modelling algorithms 

reduceModelcomplexity <- FALSE
reduceModelcomplexityThreshold <- 5

algorithms <- c("BRT","XGBOOST","MBOOST") #  c("BRT","MBOOST","XGBOOST,"MaxEnt","MPNN")

xgboostShrinkage <- seq(from = 0.1, to = 0.5, by = 0.1) 
xgboostGamma <- seq(from = 0, to = 5, by = 1)
xgboostDepth <- seq(from = 2, to = 3, by = 1)
xgboostRounds <- seq(from = 10, to = 100, by = 10)

mboostShrinkage <- seq(from = 0.1, to = 0.5, by = 0.1) 
mboostDF <- seq(from = 5, to = 50, by = 5)
mboostMStop <- seq(from = 10, to = 100, by = 10)

brtLearning <- c(0.1,0.01,0.001) 
brtTreeDepth <- seq(from = 2, to = 3, by = 1)
brtNTrees <- seq(100,1000,by=100)

mpnnHidden <- c(1,2,3) # 1
mpnnItereractions <- 100 # c(10,25,50,100) 

maxentBetamultiplier <- c(10,15,20)
maxentFeature <- data.frame( span.1 = c("linear=false","quadratic=false","hinge=true","threshold=false","product=false"),
                             span.2 = c("linear=false","quadratic=false","hinge=true","threshold=true","product=false"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------