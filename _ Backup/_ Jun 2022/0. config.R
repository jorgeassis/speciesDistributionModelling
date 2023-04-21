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

nCores <- 8

# ---------------------

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0.1),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_blank()
  )

mainTheme <- theme(
  text = element_text(size=12) ,
  panel.background = element_rect(fill = "#F0F0F0", colour = "#F0F0F0",size = 0, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',colour = "#EFEFEF"), 
  axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
  axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
  axis.text.x =element_text(size=11, margin = margin(t = 8, r = 0, b = 0, l = 0)),
  axis.text.y =element_text(size=11, margin = margin(t = 0, r = 8, b = 0, l = 0)))

# ---------------------

regionBuffer <- 20

autocorrelationClassDistance <- 5
autocorrelationMaxDistance <- 250
autocorrelationSubsetRecords <- 1000
autocorrelationSignif <- 0.05

# -------------------------------------------------
# Raster layers and autocorrelation

relocateType <- "distance" # distance / nonDistance
relocateSpeciesDistance <- 20
relocateSpeciesDepth <- TRUE # only if relocateType == "distance"

# -------------------------------------------------
# Modelling algorithms 

algorithms <- "BRT" # c("BRT","MBOOST")
useWeights <- TRUE

# -------------------------------------------------
# Cross-validation
# References: blockCV. Methods Ecol Evol. 2019; 10:225–232. https://doi.org/10.1111/2041-210X.13107

cvIndex <- "tss" # tss area
cvKFolds <- 6
cvType <- "randomBlocks" # "randomBlocks" "blocksLatitudinal" "blocksLongitudinal" 
# cvRemoveEdges <- TRUE # only for blocksLatitudinal blocksLongitudinal

# -------------------------------------------------
# Pseudo-absences

paRatio <- 1 # if ratio.p.pa < 1 is prevalence ( 1/ratio.p.pa ); If > 1 is absolute
paMinimum <- 1000
paType <- "random"
paSREThrehold <- 0.005 # 0.025

paBiasKernelSurface <- TRUE 
paBiasKernelSurfaceProb <- 0.01
paEnvironmentStratification <- TRUE
paMindistance <- 100

# ------------------------------------------------------------------------------------

mboostShrinkage <- seq(from = 0.25, to = 1, by = 0.25) 
mboostDF <- seq(from = 2, to = 12, by = 2)
mboostIterations <- seq(from = 50, to = 250, by = 50)

brtLearning <- c(0.005,0.001) 
brtTreeDepth <- 2:6
brtBagFraction <- 0.5 # seq(from = 0.5, to = 0.75, by = 0.25)
brtMaxTrees <- 1000

mpnnHidden <- c(1,2,3) # 1
mpnnItereractions <- 100 # c(10,25,50,100) 

maxentBetamultiplier <- c(10,15,20)
maxentFeature <- data.frame( span.1 = c("linear=false","quadratic=false","hinge=true","threshold=false","product=false"),
                             span.2 = c("linear=false","quadratic=false","hinge=true","threshold=true","product=false"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
