## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
gc(reset=TRUE)

# -----------------------------
# Boxplot performance

data <- read.csv(paste0(stackResultsFolder,"/performance.csv"), sep=",")

data <- rbind(
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "BRT" , Value = data$meanCV.sensitivity.BRT ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "BRT" , Value = data$meanCV.auc.BRT ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "BRT" , Value = rescale(data$meanCV.boyce.BRT, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "AdaBoost" , Value = data$meanCV.sensitivity.MBOOST ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "AdaBoost" , Value = data$meanCV.auc.MBOOST ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "AdaBoost" , Value = rescale(data$meanCV.boyce.MBOOST, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Cross validation", Index = "Sensitivity" , Algorithm= "XgBoost" , Value = data$meanCV.sensitivity.XGBOOST ),
  data.frame( Group= "Cross validation", Index = "AUC" , Algorithm= "XgBoost" , Value = data$meanCV.auc.XGBOOST ),
  data.frame( Group= "Cross validation", Index = "Boyce" , Algorithm= "XgBoost" , Value = rescale(data$meanCV.boyce.XGBOOST, to = c(0, 1), range = c(-1, 1)) ),
  
  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "zEnsemble" , Value = data$sensitivity.Ensemble ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "zEnsemble" , Value = data$auc.Ensemble ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "zEnsemble" , Value = rescale(data$boyce.Ensemble, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "zEnsembleReach" , Value = data$sensitivity.EnsembleReachab ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "zEnsembleReach" , Value = data$auc.EnsembleReachab ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "zEnsembleReach" , Value = rescale(data$boyce.EnsembleReachab, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "BRT" , Value = data$sensitivity.BRT ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "BRT" , Value = data$auc.BRT ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "BRT" , Value = rescale(data$boyce.BRT, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "AdaBoost" , Value = data$sensitivity.MBOOST ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "AdaBoost" , Value = data$auc.MBOOST ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "AdaBoost" , Value = rescale(data$boyce.MBOOST, to = c(0, 1), range = c(-1, 1)) ),

  data.frame( Group= "Final prediction", Index = "Sensitivity" , Algorithm= "XgBoost" , Value = data$sensitivity.XGBOOST ),
  data.frame( Group= "Final prediction", Index = "AUC" , Algorithm= "XgBoost" , Value = data$auc.XGBOOST ),
  data.frame( Group= "Final prediction", Index = "Boyce" , Algorithm= "XgBoost" , Value = rescale(data$boyce.XGBOOST, to = c(0, 1), range = c(-1, 1)) ))

data <- data[complete.cases(data),]

a <- ggplot(data, aes(x=Value, y=Algorithm, fill=Index)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_boxplot(lwd=0.25, outlier.size=0.4, outlier.colour="#949494") +
      coord_flip() +
      scale_fill_manual(name=NULL, values=c("#FFE600","#F29A42","#E67B73")) +
      xlab("Model performace") +
      facet_grid(~Group, scales = "free_x", space = "free_x", switch = "x") + 
      theme(axis.title.x = element_blank()) +
      theme(strip.placement = "outside") +
      theme(strip.background =element_rect(fill="white",colour = "white", size = 1)) +
      theme(legend.position="none") 
a

# Relative Contribution

data <- read.csv(paste0(stackResultsFolder,"/relativeContribution.csv"), sep=",")
data <- data[,which(grepl("ensemble",names(data)))]
names(data) <- gsub("\\.Ensemble","",names(data))
names(data) <- gsub("\\.ensemble","",names(data))

data <- tidyr::gather(data)
names(data) <- c("Predictor","Contribution")

data$Predictor <- with(data, reorder(Predictor , -Contribution, median , na.rm=T))

b <- ggplot(data, aes(x=Predictor, y=Contribution, fill=Contribution)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_boxplot(width=0.5, fill="#DE4D43", lwd=0.25, outlier.size=0.4, outlier.colour="#949494") +
      scale_y_continuous(limits=c(0,80), breaks=seq(0,80,20)) +
      labs(y="Relative contribution (%)")+
      theme(legend.position="none") + geom_hline(yintercept=5,linetype=2, size=0.25)


library(gridExtra)

pdf(paste0(stackResultsFolder,"./performanceContributionPlot.pdf"), width=12, height=5)
grid.arrange(a, b, nrow = 1)
dev.off()
