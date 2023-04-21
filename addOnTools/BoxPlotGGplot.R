## --------------------------------------------------
## --------------------------------------------------

closeAllConnections()
gc(reset=TRUE)

# Boxplot performance

data <- read.csv(paste0(stackResultsFolder,"/performance.csv"), sep=",")

data <- rbind(
  
  data.frame( Index = "TSS" , Algorithm= "zEnsemble" , Value = data$tss.ensemble ),
  data.frame( Index = "AUC" , Algorithm= "zEnsemble" , Value = data$auc.ensemble ),
  data.frame( Index = "DE" , Algorithm= "zEnsemble" , Value = data$deviance.ensemble ),

  data.frame( Index = "TSS" , Algorithm= "zEnsemble" , Value = data$tss.ensembleReachab ),
  data.frame( Index = "AUC" , Algorithm= "zEnsembleReach" , Value = data$auc.ensembleReachab ),
  data.frame( Index = "DE" , Algorithm= "zEnsembleReach" , Value = data$deviance.ensembleReachab ),

  data.frame( Index = "TSS" , Algorithm= "BRT" , Value = data$tss.BRT ),
  data.frame( Index = "AUC" , Algorithm= "BRT" , Value = data$auc.BRT ),
  data.frame( Index = "DE" , Algorithm= "BRT" , Value = data$deviance.BRT ),

  data.frame( Index = "TSS" , Algorithm= "AdaBoost" , Value = data$tss.MBOOST ),
  data.frame( Index = "AUC" , Algorithm= "AdaBoost" , Value = data$auc.MBOOST ),
  data.frame( Index = "DE" , Algorithm= "AdaBoost" , Value = data$deviance.MBOOST ),

  data.frame( Index = "TSS" , Algorithm= "XgBoost" , Value = data$tss.XGBOOST ),
  data.frame( Index = "AUC" , Algorithm= "XgBoost" , Value = data$auc.XGBOOST ),
  data.frame( Index = "DE" , Algorithm= "XgBoost" , Value = data$deviance.XGBOOST ))

a <- ggplot(data, aes(x=Value, y=Algorithm, fill=Index)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_boxplot(lwd=0.25, outlier.size=0.4, outlier.colour="#949494") +
      coord_flip() +
      scale_fill_manual(values=c("#FFE600","#F29A42","#E67B73")) +
      theme(legend.position="none")

# Relative Contribution

data <- read.csv(paste0(stackResultsFolder,"/relativeContribution.csv"), sep=",")
data <- data[,which(grepl("ensemble",names(data)))]
names(data) <- gsub("\\.ensemble","",names(data))
names(data)[6] <- "Wave energy"

data <- tidyr::gather(data)
names(data) <- c("Predictor","Contribution")

data$Predictor <- with(data, reorder(Predictor , -Contribution, median , na.rm=T))

b <- ggplot(data, aes(x=Predictor, y=Contribution, fill=Contribution)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_boxplot(width=0.5, fill="#DE4D43", lwd=0.25, outlier.size=0.4, outlier.colour="#949494") +
      scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
      labs(y="Relative contribution (%)")+
      theme(legend.position="none")


library(gridExtra)

pdf(paste0(stackResultsFolder,"./performanceContributionPlot.pdf"), width=12, height=5)
grid.arrange(a, b, nrow = 1)
dev.off()
