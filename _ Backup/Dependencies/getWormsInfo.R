
library(worms)
speciesListWorms <- data.frame()

for( sp in speciesPredicted) {
  
  speciesList.i <- NULL
  speciesList.j <- NULL
  
  tryCatch( speciesList.i <- wormsbynames( sp , verbose = FALSE , marine_only = "true" ) , error=function(e) { } )
  
  if( is.null(speciesList.i) ) { 
    
    tryCatch( speciesList.i <- wormsbymatchnames( sp , FALSE) , error=function(e) { } ) 

  }
  
  if( ! is.null(speciesList.i) ) { 
    
    if( sum(sapply(colnames(speciesList.i) , function(x) { taxaCrossChecker %in% speciesList.i[,x] } )) != 1 ) { stop(paste0("Error :: ",sp))  }
    if( sum(sapply(colnames(speciesList.i) , function(x) { taxaCrossChecker %in% speciesList.i[,x] } )) == 1 ) { speciesListWorms <- rbind(speciesListWorms,speciesList.i)  }
    
  }

}
