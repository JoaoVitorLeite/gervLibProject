#Change max/min of the graph scale, line 39 and 40

args <- commandArgs(trailingOnly = TRUE)

figPlotPivots <- function(fileName, pivotsIndex){
  
  library("ggplot2")
  library("gridExtra")
  
  data <- read.csv(fileName, header = F)
  colnames(data) <- c("x", "y")
  
  #index must start at 1
  for(i in 1:length(pivotsIndex)){
    pivotsIndex[i] <- pivotsIndex[i]+1
  }
  
  pivots <- c()
  points <- c()
  
  for(i in seq(1:nrow(data))){
    
    if(i %in% pivotsIndex) pivots <- c(pivots, i)
    else points <- c(points, i)
    
  }
  
  dataPoints <- data.frame(x = integer(0), y = integer(0))
  dataPivots <- data.frame(x = integer(0), y = integer(0))
  
  for(i in pivots){
    dataPivots[nrow(dataPivots)+1,] <- data[i,]
  }
  
  for(i in points){
    dataPoints[nrow(dataPoints)+1,] <- data[i,]
  }

  pltBase <- ggplot(dataPoints, aes(x=x,y=y)) + 
    geom_point(size=3) + 
    scale_x_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
    scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20))

  listOfPlots <- list()
   
  lab <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  x <- 1
  
  for(i in 1:nrow(dataPivots)){
    
    pltAux <- pltBase + 
      geom_point(data=dataPivots[1:i,], mapping = aes(x=x,y=y, color="red"), size=3) +
      scale_shape_identity() +
      theme_classic() +
      theme(text = element_text(size=12), legend.position = "none", axis.text=element_text(size=12)) +
      xlab(lab[[x]]) +
      ylab("")
    
    if(i != nrow(dataPivots)) pltAux <- pltAux + 
        geom_point(data=dataPivots[(i+1):nrow(dataPivots),], mapping = aes(x=x,y=y), size=3)
    
    listOfPlots[[paste("plt",i, sep="")]] <- pltAux
    
    x <- x + 1
    
    pltAux
    
  }

  do.call("grid.arrange", c(listOfPlots, ncol=3))

}

figPlotPivots(args[1], args[2])
