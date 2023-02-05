#Change Time converter

plotResults <- function(filePathList, nameFigList)
{
  
  library("ggplot2")
  library("forcats")
    
  datalist <- lapply(filePathList, function(x)read.csv(x, sep=',', header = T))
  df <- data.frame(index = character(0), k = integer(0), Time = integer(0), Count = integer(0), Leaf = integer(0))
  
  x <- as.integer(1)
  for(i in datalist)
  {
    
    i['index'] = nameFigList[x]
    x <- x + 1
    
    for(j in seq(5,100,by=5)){

      subsetdf <- subset(i, k == j)
      df[nrow(df)+1,] <- c(i$index[[1]], j, median(as.numeric(subsetdf$Time)), median(as.numeric(subsetdf$Count)), median(as.numeric(subsetdf$Leaf)))
      
    }
    
  }
  
  df[with(df, order(index,k)), ]
  write.csv(df, "out.csv", row.names=FALSE)
  

  
  # df <- read.csv(filePath, sep=',', header=T)
  # df$Time <- df$Time/10^6 
  # 
  # dfPlot <- data.frame(k = integer(0), Time = integer(0), Count = integer(0), Leaf = integer(0))
  # 
  # for(i in seq(5, 100, by=5)){
  #   
  #   subsetDF <- df[df$k == i,]
  #   dfPlot[nrow(dfPlot)+1,] <- c(i, median(subsetDF$Time), median(subsetDF$Count), median(subsetDF$Leaf))
  #             
  # }
  # 
  # plot_k_time <- ggplot(dfPlot, aes(k, Time)) + scale_x_continuous(breaks = seq(from = 0, to = 100, by = 5)) + geom_line() + geom_point() 
  # plot_k_count <- ggplot(dfPlot, aes(k, Count)) + scale_x_continuous(breaks = seq(from = 0, to = 100, by = 5)) + geom_line() + geom_point()
  # plot_k_leaf <- ggplot(dfPlot, aes(k, Leaf)) + scale_x_continuous(breaks = seq(from = 0, to = 100, by = 5)) + geom_line() + geom_point()
  # 
  # ggsave(paste0(nameFig,"_time.png"), plot=plot_k_time, device="png", dpi=300)
  # ggsave(paste0(nameFig,"_count.png"), plot=plot_k_count, device="png", dpi=300)
  # ggsave(paste0(nameFig,"_leaf.png"), plot=plot_k_leaf, device="png", dpi=300)
  
}