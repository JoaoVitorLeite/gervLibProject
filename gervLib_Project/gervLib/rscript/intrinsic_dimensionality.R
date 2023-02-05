args <- commandArgs(trailingOnly = TRUE)

intrinsic_dim <- function(datasetPath, separator, distanceFunction){

  library("parallelDist")
  
  df <- read.csv(datasetPath, sep = separator, header=F)

  mat <- parallelDist(as.matrix(df), method = distanceFunction, threads=12)
  
  mean_values <- mean(mat)
  std_values <- sd(mat)
  ans <- (mean_values*mean_values)/(2*std_values*std_values)

  return(ans)
  
}

intrinsic_dim(args[1], args[2], args[3])
