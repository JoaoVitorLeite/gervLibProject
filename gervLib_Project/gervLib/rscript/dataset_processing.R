args <- commandArgs(trailingOnly = TRUE)

dataset_processing <- function(fileName, separator, outputFileName){
  
  
  df <- read.csv(fileName, sep = separator, header = F)
  df <- df[!duplicated(df),]
  df <- as.data.frame(lapply(df, min_max))
  write.table(df, file=paste0(getwd(), "/", outputFileName), sep=",", row.names = F, col.names = F)
  
}

min_max <- function(x){
  if(sum(x) == 0.0)
    return(x)
  else
    return((x - min(x))/(max(x) - min(x)))
}

dataset_processing(args[1], args[2], args[3])
