args <- commandArgs(trailingOnly = TRUE)

sample_dataset <- function(filePath, separator, sizeOfSample, sampleFileName = "sample.csv", seed = 100){
  
  set.seed(100)
  df <- data.frame(read.csv(filePath, sep = separator, header=F))
  df <- df[sample(nrow(df), sizeOfSample), ]
  write.table(df, file=paste0(getwd(), "/", sampleFileName), sep=",", row.names = F, col.names = F)

}

if(length(args) == 5){
  
  sample_dataset(args[1], args[2], args[3], args[4], args[5])
  
} else if(length(args) == 4){
  
  sample_dataset(args[1], args[2], args[3], args[4])
  
} else if(length(args) == 3){
  
  sample_dataset(args[1], args[2], args[3])
  
}
