args <- commandArgs(trailingOnly = TRUE)

#Split 0.9/0.1 => Train/Test
split <- function(pathOfDataset, separator, seed = 100, trainFileName, testFileName)
{

  data <- read.csv(pathOfDataset, sep = separator, header=F)
  set.seed(seed) 
  sample <- sample.int(n = nrow(data), size = floor(.9*nrow(data)), replace = F)
  train <- data[sample, ]
  test  <- data[-sample, ]  
  
  write.table(train, file = trainFileName, sep = ",", col.names = FALSE, row.names = FALSE)
  write.table(test, file = testFileName, sep = ",", col.names = FALSE, row.names = FALSE)
  
}

if(length(args) == 5){
  
  split(args[1], args[2], args[3], args[4], args[5])
  
} else if(length(args) == 4){
  
  split(args[1], args[2], 100, args[3], args[4])
  
}