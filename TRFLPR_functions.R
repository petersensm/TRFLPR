# functions for TRFLPR
# v 0.2
# 2014_FEB_17_jwd_edit 
# Joe found that we didn't have a case written where the last bin needed to be merged with the previous one

require(reshape2)
require(plyr)

setClass("Bin", representation(fragments = "data.frame", mean = "numeric", min = "numeric", max = "numeric", count = "numeric", merged = "logical"),
         prototype(fragments = data.frame(), mean = 0, min = 0, max = 0, count = 0, merged = FALSE)) 

# function to parse into bins
bincimate <- function(dataframe, threshold) {
  result <- list()
  # sort the data by fragment size
  df <- dataframe[order(dataframe$Size),]
  # setting stuff for checking against threshold to first fragment
  center <- df[1, "Size" ]
  running.sum <- center
  count <- 1
  # stock first bin with data from first fragment
  binIndex = 1
  result[binIndex] <- new("Bin")
  result[[binIndex]]@fragments <- rbind(result[[binIndex]]@fragments, df[1, ])
  result[[binIndex]]@mean <- df[1, "Size" ]
  result[[binIndex]]@min <- df[1, "Size" ]
  result[[binIndex]]@max <- df[1, "Size" ]
  result[[binIndex]]@count <- 1
  
  for (i in 2:length(df$Size)) {
    if(abs(df[i, "Size" ] - center) > threshold) {
      # make 
      binIndex <- binIndex + 1
      result[binIndex] <- new("Bin")
      # reset stuff for calculating within threshold
      running.sum <- 0
      count <- 0
    }                                                   
    result[[binIndex]]@fragments <- rbind(result[[binIndex]]@fragments, df[i, ])
    running.sum <- df[i, "Size" ] + running.sum
    count <- count + 1
    center <- running.sum/count
    result[[binIndex]]@mean <- center
    result[[binIndex]]@count <- result[[binIndex]]@count + 1 
    if(result[[binIndex]]@min == 0 | result[[binIndex]]@min > df[i, "Size" ]) {
      result[[binIndex]]@min <- df[i, "Size" ]
    }
    if(result[[binIndex]]@max < df[i, "Size" ]) {
      result[[binIndex]]@max <- df[i, "Size" ]
    }
  }
  result
}

# function to merge bins - internal function
MergeBins <- function(bin1, bin2) {
  
  #s <- sprintf("Merging Bins %f and %f", bin1@mean, bin2@mean)
  #print(s)
  
  result <- new("Bin")
  result@fragments <- rbind(bin1@fragments, bin2@fragments) 
  result@mean <- mean(result@fragments$Size)
  result@min <- min(result@fragments$Size)
  result@max <- max(result@fragments$Size)
  result@count <- length(result@fragments$Size)
  result
}


# function to calc diff btwn bins - internal function
BinDifferences <- function(Bins){
  result <- vector()
  for(i in 2:length(Bins)) {
    d1 <- abs(Bins[[i-1]]@mean - Bins[[i]]@mean)
    result[i] <- d1
  }
  
  result
}

# internal function
Non.overlapping.bins <- function(bin1, bin2) {
  #print(intersect(unique(bin1@fragments$Sample.File.Name), unique(bin2@fragments$Sample.File.Name)))
  length(intersect(unique(bin1@fragments$Sample.File.Name), unique(bin2@fragments$Sample.File.Name))) == 0
}

findMergeCandidate <- function(Bins, threshold){
  
  BD <- BinDifferences(Bins)
  minIndex = 0
  minValue = BD[2]
  
  for (i in 2:length(BD)){
    if(BD[i] < threshold & Non.overlapping.bins(Bins[[i]], Bins[[i-1]]) ){
      if(BD[i] <= minValue){
        minIndex <- i
        minValue <- BD[i]
      }
    }
  }
  
  minIndex
}

# function to lump bins
Lump.bins <- function(Bins, threshold){
  
  newBins <- Bins
  candidateIndex <- findMergeCandidate(newBins, threshold)
  
  while(candidateIndex > 0){
   # print(candidateIndex)  
    if(candidateIndex > 2){
      if(candidateIndex < length(newBins)){
      newBins <- c(newBins[1:(candidateIndex-2)],
                   MergeBins(newBins[[candidateIndex-1]],
                             newBins[[candidateIndex]]),
                   newBins[(candidateIndex+1):length(newBins)])
      } else {
           newBins <- c(newBins[1:(candidateIndex-2)],
                   MergeBins(newBins[[candidateIndex-1]],
                             newBins[[candidateIndex]]))
      }
    } else {
      
      newBins <- c(list(),
                   MergeBins(newBins[[candidateIndex-1]],
                             newBins[[candidateIndex]]),
                   newBins[(candidateIndex+1):length(newBins)])
    }
    
    
    candidateIndex <- findMergeCandidate(newBins, threshold)
    
  }
  
  newBins
}

# function to tag bins with a sequential number
# run after Lump.bins but before converting to a dataframe

tagBins <- function(Bins) {
  result <- list()
  i <- 1
  for(b in Bins){
    b@fragments$tag <- rep(i)
    result <- c(result, b)
    i <- i + 1;
  }
  result
}

# function to spit a data frame out of object class Bins
Bins.to.data.frame <- function(Bins) {
  result <- data.frame()
  for(b in Bins){
    # print(b@fragments$tag)
    b@fragments$Bin_mean <- rep(b@mean)
    result <- rbind(result, b@fragments)
  }
  result
}
