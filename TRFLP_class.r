# defining a new class for TRFLPs

setClass("Bin", representation(fragments = "data.frame", mean = "numeric", min = "numeric", max = "numeric", count = "numeric"),
         prototype(fragments = data.frame(), mean = 0, min = 0, max = 0, count = 0))

Temp <- new("Bin") 

Temp@fragments
Temp@mean

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

Testflight <- bincimate(Blue.5, 0.25)

Testflight[2]

# function to spit a data frame out of object class Bins
Bins.to.data.frame <- function(Bins) {
  result <- data.frame()
  for(b in Bins){
    print(b@fragments$tag)
    b@fragments$Bin <- rep(b@mean)
    result <- rbind(result, b@fragments)
  }
  result
}

Testframe <- Bins.to.data.frame(Testflight)

length(unique(Testframe$Bin))

MergeBins <- function(bin1, bin2) {
  result <- new("Bin")
  result@fragments <- rbind(bin1@fragments, bin2@fragments)
  result@mean <- mean(result@fragments$Size)
  result@min <- min(result@fragments$Size)
  result@max <- max(result@fragments$Size)
  result@count <- length(result@fragments$Size)
  result
}

Testflight[1]
Testflight[2]

MergeBins(Testflight[[1]], Testflight[[2]])

BinDifferences <- function(Bins){
  result <- list()
  for(i in 2:(length(Bins)-1)) {
  d1 <- abs(Bins[[i-1]]@mean - Bins[[i]]@mean)
  d2 <- abs(Bins[[i]]@mean - Bins[[i+1]]@mean)
  v=c(d1,d2)
  result[[i]] <- v 
  }
  result
}

BD <- BinDifferences(Testflight)
BD

prevDiff <- function(ds){
  ds[1]
}

nextDiff <- function(ds){
  ds[2]
}

Criterion1 <- function(Bins, BD, threshold) {
  result <- list()
  print(length(BD))
  for(i in 2:length(BD)){
    print(i)
    if(prevDiff(BD[[i]]) < threshold & nextDiff(BD[[i]]) < threshold){
        if(prevDiff(BD[[i]]) < nextDiff(BD[[i]])){
          result <- c(result, MergeBins(Bins[[i-1]], Bins[[i]]))
        } else if(nextDiff(BD[[i+1]]) < threshold & nextDiff(BD[[i+1]]) < prevDiff(BD[[i+1]])){ 
          result <- c(result, Bins[[i-1]])
        }  
    } else if(prevDiff(BD[[i]]) < threshold){
      result <- c(result, MergeBins(Bins[[i-1]], Bins[[i]]))
    } else if(nextDiff(BD[[i]]) < threshold){
      if(i == length(BD)){
        result <- c(result, MergeBins(Bins[[i]], Bins[[i+1]]))
      }else if(nextDiff(BD[[i+1]]) < threshold & nextDiff(BD[[i+1]]) < prevDiff(BD[[i+1]])){ 
          result <- c(result, Bins[[i]])
      } else {
        #result <- c(result, Bins[[i-1]])
      }
    } else {
      if(i == 2){
        result <- c(result, Bins[[i-1]], Bins[[i]])
      }else if(i == length(BD)){
        print("i am here, sucker")
        result <- c(result, Bins[[i]], Bins[[i+1]])
      } else {
        result <- c(result, Bins[[i]])
      }
    }    
  }
  result
}

tagBins <- function(Bins){
  result <- list()
  i <- 1
  for(b in Bins){
    
    b@fragments$tag <- rep(i)
    result <- c(result, b)
    i <- i + 1;
  }
  
  result
}

tf2 <- tagBins(Testflight)
length(tf2)
c1 <- Criterion1(tf2,BD,0.5)
fuck[fuck$tag==83,]
fuck <- Bins.to.data.frame(c1)
tail(fuck)
duplicated(fuck)
setdiff(unique(row.names(Testframe)), unique(row.names(fuck)))
length(setdiff(unique(row.names(Testframe)), unique(row.names(fuck)))) # 83
2221-2140 # 81
?unique
BD[76:78]
Testframe["11752", ]
Testframe["48",]
Testframe["14967",]



Testflight[[2]]
c1[[2]]
