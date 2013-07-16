# defining a new class for TRFLPs

setClass("Bin", representation(fragments = "data.frame", mean = "numeric", min = "numeric", max = "numeric", count = "numeric", merged = "logical"),
         prototype(fragments = data.frame(), mean = 0, min = 0, max = 0, count = 0, merged = FALSE)) 


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


MergeBins <- function(bin1, bin2) {

  s <- sprintf("Merging Bins %f and %f", bin1@mean, bin2@mean)
  print(s)

  result <- new("Bin")
  result@fragments <- rbind(bin1@fragments, bin2@fragments)
  if(length(result@fragments[duplicated(result@fragments),]$Size) > 0){
      print(result@fragments[duplicated(result@fragments),])
      print("found duplicates")
  }
          
  result@mean <- mean(result@fragments$Size)
  result@min <- min(result@fragments$Size)
  result@max <- max(result@fragments$Size)
  result@count <- length(result@fragments$Size)
  result
}



BinDifferences <- function(Bins){
  result <- vector()
  for(i in 2:length(Bins)) {
  d1 <- abs(Bins[[i-1]]@mean - Bins[[i]]@mean)
  result[i] <- d1
  }
      
  result
}



findMergeCandidate <- function(Bins, threshold){

    BD <- BinDifferences(Bins)
    minIndex = 0
    minValue = BD[2]
    
    for (i in 2:length(BD)){
        if(BD[i] < threshold ){
            if(BD[i] <= minValue){
                minIndex <- i
                minValue <- BD[i]
            }
        }
    }
    
    minIndex
}



Criterion3 <- function(Bins, threshold){

    newBins <- Bins
    candidateIndex <- findMergeCandidate(newBins, threshold)
    
    while(candidateIndex > 0){

        if(candidateIndex > 2){
            newBins <- c(newBins[1:(candidateIndex-2)],
                         MergeBins(newBins[[candidateIndex-1]],
                                   newBins[[candidateIndex]]),
                         newBins[(candidateIndex+1):length(newBins)])
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



Criterion2 <- function(Bins, BD, threshold){
    result <- list ()
    for(i in 2:length(BD)){
        if(prevDiff(BD[[i]]) < threshold & nextDiff(BD[[i]]) < threshold){
            if(prevDiff(BD[[i]]) <= nextDiff(BD[[i]])){

                if(!Bins[[i-1]]@merged){
                    Bins[[i-1]]@merged <- TRUE
                    Bins[[i]]@merged <- TRUE     
                    result <- c(result, MergeBins(Bins[[i-1]], Bins[[i]]))
                } else {
                    if( nextDiff(BD[[i]]) < threshold ){
                        # continue on our way
                    } else {
                        result <- c(result, Bins[[i]])
                    }
                }
                if(i == length(BD)){
                    result <- c(result, Bins[[i+1]])
                }
                
            } else if (nextDiff(BD[[i]]) < prevDiff(BD[[i]])){

                if( i == 2 & i < length(BD)
                   & nextDiff(BD[[i+1]]) < prevDiff(BD[[i+1]])){

                    result <- c(result, Bins[[i-1]], Bins[[i]])

                } else if( i > 2 & i < length(BD)
                          & nextDiff(BD[[i+1]]) < prevDiff(BD[[i+1]])){

                    result <- c(result, Bins[[i]])

                } else if ( i == length(BD) ){

                    result <- c(result, MergeBins(Bins[[i]], Bins[[i+1]]))

                } else {

                    result <- c(result, Bins[[i]])
                }
                
            }

        } else if (prevDiff(BD[[i]]) < threshold) {

            if(!Bins[[i-1]]@merged){
                Bins[[i-1]]@merged <- TRUE
                Bins[[i]]@merged <- TRUE
                result <- c(result, MergeBins(Bins[[i-1]], Bins[[i]]))
            } else {
                if( nextDiff(BD[[i]]) < threshold ){
                                        # continue on our way
                } else {
                    result <- c(result, Bins[[i]])
                }
            }
            if( i == length(BD)){
                result <- c(result, Bins[[i+1]])
            }
        } else if (nextDiff(BD[[i]]) < threshold){
  
            result <- c(result, Bins[[i-1]])

            if(i < length(BD) && nextDiff(BD[[i+1]]) < prevDiff(BD[[i+1]])){
                result <- c(result, Bins[[i]])
            } else if (i == length(BD)){
                result <- c(result, MergeBins(Bins[[i]], Bins[[i+1]]))
            }
        } else {

            result <- c(result, Bins[[i]]);
            if( i == length(BD) ){
                result <- c(result, Bins[[i+1]])
            }
        }
        

    }

    result
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



findBin <- function(Bins, m){

    found <- NULL
    for( b in Bins){
        print(b@mean)
        if(b@mean == m){
            found <- b
            break
        }
    }
    found
}



# testing

Testflight <- bincimate(Blue.5, 0.25)
tf2 <- tagBins(Testflight)
c1 <- Criterion3(tf2,0.5)
df2 <- Bins.to.data.frame(c1)
length(df2$tag)
