# TRFLP binning program in R

# Step 1: import data from peak scanner ---------
#   remove all dyes not of interest 
#   split process for each dye

library(reshape2)

master <- read.csv ("TRFLP Master_Data.csv")

class(master)
str(master)
# Null's out columns not needed for analysis
master[17:26] <- list(NULL) 
master [3:11] <- list(NULL)

# Dataframe of Dye and Sample peak split apart
master_2 <-  colsplit(master$Dye.Sample.Peak, ",", c('Dye', 'Sample.Peak'))

# New master Dataframe with the seperated Dye and Sample Peak columns
master <- cbind(master, master_2)

master [3] <- list(NULL)

#Reorders columns for asthetic purposes
master <- master[, c(1,2,7,8,3,4,5,6)]

# master[c("Size")][is.na(master[c("Size")])] <- 0

# subset dyes
    Blue <- subset(master, Dye == "B") # 952 obs
    
    Green <- subset(master, Dye == "G")


# Step 2: set for threshold abnormally sized fragments (< 50 or > 600) --------


# save frags with NA's to another dataset for later
  Blue.missing <- Blue[is.na(Blue$Size), ] # 96 obs. # list of samples with no data for a future error/summary log
    Blue.missing$Frag.Quality <- rep("missing", length(Blue.missing$Size))
     
# keep frags w/o na's
  Blue <- Blue[complete.cases(Blue$Size), ]  # 856 obs. good samples
  

#Subset for Numeric variables but you can also use factors

# remove small frags
  Blue.toosmall <- Blue[Blue$Size < 50.0, ] # 33 obs. # list of frags too small  for a future error/summary log
    Blue.toosmall$Frag.Quality <- rep("too_small", length(Blue.toosmall$Size))
  Blue <- Blue[Blue$Size >= 50.0,] # 823 obs.
  
# remove big frags
  Blue.toobig <- Blue[Blue$Size > 600.0,] # 2 obs. # list of frags too big  for a future error/summary log
    Blue.toobig$Frag.Quality <- rep("too_big", length(Blue.toobig$Size))
  Blue <- Blue[Blue$Size <= 600.0,] # 821 obs

# drops sample names if they have no analyzable blue frags
  Blue$Sample.File.Name <- Blue$Sample.File.Name[,drop = TRUE]




# Step 3: calculate relative peak area for each frag in each sample-------


# make a function -- assumes the sample name column is Sample.File.Name
relative.abundance <- function(input) {
    datain <- input
    dataout <- rbind()
    SampleNames <- unique(datain$Sample.File.Name)
  for (i in seq(unique(datain$Sample.File.Name))) {
    x.i <- SampleNames[i]
    Sample.i <- datain[datain$Sample.File.Name == x.i,]
    sum.i <- sum(Sample.i$Area.in.Point)
    Sample.i$Relative.Area <- (Sample.i$Area.in.Point)/sum.i*100
    dataout <- rbind(dataout, Sample.i)
  }
  dataout
}

# Blue.3 is new dataframe with relative abundance column 
    Blue.3 <- relative.abundance(Blue)

    # check that rel ab adds to 100 for each sample
    tapply(Blue.3$Relative.Area, Blue.3$Sample.File.Name, sum)

  
# Step 4: removing unrepeatable peaks ----- 


# remove frags with relative area < 1%
  
  # save the stuff to be deleated for summary log
    Blue.tinyarea <- Blue.3[Blue.3$Relative.Area < 1, ] 
    Blue.tinyarea $Frag.Quality <- rep("tiny_area", length(Blue.tinyarea$Size))
    Blue.tinyarea [9] <- list(NULL)

# keep stuff with areas bigger than or = to 1%
    Blue.3.2 <- Blue.3[Blue.3$Relative.Area >= 1,] # 735 obs.

# repeat realtive area calcs with remaining fragments
    Blue.4 <- relative.abundance(Blue.3.2)
    
  # check that rel ab adds to 100 for each sample
    tapply(Blue.4$Relative.Area, Blue.4$Sample.File.Name, sum)   
    

#Creates a final dataframe of good Samples
  Blue.good <- Blue.4
    Blue.good$Frag.Quality <- rep("Good", length(Blue.good$Size))
    Blue.good [9] <- list(NULL)

# Dataframe with all samples merged back together with Quality tag at end to specify where it came from
Blue_Quality_Master <- rbind(Blue.good, Blue.missing, Blue.tinyarea, Blue.toobig, Blue.toosmall)


# Steps 5-6: binning OTUs using a .25 difference from center of bin------

# set some initial stuff

# defining a new class for TRFLPs

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



# function to spit a data frame out of object class Bins
Bins.to.data.frame <- function(Bins) {
  result <- data.frame()
  for(b in Bins){
    # print(b@fragments$tag)
    b@fragments$Bin <- rep(b@mean)
    result <- rbind(result, b@fragments)
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

# testing

Blue.bins1 <- bincimate(Blue.4, 0.25)
Blue.bins1.dataframe <- tagBins(Blue.bins1 )

Blue.bins.final <- Lump.bins(Blue.bins1, 0.5)
Blue.bins.final.tag <- tagBins(Blue.bins.final)
Blue.bins.final.dataframe <- Bins.to.data.frame(Blue.bins.final.tag)
# length(df2$tag)
require(reshape2)
#Counts <- dcast(df2, Sample.File.Name ~ Bin, length)


# Step 7: OTU sample matix for specified Dye------

Blue_Matrix <- dcast(Blue.bins.final.dataframe, Sample.File.Name ~ Bin, value.var = "Relative.Area", sum, fill = 0, drop = F, margins = T)

write.csv(Blue_Matrix, "Blue_Matrix.csv", row.names=F)  # creates new .csv file with all the data combined

# write.table(Blue_Matrix, "Blue_Matrix.txt", row.names=F)  # creates new .csv file with all the data combined




# Steps 2-7 Repeated using the Green Dye -------------

# Step 2: set for threshold abnormally sized fragments (< 50 or > 600) --------


# save frags with NA's to another dataset for later
Green.missing <- Green[is.na(Green$Size), ] # 96 obs. # list of samples with no data for a future error/summary log
Green.missing$Frag.Quality <- rep("missing", length(Green.missing$Size))

# keep frags w/o na's
Green <- Green[complete.cases(Green$Size), ]  # 856 obs. good samples


#Subset for Numeric variables but you can also use factors

# remove small frags
Green.toosmall <- Green[Green$Size < 50.0, ] # 33 obs. # list of frags too small  for a future error/summary log
Green.toosmall$Frag.Quality <- rep("too_small", length(Green.toosmall$Size))
Green <- Green[Green$Size >= 50.0,] # 823 obs.

# remove big frags
Green.toobig <- Green[Green$Size > 600.0,] # 2 obs. # list of frags too big  for a future error/summary log
Green.toobig$Frag.Quality <- rep("too_big", length(Green.toobig$Size))
Green <- Green[Green$Size <= 600.0,] # 821 obs

# drops sample names if they have no analyzable Green frags
Green$Sample.File.Name <- Green$Sample.File.Name[,drop = TRUE]




# Step 3: calculate relative peak area for each frag in each sample-------


# make a function -- assumes the sample name column is Sample.File.Name
relative.abundance <- function(input) {
  datain <- input
  dataout <- rbind()
  SampleNames <- unique(datain$Sample.File.Name)
  for (i in seq(unique(datain$Sample.File.Name))) {
    x.i <- SampleNames[i]
    Sample.i <- datain[datain$Sample.File.Name == x.i,]
    sum.i <- sum(Sample.i$Area.in.Point)
    Sample.i$Relative.Area <- (Sample.i$Area.in.Point)/sum.i*100
    dataout <- rbind(dataout, Sample.i)
  }
  dataout
}

# Green.3 is new dataframe with relative abundance column 
Green.3 <- relative.abundance(Green)

# check that rel ab adds to 100 for each sample
tapply(Green.3$Relative.Area, Green.3$Sample.File.Name, sum)


# Step 4: removing unrepeatable peaks ----- 


# remove frags with relative area < 1%

# save the stuff to be deleated for summary log
Green.tinyarea <- Green.3[Green.3$Relative.Area < 1, ] 
Green.tinyarea $Frag.Quality <- rep("tiny_area", length(Green.tinyarea$Size))
Green.tinyarea [9] <- list(NULL)

# keep stuff with areas bigger than or = to 1%
Green.3.2 <- Green.3[Green.3$Relative.Area >= 1,] # 735 obs.

# repeat realtive area calcs with remaining fragments
Green.4 <- relative.abundance(Green.3.2)

# check that rel ab adds to 100 for each sample
tapply(Green.4$Relative.Area, Green.4$Sample.File.Name, sum)   


#Creates a final dataframe of good Samples
Green.good <- Green.4
Green.good$Frag.Quality <- rep("Good", length(Green.good$Size))
Green.good [9] <- list(NULL)

# Dataframe with all samples merged back together with Quality tag at end to specify where it came from
Green_Quality_Master <- rbind(Green.good, Green.missing, Green.tinyarea, Green.toobig, Green.toosmall)


# Steps 5-6: binning OTUs using a .25 difference from center of bin------

# set some initial stuff

# defining a new class for TRFLPs

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



# function to spit a data frame out of object class Bins
Bins.to.data.frame <- function(Bins) {
  result <- data.frame()
  for(b in Bins){
    # print(b@fragments$tag)
    b@fragments$Bin <- rep(b@mean)
    result <- rbind(result, b@fragments)
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

# testing

Green.bins1 <- bincimate(Green.4, 0.25)
Green.bins1.dataframe <- tagBins(Green.bins1 )

Green.bins.final <- Lump.bins(Green.bins1, 0.5)
Green.bins.final.tag <- tagBins(Green.bins.final)
Green.bins.final.dataframe <- Bins.to.data.frame(Green.bins.final.tag)
# length(df2$tag)
require(reshape2)
#Counts <- dcast(df2, Sample.File.Name ~ Bin, length)


# Step 7: OTU sample matix for specified Dye------

Green_Matrix <- dcast(Green.bins.final.dataframe, Sample.File.Name ~ Bin, value.var = "Relative.Area", sum, fill = 0, drop = F, margins = T)

# write.table(Green_Matrix, "Green_Matrix.txt", row.names=F)  # creates new .csv file with all the data combined

write.csv(Green_Matrix, "Green_Matrix.csv", row.names=F)  # creates new .csv file with all the data combined
