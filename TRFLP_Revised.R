# TRFLP binning program in R

# sheryl's stuff for starting
# Clear the field and Check what objects are currently in memory 
ls()
rm( list=ls() )   # clear all objects 
ls()

# Step 1: import data from peak scanner and clean up ---------
#   remove all dyes not of interest 
#   split process for each dye

#   how to? 
#   check what columns are labeled from out of peak scanner
#   split column with dye/sample into 2 columns, so we can subset by dye
#   use subset command to pull out blue and green seperately
#   output: 2 dataframes: one blue and one green

library(reshape2)
setwd("C:/Users/sherry/Documents/GitHub/TRFLPR/TRFLP_project_files")
# read in data
#master <- read.csv ("CVNP_trflp_data.csv")
master <- read.csv ("TRFLP Master_Data_Colin.csv")

class(master)
str(master)

master[17:26] <- list(NULL)
master [3:11] <- list(NULL)

master_2 <-  colsplit(master$Dye.Sample.Peak, ",", c('Dye', 'Sample.Peak'))

master <- cbind(master, master_2)

master [3] <- list(NULL)

master <- master[, c(1,2,7,8,3,4,5,6)]

# so given that there are still NAs in the summs below, I think it is reading the 0 replacement as NA then? (NO)
#   realized later that it was because Sample.name is a factor and we removed some of the values without dropping their labels
#   so we needed to add a statement to drop unused levels -- I put into step 2
# we could remove NA's by asking for complete cases - but that would drop out samples 
# (sorry, this wasn't the problem, both ways of getting rid of NAs work and fail in the same ways!)

# master[c("Size")][is.na(master[c("Size")])] <- 0

str(master)
summary(master) # 363 NAs!

# decided to remove NAs after subsetting for purposes of generating a summary file with an account of data quality
# master3 <- master[complete.cases(master$Size), ] 
# summary(master3)
# # diff of 363 from original!
# 5554-5191 # 363

# subset dyes
    Blue <- subset(master, Dye == "B") # 952 obs
    
    Green <- subset(master, Dye == "G")

# make this a function....
# Step 2: set for threshold abnormally sized fragments --------
#   < 50 or > 600

#   how to?
#   use a subset fnc with logic 
#   maybe??? Data[Data$Size >= 50 & Data$Size <= 600, ]
#   output: thresholded data (one for each color)

# remove NAs
  # save frags with NA's to another dataset for later
  Blue.missing <- Blue[is.na(Blue$Size), ] # 96 obs. # list of samples with no data for a future error/summary log
    Blue.missing$Frag.Quality <- rep("missing", length(Blue.missing$Size))
     # keep frags w/o na's
  Blue <- Blue[complete.cases(Blue$Size), ]  # 856 obs. good samples
  
# Subset for Numeric variables but you can also use factors
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

# Step 3: calculate relative peak area for each frag in each sample

#   how to?
#   some kind of loop that repeats for each sample
#   for each sample, calc sum of area of fragments, 
#     use to calculate % area of frag in the sample (loop in loop?) and put in new column labled relative area or something
#   do some internal check here, even if for our own benefit, to see that the relative area of frags in a sample = 100
#     sheryl: reshape?
#   output: dataframe with new column (relative peak area)

#### this is all the stuff that we can later trash that went into getting step 3 done ---------------------
# n <- 86 # 86 ?
# skipped right to relative abundance
# corrected1percent <- rep(NA,821)

# # get sums of area under peak for all samples
#     SummedAbund <- tapply(Blue$Area.in.Point, Blue$Sample.File.Name, sum, na.rm = TRUE)
#     SummedAbund_dataframe <- as.data.frame(SummedAbund) 
# 
# # don't need now
#     # SummedAbund_dataframe <- na.omit(SummedAbund_dataframe)
#     
#     # str(Blue)
#     # 
#     # NumSamp <- length(unique(Blue$Sample.File.Name)) # 66 samples
#     # 
#     SampleNames <- unique(Blue$Sample.File.Name)
# 
# # troubleshooting phantom sample names -- we can remove this later
#     # # where are the NA's coming from?
#     # # I don't understand the tapply results!
#     # ?tapply
#     # # are they coming from samples with only one fragment? - NO
#     # Blue[Blue$Sample.File.Name == "11e1s1.fsa",] # no data AND I don't see this isn the file!
#     # Blue[Blue$Sample.File.Name == "10e2s2.fsa",] # data
#     # levels(Blue$Sample.File.Name) #78 - one is blank!
#     # unique(Blue$Sample.File.Name) #66
#     # # need to drop levels earlier?
#     # 
#     # Blue$Sample.File.Name <- Blue$Sample.File.Name[,drop = TRUE]
#     # levels(Blue$Sample.File.Name) 
#     # unique(Blue$Sample.File.Name)
#     # # now both 66
#     # # will do earlier after removing incomplete cases from master
#     # 
#     # # now just 66 observations!
#     # SummedAbund <- tapply(Blue$Area.in.Point, Blue$Sample.File.Name, sum, na.rm = TRUE)
#     # SummedAbund_dataframe <- as.data.frame(SummedAbund) 
# 
# # I used this + sequence term to   
# # for (i in 1:length(Blue$Area.in.Point)){
# #     
# #   if (Blue$Sample.File.Name[i] == "10e1s1.fsa")  Blue$corrected1percent[i] <- ((Blue$Area.in.Point[i])/(SummedAbund_dataframe$SummedAbund[1]))
# # 
# # }
# 
# 
# 
# # trying parts of function
# # sample1 <- Blue[Blue$Sample.File.Name == "10e1s1.fsa",]
# # sum1 <- sum(sample1$Area.in.Point)
# # sum1
# # sample1$Relative.Area <- sample1$Area.in.Point/sum1*100
# 
# # make new column for relative area
# 
#   Blue$Relative.Area <-rep(0,length(Blue$Area.in.Point))
# 


# # yahoo!
# output <- rbind()
# for (i in seq(unique(Blue$Sample.File.Name))) {
#   x.i <- SampleNames[i]
#   Sample.i <- Blue[Blue$Sample.File.Name == x.i,]
#   sum.i <- sum(Sample.i$Area.in.Point)
#   Sample.i$Relative.Area <- (Sample.i$Area.in.Point)/sum.i*100
#   output <- rbind(output, Sample.i)
#   }
# 
# tapply(output$Relative.Area, output$Sample.File.Name, sum) #Fine

# final product of step 3 ---------------
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

# try it out -- sweet!
    # keep labeling "Blue" for now - no (may have been problem with renaming, but haven't looked into it)
    Blue.3 <- relative.abundance(Blue)

    # check that rel ab adds to 100 for each sample
    tapply(Blue.3$Relative.Area, Blue.3$Sample.File.Name, sum)  # OK, lets not recycle Blue

  
# Step 4: removing unrepeatable peaks --- more thresholds! this time on relative area

#   how to?
#   subset to remove frags with relative area < 1%
#   maybe??? Data[Data$Rel >= 1, ]
#   then recalculate relative area to reflect change in sum of peak area
#     more loops?
#   output: dataframe with updated (or new?) relative area column (relative peak area)
#   stretch: would be nice if these steps gave reports about how much got removed

# remove frags with relative area <1
  # save the stuff to be deleated for summary log
    Blue.tinyarea <- Blue.3[Blue.3$Relative.Area < 1, ] # 86 obs. # list of frags with small relatvive areas  for a future error/summary log
    Blue.tinyarea $Frag.Quality <- rep("tiny_area", length(Blue.tinyarea$Size))
    Blue.tinyarea [9] <- list(NULL)

# keep stuff with areas bigger than or = to 1
    Blue.3.2 <- Blue.3[Blue.3$Relative.Area >= 1,] # 735 obs.

# repeat realtive area calcs with remaining fragments
    Blue.4 <- relative.abundance(Blue.3.2)
    # check that rel ab adds to 100 for each sample
    tapply(Blue.4$Relative.Area, Blue.4$Sample.File.Name, sum)  # looks good 
    # & matches excel output

    example.sample <- Blue.3[Blue.3$Sample.File.Name == "14e1s1.fsa",]
    example.sample.chopped <- Blue.4[Blue.4$Sample.File.Name == "14e1s1.fsa",]

Blue.good <- Blue.4
  Blue.good$Frag.Quality <- rep("Good", length(Blue.good$Size))
  Blue.good [9] <- list(NULL)

Blue_Quality_Master <- rbind(Blue.good, Blue.missing, Blue.tinyarea, Blue.toobig, Blue.toosmall)
# Steps 5-6: binning OTUs 

#   how to?
#   sort dataset by fragment size
#   I wonder if there isn't something that can automatically bin for us
#   then junk could maybe be all in one function or loop
#     new column for difference between one fragment size and the next
#     function to set breaks
#     new column for OTU and function (loop?) to label with OTU (would be ideal to label w/median frag size)
#   output: dataframe with new columns (diff and OTU)

# binning version 2 (.25 from center of bin) 
# set some initial stuff

# bincimate takes product blue.4

# functions put into package so far ---------------------
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

# function to check criteria for merging - internal function
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

# testing functions ---------------

Blue.bins1 <- bincimate(Blue.4, 0.25)
Blue.bins1.dataframe <- tagBins(Blue.bins1 )

Blue.bins.final <- Lump.bins(Blue.bins1, 0.5)
Blue.bins.final.tag <- tagBins(Blue.bins.final)
Blue.bins.final.dataframe <- Bins.to.data.frame(Blue.bins.final.tag)
# length(df2$tag)
require(reshape2)
#Counts <- dcast(df2, Sample.File.Name ~ Bin, length)


# Step 7: OTU sample matix

#   how to?
#   use reshape to take sum of relative area by sample and OTU
#     something like: Blue_Matrix <- dcast(Data, sample ~ OTU, sum, fill = 0, drop = F)
#   output: OTU matrix ready for community analysis

Blue_Matrix <- dcast(Blue.bins.final.dataframe, Sample.File.Name ~ Bin, value.var = "Relative.Area", sum, fill = 0, drop = F, margins = T)
