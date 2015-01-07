# TRFLP binning program in R - Demo with slightly tiddier 
# Joe informed me that this isn't really a demo...anyway working to replace with a "how to"
# FYI -- see TRFLP-functions 2014_FEB_17 for fixed lump bins function - KK's data helped us find a scenario we missed

# sheryl's stuff for starting ----
# Clear the field and Check what objects are currently in memory 
ls()
rm( list=ls() )   # clear all objects 
ls()

# set working directory
setwd("C:/Users/sherry/Documents/GitHub/TRFLPR/TRFLP_project_files")

# required packages
require(reshape2)

# Step 1: import data from peak scanner and clean up ---------
#   remove all dyes not of interest 
#   split out label for dye
# TODO: generalize this stuff

# read in data already analyzed in peak scanner
  #master <- read.csv ("CVNP_trflp_data.csv")
  master <- read.csv ("TRFLP Master_Data_Colin.csv")

# remove columns not needed
#   names(master)
#   master[17:26] <- list(NULL)
#   master [3:11] <- list(NULL)
  master.1 <- master[,c("Status", "Sample.Name", "Dye.Sample.Peak",  
                      "Sample.File.Name", "Size","Height", "Area.in.Point")]

# split out dye label
  master.1.2 <-  colsplit(master.1$Dye.Sample.Peak, ",", c('Dye', 'Sample.Peak'))
  master.2 <- cbind(master.1, master.1.2)
  master.2 [3] <- list(NULL) # purely asethetic

#reorder the columns for aesthetic purposes
  master.2 <- master.2[, c(1,2,7,8,3,4,5,6)]

#   str(master)
#   summary(master.2) 
  # fyi decided to remove NAs after subsetting 
  # for purposes of generating a summary file with an account of data quality

# subset dyes
    Blue <- subset(master.2, Dye == "B") 
   # Green <- subset(master, Dye == "G")

# Step 2: remove NAs and fragmements outside threshold --------
# threshold  < 50 or > 600
# save removed fragments for a quality control file
# fyi prob bad idea to keep recycling name!
# TODO: generalize/functionize this stuff

# remove NAs
  # save frags with NA's to another dataset for later
  Blue.missing <- Blue[is.na(Blue$Size), ] # 96 obs. # list of samples with no data for a future error/summary log
    Blue.missing$Frag.Quality <- rep("missing", length(Blue.missing$Size))
     # keep frags w/o na's
  Blue <- Blue[complete.cases(Blue$Size), ]  #good samples
  
# remove small frags
  Blue.toosmall <- Blue[Blue$Size < 50.0, ] #  obs. # list of frags too small  for a future error/summary log
    Blue.toosmall$Frag.Quality <- rep("too_small", length(Blue.toosmall$Size))
  Blue <- Blue[Blue$Size >= 50.0,] #  obs.
  
# remove big frags
  Blue.toobig <- Blue[Blue$Size > 600.0,] # obs. # list of frags too big  for a future error/summary log
    Blue.toobig$Frag.Quality <- rep("too_big", length(Blue.toobig$Size))
  Blue <- Blue[Blue$Size <= 600.0,] #  obs

# drops sample names if they have no analyzable blue frags
  Blue$Sample.File.Name <- Blue$Sample.File.Name[,drop = TRUE]

# Step 3: calculate relative peak area for each frag in each sample -------
#   for each sample, calc sum of area of fragments, 
#     use to calculate % area of frag in the sample (loop in loop?) and put in new column labled relative area or something
#   output: dataframe with new column (relative peak area)
# fyi - 11/21/2013 learned you could probably do this with plyr

# relative abundance fnxn assumes the sample name column is Sample.File.Name 
# and the peak area column is Area.in.Point
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

# fyi - 11/21/2013 learned you could probably do this with plyr
    # same deal via plyr TODO - double check and consider pros and cons
    Blue.3ply <- ddply(Blue, .(Sample.File.Name), 
                            mutate, 
                       Relative.Area = Area.in.Point/sum(Area.in.Point)*100)

# run relative abundance fnxn
    Blue.3 <- relative.abundance(Blue)
    # check that rel ab adds to 100 for each sample
    #tapply(Blue.3$Relative.Area, Blue.3$Sample.File.Name, sum)
    ddply(Blue.3, .(Sample.File.Name), summarize, Total.Rel.Area = sum(Relative.Area))

# Step 4: removing unrepeatable peaks --- peaks with relative area below a threshold of 1% -----
#   subset to remove frags with relative area < 1% -- save to a quality control file
#   then recalculate relative area to reflect change in sum of peak area

# remove frags with relative area <1
  # save the stuff to be deleated for summary log
    Blue.tinyarea <- Blue.3[Blue.3$Relative.Area < 1, ] # obs. # list of frags with small relatvive areas  for a future error/summary log
    Blue.tinyarea$Frag.Quality <- rep("tiny_area", length(Blue.tinyarea$Size))
    #Blue.tinyarea[9] <- list(NULL)
    Blue.tinyarea$Relative.Area <- NULL

# keep stuff with areas bigger than or = to 1
    Blue.3.2 <- Blue.3[Blue.3$Relative.Area >= 1,] # 735 obs.

# repeat realtive area calcs with remaining fragments
    Blue.4 <- relative.abundance(Blue.3.2)
    # check that rel ab adds to 100 for each sample
    #tapply(Blue.4$Relative.Area, Blue.4$Sample.File.Name, sum)  # looks good
    ddply(Blue.4, .(Sample.File.Name), summarize, Total.Rel.Area = sum(Relative.Area))
    # & matches excel output

#     example.sample <- Blue.3[Blue.3$Sample.File.Name == "14e1s1.fsa",]
#     example.sample.chopped <- Blue.4[Blue.4$Sample.File.Name == "14e1s1.fsa",]

# output "good" fragments that will be processed futher
  Blue.good <- Blue.4
  Blue.good$Frag.Quality <- rep("Good", length(Blue.good$Size))
  #Blue.good [9] <- list(NULL)
  Blue.good$Relative.Area <- NULL
  # might go back and add a blank relative area column to others!
 
# wrap up fragment quality report -- will need some work
  Blue_Quality_Master <- rbind(Blue.good, Blue.missing, Blue.tinyarea, Blue.toobig, Blue.toosmall)
  count(Blue_Quality_Master, "Frag.Quality")

# Steps 5-6: binning OTUs ----------------

# binning version 2 (.25 from center of bin) --  essentially a k means clustering
# todo - kmeans etc. other clustering considerations
# bincimate function currenlty takes product Blue.4 

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

# Run binning on blue fragments functions ---------------
# cluster fragments into bins separated by a minimum of 0.25 base pairs between outer fragments
Blue.bins1 <- bincimate(Blue.4, 0.25)
  # look at some bins for demo - in class bin form
  Blue.bins1[2]
  # in a dataframe form
  Blue.bins1.dataframe <- Bins.to.data.frame(tagBins(Blue.bins1))
  head(Blue.bins1.dataframe)
  tail(Blue.bins1.dataframe)

# consolidate bins whose fragment sizes are within .5 bp of center
Blue.bins.final <- Lump.bins(Blue.bins1, 0.5)

# label bins with an easy tag
Blue.bins.final.tag <- tagBins(Blue.bins.final)
# convert data from class bins to a dataframe
Blue.bins.final.dataframe <- Bins.to.data.frame(Blue.bins.final.tag)
  # take a look
  Blue.bins.final[26]
  head(Blue.bins.final.dataframe)
  tail(Blue.bins.final.dataframe)

# TODO: write fnxn to somehow show the merged bins or spit out a dataframe of the ones merged w/before and after lables

# Step 7: OTU sample matix
#   use reshape to take sum of relative area by sample and OTU
#   output: OTU matrix ready for community analysis

# TODO: round the Bins and relative areas before chucking into the matrix
# sample x bin mean
Blue_Matrix <- dcast(Blue.bins.final.dataframe, Sample.File.Name ~ Bin, value.var = "Relative.Area", sum, fill = 0, drop = F, margins = T)
# sample x arbitrary tag
Blue_Matrix <- dcast(Blue.bins.final.dataframe, Sample.File.Name ~ tag, value.var = "Relative.Area", sum, fill = 0, drop = F, margins = T)
