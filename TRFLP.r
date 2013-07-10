# TRFLP binning program in R

# Step 1: import data from peak scanner ---------
#   remove all dyes not of interest 
#   split process for each dye

#   how to? 
#   check what columns are labeled from out of peak scanner
#   split column with dye/sample into 2 columns, so we can subset by dye
#   use subset command to pull out blue and green seperately
#   output: 2 dataframes: one blue and one green

library(reshape2)

master <- read.csv ("CVNP_trflp_data.csv")

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
  # keep stuff with areas bigger than or = to 1
    Blue.3.2 <- Blue.3[Blue.3$Relative.Area >= 1,] # 735 obs.

# repeat realtive area calcs with remaining fragments
    Blue.4 <- relative.abundance(Blue.3.2)
    # check that rel ab adds to 100 for each sample
    tapply(Blue.4$Relative.Area, Blue.4$Sample.File.Name, sum)  # looks good 
    # & matches excel output

    example.sample <- Blue.3[Blue.3$Sample.File.Name == "14e1s1.fsa",]
    example.sample.chopped <- Blue.4[Blue.4$Sample.File.Name == "14e1s1.fsa",]


# Step 5: binning OTUs Part I

#   how to?
#   sort dataset by fragment size
#   I wonder if there isn't something that can automatically bin for us
#   then junk could maybe be all in one function or loop
#     new column for difference between one fragment size and the next
#     function to set breaks
#     new column for OTU and function (loop?) to label with OTU (would be ideal to label w/median frag size)
#   output: dataframe with new columns (diff and OTU)

Blue.4[1:10,]
# sort the data
Sorted.data <- Blue.4[order(Blue.4$Size),]
Sorted.data[1:10,]

# 1 gives first row in sorted data
Sorted.data[1,]
# make bin column
Sorted.data$Bin <- rep("NA")

# set bin to 1
bin <- 1
# set first bin to bin
Sorted.data[1, "Bin" ] = bin
# check
Sorted.data[1, "Bin" ]

# set size in first row to previous
previous <- Sorted.data[1, "Size" ]
previous

Sorted.data[1, "Size" ]

row.names(Sorted.data)
seq(row.names(Sorted.data))

# not sure how to start at second row....ignore for now!
for (i in row.names(Sorted.data)) {
  if(abs(Sorted.data[i, "Size" ]) - previous > 0.25) {
    Sorted.data[i, "Bin" ] = bin + 1
                                                     } else {
                                                       Sorted.data[i, "Bin" ] = bin
                                                     }
previous = Sorted.data[i, "Size" ]
}

# Q how to get bin to increment
abs(Sorted.data[i, "Size" ]) - previous


# Step 6: checking binning
# will need some summary tables of OTU's and either a way to manually fix problems
# count of each OTU per sample

#   how to?
#   count of otu per sample -- could be done with reshape using cast formula 
#     something like: counts <- dcast(long_data, sample ~ OTU) # require(reshape2), uses count by default
#   some flag to check if the counts are > 1 for an OTU
#     either a loop that checks using logic if the contents of each OTU column are >=2
#     could convert to long format (melt using reshape), use logical subset to find OTU's with >=2
#   if there is something wrong, then what???? -- might have to troubleshoot manually ---
#     could use subset to pull out ptoblem OTUs and either maunally or via R go through items in step6b.
#     replace problem OTU(s) with newly munged OTU(s)
#   output: hopefully a fixed dataframe with only one of each OTU per sample

# Step 7: check of binning against median to see if there are some long strung out OTUs
# not sure what to do if there are long strung out OTUs b/c there aren't natural breaks

#   how to?
#   calculate median frag size for each OTU
#   pull out min and max
#   calc diff from median for min and max frag in otu 
#   set some flag for otus that are >< .25 (maybe even more) from median
#   might just have something that pulls these out for the user to look at
#   output ?????

# Step 8: OTU sample matix

#   how to?
#   use reshape to take sum of relative area by sample and OTU
#     something like: Blue_Matrix <- dcast(Data, sample ~ OTU, sum, fill = 0, drop = F)
#   output: OTU matrix ready for community analysis


