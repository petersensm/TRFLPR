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

master[c("Size")][is.na(master[c("Size")])] <- 0

Blue <- subset(master, Dye == "B")

Green <- subset(master, Dye == "G")



# Step 2: set for threshold abnormally sized fragments
#   < 50 or > 600

#   how to?
#   use a subset fnc with logic 
#   maybe??? Data[Data$Size >= 50 & Data$Size <= 600, ]
#   output: thresholded data (one for each color)

#Subset for Numeric variables but you can also use factors
  Blue <- Blue[Blue$Size >= 50.0,]
  Blue <- Blue[Blue$Size <= 600.0,]


# Step 3: calculate relative peak area for each frag in each sample

#   how to?
#   some kind of loop that repeats for each sample
#   for each sample, calc sum of area of fragments, 
#     use to calculate % area of frag in the sample (loop in loop?) and put in new column labled relative area or something
#   do some internal check here, even if for our own benefit, to see that the relative area of frags in a sample = 100
#     sheryl: reshape?
#   output: dataframe with new column (relative peak area)

n <- 86
corrected1percent <- rep(NA,821)

SummedAbund <- tapply(Blue$Area.in.Point, Blue$Sample.File.Name, sum, na.rm = TRUE)
SummedAbund_dataframe <- as.data.frame(SummedAbund)
SummedAbund_dataframe <- na.omit(SummedAbund_dataframe)

NumSamp <- length(unique(Blue$Sample.File.Name))

SampleName <- unique(Blue$Sample.File.Name)

for (k in 1:NumSamp) {
  if SampleName[k]
}

  
for (i in 1:length(Blue$Area.in.Point)){
    
  if (Blue$Sample.File.Name[i] == "10e1s1.fsa")  Blue$corrected1percent[i] <- ((Blue$Area.in.Point[i])/(SummedAbund_dataframe$SummedAbund[1]))

}


# Step 4: removing unrepeatable peaks --- more thresholds! this time on relative area

#   how to?
#   subset to remove frags with relative area < 1%
#   maybe??? Data[Data$Rel >= 1, ]
#   then recalculate relative area to reflect change in sum of peak area
#     more loops?
#   output: dataframe with updated (or new?) relative area column (relative peak area)
#   stretch: would be nice if these steps gave reports about how much got removed

# Step 5: binning OTUs Part I

#   how to?
#   sort dataset by fragment size
#   then junk could maybe be all in one function or loop
#     new column for difference between one fragment size and the next
#     function to set breaks
#     new column for OTU and function (loop?) to label with OTU
#   output: dataframe with new columns (diff and OTU)

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


