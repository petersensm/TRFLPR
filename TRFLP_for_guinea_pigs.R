# TRFLP R code for working from the functions in the package
# v 0.3
# mod smp feb 21 2014
# last mod smp mar 14 2014

# blue

# Setup  -----------
# install and load TRFLPR package if not already done
    # 1. put zip file of package on computer
    # 2. then run this line of code:
      # fill the ...with the details for the path to the file
      install.packages("~/TRFLP/TRFLPR_0.1.zip", repos = NULL)
# try this
      install.packages("C:/.../TRFLP/TRFLPR_0.1.zip", repos = NULL)
    # 2.v2. OR click this stuff in RStudio:
      # click on tools menu above
      # select "install packages..."
      # for Install from: select "Package Archive file (.zip; .tar.gz)
      # for Package archive: select "browse.." and navigate to the zip file

# load the TRFLPR package
    # via code:
    library(TRFLPR)
    # fyi, it will also load reshape2
    # OR click packages tab in lower left pane of RStudio, scroll down to TRFLPR, click the check box

    # Since "TRFLPR" did not load in our first roll out to guinea pigs, source in the functions
    # this wouldn't work on CH's computer
    # fill the ...with the details for the path to the file
    source("C:/.../TRFLPR_functions.r") 
    source("C:/Users/sherry/Documents/GitHub/TRFLPR/TRFLPR_functions.r")

# install and load required packages, if not already done
install.packages(reshape2)
install.packages(plyr)
require(reshape2)
require(plyr)

# Step 1: import data from peak scanner and clean up ---------
#   remove all dyes not of interest 
#   split out label for dye

    # read in data already analyzed in peak scanner
    master <- read.csv ("C:/Users/sherry/Documents/GitHub/TRFLPR/TRFLP_project_files/trflp_data_KK_140213.csv")
    
    # remove columns not needed, but keep these
    master.1 <- master[,c("Status", "Sample.Name", "Dye.Sample.Peak",  
                          "Sample.File.Name", "Size","Height", "Area.in.Point")]
    
    # split out dye label, requires reshape2 package
    master.1.2 <-  colsplit(master.1$Dye.Sample.Peak, ",", c('Dye', 'Sample.Peak'))
    master.2 <- cbind(master.1, master.1.2)
    master.2 [3] <- list(NULL) # purely asethetic
    
    #reorder the columns for aesthetic purposes
    master.2 <- master.2[, c(1,2,7,8,3,4,5,6)]
    
    # subset dyes
    Blue <- subset(master.2, Dye == "B") 
    #Green <- subset(master.2, Dye == "G")

# Step 2: remove NAs and fragmements outside threshold --------
# threshold  < 50 or > 600
# save removed fragments for a quality control file

    # remove NAs
    # save frags with NA's to another dataset for later
    Blue.missing <- Blue[is.na(Blue$Size), ] # 96 obs. # list of samples with no data for a future error/summary log
    Blue.missing$Frag.Quality <- rep("missing", length(Blue.missing$Size))
    Blue.missing$Relative.Area <- rep(NA)
    # keep frags w/o na's
    Blue <- Blue[complete.cases(Blue$Size), ]  #good samples
    # levels(Blue$Sample.File.Name)
    
    # remove small frags
    Blue.toosmall <- Blue[Blue$Size < 50.0, ] #  obs. # list of frags too small  for a future error/summary log
    Blue.toosmall$Frag.Quality <- rep("too_small", length(Blue.toosmall$Size))
    Blue.toosmall$Relative.Area <- rep(NA)
    Blue <- Blue[Blue$Size >= 50.0,] #  obs.
    
    # remove big frags
    Blue.toobig <- Blue[Blue$Size > 600.0,] # obs. # list of frags too big  for a future error/summary log
    Blue.toobig$Frag.Quality <- rep("too_big", length(Blue.toobig$Size))
    # no too big frags for KKs data so tried to work around to get 10 variables
    #Blue.toobig$Relative.Area <- rep(NA)
    Blue.toobig$Relative.Area <- Blue.toobig$Frag.Quality
    Blue <- Blue[Blue$Size <= 600.0,] #  obs

    # other stuff we wanna remove
    Blue.removed <- Blue[Blue$Sample.File.Name == "blank1.fsa" | Blue$Sample.File.Name == "blank3.fsa",]
    Blue.removed$Frag.Quality <- rep("manually removed", length(Blue.removed$Size))
    Blue.removed$Relative.Area <- rep(NA)
    Blue <- Blue[Blue$Sample.File.Name != "blank1.fsa" & Blue$Sample.File.Name != "blank3.fsa",]
    
    # drops sample names if they have no analyzable blue frags
    Blue$Sample.File.Name <- Blue$Sample.File.Name[,drop = TRUE]
    # check
#     levels(Blue$Sample.File.Name)

# Step 3: calculate relative peak area for each frag in each sample -------
#   for each sample, calc sum of area of fragments, 
#     use to calculate % area of frag in the sample (loop in loop?) and put in new column labled relative area or something
#   output: dataframe with new column (relative peak area)

    # call relative abundance fnxn 
    # assumes the sample name column is Sample.File.Name 
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

#     # alt relative abund via plyr 
#     TODO - double check and consider pros and cons
#     Blue.3ply <- ddply(Blue, .(Sample.File.Name), 
#                        mutate, 
#                        Relative.Area = Area.in.Point/sum(Area.in.Point)*100)

    # run relative abundance fnxn
    Blue.3 <- relative.abundance(Blue)
    # check that rel ab adds to 100 for each sample
    ddply(Blue.3, .(Sample.File.Name), summarize, Total.Rel.Area = sum(Relative.Area))
    
#### we're getting two out of 5 blanks -- found each had 1 blue peak : ()
## we tried running with, but let's remove them in step 2 and rerun stuff above

# Step 4: removing unrepeatable peaks --- peaks with relative area below a threshold of 1% -----
#   subset to remove frags with relative area < 1% -- save to a quality control file
#   then recalculate relative area to reflect change in sum of peak area

    # remove frags with relative area <1
    # save the stuff to be deleated for summary log
    Blue.tinyarea <- Blue.3[Blue.3$Relative.Area < 1, ] # obs. # list of frags with small relatvive areas  for a future error/summary log
    Blue.tinyarea$Frag.Quality <- rep("tiny_area", length(Blue.tinyarea$Size))
    #Blue.tinyarea[9] <- list(NULL)

    
    # keep stuff with areas bigger than or = to 1
    Blue.3.2 <- Blue.3[Blue.3$Relative.Area >= 1,] # 
    
    # repeat realtive area calcs with remaining fragments
    Blue.4 <- relative.abundance(Blue.3.2)
    # check that rel ab adds to 100 for each sample
    ddply(Blue.4, .(Sample.File.Name), summarize, Total.Rel.Area = sum(Relative.Area))

    # output "good" fragments that will be processed futher
    Blue.good <- Blue.4
    Blue.good$Frag.Quality <- rep("Good", length(Blue.good$Size))

### check here
    # wrap up fragment quality report
    Blue_Quality_Master <- rbind(Blue.good, Blue.missing, Blue.toosmall, Blue.toobig, Blue.removed, Blue.tinyarea)
    Blue_Fragment_Quality_summary <- dcast(Blue_Quality_Master, 
                           Sample.File.Name ~ Frag.Quality, length, drop = F, margins = T) 
    # TODO would like to reorder and if a category has nothing would like to repeat zeroes
    # intersting, so b/c i didn't drop until the end, the samples in the data set that don't generate blue frags at all are here too?
    # 3 blanks and two samples all showjust zeroes

# Steps 5-6: binning OTUs ----------------
# binning (.25 from center of bin) --  essentially a k means clustering
# bincimate function currenlty takes product Blue.4 

    # cluster fragments into bins separated by a minimum of 0.25 base pairs between outer fragments
    Blue.bins1 <- bincimate(Blue.4, 0.25)
    # consolidate bins whose fragment sizes are within .5 bp of center
    Blue.bins.final <- Lump.bins(Blue.bins1, 0.5) 

    # label bins with an easy tag
    Blue.bins.final.tag <- tagBins(Blue.bins.final)
    # convert data from class bins to a dataframe
    Blue.bins.final.dataframe <- Bins.to.data.frame(Blue.bins.final.tag)
    
    # create a bin label that inculdes the dye name
    # to do: concatenate labels -- easy fix for now, but make part of the tagging or converting to dataframe, might want some other label!
        # rounds the bin mean to 3 digits
        # see ?paste if you want something other than spaces btwn parts of the label 
    Blue.bins.final.dataframe$Bin_label <- paste(Blue.bins.final.dataframe$Dye, 
                                                 Blue.bins.final.dataframe$tag,
                                                 round(Blue.bins.final.dataframe$Bin_mean, 3))

# Step 7: OTU sample matix
#   use reshape to take sum of relative area by sample and OTU
#   output: OTU matrix ready for community analysis

# these two versions of the same output -- choose the one you want and comment out the rest!
    # sample x bin  (with margin totals)
    Blue_Matrix <- dcast(Blue.bins.final.dataframe, 
                         Sample.File.Name ~ Bin_label, 
                         value.var = "Relative.Area", sum, 
                         fill = 0, drop = F, margins = T)
    
    # sample x  (without margin totals)
    Blue_Matrix <- dcast(Blue.bins.final.dataframe, 
                         Sample.File.Name ~ Bin_label, 
                         value.var = "Relative.Area", sum, 
                         fill = 0, drop = F, margins = F)

    # to get a count of samples in each bin:
    # you can use just to confirm that there are no duplicate fragments in a bin.  Should only be 1s and 0s!!!
    Blue_MatrixCounts <- dcast(Blue.bins.final.dataframe, 
                     Sample.File.Name ~ Bin_label, 
                     value.var = "Relative.Area", length, 
                     fill = 0, drop = F, margins = F)

#write to csv for export
write.csv(Blue_Matrix,file="TRFLPblue140122.csv", row.names=F)

#### For combining matricies of two TRF dyes or restriction enzymes on same samples into one dataset:

###run to compare two dataframes for "incomparables"
length(intersect(unique(Blue_Matrix$Sample.File.Name), unique(Green_Matrix$Sample.File.Name)))


#Merge datasets, the all.x command keeps all rows instead of deleteing incomparable.
MergeTrials<-merge(Blue_Matrix, Green_Matrix)


write.csv(MergeTrials, file= "TRFLPcombinedBlue_Green140314.csv")

