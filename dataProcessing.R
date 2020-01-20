library(reshape2)
# plyr does not like dplyr!
library(plyr)
# ============================================================================
#
#  DESCRIPTION: Data analysis for Fíji pHlorin workflow
#              
#       AUTHOR: Christopher Schmied, 
#      CONTACT: schmied@dzne.de
#     INSITUTE: Leibniz-Forschungsinstitut f r Molekulare Pharmakologie (FMP)
#               Cellular Imaging - Core facility
#               Campus Berlin-Buch
#               Robert-Roessle-Str. 10
#               13125 Berlin, Germany
#
#         BUGS:
#        NOTES: 
# DEPENDENCIES: reshape2: install.packages("reshape2")
#               plyr: install.packages("plyr")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2018-08-07
#
# ============================================================================
# function that adds metadata as a column
read_table_filename <- function(filename){
  ret <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
  # extracts from filename the metadata and adds as column
  ret$name <- regmatches(basename(filename), regexpr('.+?(?=_Roi-)', basename(filename), perl=TRUE))
  ret$roi <- regmatches(basename(filename), regexpr('(?<=_Roi-)\\d+\\d+\\d+(?=_)', basename(filename), perl=TRUE))
  ret
}

# collect all lists with specified labels converts framerate to time with timeresolution
collectList <- function(inputdir, label, timeRes){
  
  match <- paste0("(\\.*)Roi-\\d+\\d+\\d+_",label,".csv")
  pHlorin.list <- list.files(inputdir, recursive=TRUE, pattern = match, full.names = TRUE )
  
  # lapply is executing a function over a vector
  # in this case reads all the files in the area.files vector
  # and writes them into a table
  pHlorin.table <- llply(pHlorin.list, read_table_filename)
  
  # now rbind is combining them all into one list
  combine.table1 <- do.call("rbind", pHlorin.table)

  
  # rorders columns
  combine.table2 <- combine.table1[,c(11,1,12,2,3,4,5,6,7,8,9,10)]
  
  # column headers to lower case
  combine.table2 <- rename(combine.table2, c("X" = "frame", 
                           "Area1" = "area", 
                           "Mean1" = "mean",
                           "StdDev1" = "stdDev",
                           "Mode1" = "mode",
                           "Min1" = "min",
                           "Max1" = "max",
                           "IntDen1" = "intDen",
                           "Median1" = "median",
                           "RawIntDen1" = "rawIntDen"))
  
  
  # converts string to numeric and converts frames to seconds
  combine.table2$frame <- as.numeric(as.character(combine.table2$frame))
  combine.table2$time <- combine.table2$frame - 1
  combine.table2$time <- combine.table2$time * timeRes
  

  # create long table
  melt.table <- melt(combine.table2, id.vars = c("name","frame","time", "roi"))
  
  return (melt.table)
  
}

# calculates average mean intensity per frame
calcMean <- function(inputTable) {
  
  mean.inputTable <- subset(inputTable, variable == "mean" )
  average.outputTable <- ddply(mean.inputTable, as.quoted(c("name", "frame", "time")), summarise, mean=mean(value), sd = sd(value))
  return (average.outputTable )
  
}

# calculates the number of boutons per movie
calcBoutonNumber <- function(inputTable){
  
  # calcluates boutonsNumber
  count <- as.data.frame(table(inputTable$name))
  
  boutons.list <- list()
  
  for (names in count$Var1){
    
    roi.table <-subset(inputTable, name == names & frame == 1 & variable == "area")
    roiNumber <- nrow(roi.table)
    boutons.list[[names]] <- roiNumber
    
  }
  
  # combine list of dataframes into one dataframe
  boutonsNumbers <- ldply(boutons.list, data.frame)
  
  # renames column x into time
  boutonsNumbers <- rename(boutonsNumbers, c(".id"="name"))
  boutonsNumbers <- rename(boutonsNumbers, c("X..i.."="number_segments"))
  
  return (boutonsNumbers)
  
}

# subtracts background mean intensity from signal mean intensity
backgroundSubtraction <- function(signalTable, backgroundTable){
  
  forBackgroundSubtraction <- merge(signalTable, backgroundTable, by=c("name", "frame", "time"), suffixes=c(".sig",".back"))
  
  # normalize mean signal with mean background intensity
  forBackgroundSubtraction$mean.corr <- forBackgroundSubtraction$mean.sig - forBackgroundSubtraction$mean.back
  
  return (forBackgroundSubtraction)
  
}

# normalizes the data to the initial signal before the stimulation
# normalizes the surface to 0
surfaceNormalisation <- function(inputTable, frameStim) {
  
  # surface normalization
  count <- as.data.frame(table(inputTable$name))
  
  surf.list <- list()
  
  for (names in count$Var1){
    
    newdata <- subset(inputTable, name == names)
    
    # get basedata > data before stiumlation
    basedata <- subset(inputTable, frame <= frameStim & name == names)
    
    # calculates mean intensities of basedata
    baseline <- colMeans( basedata[c("mean.corr")], na.rm = FALSE, dims = 1)
    
    # devides the data by the baseline
    newdata$surf_norm <- sapply(newdata$mean.corr, function(x){x / baseline})
    surf.list[[names]] <- newdata
    
  }
  
  # combine list of dataframes into one dataframe
  surf.table <- ldply(surf.list, data.frame)
  
  return (surf.table)
  
}


# peak normalization normalises the peak to 1
peakNormalisation <- function(inputTable) {
  
  count <- as.data.frame(table(inputTable$name))
  
  max.list <- list()
  
  for (names in count$Var1){
    
    maxdata <- subset(inputTable, name == names)
    
    # gets the maxima fo the data
    max = max(maxdata$surf_norm, na.rm = FALSE)
    
    # gets the normalisation factor
    max_norm = max - 1
    
    # divides the data by the max normalisation factor to set max to 1
    maxdata$peak_norm <- sapply(maxdata$surf_norm, function(x){(x - 1) /  max_norm })
    max.list[[names]] <- maxdata
    
  }
  
  # combine list of dataframes into one dataframe
  peak.table <- ldply(max.list, data.frame)
  
  return (peak.table)
  
}

# sort data per frame
sortFrames <- function(inputTable) {
  
  count <- as.data.frame(table(inputTable$name))
  
  sorted.list <- list()
  
  for (names in count$Var1){
    
    pername <- subset(inputTable, name == names)
    
    # sorts frames
    sorted <- pername[order(pername$frame),]
    
    sorted.list[[names]] <- sorted 
    
  }
  
  # combine list of dataframes into one dataframe
  sortedTable <- ldply(sorted.list, data.frame)
  
  return (sortedTable)
}

# ============================================================================
# generate table.signal and table.background via:
# further settings
# labelSignal = "signal"
# labelBackground = "background"
# table.signal <- collectList(indir, labelSignal, timeResolution)
# table.background <- collectList(indir, labelBackground, timeResolution)
# avg.signal <- calcMean(table.signal)
# avg.background <- calcMean(table.background)

processData <- function(inputdir, frameStimulation, avg.signal, avg.background ){
  
  # subtracts background mean intensity from signal mean intensity
  subtractedBackground <- backgroundSubtraction(avg.signal, avg.background)
  
  # normalizes the data to the initial signal before the stimulation
  # takes background subtracted data and the stimulation frame
  surfaceNormalized <- surfaceNormalisation(subtractedBackground, frameStimulation)
  
  # peak normalization
  peakNormalized <- peakNormalisation(surfaceNormalized)

  
  finalTable <- sortFrames(peakNormalized)
  
  return (finalTable)
  
}
