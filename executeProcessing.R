setwd("/data1/FMP_Docs/Repositories/scripts_FMP/pHluorinShiny/")
library(gridExtra)
source("dataProcessing.R")
source("saveData.R")
source("plotData.R")
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
# DEPENDENCIES: ggplot2: install.packages("ggplot2")
#               xlsx: install.packages("gxlsx")
#               reshape2: install.packages("reshape2")
#               plyr: install.packages("plyr")
#               gridExtra: install.packages("gridExtra")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2020-01-22
#
# ============================================================================
# where to get the files
indir = "/data1/FMP_Docs/Repositories/minimal-datasets_FMP/pHlorin_TS/TestOne/"

# where to save the data
outdir = "/data1/FMP_Docs/Repositories/minimal-datasets_FMP/pHlorin_TS/ROutput/"

# ============================================================================
# resultsname
resultname = "Test"

# Time resolution in seconds
timeResolution = 2

# when stimulation happens
# these frames are used for calcuating the mean intensity
# then this value is used for the surface normalization
frameStimulation = 5

# further settings
labelSignal = "Spot"
labelBackground = "background"

# get raw data
table.signal <- collectList(indir, labelSignal, timeResolution)
table.background <- collectList(indir, labelBackground, timeResolution)

# calculates average mean intensity per frame
avg.signal <- calcMean(table.signal)
avg.background <- calcMean(table.background)

# generate final table
finalTable <- processData(indir, frameStimulation, avg.signal, avg.background)

# save files
writeToCsv(outdir, resultname, table.signal, table.background, finalTable)
#writeToXlsx(outdir, resultname, table.signal, table.background, finalTable)

# ============================================================================
# Plotting
plot.list <- plotMeans(avg.signal, avg.background, finalTable)

grid.arrange(grobs = plot.list, ncol = 1, top = "Processing results")


# ============================================================================
# create result plots
# here you may change the titles that are used for the plots
# the titles are the strings between ""
# saving is manual
plotRawMean(table.signal, "Raw grey values of active boutons")
plotRawArea(table.signal, "active boutons")

# plots Raw grey values and area of background
plotRawMean(table.background, "Raw grey values of background")
plotRawArea(table.background, "background")

plotMean(finalTable, "surf_norm", "- background and surface normalized", FALSE)
plotMean(finalTable, "peak_norm", "- background, surface and peak normalized", FALSE)

# ===========================================================================

