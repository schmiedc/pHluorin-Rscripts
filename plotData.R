library(ggplot2)
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
#               gridExtra: install.packages("gridExtra")
#
#      VERSION: 1.0.0
#      CREATED: 2018-05-24
#     REVISION: 2018-08-07
#
# ============================================================================
# plotting the raw values of all the experiments
plotRawMean <- function(dataTable, pagetitle){
  
  mean.signal <- subset(dataTable, variable == "mean" )
  
  count <- as.data.frame(table(mean.signal$name))
  
  mean.list <- list()
  plots <- list()
  
  for (names in count$Var1){
    
    one.table <- subset(mean.signal, name == names )
    
    mean.list[[names]] <- one.table
    
    plots[[names]] <- ggplot(data=one.table, aes(x=time, y=value, colour=variable, group=roi)) +
      geom_line() + 
      guides(colour=FALSE)  + 
      theme_light() +
      xlab("time (s)") + 
      ylab("Fluorescence intensity (a.u.)") + 
      #ylim(0, 800) +
      ggtitle(paste0("Raw data ", names)) 
  }
  plots <- grid.arrange(grobs = plots, ncol = 2, top = pagetitle)
  
  return (plots)
  
}

plotRawArea <- function(dataTable, label) {
  
  # create box plot the show the size distribution of the segmented objects
  area.signal <- subset(dataTable, variable == "area" )
  
  area.signal <- subset(area.signal, frame == "1" )
  
  areaplot <- ggplot(area.signal, aes(x=name, y=value)) +
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
    xlab("Movie name") + 
    ylab("Area (square micrometer)") +
    ggtitle(paste0("Area of ", label, " ROIs"))
  
  areaplot  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
}

# plotting the raw values of all the experiments
plotMean <- function(dataTable, column, title, sd){
  
  if (sd == TRUE){
    
    dataTable$high <- with(dataTable, dataTable$mean + dataTable$sd)
    dataTable$low <-  with(dataTable, dataTable$mean - dataTable$sd)
    plot <- ggplot(data=dataTable, aes(x=time, y=mean, colour=name, group=name, fill = name)) +
      geom_ribbon(aes(ymin = low, ymax = high, colour=name, group=name, fill = name ), alpha=.3) +
      geom_line(colour = "black") + 
      guides(colour=FALSE)  + 
      theme_light() +
      xlab("time (s)") + 
      ylab("Avg. Fluorescence intensity (a.u.)") + 
      #ylim(0, 450) +
      ggtitle(paste0("Average Grey value of ", title)) 
    
  } else if (sd == FALSE) {
    
    varname <- c(column)
    plot <- ggplot(data=dataTable, aes(x=time, y=get(varname[1]), colour=name, group=name, fill = name)) +
      geom_ribbon(aes(ymin = get(varname[1]), ymax = get(varname[1]), colour=name, group=name, fill = name ), alpha=.3) +
      geom_line() +
      guides(colour=FALSE)  + 
      theme_light() +
      xlab("time (s)") + 
      ylab("Avg. Fluorescence intensity (a.u.)") + 
      ggtitle(paste0("Average Grey value of active boutons ", title)) 
    
  }
  
  return (plot)
  
}

plotMeans <- function(avg.signal, avg.background, finalTable) {
  
  plot.list  <- list()
  plot.list[["Signal"]] <- plotMean(avg.signal, mean, "Signal", TRUE)
  plot.list[["Background"]] <- plotMean(avg.background, mean, "Background", TRUE)
  #plot.list[["Corrected"]] <- plotMean(finalTable, "mean.corr", "background subtracted", FALSE)
  #plot.list[["Surface Norm"]] <- plotMean(finalTable, "surf_norm", "surface normalized", FALSE)
  plot.list[["Peak Norm"]] <- plotMean(finalTable, "peak_norm", "peak normalized", FALSE)
  
  return(plot.list)
  
}
