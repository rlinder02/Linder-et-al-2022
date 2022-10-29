##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2019
# Title:       Plotting raw haplotype frequencies 
# Description: Plotting haplotype frequencies of the base population, along with the treatment replicates and ypd-only replicates; also making a table of raw haplotype frequencies and heterozygosity 

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory
library(tictoc)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(readxl)
library(gprofiler2)
library(gridExtra)
library(egg)
library(scales)
library(dplyr)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggcorrplot)
library(ggbeeswarm)
library(lemon)

source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')

# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
options("scipen"=100, "digits"=4)
# ============================================================================
# Custom functions

add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha)) 
}

# ============================================================================
# Load data

plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
colPal <- c("dodgerblue2", "blue2", "skyblue3", "darkturquoise",
            "brown", "olivedrab", "tomato",
            "chocolate", "burlywood4",
            "saddlebrown",
            "magenta4", "maroon2",
            "mediumpurple",
            "darkcyan",
            "cornsilk4",
            "mediumvioletred")
plottingColors <- colPal
colorDT <- data.table(Chemical = plottingFactors2, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
names(colPal) <- plottingFactors2
names(colPal)[11] <- 'ypd'


# ============================================================================
# Panel plot of chemicals on y-axis and days/generations on the x-axis

# first Tuesday and Wednesday, then cycle on Thursday through each week
#myPal <- RColorBrewer::brewer.pal(7, "Dark2")
#myAlphaPal <- add.alpha(colPal, 0.25)
#alphaDT <- data.table(Chemical = colorDT$Chemical, colors = myAlphaPal)
chemicals <- c("cadmium chloride", "chlorpromazine", "diamide", "glacial acetic acid", "sodium chloride", "urea", "ypd")
spoM <- "lavender"
chemCode <- 1:7
nChems <- length(chemicals)
days <- c(1,1, rep(c(1, 3, 1, 1, 1), 12))
allDays <- sort(rep(cumsum(days), nChems))
intervals <- rep(c(2, 2, 1, 1, 1), 12)
intervalSum <- 2.5 + cumsum(intervals)
intervalSum <- intervalSum
centers <- sort(rep(c(0.5, 1.5, 2.5, intervalSum), nChems))
eachCenter <- length(unique(centers))
allChemicals <- rep(chemCode, eachCenter)
allChemicals2 <- rep(chemicals, eachCenter)
widths <- rep(c(1, 1, rep(c(1, 3, 1, 1, 1), 12), 1), each = nChems)


dt <- data.table(chemical = allChemicals, Chemical = allChemicals2, center = centers, widths = widths)
colorDT <- data.table(Chemical = names(colPal), colors = colPal)
dt2 <- merge(dt, colorDT, by = "Chemical")
dt2[, colors2 := rev(colors)]
dt2[widths == 3, colors2 := spoM]
findingYPD <- diff(dt2$widths)
ypdDays <- which(findingYPD == -2) + 1
dt2[ypdDays, colors2 := "magenta4"]

gplot <- qplot(center, chemical, fill = I(dt2$colors2), colour = I("black"), data = dt2, geom = "tile", width = widths) + theme_bw() + theme(panel.grid = element_blank(), legend.position = "none", axis.line = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank()) + ylab("") + xlab("") + scale_y_continuous(breaks = c(1,2,3,4,5,6,7), labels = rev(c(chemicals))) + theme(axis.text.y = element_text(hjust = 1))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Figure1.pdf"), gplot, width = 8.5, height = 9, units = "in")

## may have to do this manually in bash:  pdfcrop input.pdf output.pdf is the format; works really well!
system2(command = "pdfcrop", args = c("/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Figure1.pdf", "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Figure1_cropped.pdf"))