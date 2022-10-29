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
library(tidyverse)
library(DescTools)
library(lemon)
library(abind)
library(cowplot)


source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')

# ============================================================================
# Set global options

defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
options("scipen"=999, "digits" = 4)
# ============================================================================
# Custom functions

find.offending.character <- function(x, maxStringLength=256){  
    print(x)
    for (c in 1:maxStringLength){
        offendingChar <- substr(x,c,c)
        #print(offendingChar) #uncomment if you want the indiv characters printed
        #the next character is the offending multibyte Character
    }    
}


# ============================================================================
# Load data

setwd("/Users/robertlinder/Dropbox/Long_lab/SEE01/Daily_bottlenecks_experiment/")
#file <- read_excel('Sexual_raw_data.xlsm', 1, col_names = FALSE)
file <- fread("Daily_bottleneck_spreadsheet.txt")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
treatKeyDT <- fread("treatment_key.txt", header = T)
treatKeyDT$chem <- gsub(".+(?=[A-Z]02)", "", treatKeyDT$Abr, perl = T)
plottingFactors <- readLines("all_chems_list.txt")
plottingFactorsDT <- data.table(plottingFactors[2:length(plottingFactors)])
setnames(plottingFactorsDT, old = "V1", new = "Chemical")
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
colorDF <- data.frame(factors = plottingFactorsDT$Chemical, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
names(colPal) <- plottingFactorsDT$Chemical

# ============================================================================
# Wrangle data to plot OD630 vs cell counts on Thursday at week 12 for 4 YPD-only controls

file <- na.omit(file)
file[Day == "Monday_a", Day := "Monday\n post sporulation"]
file[Day == "Monday_d", Day := "Monday\n pre mating"]
file[Day == "Monday_e", Day := "Monday\n post mating"]
file <- file[Day != "Monday_c"]
file$Day <- ordered(file$Day, levels=c("Monday\n post sporulation", "Monday\n pre mating", "Monday\n post mating", "Tuesday", "Wednesday", "Thursday", "Friday"))
data <- file[Sample == "YPD_R10"]
monTues <- log(data$Total_cell_number[5]/data$Total_cell_number[4])/log(2)
intoWed <- (140/500)*data$Total_cell_number[5]
tuesWed<- log(data$Total_cell_number[6]/intoWed)/log(2)
intoThurs <- data$Total_cell_number[6]/10
wedThurs <- log(data$Total_cell_number[1]/intoThurs)/log(2)
intoFri <- data$Total_cell_number[1]/10
thursFri <- log(data$Total_cell_number[7]/intoFri)/log(2)
doublingsDT <- data.table(Interval = c("M-T", "T-W", "W-Th", "Th-F", "F-M", "F-M"), Divisions = c(monTues, tuesWed, wedThurs, thursFri, 1, 1), Type = c(rep("mitosis", 5), "meiosis"))
doublingsDT$Interval <- ordered(doublingsDT$Interval, levels=c("M-T", "T-W", "W-Th", "Th-F", "F-M"))
doublingsDT$Total_cell_number <- data$Total_cell_number[c(5, 6, 1, 7, 2, 2)]
doublingsDT <- doublingsDT[-6][Interval == "F-M", Divisions := 2]
monday <- data[grepl("Monday", Day)][, Day := c("Tetrads", "Spores\n post disruption", "Mated diploids")]
monday$Day <- ordered(monday$Day, levels=c("Tetrads", "Spores\n post disruption", "Mated diploids"))

# ============================================================================
# Plot estimates of total cell numbers each day (panel B) with number of doublings between M-F (panel A)

gplotA <- ggplot(data=doublingsDT, aes(x=Interval)) + geom_bar(aes(y=Total_cell_number), stat="identity", alpha=0.5, fill = "steelblue") + theme_bw(base_size = 12) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + ylab("Cells in culture") + theme(panel.grid = element_blank()) + xlab("")

gplotB <- ggplot(data=doublingsDT, aes(x=Interval)) + geom_bar(aes(y=Divisions), stat="identity", alpha=0.5, fill = "steelblue") + scale_y_continuous(breaks = 0:5, labels = 0:5) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + xlab("")

gplotC <- ggplot(data=monday, aes(Day, Total_cell_number)) + geom_bar(stat="identity", fill="steelblue", alpha = 0.5) + scale_y_log10(limits = c(1, 10^8), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + ylab("Cells in culture") + xlab("") + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.7 )) + annotate("text", x=c(1,3), y=70000000, label= c("Monday start", "Monday end")) 

## for this plot, show on F-M (1 mitotic and 1 meiotic generation, different colors)

savePlot <- plot_grid(gplotA, NULL, gplotB, NULL, gplotC, align = 'hv', rel_heights = c(1, -.25, 1, -.25, 1), labels = c("A", "A", "B", "B", "C"), ncol = 1, hjust = -1)
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS1_daily_bottlenecks_wk11_v2.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS1_daily_bottlenecks_wk11_v2.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350) 
# ============================================================================
# Trouble-shooting