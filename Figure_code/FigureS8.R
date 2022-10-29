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
library(cowplot)
library(ggpmisc)

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




# ============================================================================
# Load data

setwd(paste0(projectDir, "Data_Storage/Phenotyping_data/"))
file <- read_excel('Sexual_raw_data.xlsm', 1, col_names = FALSE)
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
treatKeyDT <- fread("treatment_key.txt", header = T)
treatKeyDT$chem <- gsub(".+(?=[A-Z]02)", "", treatKeyDT$Abr, perl = T)
indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")

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
colorDF <- data.frame(factors = plottingFactors2, colors = plottingColors, stringsAsFactors = FALSE)
myColors <- colorDF$colors
names(colPal) <- plottingFactors2

# ============================================================================
# Wrangle data to plot number of replicates vs mean genomewide LOD score 

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/")
popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
avgLODs <- popDT[, .(avgLOD = mean(LOD)), by = "Chemical"]
chemRepDT <- data.table(Chemical = unique(popDT$Chemical), replicates = c(7,8,10,9,10,6,5))
repDT <- merge(avgLODs, chemRepDT, by = "Chemical")
avgK <- popDT[, .(avgK = mean(k)), by = "Chemical"]
avgDT <- merge(repDT, avgK, by = "Chemical")

avgK2 <- popDT[, .(LOD= mean(LOD)), by = k][order(k)]
# ============================================================================
# Wrangle data to plot number of haplotypes contributing to a test vs per-site LOD score 



# ============================================================================
# Correlate mean OD630 with mean per-site heterozygosity for each treatment

gplot <- ggplot(data=repDT, aes(replicates, avgLOD)) + geom_smooth(method = "lm", se=FALSE) + stat_poly_eq(formula = y ~ x, aes(label = ..p.value.label..), parse = TRUE,label.x.npc = "left") + geom_point(colour = "black", size = 2) + xlab("number of replicates") + ylab("average -log10(p) score") + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), legend.position = "none") + scale_colour_manual(values = colPal)

gplot2 <- ggplot(data=avgK2, aes(k, LOD)) + geom_smooth(method = "lm", se=FALSE) + stat_poly_eq(formula = y ~ x, aes(label = ..p.value.label..), parse = TRUE,label.x.npc = "left") + geom_point(colour = "black", size = 2) + xlab("number of haplotypes") + ylab("average -log10(p) score") + theme_bw(base_size = 12) + theme(panel.grid = element_blank())

savePlot <- plot_grid(gplot, gplot2, align = 'vh', labels = c("A", "B"), ncol = 1, hjust = -1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS7.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS7.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350)


# ============================================================================
# Trouble-shooting
stat_regline_equation(aes(label = ..rr.label..))