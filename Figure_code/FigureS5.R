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

setwd("/Volumes/GoogleDrive/My\ Drive/Papers/Sexual_Evolution_Paper_1/Supplement_figures_tables/")
tbleS1 <- fread('TableS2.txt', header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Phenotyping_data/"))
file <- read_excel('Sexual_raw_data.xlsm', 1, col_names = FALSE)
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))
allBind <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_All_Reps_Hets.txt"))
treatKeyDT <- fread("treatment_key.txt", header = T)
treatKeyDT$chem <- gsub(".+(?=[A-Z]02)", "", treatKeyDT$Abr, perl = T)


stand_curve <- read.table("/Users/robertlinder/Dropbox/Long_lab/AEE01/Standard_Curves_Revised_140uL.txt", header = T, as.is = T, sep = "\t")
AOC_sc <- stand_curve[c(32, 1:8),]
AOC_sc_model <- lm(Cell.Number ~ poly(OD630, 3, raw = TRUE), data = AOC_sc)

colPal <- c("dodgerblue2", "blue2", "skyblue3", "darkturquoise",
            "brown", "olivedrab", "tomato",
            "chocolate", "burlywood4",
            "saddlebrown",
            "magenta4", "maroon2",
            "mediumpurple",
            "darkcyan",
            "cornsilk4",
            "mediumvioletred")

names(colPal) <- plottingFactors
treatmentAbbrv <- fread("chemAbbrev.txt")
names(colPal) <- treatmentAbbrv$abbr
colorPal <- colPal[order(names(colPal))]
plottingFactors <- read.table("all_chems_list.txt", header = T, sep = "\t")
plottingFactors2 <- plottingFactors$Chemical
plottingFactors3 <- plottingFactors2[c(2,1,4,3,5,6,7,8,9,12,10,11,13,14,15,16)]
names(colorPal) <- plottingFactors3
# ============================================================================
# Wrangle data to plot OD630 vs heterozygosity at week 12 

file <- as.data.table(file)
#data <- file[1:(grep("wk12", ...1)[1]-1)]
file <- file[,-c(1,14,15,28:30)]
file <- file[-seq(10,nrow(file),10),]
file <- file[-seq(1,nrow(file),9),]
colnames(file) <- as.character(1:ncol(file))
odvals <- file[,1:12]
#odvals <- odvals[-seq(1,nrow(odvals),9),]
conds <- file[,13:24]
#conds <- conds[-seq(1,nrow(conds),9),]

data <- file[, .(OD630 = as.numeric(unlist(odvals)), id = unlist(conds))][, c("wk", "chem", "rep") := .(as.numeric(substr(id, 4, 5)), substr(id, 6, 10), substr(id, 14, 16))][grepl("SEE[0-9][0-9]B02", id)]
data[, chemical := lapply(.SD, function(x) {treatKeyDT$Chemical[which(treatKeyDT$chem == x)[1]] } ), .SDcols = "chem", by = seq_len(nrow(data))]
dataPlot <- data[wk <= 12]
dataPlot <- dataPlot[id != "SEE06B02CI600R06"]
dataPlot <- na.omit(dataPlot)
dataPlotCpy <- copy(dataPlot)
dataPlot2 <- dataPlot[, chemical := gsub(".+(?=[A-Z][A-Z])", "", chem, perl = T)]
avgData <- dataPlotCpy[, .(avgs = mean(OD630, na.rm = T)), by = .(chemical, wk)]
setnames(avgData, "avgs", "OD630")
avgTotal <- avgData[, .(OD630 = mean(OD630, na.rm = T)), by = chemical][, chemical := gsub("18way_", "", chemical)]
avgReps <- dataPlotCpy[, .(OD630 = mean(OD630, na.rm = T)), by = .(chemical, rep)]
avgReps[, chemical := gsub("18way_", "", chemical)]
allBind <- allBind[chemical != "18way_base_2"][order(chemical)]
allBind[, chemical := gsub("18way_|_12", "", chemical)]
setnames(allBind, old = "Replicate", new = "rep")
avgHets <- allBind[,.(avgHet = mean(het)), by = chemical]

odVsHets <- merge(avgReps, allBind[, c("het", "chemical", "rep")], by = c("chemical", "rep"))
odVsHetsAvg <- merge(avgHets, avgTotal, by = "chemical")
odVsHetsAvg[, chemical := gsub("_", " ", chemical)][chemical == "YPD", chemical := "ypd"]
reorder <- odVsHetsAvg$chemical[order(odVsHetsAvg$chemical)]
odVsHetsAvg[, chemical := factor(chemical, levels = reorder)]

names(tbleS1)[c(2:4)] <- paste0(names(tbleS1)[c(1:3)], "_est")
fracSexual <- tbleS1[, c("Drug", "Outbred sexual")]
fracSexual[, Drug := gsub('\\n.*', "", Drug)]
fracSexual[, Drug := gsub(' ', "_", Drug)]
setnames(fracSexual, old = "Drug", new = "chemical")
fracSexual[, chemical := tolower(chemical)]
fracSexual[chemical == "flucanozole", chemical := "fluconazole"][chemical == "ypd", chemical := "YPD"]
odVsClass <- merge(fracSexual, avgTotal, by = "chemical")
setnames(odVsClass, old = "Outbred sexual", new = "Outbred_sexual")
odVsClass[, chemical := gsub("_", " ", chemical)][chemical == "YPD", chemical := "ypd"]
odVsClass[, chemical := gsub("_", " ", chemical)]
odVsClass[, chemical := factor(chemical, levels = reorder)]


plottingFactors <- levels(odVsClass$chemical)

chemDT <- data.table(Chemical = plottingFactors)
fwrite(chemDT, "all_chems_list.txt", col.names = T)




# ============================================================================
# Correlate mean OD630 with mean per-site heterozygosity for each treatment

gplot <- ggplot(data=odVsHetsAvg, aes(OD630, avgHet)) + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(aes(label = ..rr.label..)) + geom_point(aes(colour = chemical), size = 2) + xlab("average within-treatment OD630") + ylab("average within-treatment heterozygosity") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + scale_colour_manual(values = colorPal)

gplot2 <- ggplot(data=odVsClass, aes(OD630, Outbred_sexual)) + geom_smooth(method = "lm", se=FALSE) + stat_regline_equation(aes(label = ..rr.label..)) + geom_point(aes(colour = chemical), size = 2) + xlab("average within-treatment OD630") + ylab("fraction classified as outbred sexual") + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + scale_colour_manual(values = colorPal)

savePlot <- ggarrange(gplot, gplot2, labels = c("A", "B"), align = "hv", ncol = 1, common.legend = TRUE, legend = "right")

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS5_07-22-2021.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS5_07-22-2021.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350)


# ============================================================================
# Trouble-shooting
align = 'hv'
