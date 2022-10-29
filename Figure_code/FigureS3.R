##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2020
# Title:       Genomewide coverage
# Description: Calculates the mean genomewide coverage per sample, excluding the mitochondria and repetitive regions, and stores this as a table. Also makes separate plots for each sample. Plots are iterated through chromosome by chromosome in a specified window.

##############################################################################

# ============================================================================
# Load packages and sourced files
## Sourced files are kept in the default working directory	

library(data.table)
library(tictoc)
library(tidyverse)
library(scales)
source('plotting/Plotting_genomewide_statistics.R')
source('seq/Position_offsetter.R')
source('seq/Coverage_calculator.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')

# ============================================================================
# Set global options
defDir <- getwd()
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

options(digits = 10)

# ============================================================================
# Custom functions

# ============================================================================
# Run analyses
## First calculate coverages genomewide and write the tables to disk.
indHapDiffsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12B02")
popDT <- do.call(rbind, indHapDiffsDTs)
freqs <- popDT[, tstrsplit(collapsedFreqs, split = ";", type.convert = TRUE, fixed = TRUE)]
collHaps <- popDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
allHaps <- cbind(popDT$chr, popDT$pos, popDT$gp, popDT$Chemical, collHaps, freqs)
freqsDF <- as.data.frame(freqs)
newFreqs <- lapply(freqsDF, function(y) gsub("^-\\d.*|^\\d.*", 1, y))
newFreqsDF <- as.data.frame(do.call(cbind, newFreqs))
newFreqsDF2 <- sapply(newFreqsDF, as.numeric)
allHaps$reps <- rowSums(newFreqsDF2, na.rm = T)
chroms <- rep(unlist(allHaps[,1]), allHaps$reps) 
positions <- rep(unlist(allHaps[,2]), allHaps$reps) 
gpositions <- rep(unlist(allHaps[,3]), allHaps$reps)
chems <- rep(unlist(allHaps[,4]), allHaps$reps)
allfreqs <- as.numeric(unlist(t(freqs)))
allfreqs2 <- na.omit(allfreqs)
hapVec <- as.character(unlist(t(collHaps)))
hapVec2 <- na.omit(hapVec)
idxDT <- data.table(chr = chroms, pos = positions, gp = gpositions, Chemical = chems, haps = hapVec2, freqs = allfreqs2)
ab3Freqs <- idxDT[haps == "AB3"]
ab3FreqDifs <- ab3Freqs[, diffs := freqs - shift(freqs), by = Chemical]
allAB3Diffs <- na.omit(ab3FreqDifs$diffs)
avgDiffs <- mean(abs(allAB3Diffs))

# ============================================================================
# Histogram of coverages genomewide

diffPlot <- ggplot() + geom_histogram(data=ab3FreqDifs , aes(x = diffs), bins = 500) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_cartesian(xlim = c(-0.1, 0.1)) + scale_y_continuous(labels = comma)
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS3_AB3_diffs_hist.pdf"), diffPlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS3_AB3_diffs_hist.png"), diffPlot, width = 8, height = 9, units = "in", dpi = 350) 


# ============================================================================
# Trouble-shooting
