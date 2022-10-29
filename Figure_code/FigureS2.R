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
library(cowplot)

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

offsets <- fread("newoffsets.txt", header = T)
g_l <- posoff_chr.bounds(offsets)   
snpCovDT <- fread("SNPtable.April23.sort_coverage_restructured.txt", header = TRUE)
repeats <- fread(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/repeatpos2.txt"), header = T)
repeats[,cp := paste(C, p, sep = "_")]
plottingBaseCovDT <- covcalc_create.chr.cov.dts(window = 2000, samplePattern = "^BAS02$")
generateBasePlottingCovsDT <- plottingBaseCovDT(as.data.frame(snpCovDT))
plottingCovsDTs <- covcalc_create.chr.cov.dts(window = 2000, samplePattern = "^SEE")
generatePlottingCovsDTs <- plottingCovsDTs(as.data.frame(snpCovDT))
savingPlottingCovsDTs <- writing_write.lists.of.dts.to.txt.files(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Genomewide_plotting_coverages_2kb_window/", analysisType = "plotting_coverage", idCol = "id")
writePlottingBaseCovFile <- savingPlottingCovsDTs(generateBasePlottingCovsDT)
writePlottingCovsFiles <- savingPlottingCovsDTs(generatePlottingCovsDTs)

gwideCovsDT <- covcalc_create.gwide.cov.dts(snpCovDT, samplePattern = "^SEE")
savingGwideCovDT <- writing_write.file(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Genomewide_avg_coverages/", fileName = "avg_gwide_coverage")
writeGwideCovsFiles <- savingGwideCovDT(gwideCovsDT)

## Generate a look-up table of low-coverage samples to filter for downstream analyses; also generate table of coverages for samples using 
indHapDiffsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12B02")
popDT <- do.call(rbind, indHapDiffsDTs)
gwideCovDTs <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Genomewide_avg_coverages/avg_gwide_coverage_DT.txt"), header = TRUE)
samplesUsing55 <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)

covCutoff <- covcalc_filter.low.cov.samples(covCutoff = 5)
lowCovSamples <- covCutoff(gwideCovDTs)

samplesUsing <- gwideCovDTs[!id %in% lowCovSamples$id & grepl("^SEE12B02", id)]
meanCov <- samplesUsing[, .(meanCov = mean(avgCov))]
rangeCov <- samplesUsing[, .(rangeCov = range(avgCov))]

samples55 <- gwideCovDTs[id %in% samplesUsing55$id]
meanCov55 <- samples55[, .(meanCov = mean(avgCov))]
rangeCov55 <- samples55[, .(rangeCov = range(avgCov))]

#meanCov <- samplesUsing[, .(meanCov = mean(avgCov))]


# ============================================================================
# Histogram of coverages genomewide

covPlotAll <- ggplot() + geom_histogram(data=samplesUsing, aes(x = avgCov), binwidth = 10) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_trans(x="log2") + scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50, 100, 200, 300, 600, 900), limits = c(5, 900)) + scale_y_continuous(name = "counts of all populations", limits = c(0, 50))
covPlot55 <- ggplot() + geom_histogram(data=samples55, aes(x = avgCov), binwidth = 10) + theme_bw(base_size = 10) + theme(panel.grid = element_blank()) + coord_trans(x="log2") + scale_x_continuous(breaks = c(5, 10, 20, 30, 40, 50, 100, 200, 300, 600, 900), limits = c(5, 900)) + scale_y_continuous(name = "counts of the 55 selected populations", limits = c(0, 50))


savePlot <- ggarrange(covPlotAll, covPlot55, align = "hv", labels = c("A", "B"), ncol = 1)

ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS2_gwide_cov_histAB.pdf"), savePlot, width = 8, height = 9, units = "in")
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/FigureS2_gwide_cov_histAB.png"), savePlot, width = 8, height = 9, units = "in", dpi = 350) 


# ============================================================================
# Trouble-shooting
